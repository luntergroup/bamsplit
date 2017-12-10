"""Split a BAM file by read supporting variants in a VCF file
"""

import argparse
import os
import operator
import pysam as ps
import ssw

def is_empty(values):
    return len(values) == 0

def parse_region(region_str):
    if ':' not in region_str or region_str[-1] == ':':
        return region_str
    else:
        contig, rest = region_str.rsplit(':', 1)
        begin, end = rest.replace(',', '').split('-')
        return contig, int(begin), int(end)

def make_bam(path, bam_in):
    return ps.AlignmentFile(path, 'wb', template=bam_in)

def make_out_bams(out_dir, bam_in, ploidy):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if len(out_dir) > 0 and out_dir[-1] != '/':
        out_dir += '/'
    result = [make_bam(out_dir + "support_" + str(i) + ".bam", bam_in) for i in range(ploidy)]
    result.append(make_bam(out_dir + "unassigned.bam", bam_in))
    return result

def fetch_reference(ref, region):
    return ref.fetch(region[0], region[1], region[2])

def fetch_tabix(tabix, region):
    if isinstance(region, tuple):
        return tabix.fetch(region[0], region[1], region[2])
    else:
        return tabix.fetch(region)

def fetch_reads(bam, region):
    return fetch_tabix(bam, region)

def fetch_variants(vcf, region):
    return fetch_tabix(vcf, region)

def write_reads(reads, bam):
    for read in reads:
        bam.write(read)

def is_phased(site, sample):
    return site.samples[sample].phased

def get_phase_region(site, sample):
    phase_set_key = "PS"
    begin = site.pos - 1
    if phase_set_key in site.samples[sample]:
        phase_set = site.samples[sample][phase_set_key]
        if not phase_set == '.': # not missing
            begin = int(phase_set.split('_')[0]) # Some PS definitions use POS_REF_ALT format
    return site.contig, begin, begin + len(site.ref)

def overlaps(lhs_region, rhs_region):
    if lhs_region[0] == rhs_region[0]:
        overlap_size = min(lhs_region[-1], rhs_region[-1]) - max(lhs_region[1], rhs_region[1])
        return overlap_size > 0
    else:
        return False

def are_in_phase(site1, site2, sample):
    return is_phased(site1, sample) and is_phased(site2, sample)\
           and overlaps(get_phase_region(site1, sample), get_phase_region(site2, sample))

def get_encompassing_region(sites):
    return sites[0].contig, sites[0].pos - 1, sites[-1].pos - 1 + len(sites[-1].ref)

def calculate_min_pad(reads, site_region):
    if sum(0 if read.is_unmapped else 1 for read in reads) > 0:
        min_read_begin = min(read.reference_start for read in reads if not read.is_unmapped)
        max_read_end   = max(read.reference_end for read in reads if not read.is_unmapped)
        min_lhs_pad = site_region[1] - min_read_begin if site_region[1] > min_read_begin else 0
        min_rhs_pad = site_region[2] - max_read_end if site_region[2] > max_read_end else 0
        return max(min_lhs_pad, min_rhs_pad)
    else:
        return 0

def indel_size(site):
    return max(abs(len(site.ref) - len(alt)) for alt in site.alts)

def max_indel_size(sites):
    return max(indel_size(site) for site in sites)

def expand(region, n):
    return region[0], max(region[1] - n, 0), region[2] + n

def is_before(read, region):
    return read.reference_end <= region[1]

def is_after(read, region):
    return read.reference_start >= region[2]

def mapped_region(read):
    return read.reference_name, read.reference_start, read.reference_end

def get_genotype_region(phased_sites, pad):
    begin = max(phased_sites[0].pos - 1 - pad, 0)
    end = phased_sites[-1].pos - 1 + len(phased_sites[-1].ref) + pad
    return phased_sites[0].contig, begin, end

def make_haplotype(alleles, region, ref_seq):
    result = ""
    ref_idx = 0
    for allele in alleles:
        next_ref_idx = allele[0] - region[1]
        result += ref_seq[ref_idx : next_ref_idx] # ref sub-sequence before allele
        result += allele[-1] # allele sequence
        ref_idx = next_ref_idx + allele[1]
    result += ref_seq[ref_idx:] # ref sub-sequence after last allele
    return result

def make_genotype(phased_sites, sample, ref, pad = 0):
    region = get_genotype_region(phased_sites, pad)
    alleles = [(site.pos - 1, len(site.ref), site.samples[sample].alleles) for site in phased_sites]
    ref_seq = fetch_reference(ref, region)
    ploidy = min(len(a[2]) for a in alleles)
    genotype = []
    for i in range(ploidy):
        genotype.append(make_haplotype([(a[0], a[1], a[2][i]) for a in alleles], region, ref_seq))
    return genotype

def calculate_alignment_score(read, haplotype):
    aligner = ssw.Aligner()
    alignment = aligner.align(reference=haplotype, query=read)
    return alignment.score

def calculate_alignment_scores(read, genotype):
    return [calculate_alignment_score(read, haplotype) for haplotype in genotype]

def split(phased_sites, ref, bam_iter, sample, bams_out, last_read=None):
    if not is_empty(phased_sites):
        phase_region = get_encompassing_region(phased_sites)
        support_region = expand(phase_region, max_indel_size(phased_sites))
        if last_read and is_before(last_read, support_region):
            bams_out[-1].write(last_read)
            last_read = None
        if not last_read:
            for read in bam_iter:
                if read.is_unmapped or is_before(read, support_region):
                    bams_out[-1].write(read)
                else:
                    last_read = read
                    break
        if last_read and not is_after(last_read, support_region):
            reads = [last_read]
            last_read = None
            for read in bam_iter:
                if read.is_unmapped or not is_after(read, support_region):
                    reads.append(read)
                else:
                    last_read = read
                    break
            genotype = make_genotype(phased_sites, sample, ref, calculate_min_pad(reads, phase_region))
            if len(genotype) + 1 > len(bams_out):
                raise ValueError("Found genotype with copy number "
                                 + str(len(genotype)) + " but the given sample ploidy is "
                                 + str(len(bams_out) - 1))
            for read in reads:
                if read.is_unmapped or len(genotype) == 1:
                    bams_out[-1].write(read)
                else:
                    scores = calculate_alignment_scores(read.query_sequence, genotype)
                    max_score_idx, max_score = min(enumerate(scores), key=operator.itemgetter(1))
                    if scores.count(max_score) > 1:
                        bams_out[-1].write(read)
                    else:
                        bams_out[max_score_idx].write(read)
    return last_read

def split_contig(region, ref, bam_in, vcf, sample, bams_out):
    bam_iter = fetch_reads(bam_in, region)
    phased_sites = []
    last_read = None
    for site in fetch_variants(vcf, region):
        if not (is_empty(phased_sites) or are_in_phase(phased_sites[0], site, sample)):
            last_read = split(phased_sites, ref, bam_iter, sample, bams_out, last_read)
            phased_sites.clear()
        phased_sites.append(site)
    last_read = split(phased_sites, ref, bam_iter, sample, bams_out, last_read)
    if last_read:
        bams_out[-1].write(last_read)
    for read in bam_iter:
        bams_out[-1].write(read)

def run_bamsplit(ref, bam_in, vcf, sample, bams_out, region=None):
    if region:
        split_contig(region, ref, bam_in, vcf, sample, bams_out)
    else:
        for contig in vcf.header.contigs:
            split_contig(contig, ref, bam_in, vcf, sample, bams_out)

def main(options):
    if options.ploidy > 1:
        ref      = ps.FastaFile(options.ref)
        bam_in   = ps.AlignmentFile(options.reads)
        vcf      = ps.VariantFile(options.variants)
        bams_out = make_out_bams(options.out_dir, bam_in, options.ploidy)
        region = parse_region(options.region) if options.region else None
        if len(vcf.header.samples) == 1:
            run_bamsplit(ref, bam_in, vcf, vcf.header.samples[0], bams_out, region)
            for bam in bams_out:
                bam.close()
                ps.index(bam.filename)
        else:
            print("ERROR: Input VCF file must contain one sample")
    else:
        print("Nothing to do for sample ploidy < 2")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--ref', type=str,
                        required=True,
                        help='Reference genome used to call variants')
    parser.add_argument('-b', '--reads', type=str,
                        required=True,
                        help='')
    parser.add_argument('-v', '--variants', type=str,
                        required=True,
                        help='Phased variant calls to split reads')
    parser.add_argument('-o', '--out_dir', type=str,
                        required=True,
                        help='Directory to output split BAMs')
    parser.add_argument('--region', type=str,
                        required=False,
                        help='Region to split')
    parser.add_argument('--ploidy', type=int,
                        required=False,
                        default=2,
                        help='Ploidy of the sample')
    parsed, unparsed = parser.parse_known_args()
    main(parsed)
