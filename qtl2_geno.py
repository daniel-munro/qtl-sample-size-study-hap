import argparse
from bisect import bisect_left
import pandas as pd
from pathlib import Path
import pysam


def genotype_code(gt: tuple, founder: bool = False) -> str:
    if gt == (None, None):
        return '-'
    elif gt == (0, 0):
        return 'A'
    elif (gt[0] == 0 and gt[1] > 0) or (gt[0] > 0 and gt[1] == 0):
        return '-' if founder else 'H'
    elif gt[0] > 0 and gt[1] > 0:
        return 'B'
    else:
        raise ValueError(f'GT not recognized: {gt}')


def genetic_pos(chrmap: pd.DataFrame, pos: int) -> float:
    r = bisect_left(chrmap['pos'], pos)
    if r == len(chrmap['pos']):
        return chrmap['cm'][r - 1]
    elif chrmap['pos'][r] == pos or r == 0:
        return chrmap['cm'][r]
    else:
        # Interpolate the genetic position.
        p_lo = chrmap['pos'][r - 1]
        p_hi = chrmap['pos'][r]
        g_lo = chrmap['cm'][r - 1]
        g_hi = chrmap['cm'][r]
        rel = (pos - p_lo) / (p_hi - p_lo)
        return g_lo + rel * (g_hi - g_lo)


parser = argparse.ArgumentParser(description='Prepare R/qtl2 input files.')
parser.add_argument('--vcf', help='VCF file for the cohort.')
parser.add_argument('--founder-vcf', help='VCF file for the founder strains.')
parser.add_argument('--snps', help='File containing list of SNPs to keep, e.g. the GBS observed SNPs.')
parser.add_argument('--gmap-dir', help='Directory containing genetic map files.')
parser.add_argument('--geno-out', help='Output file name for cohort genotype table.')
parser.add_argument('--founder-out', help='Output file name for founder genotype table.')
parser.add_argument('--pmap-out', help='Output file name for physical map.')
parser.add_argument('--gmap-out', help='Output file name for genetic map.')
args = parser.parse_args()

ID_list = open(args.snps, 'r').read().splitlines()
IDs = set(ID_list)

maps = {}
for chrom in range(1, 21):
    filename = Path(args.gmap_dir) / f'MAP4chr{chrom}.txt'
    maps[chrom] = pd.read_table(filename, sep=' ', names=['pos', 'ratio', 'cm'])

vcf = pysam.VariantFile(args.vcf)
samples = list(vcf.header.samples)
genos = {}
refs = {}
for rec in vcf.fetch():
    ID = rec.id if rec.id is not None else f'{rec.contig}:{rec.pos}'
    if IDs is None or ID in IDs:
        gt = [s['GT'] for s in rec.samples.values()]
        labels = [genotype_code(g, founder=False) for g in gt]
        genos[ID] = labels
        refs[ID] = rec.ref

ID_list = [x for x in ID_list if x in genos.keys()]
IDs = set(ID_list)

vcf = pysam.VariantFile(args.founder_vcf)
strains = list(vcf.header.samples)
founder_genos = {}
ref_mismatch = 0
remove = set()
for rec in vcf.fetch():
    ID = rec.id if rec.id is not None else f'{rec.contig}:{rec.pos}'
    if IDs is None or ID in IDs:
        gt = [s['GT'] for s in rec.samples.values()]
        labels = [genotype_code(g, founder=True) for g in gt]
        founder_genos[ID] = labels
        # assert rec.ref == refs[ID]
        if rec.ref != refs[ID]:
            # print(f'Reference mismatch: {ID} {rec.ref} vs. {refs[ID]}')
            ref_mismatch += 1
            remove.add(ID)
            del genos[ID]
            del founder_genos[ID]

if ref_mismatch > 0:
    print(f'{ref_mismatch} SNPs removed due to reference mismatch.')
    ID_list = [ID for ID in ID_list if ID not in remove]
    print(f'{len(ID_list)} shared SNPs remaining.')

with open(args.geno_out, 'w') as out:
    out.write(f'id,{",".join(samples)}\n')
    for ID in ID_list:
        out.write(f'{ID},{",".join(genos[ID])}\n')

with open(args.founder_out, 'w') as out:
    out.write(f'id,{",".join(strains)}\n')
    for ID in ID_list:
        out.write(f'{ID},{",".join(founder_genos[ID])}\n')

with open(args.pmap_out, 'w') as out:
    out.write('marker,chr,pos\n')
    for ID in ID_list:
        chrom, pos = tuple(ID.replace('chr', '').split(':'))
        pos = int(pos) / 1e6  # Units are Mb.
        out.write(f'{ID},{chrom},{pos}\n')

with open(args.gmap_out, 'w') as out:
    out.write('marker,chr,pos\n')
    for ID in ID_list:
        chrom, pos = tuple(ID.replace('chr', '').split(':'))
        gpos = genetic_pos(maps[int(chrom)], int(pos))
        out.write(f'{ID},{chrom},{round(gpos, 6)}\n')
