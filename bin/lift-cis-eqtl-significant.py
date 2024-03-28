#!/usr/bin/env python
# coding: utf-8


import sys
import pandas as pd
import numpy as np
import variants
from Bio import SeqIO
import sys

DIRECT = sys.argv[1]
HG38_FASTA = sys.argv[2]
HG19_FASTA = sys.argv[3]
CHAIN = sys.argv[4]
PICARD_JAR = '/sw/picard/picard.jar'

# DIRECT = '/net/topmed11/working/porchard/direct-preprocessing/data/direct/Table-S1.txt'
# HG38_FASTA = '/net/topmed10/working/porchard/rnaseq/data/testfasta/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta'
# HG19_FASTA = '/net/topmed10/working/porchard/rnaseq/data/fasta-hg19/hg19.fa'
# CHAIN = '/net/topmed10/working/porchard/rnaseq/data/chain/hg19ToHg38.over.chain.gz'
# PICARD_JAR = '/net/snowwhite/home/porchard/sw/picard/picard.jar'


direct = pd.read_csv(DIRECT, sep='\t')
direct = direct[[i for i in direct.columns if 'caveman' not in i.lower()]]
#tmp = direct[direct.A2_FREQ!='ME']
#(tmp.A1_FREQ + tmp.A2_FREQ.astype(float)).min()
direct['MAF'] = np.minimum(direct.A1_FREQ, 1-direct.A1_FREQ)
direct = direct[['GeneID', 'chrSNP', 'SNPpos', 'SNPid', 'REF', 'ALT', 'MAF', 'Nominal_Pval', 'Slope', 'EmpiricalAdjustedPval', 'BetaAdjustedPval', 'DiscoveryOrder', 'PvalueOrder']]
direct.chrSNP = direct.chrSNP.astype(str).map(lambda x: 'X' if x == '23' else x)
direct.SNPid = direct.SNPid.str.replace("^23:", "X:", regex=True)

to_lift = direct[['chrSNP', 'SNPpos', 'SNPid', 'REF', 'ALT']].drop_duplicates()
to_lift.columns = ['chrom', 'pos', 'id', 'ref', 'alt']
to_lift.chrom = 'chr' + to_lift.chrom.astype(str)
assert(to_lift.id.value_counts().max() == 1)



hg19 = dict()

with open(HG19_FASTA, 'r') as f:
    for record in SeqIO.parse(f, 'fasta'):
        sys.stderr.write('Loaded sequence {}\n'.format(record.id))
        hg19[record.id] = record

sys.stderr.write('Finished loading sequences\n')

to_lift['hg19_ref'] = [hg19[chrom][pos-1] for chrom, pos in zip(to_lift.chrom, to_lift.pos)]
to_lift['hg19_ref'] = to_lift['hg19_ref'].str.upper()

snps = to_lift[(to_lift.ref.str.len() == 1) & (to_lift.alt.str.len() == 1)]
#snps[snps.ref != snps.hg19_ref].head()
#assert(all(snps.ref == snps.hg19_ref))
assert((snps.ref == snps.hg19_ref).mean() >= 0.999)

to_lift = to_lift[['chrom', 'pos', 'id', 'ref', 'alt']]


lifted, rejected = variants.lift_variants(to_lift, HG38_FASTA, CHAIN, PICARD_JAR)

lifted['topmed_id'] = lifted.chrom + '_' + lifted.pos.astype(str) + '_' + lifted.ref + '_' + lifted.alt

swapped_alleles = dict(zip(lifted.id, lifted.swapped_alleles))

old_id_to_new_id = dict(zip(lifted.id, lifted.topmed_id))
topmed_id_to_chrom = dict(zip(lifted.topmed_id, lifted.chrom))
topmed_id_to_pos = dict(zip(lifted.topmed_id, lifted.pos.astype(int)))

direct_hg38 = direct
direct_hg38 = direct_hg38[~direct_hg38.SNPid.isin(rejected.id)]
direct_hg38['Slope'] = direct_hg38.Slope * np.where(direct_hg38.SNPid.map(swapped_alleles), -1, 1)
direct_hg38.SNPid = direct_hg38.SNPid.map(old_id_to_new_id)

COLUMNS = direct_hg38.columns.to_list()
direct_hg38['#chrom'] = direct_hg38.SNPid.map(topmed_id_to_chrom)
direct_hg38['start'] = (direct_hg38.SNPid.map(topmed_id_to_pos) - 1).astype(str)
direct_hg38['end'] = direct_hg38.SNPid.map(topmed_id_to_pos).astype(str)
direct_hg38 = direct_hg38[['#chrom', 'start', 'end'] + COLUMNS]
direct_hg38.to_csv(sys.stdout, sep='\t', index=False)