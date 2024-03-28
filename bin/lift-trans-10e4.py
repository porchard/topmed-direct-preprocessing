#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
import numpy as np
import variants
import sys

DIRECT = sys.argv[1]
HG38_FASTA = sys.argv[2]
CHAIN = sys.argv[3]
PICARD_JAR = '/sw/picard/picard.jar'

# DIRECT = '/net/topmed11/working/porchard/direct-preprocessing/data/direct/Pvalues_nominal_trans_eQTLs_10e4_Genes_DIRECT.txt.gz'
# HG38_FASTA = '/net/topmed10/working/porchard/rnaseq/data/testfasta/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta'
# CHAIN = '/net/topmed10/working/porchard/rnaseq/data/chain/hg19ToHg38.over.chain.gz'
# PICARD_JAR = '/net/snowwhite/home/porchard/sw/picard/picard.jar'

direct = pd.read_csv(DIRECT, sep='\t')
direct.GeneChr = direct.GeneChr.astype(str).map(lambda x: 'X' if x == '23' else x)
direct.SNPchr = direct.SNPchr.astype(str).map(lambda x: 'X' if x == '23' else x)
direct.SNPid = direct.SNPid.str.replace("^23:", "X:", regex=True)

to_lift = direct[['SNPchr', 'SNPposition', 'SNPid', 'REF', 'ALT']].drop_duplicates()
to_lift.columns = ['chrom', 'pos', 'id', 'ref', 'alt']
to_lift.chrom = 'chr' + to_lift.chrom.astype(str)
assert(to_lift.id.value_counts().max() == 1)


lifted, rejected = variants.lift_variants(to_lift, HG38_FASTA, CHAIN, PICARD_JAR)

lifted['topmed_id'] = lifted.chrom + '_' + lifted.pos.astype(str) + '_' + lifted.ref + '_' + lifted.alt

swapped_alleles = dict(zip(lifted.id, lifted.swapped_alleles))

old_id_to_new_id = dict(zip(lifted.id, lifted.topmed_id))
topmed_id_to_chrom = dict(zip(lifted.topmed_id, lifted.chrom))
topmed_id_to_pos = dict(zip(lifted.topmed_id, lifted.pos.astype(int)))

direct['MAF'] = np.minimum(direct.FreqREF, direct.FreqALT)
direct_hg38 = direct[['GeneID', 'SNPid', 'MAF', 'Pvalue', 'Slope', 'SE', 'AdjustedPvalue']]
direct_hg38 = direct_hg38[~direct_hg38.SNPid.isin(rejected.id)]
direct_hg38['Slope'] = direct_hg38.Slope * np.where(direct_hg38.SNPid.map(swapped_alleles), -1, 1)
direct_hg38.SNPid = direct_hg38.SNPid.map(old_id_to_new_id)

COLUMNS = direct_hg38.columns.to_list()
direct_hg38['#chrom'] = direct_hg38.SNPid.map(topmed_id_to_chrom)
direct_hg38['start'] = (direct_hg38.SNPid.map(topmed_id_to_pos) - 1).astype(str)
direct_hg38['end'] = direct_hg38.SNPid.map(topmed_id_to_pos).astype(str)
direct_hg38 = direct_hg38[['#chrom', 'start', 'end'] + COLUMNS]
direct_hg38.to_csv(sys.stdout, sep='\t', index=False)
