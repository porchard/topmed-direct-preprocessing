ROOT=/net/topmed11/working/porchard/direct-preprocessing
DATA=$(ROOT)/data
WORK=$(ROOT)/work
BIN=$(ROOT)/bin

ANALYSIS=$(WORK)/$@
SIF=singularity exec --bind /net $(ROOT)/general.sif

.PHONY: all

define NL


endef

data: singularity fasta-hg38 fasta-hg19 chain direct

singularity:
	singularity pull general.sif docker://porchard/general:20220406125608
	
fasta-hg38:
	mkdir -p $(DATA)/fasta/hg38
	cp /net/topmed10/working/porchard/rnaseq/data/fasta/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta $(DATA)/fasta/hg38/

fasta-hg19:
	mkdir -p $(DATA)/fasta/hg19
	cd $(DATA)/fasta/hg19 && wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz && zcat hg19.fa.gz > hg19.unsorted.fa
	cd $(DATA)/fasta/hg19 && $(SIF) python $(BIN)/sort-fasta.py hg19.unsorted.fa > hg19.fa && rm hg19.unsorted.fa && $(SIF) samtools faidx hg19.fa

direct:
	mkdir -p $(DATA)/$@
	# description of tables: https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-40569-3/MediaObjects/41467_2023_40569_MOESM3_ESM.pdf
	cd $(DATA)/$@ && wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-40569-3/MediaObjects/41467_2023_40569_MOESM4_ESM.txt -O Table-S1.txt
	cd $(DATA)/$@ && wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-40569-3/MediaObjects/41467_2023_40569_MOESM6_ESM.txt -O Table-S3.txt
	cd $(DATA)/$@ && wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-40569-3/MediaObjects/41467_2023_40569_MOESM10_ESM.txt -O Table-S7.txt
	cd $(DATA)/$@ && wget https://zenodo.org/records/7521410/files/Pvalues_nominal_trans_eQTLs_10e4_Genes_DIRECT.txt.gz
	
chain:
	mkdir -p $(DATA)/$@
	wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz --directory-prefix $(DATA)/$@/

lift-and-tabix-trans-10e4: ANALYSIS=$(WORK)/lift-and-tabix/trans-10e4
lift-and-tabix-trans-10e4:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --direct $(DATA)/direct/Pvalues_nominal_trans_eQTLs_10e4_Genes_DIRECT.txt.gz --hg38_fasta $(DATA)/fasta/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta --chain $(DATA)/chain/hg19ToHg38.over.chain.gz --results $(ANALYSIS)/results $(ROOT)/lift-and-tabix-trans-10e4.nf &

lift-and-tabix-trans-significant: ANALYSIS=$(WORK)/lift-and-tabix/trans-significant
lift-and-tabix-trans-significant:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --direct $(DATA)/direct/Table-S7.txt --hg19_fasta $(DATA)/fasta/hg19/hg19.fa --hg38_fasta $(DATA)/fasta/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta --chain $(DATA)/chain/hg19ToHg38.over.chain.gz --results $(ANALYSIS)/results $(ROOT)/lift-and-tabix-trans-significant.nf &

lift-and-tabix-cis-eqtl-significant: ANALYSIS=$(WORK)/lift-and-tabix/cis-eqtl-significant
lift-and-tabix-cis-eqtl-significant:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --direct $(DATA)/direct/Table-S1.txt --hg19_fasta $(DATA)/fasta/hg19/hg19.fa --hg38_fasta $(DATA)/fasta/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta --chain $(DATA)/chain/hg19ToHg38.over.chain.gz --results $(ANALYSIS)/results $(ROOT)/lift-and-tabix-cis-eqtl-significant.nf &
