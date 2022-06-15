#!/bin/bash
aptPath="refs/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/"
binPath="resources/apt_2.11.4_linux_64_bit_x86_binaries/bin/"
outDir="data/logR_BAF/"
CELlist="data/cel_files"
ct=$( basename $CELlist _cel_files )

NOW=$(date +%F-%T)
echo -e "${NOW}"

if ! [ -d "${outDir}" ]; then
	mkdir -p ${outDir}
fi

echo -e "Generate genotype calls"
# generate genotyping calls
${binPath}apt-probeset-genotype -c ${aptPath}GenomeWideSNP_6.cdf \
        -a birdseed \
        --read-models-birdseed ${aptPath}GenomeWideSNP_6.birdseed.models \
        --special-snps ${aptPath}GenomeWideSNP_6.specialSNPs \
        --out-dir ${outDir}raw/ \
        --cel-files ${CELlist}

NOW=$(date +%F-%T)
echo -e "${NOW}"

echo -e "Allele-specific signal extraction"
# allele specific signal extraction
${binPath}apt-probeset-summarize --cdf-file ${aptPath}GenomeWideSNP_6.cdf \
	--analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true \
	--target-sketch ${aptPath}hapmap.quant-norm.normalization-target.txt \
	--out-dir ${outDir}raw/ \
	--cel-files ${CELlist}

NOW=$(date +%F-%T)
echo -e "${NOW}"

echo -e "LLR and BAF"
# LRR and BAF calculation
${binPath}normalize_affy_geno_cluster.pl ${aptPath}hapmap.genocluster ${outDir}raw/quant-norm.pm-only.med-polish.expr.summary.txt \
	-locfile ${aptPath}affygw6.hg19.pfb \
	-out ${outDir}raw/lrr_baf.txt

NOW=$(date +%F-%T)
echo -e "${NOW}"

echo -e "Split samples"
${binPath}kcolumn.pl ${outDir}raw/lrr_baf.txt split 2 -tab -head 3 -name -out gw6

if ! [ -d "${outDir}sample_level_affy_output/" ]; then
        mkdir -p ${outDir}sample_level_affy_output/
fi

mv gw6* ${outDir}sample_level_affy_output/

NOW=$(date +%F-%T)
echo -e "${NOW}"

echo -e "Make ASCAT-ready input"
if ! [ -d "${outDir}sample_files/" ]; then
        mkdir -p ${outDir}sample_files/
fi
Rscript scripts/affy_to_ascat_format.R ${outDir}

echo -e "fileid\tlogr\tbaf\tpenalty" > data/cellLine_CEL_file_mapping.tsv
paste <(basename -a $(ls data/logR_BAF/sample_files/*_LogR.txt | sed 's/_tumour.*$//')) <(ls data/logR_BAF/sample_files/*_LogR.txt) | \
paste - <(ls data/logR_BAF/sample_files/*_BAF.txt) | sed 's/\(\S\+\)\t\(\S\+\)\t\(\S\+\)/\1\t\2\t\3\t70/' >> data/cellLine_CEL_file_mapping.tsv

NOW=$(date +%F-%T)
echo -e "${NOW}"
echo -e "Done"
