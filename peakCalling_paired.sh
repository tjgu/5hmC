#!/bin/bash

one=$1 # input fastq
ID=$2 # input fastq file ID
dirres=$3 # directory for saving the results
dirgenome=$4 # reference genome
two=$5 # for the second end of the paired-end samples

echo -e "inputFile=$one"
echo -e "ID=$ID"
echo -e "dirres=$dirres"
echo -e "dirgenome=$dirgenome"

echo -e "in the third submit3:\none=${one}"
echo -e "python version: "
python -V

module load samtools/1.10 bedtools/2.29.2

dirref=/ufrc/data/reference/icbr/human/ensembl_GRCh38/bwa/Homo_sapiens_assembly38.fasta
ctrlBam=/orange/icbrbi/tgu/Jinying_Zhao/control_merged_NovogeneInputs.bam.sor.bam
blist=/ufrc/data/reference/icbr/human/ensembl_GRCh38/bowtie2/hg38-blacklist.v2.bed.gz
if [[ ${one} == *Yerkes* ]]; then
	echo "escape Yerkes"
elif [[ ${one} == *Admera* ]]; then
	inputBam=${dirres}/${ID}_Admera/${ID}_b2_rmMultiQual.bam
        samtools sort -@ 4 -T ./temp $inputBam --output-fmt BAM -o ${inputBam}.sor.bam
        samtools index ${inputBam}.sor.bam
        bedtools intersect -v -abam ${inputBam}.sor.bam -b ${blist} >${inputBam}.sor.bam.fil.bam

        bedtools bamtobed -i ${inputBam}.sor.bam.fil.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/(m2+0.1)}' >${dirres}/${ID}_Admera/${ID}_b2_rmMultiQual.bam.sor.bam.fil.bam.pbc

        rawPath=$PATH
        export PATH=/home/tgu/.conda/envs/run_spp/bin:$rawPath
        bedtools bamtobed -i ${inputBam}.sor.bam.fil.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -nc >${inputBam}.sor.bam.fil.bam_bed.gz
        READ_LEN=$(head -n 100 <(zcat ${inputBam}.sor.bam.fil.bam_bed.gz) | awk 'function abs(v) {{return v < 0 ? -v : v}} BEGIN{{sum=0}} {{sum+=abs($3-$2)}} END{{print int(sum/NR)}}')

        EXCLUSION_RANGE_MIN=-500
        #EXCLUSION_RANGE_MAX=`expr ${READ_LEN} + 10`
        EXCLUSION_RANGE_MAX=$((READ_LEN + 10))
	echo -e "READ_LEN=$READ_LEN \t EXCLUSION_RANGE_MIN=$EXCLUSION_RANGE_MIN \t EXCLUSION_RANGE_MAX=$EXCLUSION_RANGE_MAX"
        CC_SCORES_FILE=${inputBam}.sor.bam.fil.bam.cc.qc
        CC_PLOT_FILE=${inputBam}.sor.bam.fil.bam.cc.plot.pdf
        #module load conda
        #source activate run_spp
	if [ -r $CC_SCORES_FILE ]; then
		rm $CC_SCORES_FILE $CC_PLOT_FILE -rf
	fi
	
        run_spp.R -c=${inputBam}.sor.bam.fil.bam -p=4 -filtchr=chrM -savp=${CC_PLOT_FILE} -out=${CC_SCORES_FILE} -x=${EXCLUSION_RANGE_MIN}:${EXCLUSION_RANGE_MAX}
	#conda deactivate
        export PATH=$rawPath
        fragSize=`awk '{print $3}' ${CC_SCORES_FILE}  | perl -ne '@line=split(/,/); print $line[0];'`
	
	rm ${inputBam}.sor.bam.fil.bam_bed.gz -rf

	module load macs/2.2.7.1
	mkdir -p ${dirres}/${ID}_Admera/macs2_broad_ctrl_pe_fil_dn_lam
	macs2 callpeak -t ${inputBam}.sor.bam.fil.bam -c $ctrlBam \
	--scale-to small \
	--nomodel  \
	--broad \
	--extsize ${fragSize}\
	--shift 0 \
 	-f BAMPE -g 2913022398 \
	-n ${ID}_q0.01 \
	--outdir ${dirres}/${ID}_Admera/macs2_broad_ctrl_pe_fil_dn_lam \
	--keep-dup all \
	-B \
	--SPMR \
	-q 1e-02

	macs2 callpeak -t ${inputBam}.sor.bam.fil.bam -c $ctrlBam \
	--scale-to small \
	--nomodel  \
	--broad \
	--extsize ${fragSize}\
	--shift 0 \
 	-f BAMPE -g 2913022398 \
	-n ${ID}_p0.05 \
	--outdir ${dirres}/${ID}_Admera/macs2_broad_ctrl_pe_fil_dn_lam \
	--keep-dup all \
	-B \
	--SPMR \
	-p 5e-02

elif [[ ${one} == *Novogene* ]]; then
	inputBam=${dirres}/${ID}_Novogene/${ID}_b2_rmMultiQual.bam
        samtools sort -@ 4 -T ./temp $inputBam --output-fmt BAM -o ${inputBam}.sor.bam
        samtools index ${inputBam}.sor.bam
        bedtools intersect -v -abam ${inputBam}.sor.bam -b ${blist} >${inputBam}.sor.bam.fil.bam

        bedtools bamtobed -i ${inputBam}.sor.bam.fil.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/(m2+0.1)}' >${dirres}/${ID}_Novogene/${ID}_b2_rmMultiQual.bam.sor.bam.fil.bam.pbc

        bedtools bamtobed -i ${inputBam}.sor.bam.fil.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -nc >${inputBam}.sor.bam.fil.bam_bed.gz
        READ_LEN=$(head -n 100 <(zcat ${inputBam}.sor.bam.fil.bam_bed.gz) | awk 'function abs(v) {{return v < 0 ? -v : v}} BEGIN{{sum=0}} {{sum+=abs($3-$2)}} END{{print int(sum/NR)}}')

        EXCLUSION_RANGE_MIN=-500
        #EXCLUSION_RANGE_MAX=`expr ${READ_LEN} + 10`
        EXCLUSION_RANGE_MAX=$((READ_LEN + 10))
        CC_SCORES_FILE=${inputBam}.sor.bam.fil.bam.cc.qc
        CC_PLOT_FILE=${inputBam}.sor.bam.fil.bam.cc.plot.pdf
        rawPath=$PATH
        export PATH=/home/tgu/.conda/envs/run_spp/bin:$rawPath
        #module load conda
        #source activate run_spp
	if [ -r $CC_SCORES_FILE ]; then
		rm $CC_SCORES_FILE $CC_PLOT_FILE -rf
	fi
        run_spp.R -c=${inputBam}.sor.bam.fil.bam -p=4 -filtchr=chrM -savp=${CC_PLOT_FILE} -out=${CC_SCORES_FILE} -x=${EXCLUSION_RANGE_MIN}:${EXCLUSION_RANGE_MAX}
	#conda deactivate
        export PATH=$rawPath
        fragSize=`awk '{print $3}' ${CC_SCORES_FILE}  | perl -ne '@line=split(/,/); print $line[0];'`

	rm ${inputBam}.sor.bam.fil.bam_bed.gz -rf	

	module load macs/2.2.7.1
	mkdir -p ${dirres}/${ID}_Novogene/macs2_broad_ctrl_pe_fil_dn_lam
	macs2 callpeak -t ${inputBam}.sor.bam.fil.bam -c $ctrlBam \
	--scale-to small \
	--nomodel  \
	--broad \
	--extsize ${fragSize}\
	--shift 0 \
 	-f BAMPE -g 2913022398 \
	-n ${ID}_q0.01 \
	--outdir ${dirres}/${ID}_Novogene/macs2_broad_ctrl_pe_fil_dn_lam \
	--keep-dup all \
	-B \
	--SPMR \
	-q 1e-02

	macs2 callpeak -t ${inputBam}.sor.bam.fil.bam -c $ctrlBam \
	--scale-to small \
	--nomodel  \
	--broad \
	--extsize ${fragSize}\
	--shift 0 \
 	-f BAMPE -g 2913022398 \
	-n ${ID}_p0.05 \
	--outdir ${dirres}/${ID}_Novogene/macs2_broad_ctrl_pe_fil_dn_lam \
	--keep-dup all \
	-B \
	--SPMR \
	-p 5e-02

fi
