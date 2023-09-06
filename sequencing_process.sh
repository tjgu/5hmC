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

dirref=/ufrc/data/reference/icbr/human/ensembl_GRCh38/bwa/Homo_sapiens_assembly38.fasta

echo -e "in the third submit3:\none=${one}"
if [[ ${one} == *Yerkes* ]]; then
	echo -e "start analyzing Yerkes files"
	mkdir -p ${dirres}/${ID}_Yerkes
	fastqc -t 8 ${one} -o ${dirres}/${ID}_Yerkes
	perl -ne 'if(/Filename<\/td><td>(.*)<\/td><\/tr><tr><td>File type.*Encoding<\/td><td>(.*)<\/td><\/tr><tr><td>Total Sequences<\/td><td>(\d+)<\/td><\/tr><tr><td>Sequences flagged as poor quality<\/td><td>(.*)<\/td><\/tr><tr><td>Sequence length<\/td><td>(\d+)<\/td><\/tr><tr><td>%GC<\/td><td>(\d+)/){print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t";} if(/alt=\"\[(\w+)\]\"\/><a href=\"\#M10\">Adapter Content<\/a>/){print $1,"\n"}' ${dirres}/${ID}_Yerkes/*001_fastqc.html >>$dirres/quality_matrix_sequencing_Yerkes.txt

	trimmomatic SE -threads 8 -phred33 ${one} ${dirres}/${ID}_Yerkes/${ID}_R1_trim.fq.gz ILLUMINACLIP:${dirres}/adapters.fa:1:22:11 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40 &>${dirres}/${ID}_Yerkes/${ID}_trimlog_Yerkes.txt

	fastqc -t 8 ${dirres}/${ID}_Yerkes/${ID}_R1_trim.fq.gz
	perl -ne 'if(/Filename<\/td><td>(.*)<\/td><\/tr><tr><td>File type.*Encoding<\/td><td>(.*)<\/td><\/tr><tr><td>Total Sequences<\/td><td>(\d+)<\/td><\/tr><tr><td>Sequences flagged as poor quality<\/td><td>(.*)<\/td><\/tr><tr><td>Sequence length<\/td><td>(.*)<\/td><\/tr><tr><td>%GC<\/td><td>(\d+)/){print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t";} if(/alt=\"\[(\w+)\]\"\/><a href=\"\#M10\">Adapter Content<\/a>/){print $1,"\n"}' ${dirres}/${ID}_Yerkes/*trim_fastqc.html >>$dirres/quality_matrix_sequencing_Yerkes.txt
	perl -ne 'if(/share\/Gu5hmC\/bams\/(\S+)_Yerkes\/\S+.gz /){print $1,"\t";} if(/Input Read Pairs: (\d+) Both Surviving: (\d+) \((.*)\%\) Forward Only Surviving: /){print $1, "\t", $2, "\t", $3, "\n"}elsif(/Input Reads: (\d+) Surviving: (\d+) \((.*)\%\) Dropped:/){print $1, "\t", $2, "\t", $3, "\n"}' ${dirres}/${ID}_Yerkes/${ID}_trimlog_Yerkes.txt >>$dirres/quality_matrix_trim_Yerkes.txt

	bowtie2 -X2000 --mm --threads 8 -x $dirgenome $dirres/${ID}_Yerkes/${ID}_R1_trim.fq.gz 2>${dirres}/${ID}_Yerkes/${ID}_b2log_Yerkes.txt | samtools view -bt ${dirref}.fai -o $dirres/${ID}_Yerkes/${ID}_b2.bam -

	echo -ne "${ID}\t" >>$dirres/quality_matrix_b2_Yerkes.txt
	perl -ne 'if(/^(\d+) reads; of these/){print $1,"\t";} if(/\((.*%)\) aligned 0 times/){print $1,"\t";} if(/\((.*%)\) aligned exactly 1 time/){print $1,"\t";} if(/\((.*%)\) aligned >1 times/){print $1,"\t";} if(/^(.*%) overall alignment rate/){print $1,"\n";}' ${dirres}/${ID}_Yerkes/${ID}_b2log_Yerkes.txt >>$dirres/quality_matrix_b2_Yerkes.txt

	java -Djava.io.tmpdir=`pwd`/temp -jar  $HPC_PICARD_DIR/picard.jar SortSam I=$dirres/${ID}_Yerkes/${ID}_b2.bam O=$dirres/${ID}_Yerkes/${ID}_b2.sor.bam SO=queryname
        java -jar $HPC_PICARD_DIR/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$dirres/${ID}_Yerkes/${ID}_b2.sor.bam O=$dirres/${ID}_Yerkes/${ID}_b2.dedup.sor.bam METRICS_FILE=$dirres/${ID}_Yerkes/${ID}_b2.mkdup.matrics.txt
	samtools view -b -q 10 ${dirres}/${ID}_Yerkes/${ID}_b2.dedup.sor.bam >${dirres}/${ID}_Yerkes/${ID}_b2_rmMultiQual.bam
	module load rseqc/3.0.0
	/apps/rseqc/3.0.0/bin/bam_stat.py -i ${dirres}/${ID}_Yerkes/${ID}_b2_rmMultiQual.bam -q 4 &>${dirres}/${ID}_Yerkes/${ID}_b2_rmMultiQual.bam.stat.txt
	/apps/rseqc/3.0.0/bin/bam_stat.py -i ${dirres}/${ID}_Yerkes/${ID}_b2.bam -q 4 &>${dirres}/${ID}_Yerkes/${ID}_b2.bam.stat.txt
	echo -ne "${ID}_rawBam\t" >>$dirres/quality_matrix_b2_bamStat_Yerkes.txt
	perl -e '$TotalR=0; $Unique=0; $read1=0; $read2=0; $strand1=0; $strand2=0; $nonSplice=0; $mappedPair=0; $unMapped=0; $nonUnique=0; $PairedChr=0; $nonPri=0; $splice=0; $dups=0; while(<>){if(/Total records:\s+(\d+)/){$TotalR=$1;} if(/mapq >= mapq_cut \(unique\):\s+(\d+)/){$Unique=$1;} if(/Read-1:\s+(\d+)/){$read1=$1} if(/Read-2:\s+(\d+)/){$read2=$1} if(/Reads map to .\+.:\s+(\d+)/){$strand1=$1} if(/Reads map to .\-.:\s+(\d+)/){$strand2=$1;} if(/Non-splice reads:\s+(\d+)/){$nonSplice=$1} if(/Reads mapped in proper pairs:\s+(\d+)/){$mappedPair=$1} if(/Unmapped reads:\s+(\d+)/){$unMapped=$1} if(/mapq < mapq_cut \(non-unique\):\s+(\d+)/){$nonUnique=$1} if(/Proper-paired reads map to different chrom:(\d+)/){$PairedChr=$1} if(/Non primary hits\s+(\d+)/){$nonPri=$1} if(/Splice reads:\s+(\d+)/){$splice=$1} if(/Optical\/PCR duplicate:\s+(\d+)/){$dups=$1}} print $TotalR,"\t", $Unique,"\t", $read1,"\t", $read2,"\t", $strand1,"\t", $strand2,"\t", $nonSplice,"\t", $mappedPair,"\t", $unMapped,"\t", $nonUnique,"\t", $PairedChr,"\t", $nonPri,"\t", $splice,"\t", $dups,"\n"' ${dirres}/${ID}_Yerkes/${ID}_b2.bam.stat.txt >>$dirres/quality_matrix_b2_bamStat_Yerkes.txt
	echo -ne "${ID}_filBam\t" >>$dirres/quality_matrix_b2_bamStat_Yerkes.txt
	perl -e '$TotalR=0; $Unique=0; $read1=0; $read2=0; $strand1=0; $strand2=0; $nonSplice=0; $mappedPair=0; $unMapped=0; $nonUnique=0; $PairedChr=0; $nonPri=0; $splice=0; $dups=0; while(<>){if(/Total records:\s+(\d+)/){$TotalR=$1;} if(/mapq >= mapq_cut \(unique\):\s+(\d+)/){$Unique=$1;} if(/Read-1:\s+(\d+)/){$read1=$1} if(/Read-2:\s+(\d+)/){$read2=$1} if(/Reads map to .\+.:\s+(\d+)/){$strand1=$1} if(/Reads map to .\-.:\s+(\d+)/){$strand2=$1;} if(/Non-splice reads:\s+(\d+)/){$nonSplice=$1} if(/Reads mapped in proper pairs:\s+(\d+)/){$mappedPair=$1} if(/Unmapped reads:\s+(\d+)/){$unMapped=$1} if(/mapq < mapq_cut \(non-unique\):\s+(\d+)/){$nonUnique=$1} if(/Proper-paired reads map to different chrom:(\d+)/){$PairedChr=$1} if(/Non primary hits\s+(\d+)/){$nonPri=$1} if(/Splice reads:\s+(\d+)/){$splice=$1} if(/Optical\/PCR duplicate:\s+(\d+)/){$dups=$1}} print $TotalR,"\t", $Unique,"\t", $read1,"\t", $read2,"\t", $strand1,"\t", $strand2,"\t", $nonSplice,"\t", $mappedPair,"\t", $unMapped,"\t", $nonUnique,"\t", $PairedChr,"\t", $nonPri,"\t", $splice,"\t", $dups,"\n"' ${dirres}/${ID}_Yerkes/${ID}_b2_rmMultiQual.bam.stat.txt >>$dirres/quality_matrix_b2_bamStat_Yerkes.txt
	rm $dirres/${ID}_Yerkes/${ID}_b2.sor.bam
	rm $dirres/${ID}_Yerkes/${ID}*gz

elif [[ ${one} == *Novogene* ]]; then
	mkdir -p ${dirres}/${ID}_Novogene
	fastqc -t 8 ${one} -o ${dirres}/${ID}_Novogene
	fastqc -t 8 ${two} -o ${dirres}/${ID}_Novogene
	perl -ne 'if(/Filename<\/td><td>(.*)<\/td><\/tr><tr><td>File type.*Encoding<\/td><td>(.*)<\/td><\/tr><tr><td>Total Sequences<\/td><td>(\d+)<\/td><\/tr><tr><td>Sequences flagged as poor quality<\/td><td>(.*)<\/td><\/tr><tr><td>Sequence length<\/td><td>(\d+)<\/td><\/tr><tr><td>%GC<\/td><td>(\d+)/){print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t";} if(/alt=\"\[(\w+)\]\"\/><a href=\"\#M10\">Adapter Content<\/a>/){print $1,"\n"}' ${dirres}/${ID}_Novogene/*_?_fastqc.html >>$dirres/quality_matrix_sequencing_Novogene.txt

	trimmomatic PE -threads 8 -phred33 ${one} ${two} ${dirres}/${ID}_Novogene/${ID}_paired_1.fq.gz ${dirres}/${ID}_Novogene/${ID}_unpaired_1.fq.gz ${dirres}/${ID}_Novogene/${ID}_paired_2.fq.gz ${dirres}/${ID}_Novogene/${ID}_unpaired_2.fq.gz ILLUMINACLIP:${dirres}/adapters.fa:1:22:11 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40 &>${dirres}/${ID}_Novogene/${ID}_trimlog_Novogene.txt
	fastqc -t 8 ${dirres}/${ID}_Novogene/${ID}_paired_1.fq.gz
	fastqc -t 8 ${dirres}/${ID}_Novogene/${ID}_paired_2.fq.gz
	perl -ne 'if(/Filename<\/td><td>(.*)<\/td><\/tr><tr><td>File type.*Encoding<\/td><td>(.*)<\/td><\/tr><tr><td>Total Sequences<\/td><td>(\d+)<\/td><\/tr><tr><td>Sequences flagged as poor quality<\/td><td>(.*)<\/td><\/tr><tr><td>Sequence length<\/td><td>(.*)<\/td><\/tr><tr><td>%GC<\/td><td>(\d+)/){print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t";} if(/alt=\"\[(\w+)\]\"\/><a href=\"\#M10\">Adapter Content<\/a>/){print $1,"\n"}' ${dirres}/${ID}_Novogene/*paired_?_fastqc.html >>$dirres/quality_matrix_sequencing_Novogene.txt
	perl -ne 'if(/share\/Gu5hmC\/bams\/(\S+)_Novogene\/\S+.gz /){print $1,"\t";} if(/Input Read Pairs: (\d+) Both Surviving: (\d+) \((.*)\%\) Forward Only Surviving: /){print $1, "\t", $2, "\t", $3, "\n"}elsif(/Input Reads: (\d+) Surviving: (\d+) \((.*)\%\) Dropped:/){print $1, "\t", $2, "\t", $3, "\n"}' ${dirres}/${ID}_Novogene/${ID}_trimlog_Novogene.txt >>$dirres/quality_matrix_trim_Novogene.txt

	bowtie2 -X2000 --mm --threads 8 -x $dirgenome -1 $dirres/${ID}_Novogene/${ID}_paired_1.fq.gz -2 $dirres/${ID}_Novogene/${ID}_paired_2.fq.gz 2>${dirres}/${ID}_Novogene/${ID}_b2log_Novogene.txt | samtools view -bt ${dirref}.fai -o $dirres/${ID}_Novogene/${ID}_b2.bam -
	echo -ne "${ID}\t" >>$dirres/quality_matrix_b2_Novogene.txt
	perl -ne 'if(/^(\d+) reads; of these/){print $1,"\t";} if(/\((.*%)\) aligned concordantly 0 times/){print $1,"\t";} if(/\((.*%)\) aligned concordantly exactly 1 time/){print $1,"\t";} if(/\((.*%)\) aligned concordantly >1 times/){print $1,"\t";} if(/\((.*%)\) aligned discordantly 1 time/){print $1,"\t";} if(/\((.*%)\) aligned exactly 1 time/){print $1,"\t";} if(/\((.*%)\) aligned >1 times/){print $1,"\t";} if(/^(.*%) overall alignment rate/){print $1,"\n";}' ${dirres}/${ID}_Novogene/${ID}_b2log_Novogene.txt >>$dirres/quality_matrix_b2_Novogene.txt

	java -Djava.io.tmpdir=`pwd`/temp -jar  $HPC_PICARD_DIR/picard.jar SortSam I=$dirres/${ID}_Novogene/${ID}_b2.bam O=$dirres/${ID}_Novogene/${ID}_b2.sor.bam SO=queryname
        java -jar $HPC_PICARD_DIR/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$dirres/${ID}_Novogene/${ID}_b2.sor.bam O=$dirres/${ID}_Novogene/${ID}_b2.dedup.sor.bam METRICS_FILE=$dirres/${ID}_Novogene/${ID}_b2.mkdup.matrics.txt

	samtools view -b -q 10 ${dirres}/${ID}_Novogene/${ID}_b2.dedup.sor.bam >${dirres}/${ID}_Novogene/${ID}_b2_rmMultiQual.bam
	module load rseqc/3.0.0
	/apps/rseqc/3.0.0/bin/bam_stat.py -i ${dirres}/${ID}_Novogene/${ID}_b2_rmMultiQual.bam -q 4 &>${dirres}/${ID}_Novogene/${ID}_b2_rmMultiQual.bam.stat.txt
	/apps/rseqc/3.0.0/bin/bam_stat.py -i ${dirres}/${ID}_Novogene/${ID}_b2.bam -q 4 &>${dirres}/${ID}_Novogene/${ID}_b2.bam.stat.txt
	echo -ne "${ID}_rawBam\t" >>$dirres/quality_matrix_b2_bamStat_Novogene.txt
	perl -e '$TotalR=0; $Unique=0; $read1=0; $read2=0; $strand1=0; $strand2=0; $nonSplice=0; $mappedPair=0; $unMapped=0; $nonUnique=0; $PairedChr=0; $nonPri=0; $splice=0; $dups=0; while(<>){if(/Total records:\s+(\d+)/){$TotalR=$1;} if(/mapq >= mapq_cut \(unique\):\s+(\d+)/){$Unique=$1;} if(/Read-1:\s+(\d+)/){$read1=$1} if(/Read-2:\s+(\d+)/){$read2=$1} if(/Reads map to .\+.:\s+(\d+)/){$strand1=$1} if(/Reads map to .\-.:\s+(\d+)/){$strand2=$1;} if(/Non-splice reads:\s+(\d+)/){$nonSplice=$1} if(/Reads mapped in proper pairs:\s+(\d+)/){$mappedPair=$1} if(/Unmapped reads:\s+(\d+)/){$unMapped=$1} if(/mapq < mapq_cut \(non-unique\):\s+(\d+)/){$nonUnique=$1} if(/Proper-paired reads map to different chrom:(\d+)/){$PairedChr=$1} if(/Non primary hits\s+(\d+)/){$nonPri=$1} if(/Splice reads:\s+(\d+)/){$splice=$1} if(/Optical\/PCR duplicate:\s+(\d+)/){$dups=$1}} print $TotalR,"\t", $Unique,"\t", $read1,"\t", $read2,"\t", $strand1,"\t", $strand2,"\t", $nonSplice,"\t", $mappedPair,"\t", $unMapped,"\t", $nonUnique,"\t", $PairedChr,"\t", $nonPri,"\t", $splice,"\t", $dups,"\n"' ${dirres}/${ID}_Novogene/${ID}_b2.bam.stat.txt >>$dirres/quality_matrix_b2_bamStat_Novogene.txt
	echo -ne "${ID}_filBam\t" >>$dirres/quality_matrix_b2_bamStat_Novogene.txt
	perl -e '$TotalR=0; $Unique=0; $read1=0; $read2=0; $strand1=0; $strand2=0; $nonSplice=0; $mappedPair=0; $unMapped=0; $nonUnique=0; $PairedChr=0; $nonPri=0; $splice=0; $dups=0; while(<>){if(/Total records:\s+(\d+)/){$TotalR=$1;} if(/mapq >= mapq_cut \(unique\):\s+(\d+)/){$Unique=$1;} if(/Read-1:\s+(\d+)/){$read1=$1} if(/Read-2:\s+(\d+)/){$read2=$1} if(/Reads map to .\+.:\s+(\d+)/){$strand1=$1} if(/Reads map to .\-.:\s+(\d+)/){$strand2=$1;} if(/Non-splice reads:\s+(\d+)/){$nonSplice=$1} if(/Reads mapped in proper pairs:\s+(\d+)/){$mappedPair=$1} if(/Unmapped reads:\s+(\d+)/){$unMapped=$1} if(/mapq < mapq_cut \(non-unique\):\s+(\d+)/){$nonUnique=$1} if(/Proper-paired reads map to different chrom:(\d+)/){$PairedChr=$1} if(/Non primary hits\s+(\d+)/){$nonPri=$1} if(/Splice reads:\s+(\d+)/){$splice=$1} if(/Optical\/PCR duplicate:\s+(\d+)/){$dups=$1}} print $TotalR,"\t", $Unique,"\t", $read1,"\t", $read2,"\t", $strand1,"\t", $strand2,"\t", $nonSplice,"\t", $mappedPair,"\t", $unMapped,"\t", $nonUnique,"\t", $PairedChr,"\t", $nonPri,"\t", $splice,"\t", $dups,"\n"' ${dirres}/${ID}_Novogene/${ID}_b2_rmMultiQual.bam.stat.txt >>$dirres/quality_matrix_b2_bamStat_Novogene.txt
	rm $dirres/${ID}_Novogene/${ID}_b2.sor.bam
	rm $dirres/${ID}_Novogene/${ID}*gz

elif [[ ${one} == *Admera* ]];then
	mkdir -p ${dirres}/${ID}_Admera
	fastqc -t 8 ${one} -o ${dirres}/${ID}_Admera
	fastqc -t 8 ${two} -o ${dirres}/${ID}_Admera
	perl -ne 'if(/Filename<\/td><td>(.*)<\/td><\/tr><tr><td>File type.*Encoding<\/td><td>(.*)<\/td><\/tr><tr><td>Total Sequences<\/td><td>(\d+)<\/td><\/tr><tr><td>Sequences flagged as poor quality<\/td><td>(.*)<\/td><\/tr><tr><td>Sequence length<\/td><td>(\d+)<\/td><\/tr><tr><td>%GC<\/td><td>(\d+)/){print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t";} if(/alt=\"\[(\w+)\]\"\/><a href=\"\#M10\">Adapter Content<\/a>/){print $1,"\n"}' ${dirres}/${ID}_Admera/*_001_fastqc.html >>$dirres/quality_matrix_sequencing_Admera.txt

	trimmomatic PE -threads 8 -phred33 ${one} ${two} ${dirres}/${ID}_Admera/${ID}_paired_1.fq.gz ${dirres}/${ID}_Admera/${ID}_unpaired_1.fq.gz ${dirres}/${ID}_Admera/${ID}_paired_2.fq.gz ${dirres}/${ID}_Admera/${ID}_unpaired_2.fq.gz ILLUMINACLIP:${dirres}/adapters.fa:1:22:11 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40 &>${dirres}/${ID}_Admera/${ID}_trimlog_Admera.txt

	fastqc -t 8 ${dirres}/${ID}_Admera/${ID}_paired_1.fq.gz
	fastqc -t 8 ${dirres}/${ID}_Admera/${ID}_paired_2.fq.gz
	perl -ne 'if(/Filename<\/td><td>(.*)<\/td><\/tr><tr><td>File type.*Encoding<\/td><td>(.*)<\/td><\/tr><tr><td>Total Sequences<\/td><td>(\d+)<\/td><\/tr><tr><td>Sequences flagged as poor quality<\/td><td>(.*)<\/td><\/tr><tr><td>Sequence length<\/td><td>(.*)<\/td><\/tr><tr><td>%GC<\/td><td>(\d+)/){print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t";} if(/alt=\"\[(\w+)\]\"\/><a href=\"\#M10\">Adapter Content<\/a>/){print $1,"\n"}' ${dirres}/${ID}_Admera/*paired_?_fastqc.html >>$dirres/quality_matrix_sequencing_Admera.txt
	perl -ne 'if(/share\/Gu5hmC\/bams\/(\S+)_Admera\/\S+.gz /){print $1,"\t";} if(/Input Read Pairs: (\d+) Both Surviving: (\d+) \((.*)\%\) Forward Only Surviving: /){print $1, "\t", $2, "\t", $3, "\n"}elsif(/Input Reads: (\d+) Surviving: (\d+) \((.*)\%\) Dropped:/){print $1, "\t", $2, "\t", $3, "\n"}' ${dirres}/${ID}_Admera/${ID}_trimlog_Admera.txt >>$dirres/quality_matrix_trim_Admera.txt

	bowtie2 -X2000 --mm --threads 8 -x $dirgenome -1 $dirres/${ID}_Admera/${ID}_paired_1.fq.gz -2 $dirres/${ID}_Admera/${ID}_paired_2.fq.gz 2>${dirres}/${ID}_Admera/${ID}_b2log_Admera.txt | samtools view -bt ${dirref}.fai -o $dirres/${ID}_Admera/${ID}_b2.bam -
	echo -ne "${ID}\t" >>$dirres/quality_matrix_b2_Admera.txt
	perl -ne 'if(/^(\d+) reads; of these/){print $1,"\t";} if(/\((.*%)\) aligned concordantly 0 times/){print $1,"\t";} if(/\((.*%)\) aligned concordantly exactly 1 time/){print $1,"\t";} if(/\((.*%)\) aligned concordantly >1 times/){print $1,"\t";} if(/\((.*%)\) aligned discordantly 1 time/){print $1,"\t";} if(/\((.*%)\) aligned exactly 1 time/){print $1,"\t";} if(/\((.*%)\) aligned >1 times/){print $1,"\t";} if(/^(.*%) overall alignment rate/){print $1,"\n";}' ${dirres}/${ID}_Admera/${ID}_b2log_Admera.txt >>$dirres/quality_matrix_b2_Admera.txt

	java -Djava.io.tmpdir=`pwd`/temp -jar  $HPC_PICARD_DIR/picard.jar SortSam I=$dirres/${ID}_Admera/${ID}_b2.bam O=$dirres/${ID}_Admera/${ID}_b2.sor.bam SO=queryname
	java -jar $HPC_PICARD_DIR/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$dirres/${ID}_Admera/${ID}_b2.sor.bam O=$dirres/${ID}_Admera/${ID}_b2.dedup.sor.bam METRICS_FILE=$dirres/${ID}_Admera/${ID}_b2.mkdup.matrics.txt

	samtools view -b -q 10 ${dirres}/${ID}_Admera/${ID}_b2.dedup.sor.bam >${dirres}/${ID}_Admera/${ID}_b2_rmMultiQual.bam
	module load rseqc/3.0.0
	/apps/rseqc/3.0.0/bin/bam_stat.py -i ${dirres}/${ID}_Admera/${ID}_b2_rmMultiQual.bam -q 4 &>${dirres}/${ID}_Admera/${ID}_b2_rmMultiQual.bam.stat.txt
	/apps/rseqc/3.0.0/bin/bam_stat.py -i ${dirres}/${ID}_Admera/${ID}_b2.bam -q 4 &>${dirres}/${ID}_Admera/${ID}_b2.bam.stat.txt
	echo -ne "${ID}_rawBam\t" >>$dirres/quality_matrix_b2_bamStat_Admera.txt
	perl -e '$TotalR=0; $Unique=0; $read1=0; $read2=0; $strand1=0; $strand2=0; $nonSplice=0; $mappedPair=0; $unMapped=0; $nonUnique=0; $PairedChr=0; $nonPri=0; $splice=0; $dups=0; while(<>){if(/Total records:\s+(\d+)/){$TotalR=$1;} if(/mapq >= mapq_cut \(unique\):\s+(\d+)/){$Unique=$1;} if(/Read-1:\s+(\d+)/){$read1=$1} if(/Read-2:\s+(\d+)/){$read2=$1} if(/Reads map to .\+.:\s+(\d+)/){$strand1=$1} if(/Reads map to .\-.:\s+(\d+)/){$strand2=$1;} if(/Non-splice reads:\s+(\d+)/){$nonSplice=$1} if(/Reads mapped in proper pairs:\s+(\d+)/){$mappedPair=$1} if(/Unmapped reads:\s+(\d+)/){$unMapped=$1} if(/mapq < mapq_cut \(non-unique\):\s+(\d+)/){$nonUnique=$1} if(/Proper-paired reads map to different chrom:(\d+)/){$PairedChr=$1} if(/Non primary hits\s+(\d+)/){$nonPri=$1} if(/Splice reads:\s+(\d+)/){$splice=$1} if(/Optical\/PCR duplicate:\s+(\d+)/){$dups=$1}} print $TotalR,"\t", $Unique,"\t", $read1,"\t", $read2,"\t", $strand1,"\t", $strand2,"\t", $nonSplice,"\t", $mappedPair,"\t", $unMapped,"\t", $nonUnique,"\t", $PairedChr,"\t", $nonPri,"\t", $splice,"\t", $dups,"\n"' ${dirres}/${ID}_Admera/${ID}_b2.bam.stat.txt >>$dirres/quality_matrix_b2_bamStat_Admera.txt
	echo -ne "${ID}_filBam\t" >>$dirres/quality_matrix_b2_bamStat_Admera.txt
	perl -e '$TotalR=0; $Unique=0; $read1=0; $read2=0; $strand1=0; $strand2=0; $nonSplice=0; $mappedPair=0; $unMapped=0; $nonUnique=0; $PairedChr=0; $nonPri=0; $splice=0; $dups=0; while(<>){if(/Total records:\s+(\d+)/){$TotalR=$1;} if(/mapq >= mapq_cut \(unique\):\s+(\d+)/){$Unique=$1;} if(/Read-1:\s+(\d+)/){$read1=$1} if(/Read-2:\s+(\d+)/){$read2=$1} if(/Reads map to .\+.:\s+(\d+)/){$strand1=$1} if(/Reads map to .\-.:\s+(\d+)/){$strand2=$1;} if(/Non-splice reads:\s+(\d+)/){$nonSplice=$1} if(/Reads mapped in proper pairs:\s+(\d+)/){$mappedPair=$1} if(/Unmapped reads:\s+(\d+)/){$unMapped=$1} if(/mapq < mapq_cut \(non-unique\):\s+(\d+)/){$nonUnique=$1} if(/Proper-paired reads map to different chrom:(\d+)/){$PairedChr=$1} if(/Non primary hits\s+(\d+)/){$nonPri=$1} if(/Splice reads:\s+(\d+)/){$splice=$1} if(/Optical\/PCR duplicate:\s+(\d+)/){$dups=$1}} print $TotalR,"\t", $Unique,"\t", $read1,"\t", $read2,"\t", $strand1,"\t", $strand2,"\t", $nonSplice,"\t", $mappedPair,"\t", $unMapped,"\t", $nonUnique,"\t", $PairedChr,"\t", $nonPri,"\t", $splice,"\t", $dups,"\n"' ${dirres}/${ID}_Admera/${ID}_b2_rmMultiQual.bam.stat.txt >>$dirres/quality_matrix_b2_bamStat_Admera.txt
	rm $dirres/${ID}_Admera/${ID}_b2.sor.bam
	rm $dirres/${ID}_Admera/${ID}*gz

fi


echo "***** finish fastqc trimming and alignment for ${ID} *****"

date

