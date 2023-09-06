# 5hmC
The codes shown here were used to perform the alignment of the fastq files to human genome (hg38), peak calling, elastic net selection, and regression for identifying AD related 5hmCs.

## analysis process
Alignment:

    sequencing_process.sh


Peak calling:

    peakCalling_paired.sh
  
    peakCalling_single.sh
  
    combine_peaks.R


Elastic Net selection:

    en_sel.R


Regression:

    glm.R
