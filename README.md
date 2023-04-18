# GLONG mouse TEspeX  

This is the analysis of the TE families expressed in the GLONG mouse procssed as described in the mouse_TEspeX Github repository  
The mouse_te_v0.01.sh for PE reads was used for the configuration  

1- The mouse TE counting is started with the mouse_te_v0.01.pbs command after putting the different fastq files in a folder  
2- The RNA-Seq_mouse pipeline is usually ran before the TE counting, allowing the merging of the two types of datasets  
3- The STAR counts from the RNA-Seq_mouse pipeline is placed inside the same folder as in 1) for differentital analysis  

