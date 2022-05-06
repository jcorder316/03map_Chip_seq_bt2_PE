#!/bin/bash
source ~/.bash_profile
p=..
p1=$p/01RAW/
tri=~/001_sofwares/Trimmomatic-0.36
BA=$p/02BAMS
## path to Bt2 index
geno10=~/002_Genomes/01_Mm/mm10/bt2/mm10
## path to Trimmomatic Bt2 index
ADAPTERS=~/001_sofwares/Trimmomatic-0.36/adapters/TruSeq3-SE.fa	
igv=$p/03IGV
FA=$p/00QC
### FIRST prepare the name of the files PROPERLY!! a tip is to add 01, 02 at the beginning of the name of the sample, to make downstream analysis easier
######### second run for one sample fastqc to see how to trim the samples, USUALLY all the samples from the run will require similar trimming settings ####
		## open terminal in the 01RAW folder and run this for one sample fastqc sample_1.fastq.gz then see the result in the browser to see what trim1 and trim2 to take
## fastqc sample_1.fastq.gz

######### COMMEND:USUALLY all the samples from the run will require similar trimming settings ####

trim1="10"
trim2="70"
mini="15"

##### WHAT THIS SCRIPT WILL PRODUCE? mapped bam files (removing PCR duplicates) and BigWig files
##IMPORTANT- BEFORE running the script "run the echo line31 and 32 to see if the script read the Raw fastq file properly!..
mkdir ../00QC ../02BAMS ../03IGV
for R1 in $p1/*_1.fq.gz
 do
   R1_trim=${R1//.fastq.gz/_trim.fq.gz}
   yolo=$(echo "$R1" | rev | cut -c 9- | rev)
echo "$R1 and ${yolo}_trim.fq.gz"
   echo " java -jar $tri/trimmomatic-0.36.jar PE -threads 8 -phred33 ${samp} ${yolo}_2.fq.gz ${yolo}_fw_pair.fq ${yolo}_fw_unpair.fq.gz ${yolo}_rw_pair.fq ${yolo}_rw_unpair.fq.gz ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$mini CROP:$trim2 HEADCROP:$trim1 >> trimmomatic.cmds"
#java -jar $tri/trimmomatic-0.36.jar SE -threads 4 -phred33 $R1 ${yolo}_trim.fq.gz ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$mini CROP:$trim2 HEADCROP:$trim1
java -jar $tri/trimmomatic-0.36.jar PE -threads 8 -phred33 ${samp} ${yolo}_2.fq.gz ${yolo}_fw_pair.fq ${yolo}_fw_unpair.fq.gz ${yolo}_rw_pair.fq ${yolo}_rw_unpair.fq.gz ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$mini CROP:$trim2 HEADCROP:$trim1 >> trimmomatic.cmds

		bowtie2 -p8 -x $geno10 -1 ${yolo}_fw_pair.fq -2 ${yolo}_rw_pair.fq -S ${yolo}_mm10.sam
  
			samtools view -Sb -u ${yolo}_mm10.sam | samtools sort - -o ${yolo}_s.bam
				
			java -jar ~/001_sofwares/picard-tools-1.119/MarkDuplicates.jar INPUT=${yolo}_s.bam OUTPUT=${yolo}_s_rmd.bam METRICS_FILE=${yolo}_rmd_dupl_INFO.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true
#			 	
			bamToBed -i ${yolo}_s_rmd.bam  > ${yolo}_s_rmd.bed
						
			makeTagDirectory ${yolo}_ucsc ${yolo}_s_rmd.bed
			
			samtools index ${yolo}_s_rmd.bam	
      
      bamCoverage -b ${yolo}_s_rmd.bam -bs 20 --smoothLength 40 -p max  --normalizeUsing RPKM -e 150 -o ${yolo}_RPKM_mm10.bw
					
mv ${yolo}_RPKM_mm10.bw $igv
### this files could be removed , Because they are heavy and other files contains similar infomation
rm ${yolo}_s_rmd.bed 
rm ${yolo}_mm10.sam
rm ${yolo}_s.bam
### FQ for all RAW samples#### It could also run for TRIM
fastqc $R1
#fastqc ${yolo}_trim.fq.gz
rm ${yolo}_trim.fq.gz
mv ${yolo}_fastqc.html $FA
mv ${yolo}_fastqc.zip $FA
mv ${yolo}_rmd_dupl_INFO.txt $FA
mv ${yolo}_s_rmd.bam $BA
mv ${yolo}_ucsc $BA
mv ${yolo}_s_rmd.bam.bai $BA
done

## this work if you multiqc is installed in the same enviroment where you run the script. otherwide could be run separate it.
cd ../00QC/multiqc * -n 001NAME_of_QC_FILE

