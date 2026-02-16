# This document contains all the code used to analyze 2bRAD sequencing data for Jaskiel et al., 2026
# Author: Jacob Jaskiel
# I have adapted much of what is here at the beginning from Misha Matz's readme file, found here: https://github.com/z0on/2bRAD_denovo/blob/master/2bRAD_README.sh

#================================== Getting Set Up ==========================================
# (Download scripts from Misha Matz's github: https://github.com/z0on/2bRAD_denovo)
# Clone github into working directory:
git clone https://github.com/z0on/2bRAD_denovo.git

# Make everything executable (chmod) -> This way we can access everything. Just copy and paste these commands into the terminal, they will automatically run
chmod +x *.pl
chmod +x *.py
chmod +x *.R
chmod +x *.txt

# Install modules - we will need perl, bowtie2, samtools, and picard. They are pre-installed as modules on TACC and they are already installed on the SCC
# Copy and paste into terminal
module load perl
module load bowtie2
module load samtools
module load picard

# Reference Genome found here: https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_914725855.1/
# SSID fasta file
GCF_914725855.1_fThuAlb1.1_genomic.fasta

# Index transcriptome for bowtie2 mapper (submit as a job as this can take a while)
>indexing
nano indexing 
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N indexing # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
module load bowtie2
export GENOME_FASTA=GCF_914725855.1_fThuAlb1.1_genomic.fasta
bowtie2-build $GENOME_FASTA $GENOME_FASTA 

qsub indexing

>indexing2
nano indexing2
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N indexing2 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
module load samtools
samtools faidx $GENOME_FASTA

qsub indexing2

#================================== Demultiplexing and Trimming Barcodes/Poor Quality End Reads ==========================================
# Step 1: Splitting by in-read barcode, deduplicating and quality-filtering the reads
pwd
/projectnb/mullenl/jaskiel/2bRAD_fastq/2bRAD_2024
gunzip *.gz

# Creating a file of commands to run (assuming reads are in fastq files, one file per sample.)
/projectnb/mullenl/jaskiel/2bRAD_fastq/2bRAD_2024/2bRAD_trim_launch_dedup_old.pl fastq > trims   #Use the _old script if using oligos purchased before March 2022 like me

# Modify this trims file to be able to submit it as a job on the SCC, designate where the perl script is, and add '&'s to the end of the line
# It looks like this:

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N trims # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m be
module load perl
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_15_S15_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_7_S7_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_14_S14_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_23_S7_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_5_S5_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_16_S16_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_18_S2_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_13_S13_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_8_S8_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_4_S4_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_21_S5_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_11_S11_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_6_S6_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_10_S10_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_20_S4_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_24_S8_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_9_S9_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_1_S1_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_17_S1_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_3_S3_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_12_S12_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_22_S6_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_19_S3_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_2_S2_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
wait

qsub trims

# Check to see if you have expected number of *.tr0 files
ls -l *.tr0 | wc -l

# Cutadapt info:
# -q is a quality filter, used to trim low-quality ends from reads.
# The comma separated argument here trims the 5' end with a cutoff of 15 and the 3' end with a cutoff of 15.
# -m is a minimum length filter, it discards processed reads that are shorter than LENGTH provided (here, 25).
# For reference-based analysis: trimming poor quality bases off ends:
>trimse
for file in *.tr0; do
echo "cutadapt -q 15,15 -m 25 -o ${file/.tr0/}.trim $file > ${file}_trimlog.txt" >> trimse;  #This is run outside of the trimse job, and fills in trimse with the proper commands but you still need to add the header and load the relevant modules
done

nano trimse
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N trimse # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -j y # Join standard output and error to a single file
#$ -o trimse.qlog # Name the file where to redirect standard output and error
module load cutadapt 
cutadapt -q 15,15 -m 25 -o Pool_10_S10_L002_R1_001_ACCA.trim Pool_10_S10_L002_R1_001_ACCA.tr0 > Pool_10_S10_L002_R1_001_ACCA.tr0_trimlog.txt
cutadapt -q 15,15 -m 25 -o Pool_10_S10_L002_R1_001_AGAC.trim Pool_10_S10_L002_R1_001_AGAC.tr0 > Pool_10_S10_L002_R1_001_AGAC.tr0_trimlog.txt
cutadapt -q 15,15 -m 25 -o Pool_10_S10_L002_R1_001_AGTG.trim Pool_10_S10_L002_R1_001_AGTG.tr0 > Pool_10_S10_L002_R1_001_AGTG.tr0_trimlog.txt
cutadapt -q 15,15 -m 25 -o Pool_10_S10_L002_R1_001_CATC.trim Pool_10_S10_L002_R1_001_CATC.tr0 > Pool_10_S10_L002_R1_001_CATC.tr0_trimlog.txt
cutadapt -q 15,15 -m 25 -o Pool_10_S10_L002_R1_001_CTAC.trim Pool_10_S10_L002_R1_001_CTAC.tr0 > Pool_10_S10_L002_R1_001_CTAC.tr0_trimlog.txt
cutadapt -q 15,15 -m 25 -o Pool_10_S10_L002_R1_001_GACT.trim Pool_10_S10_L002_R1_001_GACT.tr0 > Pool_10_S10_L002_R1_001_GACT.tr0_trimlog.txt
cutadapt -q 15,15 -m 25 -o Pool_10_S10_L002_R1_001_GCTT.trim Pool_10_S10_L002_R1_001_GCTT.tr0 > Pool_10_S10_L002_R1_001_GCTT.tr0_trimlog.txt
cutadapt -q 15,15 -m 25 -o Pool_10_S10_L002_R1_001_GTGA.trim Pool_10_S10_L002_R1_001_GTGA.tr0 > Pool_10_S10_L002_R1_001_GTGA.tr0_trimlog.txt
#only first 8 rows shown

# Submit the trimse job:
qsub trimse

# Do we have expected number of *.trim files created?
ls -l *.trim | wc -l #yes, proceed!

#================================== Mapping to Reference Genome ==========================================

export GENOME_FASTA=GCF_914725855.1_fThuAlb1.1_genomic.fasta

# Create a file called 'maps'
>maps
../2bRAD_2024/2bRAD_bowtie2_launch.pl '\.trim$' $GENOME_FASTA > maps

# Same as you did above for trimse, nano into maps once you created it and add your header and your 'module load line'.
nano maps
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N maps # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o maps.qlog # Name the file where to redirect standard output and error
module load bowtie2
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x GCF_914725855.1_fThuAlb1.1_genomic.fasta -U Pool_17_S1_L001_R1_001_AGAC.trim -S Pool_17_S1_L001_R1_001_AGAC.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x GCF_914725855.1_fThuAlb1.1_genomic.fasta -U Pool_11_S11_L002_R1_001_TGTC.trim -S Pool_11_S11_L002_R1_001_TGTC.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x GCF_914725855.1_fThuAlb1.1_genomic.fasta -U Pool_17_S1_L001_R1_001_GCTT.trim -S Pool_17_S1_L001_R1_001_GCTT.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x GCF_914725855.1_fThuAlb1.1_genomic.fasta -U Pool_15_S15_L002_R1_001_TGGT.trim -S Pool_15_S15_L002_R1_001_TGGT.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x GCF_914725855.1_fThuAlb1.1_genomic.fasta -U Pool_5_S5_L001_R1_001_GACT.trim -S Pool_5_S5_L001_R1_001_GACT.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x GCF_914725855.1_fThuAlb1.1_genomic.fasta -U Pool_21_S5_L001_R1_001_ACCA.trim -S Pool_21_S5_L001_R1_001_ACCA.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x GCF_914725855.1_fThuAlb1.1_genomic.fasta -U Pool_14_S14_L002_R1_001_TCAC.trim -S Pool_14_S14_L002_R1_001_TCAC.trim.bt2.sam

# Setting minimum score that a read has to achieve if the read is 33 base pairs in length (check this) so minimum score has to be 49.
# L,16,1. L=length, 16=constant, 1=multiplier of L (1xL). If you have a length of 33, then say 16+(33x1)=49.
# Matched bases=2 points. Mismatched bases=-6 points, and 2 base pair read gaps= -11 points
# (33x2)-6-11=49, therefore our alignment threshold allows 1 base pair mismatch and 1 gap, or 2 mismatches 
# Submit the job
qsub maps

# Create a list of the .sam files you made and check how many - should equal your amount of samples
ls -l *.sam | wc -l
ls *.sam > sams
cat sams | wc -l

# Next stage is compressing, sorting and indexing the SAM files, so they become BAM files:
# Create a file called s2b
>s2b
for file in *.sam; do
echo "samtools sort -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b;
done

# Add module load samtools to top of s2b file
nano s2b
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N s2b # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o s2b.qlog # Name the file where to redirect standard output and error
module load samtools
samtools sort -O bam -o Pool_10_S10_L002_R1_001_ACCA.trim.bt2.bam Pool_10_S10_L002_R1_001_ACCA.trim.bt2.sam && samtools index Pool_10_S10_L002_R1_001_ACCA.trim.bt2.bam
samtools sort -O bam -o Pool_10_S10_L002_R1_001_AGAC.trim.bt2.bam Pool_10_S10_L002_R1_001_AGAC.trim.bt2.sam && samtools index Pool_10_S10_L002_R1_001_AGAC.trim.bt2.bam
samtools sort -O bam -o Pool_10_S10_L002_R1_001_AGTG.trim.bt2.bam Pool_10_S10_L002_R1_001_AGTG.trim.bt2.sam && samtools index Pool_10_S10_L002_R1_001_AGTG.trim.bt2.bam
samtools sort -O bam -o Pool_10_S10_L002_R1_001_CATC.trim.bt2.bam Pool_10_S10_L002_R1_001_CATC.trim.bt2.sam && samtools index Pool_10_S10_L002_R1_001_CATC.trim.bt2.bam
samtools sort -O bam -o Pool_10_S10_L002_R1_001_CTAC.trim.bt2.bam Pool_10_S10_L002_R1_001_CTAC.trim.bt2.sam && samtools index Pool_10_S10_L002_R1_001_CTAC.trim.bt2.bam
samtools sort -O bam -o Pool_10_S10_L002_R1_001_GACT.trim.bt2.bam Pool_10_S10_L002_R1_001_GACT.trim.bt2.sam && samtools index Pool_10_S10_L002_R1_001_GACT.trim.bt2.bam
samtools sort -O bam -o Pool_10_S10_L002_R1_001_GCTT.trim.bt2.bam Pool_10_S10_L002_R1_001_GCTT.trim.bt2.sam && samtools index Pool_10_S10_L002_R1_001_GCTT.trim.bt2.bam

qsub s2b

# Do we have the correct number of bam files? Should be the same number as number of trim files
ls *.bam | wc -l  
# BAM files are the input into various genotype calling / popgen programs, this is the main interim result of the analysis. Archive them.

#Checking fragment length
>frag_length_by_samp
nano frag_length_by_samp
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N frag_length_by_samp # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o frag_length_by_samp.qlog # Name the file where to redirect standard output and error
module load samtools
for bam in *.bam; do
    echo "=== $bam ==="
    samtools stats "$bam" | grep '^RL' | awk '{print $2 "\t" $3}' | sort -n
    echo ""
done > all_bams_read_lengths.txt

qsub frag_length_by_samp

#Across all
samtools cat *.bam | samtools stats - > combined_stats.txt

>frag_length_all
nano frag_length_all
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N frag_length_all # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 2
#$ -j y # Join standard output and error to a single file
#$ -o frag_length_all.qlog # Name the file where to redirect standard output and error
module load samtools
for bam in *.bam; do
    samtools view "$bam" | cut -f10 | awk '{print length($1)}'
done | sort | uniq -c | sort -n > global_length_histogram.txt

qsub frag_length_all

#===================================== ANGSD "Fuzzy" Genotyping ============================================
pwd
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis

# This first run-through of ANGSD is just to take a look at base qualities and coverage depth 
# We will run angsd again with filters informed by this first run-through
# "FUZZY genotyping" with ANGSD - without calling actual genotypes but working with genotype likelihoods at each SNP. Optimal for low-coverage data (<10x).
ls *.bam > bams
ls -l *.bam | wc -l
#Urgent: Before putting anything into R, you'll need a .csv file with the bam names and other data about each individual. If the order of bams is not the same as in the bams file on SCC, then samples will be mislabelled in R figures. It's a good idea to add any relevant metadata too

#------------------------------- Assessing base qualities and coverage depth -------------------------------------
# angsd settings:
# -minMapQ 20 : only highly unique mappings (prob of erroneous mapping =< 1%)
# -baq 1 : realign around indels (not terribly relevant for 2bRAD reads mapped with --local option) 
# -maxDepth : highest total depth (sum over all samples) to assess; set to 10x number of samples
# -minInd : the minimal number of individuals the site must be genotyped in. Reset to 50% of total N at this stage.

>angsdDD
nano angsdDD
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsdDD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsdDD.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 4660 -minInd 233"
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out dd 

qsub angsdDD

# Summarizing results (using modified script by Matteo Fumagalli)
module load R
Rscript /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/plotQC.R prefix=dd

# Proportion of sites covered at >5x:
cat quality.txt

# I've removed all low quality bams that have a proportion of sites covered at >5x under 0.4; the rest of the bams are in bams_all.csv
# scp dd.pdf to laptop to look at distribution of base quality scores, fraction of sites in each sample passing coverage thresholds, and fraction of sites passing genotyping rates cutoffs. Use these to guide choices of -minQ,  -minIndDepth and -minInd filters in subsequent ANGSD runs 

#----------------------------------- Tech Rep Detection ----------------------------------
# Looking at all samples including technical replicates to confirm accuracy of sequencing and also help determine filtering thresholds
# Detecting TRs (note: lower minind makes it easier to determine tech reps, but should be raised for subsequent analyses!)

>angsd_TR_all_mi50
nano angsd_TR_all_mi50
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_TR_all_mi50 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_TR_all_mi50.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 50 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out angsd_TR_all_mi50

qsub angsd_TR_all_mi50

NSITES=`zcat angsd_TR_all_mi50.mafs.gz | wc -l` 
echo $NSITES
#31694 sites

#Trying TR detection with lower minMaf
>angsd_TR_all_mi50_minmaf01
nano angsd_TR_all_mi50_minmaf01
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_TR_all_mi50_minmaf01 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_TR_all_mi50_minmaf01.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 50 -snp_pval 1e-5 -minMaf 0.01"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out angsd_TR_all_mi50_minmaf01

qsub angsd_TR_all_mi50_minmaf01

NSITES=`zcat angsd_TR_all_mi50_minmaf01.mafs.gz | wc -l` 
echo $NSITES
#84566 vs 31694 sites with minMaf = 0.05

>angsd_TR_all_mi50_minmaf02
nano angsd_TR_all_mi50_minmaf02
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_TR_all_mi50_minmaf02 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_TR_all_mi50_minmaf02.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 50 -snp_pval 1e-5 -minMaf 0.02"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out angsd_TR_all_mi50_minmaf02

qsub angsd_TR_all_mi50_minmaf02

NSITES=`zcat angsd_TR_all_mi50_minmaf02.mafs.gz | wc -l` 
echo $NSITES
#57285 vs 31694 sites with minMaf = 0.05

#----------------------------------- Merging Tech Reps ------------------------------------------------------ 
#Also removing 3 weird "noisy" samples, 268-G273, 281-G48, and 287-G53 that seem to be poor quality/coverage and cause issues in ordination
#Pool_4_S4_L001_R1_001_TGTC.trim.bt2.bam, Pool_2_Jaskiel_Lane1_S2_L003_R1_001_GTGA.trim.bt2.bam, Pool_12_S12_L002_R1_001_TGTC.trim.bt2.bam

#Merging bams, indexing merged bams
>merge_bams
nano merge_bams
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N merge_bams # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o merge_bams.qlog # Name the file where to redirect standard output and error
module load samtools
set -euo pipefail

while IFS=$'\t' read -r sample bam1 bam2 bam3; do
  bams=()
  bams+=("$bam1")
  bams+=("$bam2")
  [[ -n "$bam3" ]] && bams+=("$bam3")

  for f in "${bams[@]}"; do
    if [[ ! -r "$f" ]]; then
      echo "ERROR: BAM file $f for sample $sample not found or not readable" >&2
      exit 1
    fi
  done

  out="${sample}.merged.bam"
  echo "Merging ${#bams[@]} replicates for $sample → $out"
  samtools merge "$out" "${bams[@]}"
  samtools index "$out"

done < bams_to_merge.txt

qsub merge_bams

#Checking to see it worked
samtools idxstats Pool_10_Jaskiel_Lane2_S10_L004_R1_001_CTAC.trim.bt2.bam | awk '{sum+=$3+$4} END{print sum}' #740546
samtools idxstats Pool_1_Jaskiel_Lane1_S1_L003_R1_001_AGAC.trim.bt2.bam | awk '{sum+=$3+$4} END{print sum}' #441491
samtools idxstats 268-G11.merged.bam | awk '{sum+=$3+$4} END{print sum}' #1182037 
#looks good

#-------------------- ANGSD with all samples (merged bams), no Auxis --------------------
#start by filtering out sites to work with

# filtering sites to work on - use only filters that do not distort allele frequency
# set minInd to 75-90% of the total number fo individuals in the project
# if you are doing any other RAD than 2bRAD or GBS, remove '-sb_pval 1e-5' from FILTERS

>angsd_AllSites_merged
nano angsd_AllSites_merged
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_AllSites_merged # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_AllSites_merged.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 288"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 8"
angsd -b bams_all_noTR_merged -GL 1 $FILTERS $TODO -P 1 -out AllSites_merged

qsub angsd_AllSites_merged

NSITES=`zcat AllSites_merged.mafs.gz | wc -l` 
echo $NSITES
#715571 vs 704405 before merging bams from technical replicates

# Collecting and indexing filter-passing sites
zcat AllSites_merged.mafs.gz | cut -f 1,2 | tail -n +2 >AllSites_merged

>indexing_AllSites_merged
nano indexing_AllSites_merged
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N indexing_AllSites_merged # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o indexing_AllSites_merged.qlog # Name the file where to redirect standard output and error
module load angsd
angsd sites index AllSites_merged

qsub indexing_AllSites_merged

#Running ANGSD to produce PCA/dendrogram for species validation, input to admixture analyses
#minind 285 (80% of 356) minQ 25 w/ -setMinDepthInd 5
#With linked sites
>angsd_all_mid5_merged_noaux
nano angsd_all_mid5_merged_noaux
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_all_mid5_merged_noaux # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_all_mid5_merged_noaux.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 285 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_all_no_aux_merged -sites AllSites_merged -GL 1 $FILTERS $TODO -P 1 -out angsd_all_mid5_merged_noaux

qsub angsd_all_mid5_merged_noaux

NSITES=`zcat angsd_all_mid5_merged_noaux.mafs.gz | wc -l` 
echo $NSITES
#2539 vs 2314 w/ LD sites filtered

#minind 285 (80% of 356) minQ 25 w/ -setMinDepthInd 5, minmaf 0.01
#With linked sites
>angsd_all_mid5_merged_noaux_minmaf01
nano angsd_all_mid5_merged_noaux_minmaf01
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_all_mid5_merged_noaux_minmaf01 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_all_mid5_merged_noaux_minmaf01.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 285 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.01"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_all_no_aux_merged -sites AllSites_merged -GL 1 $FILTERS $TODO -P 1 -out angsd_all_mid5_merged_noaux_minmaf01

qsub angsd_all_mid5_merged_noaux_minmaf01

NSITES=`zcat angsd_all_mid5_merged_noaux_minmaf01.mafs.gz | wc -l` 
echo $NSITES
#5821 vs 5375 w/ LD sites filtered

#================================= NGSadmix and ADMIXTURE (all species) ========================================
pwd
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis

# First, install:
wget popgen.dk/software/download/NGSadmix/ngsadmix32.cpp 
g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix

# --- All taxa no technical replicates or genotyping duplicates (normal filters)
# NgsAdmix for K from 1 to 7 : do not run if the dataset contains clones or genotyping replicates!
>ngsadmix_all_mid5_merged_noaux
nano ngsadmix_all_mid5_merged_noaux
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsadmix_all_mid5_merged_noaux # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o ngsadmix_all_mid5_merged_noaux.qlog # Name the file where to redirect standard output and error
for K in `seq 1 10` ;
do
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/NGSadmix -likes angsd_all_mid5_merged_noaux.beagle.gz -K $K -P 1 -o /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_all_output_files/ngsadmix_all_mid5_merged_noaux_k${K};
done

qsub ngsadmix_all_mid5_merged_noaux

# scp *qopt files to local machine, use admixturePlotting_v5.R to plot

##Trying ADMIXTURE to find optimal K for all taxa
#BCF -> VCF
module load htslib/1.16
module load bcftools/1.16
bcftools convert -O v -o angsd_all_mid5_merged_noaux.vcf angsd_all_mid5_merged_noaux.bcf

module load admixture/1.3.0
module load plink/1.90b6.4

#all_mid5_noLD_merged_noaux
plink --vcf angsd_all_mid5_merged_noaux.vcf --make-bed --allow-extra-chr 0 --out angsd_all_mid5_merged_noaux --const-fid 0

#Trying multiple replicates per K and Evanno method for better confidence in K
>admixture_all_mid5_merged_noaux
nano admixture_all_mid5_merged_noaux
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N admixture_all_mid5_merged_noaux # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o admixture_all_mid5_merged_noaux.qlog # Name the file where to redirect standard output and error
module load admixture/1.3.0
for K in $(seq 1 10); do
  for REP in $(seq 1 10); do
    SEED=$((1000 + $RANDOM % 100000))  # generates a random 4–5 digit number
    admixture --cv=10 --seed=$SEED angsd_all_mid5_merged_noaux.bed $K \
      | tee /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_all_output_files/myresult_all_mid5_merged_noaux_K${K}_rep${REP}.out
  done
done
grep "CV error" /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_all_output_files/myresult_all_mid5_merged_noaux_K*_rep*.out > /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_all_output_files/admixture_cv_errors_all_merged_noaux.txt

qsub admixture_all_mid5_merged_noaux

# --- All taxa no technical replicates or genotyping duplicates (normal filters)
# NgsAdmix for K from 1 to 7 : do not run if the dataset contains clones or genotyping replicates!
>ngsadmix_all_mid5_merged_noaux_minmaf01
nano ngsadmix_all_mid5_merged_noaux_minmaf01
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsadmix_all_mid5_merged_noaux_minmaf01 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o ngsadmix_all_mid5_merged_noaux_minmaf01.qlog # Name the file where to redirect standard output and error
for K in `seq 1 10` ;
do
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/NGSadmix -likes angsd_all_mid5_merged_noaux_minmaf01.beagle.gz -K $K -P 1 -o /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_all_output_files/ngsadmix_all_mid5_merged_noaux_minmaf01_k${K};
done

qsub ngsadmix_all_mid5_merged_noaux_minmaf01

# scp *qopt files to local machine, use admixturePlotting_v5.R to plot

##Trying ADMIXTURE to find optimal K for all taxa
#BCF -> VCF
module load htslib/1.16
module load bcftools/1.16
bcftools convert -O v -o angsd_all_mid5_merged_noaux_minmaf01.vcf angsd_all_mid5_merged_noaux_minmaf01.bcf

module load admixture/1.3.0
module load plink/1.90b6.4

#all_mid5_noLD_merged_noaux
plink --vcf angsd_all_mid5_merged_noaux_minmaf01.vcf --make-bed --allow-extra-chr 0 --out angsd_all_mid5_merged_noaux_minmaf01 --const-fid 0

#Trying multiple replicates per K and Evanno method for better confidence in K
>admixture_all_mid5_merged_noaux_minmaf01
nano admixture_all_mid5_merged_noaux_minmaf01
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N admixture_all_mid5_merged_noaux_minmaf01 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o admixture_all_mid5_merged_noaux_minmaf01.qlog # Name the file where to redirect standard output and error
module load admixture/1.3.0
for K in $(seq 1 10); do
  for REP in $(seq 1 10); do
    SEED=$((1000 + $RANDOM % 100000))  # generates a random 4–5 digit number
    admixture --cv=10 --seed=$SEED angsd_all_mid5_merged_noaux_minmaf01.bed $K \
      | tee /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_all_output_files/myresult_all_mid5_merged_noaux_minmaf01_K${K}_rep${REP}.out
  done
done
grep "CV error" /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_all_output_files/myresult_all_mid5_merged_noaux_minmaf01_K*_rep*.out > /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_all_output_files/admixture_cv_errors_all_merged_noaux_minmaf01.txt

qsub admixture_all_mid5_merged_noaux_minmaf01

#================================= Demographic Analyses ========================================
## SFS work with individual species groups
#With Linked Sites still
#Skipjack (203 individuals, so 80% is 162)
>sfs_skj_dem_merged_LD
nano sfs_skj_dem_merged_LD
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_dem_merged_LD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_dem_merged_LD.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites AllSites_skj_merged -b bams_skj_merged -GL 1 -P 1 -minInd 162 $TODO -out skj_dem_merged_LD

qsub sfs_skj_dem_merged_LD

#Yellowfin 
>sfs_yft_dem_merged_LD
nano sfs_yft_dem_merged_LD
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_dem_merged_LD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_dem_merged_LD.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites AllSites_yft_merged -b bams_yft_merged -GL 1 -P 1 -minInd 77 $TODO -out yft_dem_merged_LD

qsub sfs_yft_dem_merged_LD

#Bigeye
>sfs_bet_dem_merged_LD
nano sfs_bet_dem_merged_LD
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_dem_merged_LD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_dem_merged_LD.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites AllSites_bet_merged -b bams_bet_merged -GL 1 -P 1 -minInd 46 $TODO -out bet_dem_merged_LD

qsub sfs_bet_dem_merged_LD

# writing down 2d-SFS priors - SKJ vs YFT 
>sfs_skj_yft_merged_dem_folded_LD
nano sfs_skj_yft_merged_dem_folded_LD
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_yft_merged_dem_folded_LD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_yft_merged_dem_folded_LD.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS skj_dem_merged_LD.saf.idx yft_dem_merged_LD.saf.idx -P 1 -fold 1 > skj_yft_merged_dem_folded_LD.sfs ; realSFS fst index skj_dem_merged_LD.saf.idx yft_dem_merged_LD.saf.idx -sfs skj_yft_merged_dem_folded_LD.sfs -fstout skj_yft_merged_dem_folded_LD

qsub sfs_skj_yft_merged_dem_folded_LD

# global Fst between populations (after LD pruning and folding)
realSFS fst stats skj_yft_merged_dem_folded_LD.fst.idx
#output:
FST.Unweight[nObs:655893]:0.004028 Fst.Weight:0.249107

# The difference between unweighted and weighted values is averaging versus ratio of sums method. 
# The weighted value is derived from the ratio of the separate per-locus sums of numerator and denominator values, 
# while the unweighted value is the average of per-locus values. [If that is unclear: weighted is sum(a)/sum(a+b), while unweighted is average(a/(a+b))].

# writing down 2d-SFS priors - SKJ vs BET 
>sfs_skj_bet_merged_dem_folded_LD
nano sfs_skj_bet_merged_dem_folded_LD
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_bet_merged_dem_folded_LD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_bet_merged_dem_folded_LD.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS skj_dem_merged_LD.saf.idx bet_dem_merged_LD.saf.idx -P 1 -fold 1 > skj_bet_merged_dem_folded_LD.sfs ; realSFS fst index skj_dem_merged_LD.saf.idx bet_dem_merged_LD.saf.idx -sfs skj_bet_merged_dem_folded_LD.sfs -fstout skj_bet_merged_dem_folded_LD

qsub sfs_skj_bet_merged_dem_folded_LD

# global Fst between populations (after LD pruning and folding)
realSFS fst stats skj_bet_merged_dem_folded_LD.fst.idx
#output:
FST.Unweight[nObs:668111]:0.006065 Fst.Weight:0.278329

# writing down 2d-SFS priors - YFT vs BET 
>sfs_yft_bet_merged_dem_folded_LD
nano sfs_yft_bet_merged_dem_folded_LD
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_bet_merged_dem_folded_LD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_bet_merged_dem_folded_LD.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_dem_merged_LD.saf.idx bet_dem_merged_LD.saf.idx -P 1 -fold 1 > yft_bet_merged_dem_folded_LD.sfs ; realSFS fst index yft_dem_merged_LD.saf.idx bet_dem_merged_LD.saf.idx -sfs yft_bet_merged_dem_folded_LD.sfs -fstout yft_bet_merged_dem_folded_LD

qsub sfs_yft_bet_merged_dem_folded_LD

# global Fst between populations (after LD pruning and folding)
realSFS fst stats yft_bet_merged_dem_folded_LD.fst.idx
#output:
FST.Unweight[nObs:1651374]:0.006197 Fst.Weight:0.317648

#------------------------- Filtering for LD with ngsLD ---------------------------------
#https://github.com/fgvieira/ngsLD
#Use https://github.com/fgvieira/prune_graph since ngsLD's perl-based pruning script is SLOW and the python script isn't much better (days vs. minutes difference)

#Making AllSites files for individual species clusters before filtering for LD
#starting with skipjack mapped against yft ref genome
>angsd_skj_AllSites_merged 
nano angsd_skj_AllSites_merged
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_skj_AllSites_merged # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_skj_AllSites_merged.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 162"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 8"
angsd -b bams_skj_merged -GL 1 $FILTERS $TODO -P 1 -out AllSites_skj_merged

qsub angsd_skj_AllSites_merged

NSITES=`zcat AllSites_skj_merged.mafs.gz | wc -l` 
echo $NSITES
#990548

# Collecting and indexing filter-passing sites
zcat AllSites_skj_merged.mafs.gz | cut -f 1,2 | tail -n +2 >AllSites_skj_merged

>indexing_AllSites_skj_merged
nano indexing_AllSites_skj_merged
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N indexing_AllSites_skj_merged # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o indexing_AllSites_skj_merged.qlog # Name the file where to redirect standard output and error
module load angsd
angsd sites index AllSites_skj_merged

qsub indexing_AllSites_skj_merged

#Trying with a large maxkb distance to see the decay curve
#100kb, and minMaf = 0.01 
>ngsld_skj_unlinked_merged_0.5_100kb_minmaf01
nano ngsld_skj_unlinked_merged_0.5_100kb_minmaf01
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsld_skj_unlinked_merged_0.5_100kb_minmaf01 # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o ngsld_skj_unlinked_merged_0.5_100kb_minmaf01.qlog # Name the file where to redirect standard output and error
module load ngsld
NS=`zcat AllSites_skj_merged.geno.gz| wc -l`
NB=`cat bams_skj_merged | wc -l`
ngsLD --geno AllSites_skj_merged.geno.gz --probs 1 --n_ind $NB --n_sites $NS --max_kb_dist 100 --pos AllSites_skj_merged --out AllSites_skj_merged_0.5_100kb_minmaf01.LD --n_threads 1 --extend_out 1 --min_maf 0.01

qsub ngsld_skj_unlinked_merged_0.5_100kb_minmaf01

>LD_decay_skj
nano LD_decay_skj
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N LD_decay_skj # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o LD_decay_skj.qlog # Name the file where to redirect standard output and error
module load R
LD_IN="AllSites_skj_merged_0.5_100kb_minmaf01.LD"
SAMPLE_FRAC=0.5
LD_SUB="${LD_IN%.LD}.sub.LD"
HEADER='snp1\tsnp2\tdist\tr2_ExpG\tD\tDp\tr2\tn_samples\tmaf1\tmaf2\thap00\thap01\thap10\thap11\tchi2_1df'
awk -v p="$SAMPLE_FRAC" 'NR==1 || rand()<p' "$LD_IN" > "$LD_SUB"
sed -i "1i $HEADER" "$LD_SUB"
printf "file\n%s\n" "$LD_SUB" > skj_ld_files.tsv
Rscript --vanilla fit_LDdecay.R --ld_files skj_ld_files.tsv --header --col 3 --ld r2_ExpG --max_kb_dist 100 --fit_boot 200 --fit_bin_size 250 --fit_level 100 --plot_data --plot_scale 3 -o skj_LD_decay.pdf

qsub LD_decay_skj

#Prune with kb distance of 5kb and weight (Pearson's R^2) 0.2 based on LD plots
>prune_LD_skj_merged_0.2_5kb_minmaf01_pg
nano prune_LD_skj_merged_0.2_5kb_minmaf01_pg
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N prune_LD_skj_merged_0.2_5kb_minmaf01_pg # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o prune_LD_skj_merged_0.2_5kb_minmaf01_pg.qlog # Name the file where to redirect standard output and error
module load rust
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/prune_graph/target/release/prune_graph --in AllSites_skj_merged_0.5_100kb_minmaf01.LD --weight-field column_4 --weight-filter "column_3 <= 5000 && column_4 >= 0.2" --out 0.2_5kb_minmaf01_skj_merged_unlinked_pg

qsub prune_LD_skj_merged_0.2_5kb_minmaf01_pg

#make tab delimited and remove missing lines
sed 's/:/\t/g' 0.2_5kb_minmaf01_skj_merged_unlinked_pg > 0.2_5kb_minmaf01_skj_merged_unlinked_pg.sites
awk  '$2!=""' 0.2_5kb_minmaf01_skj_merged_unlinked_pg.sites > 0.2_5kb_minmaf01_skj_merged_unlinked_pg.sites.tmp; mv 0.2_5kb_minmaf01_skj_merged_unlinked_pg.sites.tmp 0.2_5kb_minmaf01_skj_merged_unlinked_pg.sites
sort -k1 0.2_5kb_minmaf01_skj_merged_unlinked_pg.sites > 0.2_5kb_minmaf01_skj_merged_unlinked_pg.sites_2.tmp; mv 0.2_5kb_minmaf01_skj_merged_unlinked_pg.sites_2.tmp 0.2_5kb_minmaf01_skj_merged_unlinked_pg.sites

angsd sites index 0.2_5kb_minmaf01_skj_merged_unlinked_pg.sites

>final_sites_skj_merged_0.2_5kb_minmaf01_pg
nano final_sites_skj_merged_0.2_5kb_minmaf01_pg
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N final_sites_skj_merged_0.2_5kb_minmaf01_pg # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o final_sites_skj_merged_0.2_5kb_minmaf01_pg.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -maxHetFreq 0.5 -hetbias_pval 1e-5 -minInd 162 -sb_pval 1e-5"
TODO='-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2'
angsd -sites 0.2_5kb_minmaf01_skj_merged_unlinked_pg.sites -b bams_skj_merged -GL 1 $FILTERS $TODO -P 1 -out 0.2_5kb_minmaf01_skj_merged_unlinked_pg.sites

qsub final_sites_skj_merged_0.2_5kb_minmaf01_pg

NSITES=`zcat 0.2_5kb_minmaf01_skj_merged_unlinked_pg.sites.mafs.gz | wc -l`
echo $NSITES
#28561

zcat 0.2_5kb_minmaf01_skj_merged_unlinked_pg.sites.mafs.gz | cut -f 1,2 | tail -n +2 > finalsites_skj_unlinked
angsd sites index finalsites_skj_unlinked #use this for skj analysis going forward

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Yellowfin LD filtering
>angsd_yft_AllSites_merged 
nano angsd_yft_AllSites_merged
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_yft_AllSites_merged # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_yft_AllSites_merged.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 77"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 8"
angsd -b bams_yft_merged -GL 1 $FILTERS $TODO -P 1 -out AllSites_yft_merged

qsub angsd_yft_AllSites_merged

NSITES=`zcat AllSites_yft_merged.mafs.gz | wc -l` 
echo $NSITES
#1791054

# Collecting and indexing filter-passing sites
zcat AllSites_yft_merged.mafs.gz | cut -f 1,2 | tail -n +2 >AllSites_yft_merged

>indexing_AllSites_yft_merged
nano indexing_AllSites_yft_merged
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N indexing_AllSites_yft_merged # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o indexing_AllSites_yft_merged.qlog # Name the file where to redirect standard output and error
module load angsd
angsd sites index AllSites_yft_merged

qsub indexing_AllSites_yft_merged

#Trying with a large maxkb distance to see the decay curve, then rerun based on the R^2 and elbow of the decay curve
#100kb, and minMaf = 0.01 
>ngsld_yft_unlinked_merged_0.5_100kb_minmaf01
nano ngsld_yft_unlinked_merged_0.5_100kb_minmaf01
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsld_yft_unlinked_merged_0.5_100kb_minmaf01 # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o ngsld_yft_unlinked_merged_0.5_100kb_minmaf01.qlog # Name the file where to redirect standard output and error
module load ngsld
NS=`zcat AllSites_yft_merged.geno.gz| wc -l`
NB=`cat bams_yft_merged | wc -l`
ngsLD --geno AllSites_yft_merged.geno.gz --probs 1 --n_ind $NB --n_sites $NS --max_kb_dist 100 --pos AllSites_yft_merged --out AllSites_yft_merged_0.5_100kb_minmaf01.LD --n_threads 1 --extend_out 1 --min_maf 0.01

qsub ngsld_yft_unlinked_merged_0.5_100kb_minmaf01

#Plotting LD Decay for first 100kb, subsampling half of the LD file from the previous job 
>LD_decay_yft
nano LD_decay_yft
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N LD_decay_yft # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o LD_decay_yft.qlog # Name the file where to redirect standard output and error
module load R
LD_IN="AllSites_yft_merged_0.5_100kb_minmaf01.LD"
SAMPLE_FRAC=0.5
LD_SUB="${LD_IN%.LD}.sub.LD"
HEADER='snp1\tsnp2\tdist\tr2_ExpG\tD\tDp\tr2\tn_samples\tmaf1\tmaf2\thap00\thap01\thap10\thap11\tchi2_1df'
awk -v p="$SAMPLE_FRAC" 'NR==1 || rand()<p' "$LD_IN" > "$LD_SUB"
sed -i "1i $HEADER" "$LD_SUB"
printf "file\n%s\n" "$LD_SUB" > yft_ld_files.tsv
Rscript --vanilla fit_LDdecay.R --ld_files yft_ld_files.tsv --header --col 3 --ld r2_ExpG --max_kb_dist 100 --fit_boot 200 --fit_bin_size 250 --fit_level 100 --plot_data --plot_scale 3 -o yft_LD_decay.pdf

qsub LD_decay_yft

#Prune with kb distance of 5kb and weight (R^2) 0.2 pearsons r2
>prune_LD_yft_merged_0.2_5kb_minmaf01_pg
nano prune_LD_yft_merged_0.2_5kb_minmaf01_pg
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N prune_LD_yft_merged_0.2_5kb_minmaf01_pg # job name, anything you want
#$ -l h_rt=48:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o prune_LD_yft_merged_0.2_5kb_minmaf01_pg.qlog # Name the file where to redirect standard output and error
module load rust
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/prune_graph/target/release/prune_graph --in AllSites_yft_merged_0.5_100kb_minmaf01.LD --weight-field column_4 --weight-filter "column_3 <= 5000 && column_4 >= 0.2" --out 0.2_5kb_minmaf01_yft_merged_unlinked_pg

qsub prune_LD_yft_merged_0.2_5kb_minmaf01_pg

#make tab delimited and remove missing lines
sed 's/:/\t/g' 0.2_5kb_minmaf01_yft_merged_unlinked_pg > 0.2_5kb_minmaf01_yft_merged_unlinked_pg.sites
awk  '$2!=""' 0.2_5kb_minmaf01_yft_merged_unlinked_pg.sites > 0.2_5kb_minmaf01_yft_merged_unlinked_pg.sites.tmp; mv 0.2_5kb_minmaf01_yft_merged_unlinked_pg.sites.tmp 0.2_5kb_minmaf01_yft_merged_unlinked_pg.sites
sort -k1 0.2_5kb_minmaf01_yft_merged_unlinked_pg.sites > 0.2_5kb_minmaf01_yft_merged_unlinked_pg.sites_2.tmp; mv 0.2_5kb_minmaf01_yft_merged_unlinked_pg.sites_2.tmp 0.2_5kb_minmaf01_yft_merged_unlinked_pg.sites

angsd sites index 0.2_5kb_minmaf01_yft_merged_unlinked_pg.sites

>final_sites_yft_merged_0.2_5kb_minmaf01_pg
nano final_sites_yft_merged_0.2_5kb_minmaf01_pg
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N final_sites_yft_merged_0.2_5kb_minmaf01_pg # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o final_sites_yft_merged_0.2_5kb_minmaf01_pg.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -maxHetFreq 0.5 -hetbias_pval 1e-5 -minInd 77 -sb_pval 1e-5"
TODO='-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2'
angsd -sites 0.2_5kb_minmaf01_yft_merged_unlinked_pg.sites -b bams_yft_merged -GL 1 $FILTERS $TODO -P 1 -out 0.2_5kb_minmaf01_yft_merged_unlinked_pg.sites

qsub final_sites_yft_merged_0.2_5kb_minmaf01_pg

NSITES=`zcat 0.2_5kb_minmaf01_yft_merged_unlinked_pg.sites.mafs.gz | wc -l`
echo $NSITES
#107631

zcat 0.2_5kb_minmaf01_yft_merged_unlinked_pg.sites.mafs.gz | cut -f 1,2 | tail -n +2 > finalsites_yft_unlinked
angsd sites index finalsites_yft_unlinked #use this for yft analysis going forward

#-----------------------------------------------------------------------------------------
#Bigeye LD filtering
>angsd_bet_AllSites_merged 
nano angsd_bet_AllSites_merged
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_bet_AllSites_merged # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_bet_AllSites_merged.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 46"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 8"
angsd -b bams_bet_merged -GL 1 $FILTERS $TODO -P 1 -out AllSites_bet_merged

qsub angsd_bet_AllSites_merged

NSITES=`zcat AllSites_bet_merged.mafs.gz | wc -l` 
echo $NSITES
#1902666

# Collecting and indexing filter-passing sites
zcat AllSites_bet_merged.mafs.gz | cut -f 1,2 | tail -n +2 >AllSites_bet_merged

>indexing_AllSites_bet_merged
nano indexing_AllSites_bet_merged
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N indexing_AllSites_bet_merged # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o indexing_AllSites_bet_merged.qlog # Name the file where to redirect standard output and error
module load angsd
angsd sites index AllSites_bet_merged

qsub indexing_AllSites_bet_merged

#Trying with a large maxkb distance to see the decay curve, then rerun based on the R^2 and elbow of the decay curve
#0.5 weight, 100kb, and minMaf = 0.01 
>ngsld_bet_unlinked_merged_0.5_100kb_minmaf01
nano ngsld_bet_unlinked_merged_0.5_100kb_minmaf01
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsld_bet_unlinked_merged_0.5_100kb_minmaf01 # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o ngsld_bet_unlinked_merged_0.5_100kb_minmaf01.qlog # Name the file where to redirect standard output and error
module load ngsld
NS=`zcat AllSites_bet_merged.geno.gz| wc -l`
NB=`cat bams_bet_merged | wc -l`
ngsLD --geno AllSites_bet_merged.geno.gz --probs 1 --n_ind $NB --n_sites $NS --max_kb_dist 100 --pos AllSites_bet_merged --out AllSites_bet_merged_0.5_100kb_minmaf01.LD --n_threads 1 --extend_out 1 --min_maf 0.01

qsub ngsld_bet_unlinked_merged_0.5_100kb_minmaf01

>LD_decay_bet
nano LD_decay_bet
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N LD_decay_bet # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o LD_decay_bet.qlog # Name the file where to redirect standard output and error
module load R
LD_IN="AllSites_bet_merged_0.5_100kb_minmaf01.LD"
SAMPLE_FRAC=0.5
LD_SUB="${LD_IN%.LD}.sub.LD"
HEADER='snp1\tsnp2\tdist\tr2_ExpG\tD\tDp\tr2\tn_samples\tmaf1\tmaf2\thap00\thap01\thap10\thap11\tchi2_1df'
awk -v p="$SAMPLE_FRAC" 'NR==1 || rand()<p' "$LD_IN" > "$LD_SUB"
sed -i "1i $HEADER" "$LD_SUB"
printf "file\n%s\n" "$LD_SUB" > bet_ld_files.tsv
Rscript --vanilla fit_LDdecay.R --ld_files bet_ld_files.tsv --header --col 3 --ld r2_ExpG --max_kb_dist 100 --fit_boot 200 --fit_bin_size 250 --fit_level 100 --plot_data --plot_scale 3 -o bet_LD_decay.pdf

qsub LD_decay_bet

#Prune with kb distance of 5kb and weight (R^2) 0.2 pearsons r2 
>prune_LD_bet_merged_0.2_5kb_minmaf01_pg
nano prune_LD_bet_merged_0.2_5kb_minmaf01_pg
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N prune_LD_bet_merged_0.2_5kb_minmaf01_pg # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o prune_LD_bet_merged_0.2_5kb_minmaf01_pg.qlog # Name the file where to redirect standard output and error
module load rust
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/prune_graph/target/release/prune_graph --in AllSites_bet_merged_0.5_100kb_minmaf01.LD --weight-field column_4 --weight-filter "column_3 <= 5000 && column_4 >= 0.2" --out 0.2_5kb_minmaf01_bet_merged_unlinked_pg

qsub prune_LD_bet_merged_0.2_5kb_minmaf01_pg

#make tab delimited and remove missing lines
sed 's/:/\t/g' 0.2_5kb_minmaf01_bet_merged_unlinked_pg > 0.2_5kb_minmaf01_bet_merged_unlinked_pg.sites
awk  '$2!=""' 0.2_5kb_minmaf01_bet_merged_unlinked_pg.sites > 0.2_5kb_minmaf01_bet_merged_unlinked_pg.sites.tmp; mv 0.2_5kb_minmaf01_bet_merged_unlinked_pg.sites.tmp 0.2_5kb_minmaf01_bet_merged_unlinked_pg.sites
sort -k1 0.2_5kb_minmaf01_bet_merged_unlinked_pg.sites > 0.2_5kb_minmaf01_bet_merged_unlinked_pg.sites_2.tmp; mv 0.2_5kb_minmaf01_bet_merged_unlinked_pg.sites_2.tmp 0.2_5kb_minmaf01_bet_merged_unlinked_pg.sites

angsd sites index 0.2_5kb_minmaf01_bet_merged_unlinked_pg.sites

>final_sites_bet_merged_0.2_5kb_minmaf01_pg
nano final_sites_bet_merged_0.2_5kb_minmaf01_pg
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N final_sites_bet_merged_0.2_5kb_minmaf01_pg # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o final_sites_bet_merged_0.2_5kb_minmaf01_pg.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -maxHetFreq 0.5 -hetbias_pval 1e-5 -minInd 46 -sb_pval 1e-5"
TODO='-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2'
angsd -sites 0.2_5kb_minmaf01_bet_merged_unlinked_pg.sites -b bams_bet_merged -GL 1 $FILTERS $TODO -P 1 -out 0.2_5kb_minmaf01_bet_merged_unlinked_pg.sites

qsub final_sites_bet_merged_0.2_5kb_minmaf01_pg

NSITES=`zcat 0.2_5kb_minmaf01_bet_merged_unlinked_pg.sites.mafs.gz | wc -l`
echo $NSITES
#105972

zcat 0.2_5kb_minmaf01_bet_merged_unlinked_pg.sites.mafs.gz | cut -f 1,2 | tail -n +2 > finalsites_bet_unlinked
angsd sites index finalsites_bet_unlinked #use this for bet analysis going forward

#================================= NGSRelate =============================================
#Based on https://github.com/ANGSD/NgsRelate
#NOTE: you will probably need to rerun ANGSD with -doGlf set to 3 to get a binary-formatted beagle.gz file, which is the input to NGSrelate. Although the new version says it can handle vcf/bcf format
#Run on individual species clusters because estimates of relatedness can be wonky if there's background structure
#This includes the putative divergent yellowfin cluster
#Higher values of pairwise relatedness indicate closer relationships, with R = 1 for self, 0.5 for parent/child or full siblings, 0.25 for half-siblings, and 0 for unrelated individuals (in ideal, non-inbred scenarios). 

git clone --recursive https://github.com/SAMtools/htslib
git clone https://github.com/ANGSD/ngsRelate
cd htslib/;make -j2;cd ../ngsRelate;make HTSSRC=../htslib/

#Example from git
### First we generate a file with allele frequencies (angsdput.mafs.gz) and a file with genotype likelihoods (angsdput.glf.gz).
./angsd -b filelist -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3
### Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
zcat angsdput.mafs.gz | cut -f6 |sed 1d >freq
### run NgsRelate
./ngsrelate  -g angsdput.glf.gz -n 100 -f freq  -O newres

#Skipjack (same filters as for PCA, etc)
>angsd_ngsrelate_skj
nano angsd_ngsrelate_skj
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_ngsrelate_skj # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_ngsrelate_skj.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -setMinDepthInd 5 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 162 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBCF 1 -doPost 1 -doGlf 3" #-doGlf set to 3 for binary beagle file
angsd -b bams_skj_merged -sites finalsites_skj_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_ngsrelate_skj

qsub angsd_ngsrelate_skj

#how many SNPs?
NSITES=`zcat angsd_ngsrelate_skj.mafs.gz | wc -l`
echo $NSITES
#1262

zcat angsd_ngsrelate_skj.mafs.gz | cut -f5 |sed 1d >freq_skj
gunzip -c angsd_ngsrelate_skj.glf.gz |wc -c #6143592

>ngsrelate_skj
nano ngsrelate_skj
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsrelate_skj # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 4 #ngsRelate defaults to 4 threads
#$ -j y # Join standard output and error to a single file
#$ -o ngsrelate_skj.qlog # Name the file where to redirect standard output and error
module load ngsrelate
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsRelate/ngsRelate -g angsd_ngsrelate_skj.glf.gz -n 203 -f freq_skj -z skj_samps -O skj_relate

qsub ngsrelate_skj

cp skj_relate skj_relate.tsv

#Yellowfin
>angsd_ngsrelate_yft
nano angsd_ngsrelate_yft
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_ngsrelate_yft # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_ngsrelate_yft.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -setMinDepthInd 5 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 77 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBCF 1 -doPost 1 -doGlf 3" #-doGlf set to 3 for binary beagle file
angsd -b bams_yft_merged -sites finalsites_yft_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_ngsrelate_yft

qsub angsd_ngsrelate_yft

# how many SNPs?
NSITES=`zcat angsd_ngsrelate_yft.mafs.gz | wc -l`
echo $NSITES
#5842

zcat angsd_ngsrelate_yft.mafs.gz | cut -f5 |sed 1d >freq_yft 
gunzip -c angsd_ngsrelate_yft.glf.gz |wc -c #13457664

>ngsrelate_yft
nano ngsrelate_yft
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsrelate_yft # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 4
#$ -j y # Join standard output and error to a single file
#$ -o ngsrelate_yft.qlog # Name the file where to redirect standard output and error
module load ngsrelate
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsRelate/ngsRelate -g angsd_ngsrelate_yft.glf.gz -n 96 -f freq_yft -z yft_samps -O yft_relate

qsub ngsrelate_yft

cp yft_relate yft_relate.tsv
#Inbreeding and pairwise relatedness seem very high for the divergent yft cluster. While this points to sibship, the divergence between the two groups likely inflates these estimates. 
#Splitting into main yft cluster and divergent cluster to see if/how relatedness changes

#79 individuals, 80% is 63
>angsd_ngsrelate_yft_main
nano angsd_ngsrelate_yft_main
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_ngsrelate_yft_main # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_ngsrelate_yft_main.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -setMinDepthInd 5 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 63 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBCF 1 -doPost 1 -doGlf 3" #-doGlf set to 3 for binary beagle file
angsd -b bams_yft_merged_main -sites finalsites_yft_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_ngsrelate_yft_main

qsub angsd_ngsrelate_yft_main

# how many SNPs?
NSITES=`zcat angsd_ngsrelate_yft_main.mafs.gz | wc -l`
echo $NSITES
#5271

zcat angsd_ngsrelate_yft_main.mafs.gz | cut -f5 |sed 1d >freq_yft_main
gunzip -c angsd_ngsrelate_yft_main.glf.gz |wc -c #9991920

>ngsrelate_yft_main
nano ngsrelate_yft_main
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsrelate_yft_main # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 4
#$ -j y # Join standard output and error to a single file
#$ -o ngsrelate_yft_main.qlog # Name the file where to redirect standard output and error
module load ngsrelate
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsRelate/ngsRelate -g angsd_ngsrelate_yft_main.glf.gz -n 79 -f freq_yft_main -z yft_samps_main -O yft_relate_main

qsub ngsrelate_yft_main

cp yft_relate_main yft_relate_main.tsv

#Subpop (17 individuals, 80% is 14)
>angsd_ngsrelate_yft_subpop
nano angsd_ngsrelate_yft_subpop
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_ngsrelate_yft_subpop # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_ngsrelate_yft_subpop.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -setMinDepthInd 5 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 14 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBCF 1 -doPost 1 -doGlf 3" #-doGlf set to 3 for binary beagle file
angsd -b bams_yft_merged_subpop -sites finalsites_yft_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_ngsrelate_yft_subpop

qsub angsd_ngsrelate_yft_subpop

# how many SNPs?
NSITES=`zcat angsd_ngsrelate_yft_subpop.mafs.gz | wc -l`
echo $NSITES
#4379 vs 5842 w/ all yft

zcat angsd_ngsrelate_yft_subpop.mafs.gz | cut -f5 |sed 1d >freq_yft_subpop
gunzip -c angsd_ngsrelate_yft_subpop.glf.gz |wc -c #1786224

>ngsrelate_yft_subpop
nano ngsrelate_yft_subpop
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsrelate_yft_subpop # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 4
#$ -j y # Join standard output and error to a single file
#$ -o ngsrelate_yft_subpop.qlog # Name the file where to redirect standard output and error
module load ngsrelate
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsRelate/ngsRelate -g angsd_ngsrelate_yft_subpop.glf.gz -n 17 -f freq_yft_subpop -z yft_samps_subpop -O yft_relate_subpop

qsub ngsrelate_yft_subpop

cp yft_relate_subpop yft_relate_subpop.tsv
#After running this with separately calculated allele frequencies, relatedness for both groups drops significantly, with the divergent sub pop having very low pairwise relatedness, so structure was causing the elevated values before

#Bigeye
>angsd_ngsrelate_bet
nano angsd_ngsrelate_bet
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_ngsrelate_bet # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_ngsrelate_bet.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -setMinDepthInd 5 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 46 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBCF 1 -doPost 1 -doGlf 3" #-doGlf set to 3 for binary beagle file
angsd -b bams_bet_merged -sites finalsites_bet_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_ngsrelate_bet

qsub angsd_ngsrelate_bet

# how many SNPs?
NSITES=`zcat angsd_ngsrelate_bet.mafs.gz | wc -l`
echo $NSITES
#4133

zcat angsd_ngsrelate_bet.mafs.gz | cut -f5 |sed 1d >freq_bet 
gunzip -c angsd_ngsrelate_bet.glf.gz |wc -c #5652576

>ngsrelate_bet
nano ngsrelate_bet
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsrelate_bet # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 4
#$ -j y # Join standard output and error to a single file
#$ -o ngsrelate_bet.qlog # Name the file where to redirect standard output and error
module load ngsrelate
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsRelate/ngsRelate -g angsd_ngsrelate_bet.glf.gz -n 57 -f freq_bet -z bet_samps -O bet_relate

qsub ngsrelate_bet

cp bet_relate bet_relate.tsv

# - - - - - - Removing related individuals w/ pw relatedness > 0.125 (first cousins and closer) - - - - - -
#Note: only needed for skj as yft and bet didn't have any values much greater than 0.06

#Remove the individual from each dyad with the lowest coverage >5x depth

#Pool_10_Jaskiel_Lane2_S10_L004_R1_001_GTGA.trim.bt2.bam 0.40806631193128	OR	Pool_9_Jaskiel_Lane2_S9_L004_R1_001_CTAC.trim.bt2.bam 0.464120547128452
#Pool_10_Jaskiel_Lane2_S10_L004_R1_001_TCAC.trim.bt2.bam 0.428657970232507	OR	Pool_10_Jaskiel_Lane2_S10_L004_R1_001_TGTC.trim.bt2.bam 0.446963859201678
#Pool_12_S12_L002_R1_001_CTAC.trim.bt2.bam 0.408351326671034	OR	Pool_20_S4_L001_R1_001_CATC.trim.bt2.bam 0.593753230473589
#Pool_13_Jaskiel_Lane2_S13_L004_R1_001_ACCA.trim.bt2.bam 0.408619481919043	OR	Pool_13_Jaskiel_Lane2_S13_L004_R1_001_GACT.trim.bt2.bam 0.517922860395852
#Pool_13_Jaskiel_Lane2_S13_L004_R1_001_ACCA.trim.bt2.bam 0.408619481919043	OR	Pool_20_S4_L001_R1_001_GCTT.trim.bt2.bam 0.584041855541275
#Pool_13_Jaskiel_Lane2_S13_L004_R1_001_GCTT.trim.bt2.bam 0.421166854553752	OR	Pool_8_S8_L001_R1_001_TGGT.trim.bt2.bam 0.750525547478687/Pool_9_S9_L002_R1_001_TGGT.trim.bt2.bam 0.711435166858636/281-G203.merged.bam
#Pool_17_S1_L001_R1_001_TCAG.trim.bt2.bam 0.427811543309659	OR	Pool_18_S2_L001_R1_001_TCAG.trim.bt2.bam 0.928676327925476
#Pool_19_S3_L001_R1_001_AGTG.trim.bt2.bam 0.690101875615118	OR	Pool_7_S7_L001_R1_001_CTAC.trim.bt2.bam 0.461787096218706
#Pool_5_S5_L001_R1_001_GACT.trim.bt2.bam 0.668533881808307	OR	Pool_6_S6_L001_R1_001_GACT.trim.bt2.bam 0.898437010071253/Pool_8_S8_L001_R1_001_GACT.trim.bt2.bam 0.895910531598728/281-G60.merged.bam

#remove 268-G3, 281-G8, 287-G45, 274-G5, 268-G5, 268-G96, 281-G92, 274-G51 from bams_all_merged and bams_skj_merged before running ANGSD for downstream analyses
#new bamslists are bams_all_noTR_merged_norel, bams_skj_merged_norel, and bams_all_no_aux_merged_norel

#----------------------------------- Outlier Detection with PCAdapt -------------------------------
#See https://bcm-uga.github.io/pcadapt/articles/pcadapt.html


#running angsd to create necessary files
## - - - - - - - - - - - - - - - - - - - - Just Skipjack - - - - - - - - - - - - - - - - - - -
# minind 156 (80% of 195) minQ 25 w/ -setMinDepthInd 5, LD thinned, no related individuals
>angsd_skj_mid5_noLD_merged_norel
nano angsd_skj_mid5_noLD_merged_norel
#!/bin/bash -l 
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_skj_mid5_noLD_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_skj_mid5_noLD_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -setMinDepthInd 5 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 156 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBCF 1 -doPost 1 -doGlf 2"
angsd -b bams_skj_merged_norel -sites finalsites_skj_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_skj_mid5_noLD_merged_norel

qsub angsd_skj_mid5_noLD_merged_norel

# how many SNPs?
NSITES=`zcat angsd_skj_mid5_noLD_merged_norel.mafs.gz | wc -l`
echo $NSITES
#1326 vs 1262 w/ related skj

# minind 156 (80% of 195) minQ 25 w/ -setMinDepthInd 5, LD thinned, minmaf 0.01
>angsd_skj_mid5_noLD_merged_minmaf01_norel
nano angsd_skj_mid5_noLD_merged_minmaf01_norel
#!/bin/bash -l 
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_skj_mid5_noLD_merged_minmaf01_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_skj_mid5_noLD_merged_minmaf01_norel.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -setMinDepthInd 5 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 156 -snp_pval 1e-5 -minMaf 0.01"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBCF 1 -doPost 1 -doGlf 2"
angsd -b bams_skj_merged_norel -sites finalsites_skj_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_skj_mid5_noLD_merged_minmaf01_norel

qsub angsd_skj_mid5_noLD_merged_minmaf01_norel

# how many SNPs?
NSITES=`zcat angsd_skj_mid5_noLD_merged_minmaf01_norel.mafs.gz | wc -l`
echo $NSITES
#5077 vs 4828 w/ related skj 

## - - - - - - - - - - - - - - - - - - - -  Just Yellowfin - - - - - - - - - - - - - - - - - - -
# minind 77 (80% of 96) minQ 25 w/ -setMinDepthInd 5, LD thinned, merged bams
>angsd_yft_mid5_noLD_merged
nano angsd_yft_mid5_noLD_merged
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_yft_mid5_noLD_merged # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_yft_mid5_noLD_merged.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 77 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_yft_merged -sites finalsites_yft_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_yft_mid5_noLD_merged

qsub angsd_yft_mid5_noLD_merged

NSITES=`zcat angsd_yft_mid5_noLD_merged.mafs.gz | wc -l` 
echo $NSITES
#5842

# minind 77 (80% of 96) minQ 25 w/ -setMinDepthInd 5, LD thinned, minmaf 0.01
>angsd_yft_mid5_noLD_merged_minmaf01
nano angsd_yft_mid5_noLD_merged_minmaf01
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_yft_mid5_noLD_merged_minmaf01 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_yft_mid5_noLD_merged_minmaf01.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 77 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.01"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_yft_merged -sites finalsites_yft_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_yft_mid5_noLD_merged_minmaf01

qsub angsd_yft_mid5_noLD_merged_minmaf01

NSITES=`zcat angsd_yft_mid5_noLD_merged_minmaf01.mafs.gz | wc -l` 
echo $NSITES
#17174

# - - - - - - - - - - - - - - - - - - - - Just Bigeye - - - - - - - - - - - - - - - - - - - - - - - 
# minind 46 (80% of 57) minQ 25 w/ -setMinDepthInd 5, LD thinned, merged bams
>angsd_bet_mid5_noLD_merged
nano angsd_bet_mid5_noLD_merged
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_bet_mid5_noLD_merged # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_bet_mid5_noLD_merged.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 46 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_bet_merged -sites finalsites_bet_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_bet_mid5_noLD_merged

qsub angsd_bet_mid5_noLD_merged

NSITES=`zcat angsd_bet_mid5_noLD_merged.mafs.gz | wc -l` 
echo $NSITES
#4133 

# minind 46 (80% of 57) minQ 25 w/ -setMinDepthInd 5, LD thinned, merged bams, minmaf 0.01
>angsd_bet_mid5_noLD_merged_minmaf01
nano angsd_bet_mid5_noLD_merged_minmaf01
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_bet_mid5_noLD_merged_minmaf01 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_bet_mid5_noLD_merged_minmaf01.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 46 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.01"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_bet_merged -sites finalsites_bet_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_bet_mid5_noLD_merged_minmaf01

qsub angsd_bet_mid5_noLD_merged_minmaf01

NSITES=`zcat angsd_bet_mid5_noLD_merged_minmaf01.mafs.gz | wc -l` 
echo $NSITES
#13317

#Converting bcf --> vcf --> bed, bim, fam files
#general structure
module load bcftools
bcftools convert -O v -o angsd_.vcf angsd_.bcf

#converting to bed format (also creates .fam and .bim files needed in PCAdapt R scripts)
#keep --allow-extra-chr 0 for least CV error in admixture, and remove the zero and redo for input files for PCAdapt
module load plink
module load bcftools

#all_mid5_noLD_merged_noaux_norel
plink --vcf angsd_all_mid5_noLD_merged_noaux_norel.vcf --make-bed --allow-extra-chr --out angsd_all_mid5_noLD_merged_noaux_norel --const-fid 0

#all (merged, noaux, noLD, minmaf01)
bcftools convert -O v -o angsd_all_mid5_noLD_merged_minmaf01_noaux_norel.vcf angsd_all_mid5_noLD_merged_minmaf01_noaux_norel.bcf
plink --vcf angsd_all_mid5_noLD_merged_minmaf01_noaux_norel.vcf --make-bed --allow-extra-chr --out angsd_all_mid5_noLD_merged_minmaf01_noaux_norel --const-fid 0

#skj (merged, noLD, norel)
plink --vcf angsd_skj_mid5_noLD_merged_norel.vcf --make-bed --allow-extra-chr --out angsd_skj_mid5_noLD_merged_norel --const-fid 0

#skj (merged, noLD, norel, minmaf01)
bcftools convert -O v -o angsd_skj_mid5_noLD_merged_minmaf01_norel.vcf angsd_skj_mid5_noLD_merged_minmaf01_norel.bcf
plink --vcf angsd_skj_mid5_noLD_merged_minmaf01_norel.vcf --make-bed --allow-extra-chr --out angsd_skj_mid5_noLD_merged_minmaf01_norel --const-fid 0

#yft (merged, noLD)
plink --vcf angsd_yft_mid5_noLD_merged.vcf --make-bed --allow-extra-chr --out angsd_yft_mid5_noLD_merged --const-fid 0

#yft (merged, noLD, minmaf01)
bcftools convert -O v -o angsd_yft_mid5_noLD_merged_minmaf01.vcf angsd_yft_mid5_noLD_merged_minmaf01.bcf
plink --vcf angsd_yft_mid5_noLD_merged_minmaf01.vcf --make-bed --allow-extra-chr --out angsd_yft_mid5_noLD_merged_minmaf01 --const-fid 0

#bet (merged, noLD)
plink --vcf angsd_bet_mid5_noLD_merged.vcf --make-bed --allow-extra-chr --out angsd_bet_mid5_noLD_merged --const-fid 0

#bet (merged, noLD, minmaf01)
bcftools convert -O v -o angsd_bet_mid5_noLD_merged_minmaf01.vcf angsd_bet_mid5_noLD_merged_minmaf01.bcf
plink --vcf angsd_bet_mid5_noLD_merged_minmaf01.vcf --make-bed --allow-extra-chr --out angsd_bet_mid5_noLD_merged_minmaf01 --const-fid 0

#scp bed, bim, and fam files and use PCAdapt code in R to produce txt files for outlier sites. Then use that and the list of pruned sites to make lists of neutral sites and run ANGSD for PCA etc.

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#Here I'm working with the .txt files produced in the R script for PCadapt. First, I trim them to just be the location and the length:
zgrep "^[^#]" outliers_all_mid5_noLD_merged_noaux_norel_k4.txt | awk '{print $2,$3}' >outlier_sites_all_noLD_merged_noaux_norel_k4.sites

#Then remove the quotation marks around the contig names
sed 's/"//g' outlier_sites_all_noLD_merged_noaux_norel_k4.sites > outliers_all_noLD_merged_noaux_norel_k4.sites
awk '{print $1"\t"$2}' outliers_all_noLD_merged_noaux_norel_k4.sites > outliers_all_noLD_merged_noaux_norel_k4
module load angsd 
angsd sites index outliers_all_noLD_merged_noaux_norel_k4

#Subtract outliers from LD-pruned allsites file
grep -Fv -f outliers_all_noLD_merged_noaux_norel_k4 finalsites_all_unlinked > neutral_sites_all_unlinked #check number of lines to make sure it looks right
angsd sites index neutral_sites_all_unlinked

#Repeating for all outlier lists
#All Sites, minmaf01
zgrep "^[^#]" outliers_all_mid5_noLD_merged_minmaf01_noaux_norel_k4.txt | awk '{print $2,$3}' >outlier_sites_all_noLD_merged_minmaf01_noaux_norel_k4.sites

#Then remove the quotation marks around the contig names
sed 's/"//g' outlier_sites_all_noLD_merged_minmaf01_noaux_norel_k4.sites > outliers_all_noLD_merged_minmaf01_noaux_norel_k4.sites
awk '{print $1"\t"$2}' outliers_all_noLD_merged_minmaf01_noaux_norel_k4.sites > outliers_all_noLD_merged_minmaf01_noaux_norel_k4
angsd sites index outliers_all_noLD_merged_minmaf01_noaux_norel_k4

#Subtract outliers from LD-pruned allsites file
grep -Fv -f outliers_all_noLD_merged_minmaf01_noaux_norel_k4 finalsites_all_unlinked > neutral_sites_all_unlinked_minmaf01 #check number of lines to make sure it looks right
angsd sites index neutral_sites_all_unlinked_minmaf01

####Note: no skipjack outlier sites on the ref run

#Yellowfin
zgrep "^[^#]" outliers_yft_mid5_noLD_merged_k2.txt | awk '{print $2,$3}' >outlier_sites_yft_noLD_merged_k2.sites

#Then remove the quotation marks around the contig names
sed 's/"//g' outlier_sites_yft_noLD_merged_k2.sites > outliers_yft_noLD_merged_k2.sites
awk '{print $1"\t"$2}' outliers_yft_noLD_merged_k2.sites > outliers_yft_noLD_merged_k2
angsd sites index outliers_yft_noLD_merged_k2

#Subtract outliers from LD-pruned allsites file
grep -Fv -f outliers_yft_noLD_merged_k2 finalsites_yft_unlinked > neutral_sites_yft_unlinked #check number of lines to make sure it looks right
angsd sites index neutral_sites_yft_unlinked

#Yellowfin (minmaf01)
zgrep "^[^#]" outliers_yft_mid5_noLD_merged_minmaf01_k2.txt | awk '{print $2,$3}' >outlier_sites_yft_noLD_merged_minmaf01_k2.sites

#Then remove the quotation marks around the contig names
sed 's/"//g' outlier_sites_yft_noLD_merged_minmaf01_k2.sites > outliers_yft_noLD_merged_minmaf01_k2.sites
awk '{print $1"\t"$2}' outliers_yft_noLD_merged_minmaf01_k2.sites > outliers_yft_noLD_merged_minmaf01_k2
angsd sites index outliers_yft_noLD_merged_minmaf01_k2

#Subtract outliers from LD-pruned allsites file
grep -Fv -f outliers_yft_noLD_merged_minmaf01_k2 finalsites_yft_unlinked > neutral_sites_yft_unlinked_minmaf01 #check number of lines to make sure it looks right
angsd sites index neutral_sites_yft_unlinked_minmaf01

#Bigeye
zgrep "^[^#]" outliers_bet_mid5_noLD_merged_k1.txt | awk '{print $2,$3}' >outlier_sites_bet_noLD_merged_k1.sites

#Then remove the quotation marks around the contig names
sed 's/"//g' outlier_sites_bet_noLD_merged_k1.sites > outliers_bet_noLD_merged_k1.sites
awk '{print $1"\t"$2}' outliers_bet_noLD_merged_k1.sites > outliers_bet_noLD_merged_k1
angsd sites index outliers_bet_noLD_merged_k1

#Subtract outliers from LD-pruned allsites file
grep -Fv -f outliers_bet_noLD_merged_k1 finalsites_bet_unlinked > neutral_sites_bet_unlinked #check number of lines to make sure it looks right
angsd sites index neutral_sites_bet_unlinked

#Bigeye (minmaf01)
zgrep "^[^#]" outliers_bet_mid5_noLD_merged_minmaf01_k1.txt | awk '{print $2,$3}' >outlier_sites_bet_noLD_merged_minmaf01_k1.sites

#Then remove the quotation marks around the contig names
sed 's/"//g' outlier_sites_bet_noLD_merged_minmaf01_k1.sites > outliers_bet_noLD_merged_minmaf01_k1.sites
awk '{print $1"\t"$2}' outliers_bet_noLD_merged_minmaf01_k1.sites > outliers_bet_noLD_merged_minmaf01_k1
angsd sites index outliers_bet_noLD_merged_minmaf01_k1

#Subtract outliers from LD-pruned allsites file
grep -Fv -f outliers_bet_noLD_merged_minmaf01_k1 finalsites_bet_unlinked > neutral_sites_bet_unlinked_minmaf01 #check number of lines to make sure it looks right
angsd sites index neutral_sites_bet_unlinked_minmaf01

#------------- ANGSD Runs with Neutral and Outlier Sites on Individual Species clusters --------------
#note: skipjack had no outlier sites detected so the ANGSD run(s) for that group won't have "neutral" in the name, but they are in fact neutral sites only

## - - - - - - - - - - - - - - - - - - - - Just Skipjack - - - - - - - - - - - - - - - - - - -
# minind 156 (80% of 195) minQ 25 w/ -setMinDepthInd 5, LD thinned, no related individuals
>angsd_skj_mid5_noLD_merged_norel
nano angsd_skj_mid5_noLD_merged_norel
#!/bin/bash -l 
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_skj_mid5_noLD_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_skj_mid5_noLD_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -setMinDepthInd 5 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 156 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBCF 1 -doPost 1 -doGlf 2"
angsd -b bams_skj_merged_norel -sites finalsites_skj_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_skj_mid5_noLD_merged_norel

qsub angsd_skj_mid5_noLD_merged_norel

# how many SNPs?
NSITES=`zcat angsd_skj_mid5_noLD_merged_norel.mafs.gz | wc -l`
echo $NSITES
#1326 vs 1262 w/ related skj

# minind 156 (80% of 195) minQ 25 w/ -setMinDepthInd 5, LD thinned, minmaf 0.01
>angsd_skj_mid5_noLD_merged_minmaf01_norel
nano angsd_skj_mid5_noLD_merged_minmaf01_norel
#!/bin/bash -l 
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_skj_mid5_noLD_merged_minmaf01_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_skj_mid5_noLD_merged_minmaf01_norel.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -setMinDepthInd 5 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 156 -snp_pval 1e-5 -minMaf 0.01"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBCF 1 -doPost 1 -doGlf 2"
angsd -b bams_skj_merged_norel -sites finalsites_skj_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_skj_mid5_noLD_merged_minmaf01_norel

qsub angsd_skj_mid5_noLD_merged_minmaf01_norel

# how many SNPs?
NSITES=`zcat angsd_skj_mid5_noLD_merged_minmaf01_norel.mafs.gz | wc -l`
echo $NSITES
#5077 vs 4828 w/ related skj 

#Without LD filtered sites
# minind 156 (80% of 195) minQ 25 w/ -setMinDepthInd 5, NOT LD thinned, no related individuals
>angsd_skj_mid5_merged_norel
nano angsd_skj_mid5_merged_norel
#!/bin/bash -l 
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_skj_mid5_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_skj_mid5_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -setMinDepthInd 5 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 156 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBCF 1 -doPost 1 -doGlf 2"
angsd -b bams_skj_merged_norel -sites AllSites_skj_merged -GL 1 $FILTERS $TODO -P 1 -out angsd_skj_mid5_merged_norel

qsub angsd_skj_mid5_merged_norel

# how many SNPs?
NSITES=`zcat angsd_skj_mid5_merged_norel.mafs.gz | wc -l`
echo $NSITES
#1447 vs 1326 (lost 121)

#-----------------------------------------------------------------------------------------
#Yellowfin - Neutral
# minind 77 (80% of 96) minQ 25 w/ -setMinDepthInd 5, LD thinned, merged bams
>angsd_yft_mid5_noLD_merged_neutral
nano angsd_yft_mid5_noLD_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_yft_mid5_noLD_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_yft_mid5_noLD_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 77 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_yft_merged -sites neutral_sites_yft_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_yft_mid5_noLD_merged_neutral

qsub angsd_yft_mid5_noLD_merged_neutral

NSITES=`zcat angsd_yft_mid5_noLD_merged_neutral.mafs.gz | wc -l` 
echo $NSITES
#3577 vs 5842 with all sites

# minind 77 (80% of 96) minQ 25 w/ -setMinDepthInd 5, LD thinned, minmaf 0.01
>angsd_yft_mid5_noLD_merged_minmaf01_neutral
nano angsd_yft_mid5_noLD_merged_minmaf01_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_yft_mid5_noLD_merged_minmaf01_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_yft_mid5_noLD_merged_minmaf01_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 77 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.01"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_yft_merged -sites neutral_sites_yft_unlinked_minmaf01 -GL 1 $FILTERS $TODO -P 1 -out angsd_yft_mid5_noLD_merged_minmaf01_neutral

qsub angsd_yft_mid5_noLD_merged_minmaf01_neutral

NSITES=`zcat angsd_yft_mid5_noLD_merged_minmaf01_neutral.mafs.gz | wc -l` 
echo $NSITES
#14889 vs 17174 with all sites

#Yellowfin - Outlier
# minind 77 (80% of 96) minQ 25 w/ -setMinDepthInd 5, LD thinned, merged bams
>angsd_yft_mid5_noLD_merged_outlier
nano angsd_yft_mid5_noLD_merged_outlier
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_yft_mid5_noLD_merged_outlier # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_yft_mid5_noLD_merged_outlier.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 77 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_yft_merged -sites outliers_yft_noLD_merged_k2 -GL 1 $FILTERS $TODO -P 1 -out angsd_yft_mid5_noLD_merged_outlier

qsub angsd_yft_mid5_noLD_merged_outlier

NSITES=`zcat angsd_yft_mid5_noLD_merged_outlier.mafs.gz | wc -l` 
echo $NSITES
#2266 vs 5842 with all sites

# minind 77 (80% of 96) minQ 25 w/ -setMinDepthInd 5, LD thinned, minmaf 0.01
>angsd_yft_mid5_noLD_merged_minmaf01_outlier
nano angsd_yft_mid5_noLD_merged_minmaf01_outlier
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_yft_mid5_noLD_merged_minmaf01_outlier # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_yft_mid5_noLD_merged_minmaf01_outlier.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 77 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.01"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_yft_merged -sites outliers_yft_noLD_merged_minmaf01_k2 -GL 1 $FILTERS $TODO -P 1 -out angsd_yft_mid5_noLD_merged_minmaf01_outlier

qsub angsd_yft_mid5_noLD_merged_minmaf01_outlier

NSITES=`zcat angsd_yft_mid5_noLD_merged_minmaf01_outlier.mafs.gz | wc -l` 
echo $NSITES
#2286 vs 17174 with all sites

#Without LD filters
# minind 77 (80% of 96) minQ 25 w/ -setMinDepthInd 5, NOT LD thinned, merged bams
>angsd_yft_mid5_merged_neutral
nano angsd_yft_mid5_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_yft_mid5_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_yft_mid5_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 77 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_yft_merged -sites AllSites_yft_merged -GL 1 $FILTERS $TODO -P 1 -out angsd_yft_mid5_merged_neutral

qsub angsd_yft_mid5_merged_neutral

NSITES=`zcat angsd_yft_mid5_merged_neutral.mafs.gz | wc -l` 
echo $NSITES
#6757 vs 3577 (lost 3,180)

#-----------------------------------------------------------------------------------------
#Bigeye - Neutral
# minind 46 (80% of 57) minQ 25 w/ -setMinDepthInd 5, LD thinned, merged bams
>angsd_bet_mid5_noLD_merged_neutral
nano angsd_bet_mid5_noLD_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_bet_mid5_noLD_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_bet_mid5_noLD_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 46 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_bet_merged -sites neutral_sites_bet_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_bet_mid5_noLD_merged_neutral

qsub angsd_bet_mid5_noLD_merged_neutral

NSITES=`zcat angsd_bet_mid5_noLD_merged_neutral.mafs.gz | wc -l` 
echo $NSITES
#4130 vs 4133 with all sites

# minind 46 (80% of 57) minQ 25 w/ -setMinDepthInd 5, LD thinned, merged bams, minmaf 0.01
>angsd_bet_mid5_noLD_merged_minmaf01_neutral
nano angsd_bet_mid5_noLD_merged_minmaf01_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_bet_mid5_noLD_merged_minmaf01_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_bet_mid5_noLD_merged_minmaf01_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 46 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.01"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_bet_merged -sites neutral_sites_bet_unlinked_minmaf01 -GL 1 $FILTERS $TODO -P 1 -out angsd_bet_mid5_noLD_merged_minmaf01_neutral

qsub angsd_bet_mid5_noLD_merged_minmaf01_neutral

NSITES=`zcat angsd_bet_mid5_noLD_merged_minmaf01_neutral.mafs.gz | wc -l` 
echo $NSITES
#13315 vs 13317 with all sites

#Bigeye - Outlier
# minind 46 (80% of 57) minQ 25 w/ -setMinDepthInd 5, LD thinned, merged bams
>angsd_bet_mid5_noLD_merged_outlier
nano angsd_bet_mid5_noLD_merged_outlier
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_bet_mid5_noLD_merged_outlier # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_bet_mid5_noLD_merged_outlier.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 46 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_bet_merged -sites outliers_bet_noLD_merged_k1 -GL 1 $FILTERS $TODO -P 1 -out angsd_bet_mid5_noLD_merged_outlier

qsub angsd_bet_mid5_noLD_merged_outlier

NSITES=`zcat angsd_bet_mid5_noLD_merged_outlier.mafs.gz | wc -l` 
echo $NSITES
#4 vs 4133 with all sites

# minind 46 (80% of 57) minQ 25 w/ -setMinDepthInd 5, LD thinned, merged bams, minmaf 0.01
>angsd_bet_mid5_noLD_merged_minmaf01_outlier
nano angsd_bet_mid5_noLD_merged_minmaf01_outlier
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_bet_mid5_noLD_merged_minmaf01_outlier # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_bet_mid5_noLD_merged_minmaf01_outlier.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 46 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.01"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_bet_merged -sites outliers_bet_noLD_merged_minmaf01_k1 -GL 1 $FILTERS $TODO -P 1 -out angsd_bet_mid5_noLD_merged_minmaf01_outlier

qsub angsd_bet_mid5_noLD_merged_minmaf01_outlier

NSITES=`zcat angsd_bet_mid5_noLD_merged_minmaf01_outlier.mafs.gz | wc -l` 
echo $NSITES
#3 vs 13317 with all sites

#Not LD filtered
# minind 46 (80% of 57) minQ 25 w/ -setMinDepthInd 5, NOT LD thinned, merged bams
>angsd_bet_mid5_merged_neutral
nano angsd_bet_mid5_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_bet_mid5_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_bet_mid5_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 46 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_bet_merged -sites AllSites_bet_merged -GL 1 $FILTERS $TODO -P 1 -out angsd_bet_mid5_merged_neutral

qsub angsd_bet_mid5_merged_neutral

NSITES=`zcat angsd_bet_mid5_merged_neutral.mafs.gz | wc -l` 
echo $NSITES
#4580 vs 4130 (lost 450)

#=========================================================================================
##--NGSadmix and ADMIXTURE (ind. species cluster, LD thinned, no related ind., neutral sites)
pwd
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis

# First, install:
wget popgen.dk/software/download/NGSadmix/ngsadmix32.cpp 
g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix

#Skipjack
>ngsadmix_skj_mid5_noLD_merged_norel
nano ngsadmix_skj_mid5_noLD_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsadmix_skj_mid5_noLD_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o ngsadmix_skj_mid5_noLD_merged_norel.qlog # Name the file where to redirect standard output and error
mkdir /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_skj_output_files
for K in `seq 1 5` ;
do
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/NGSadmix -likes angsd_skj_mid5_noLD_merged_norel.beagle.gz -K $K -P 1 -o /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_skj_output_files/ngsadmix_skj_mid5_noLD_merged_norel_k${K};
done

qsub ngsadmix_skj_mid5_noLD_merged_norel

# scp *qopt files to local machine, use admixturePlotting_v5.R to plot
bcftools convert -O v -o angsd_skj_mid5_noLD_merged_norel.vcf angsd_skj_mid5_noLD_merged_norel.bcf
plink --vcf angsd_skj_mid5_noLD_merged_norel.vcf --make-bed --allow-extra-chr 0 --out angsd_skj_mid5_noLD_merged_norel --const-fid 0

#Don't run this and remove once you're sure you don't need it
for K in `seq 1 5`; \
do admixture --cv angsd_skj_mid5_noLD_merged_norel.bed $K | tee myresult_skj_mid5_noLD_merged_norel_${K}.out; done
# use this to check K of least CV error:
grep -h CV myresult_skj_mid5_noLD_merged_norel_*.out

#Trying multiple replicates per K and Evanno method for better confidence in K
>admixture_skj_mid5_noLD_merged_norel
nano admixture_skj_mid5_noLD_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N admixture_skj_mid5_noLD_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o admixture_skj_mid5_noLD_merged_norel.qlog # Name the file where to redirect standard output and error
module load admixture/1.3.0
for K in $(seq 1 5); do
  for REP in $(seq 1 10); do
    SEED=$((1000 + $RANDOM % 100000))  # generates a random 4–5 digit number
    admixture --cv=10 --seed=$SEED angsd_skj_mid5_noLD_merged_norel.bed $K \
      | tee /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_skj_output_files/myresult_skj_mid5_noLD_merged_norel_K${K}_rep${REP}.out
  done
done
grep "CV error" /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_skj_output_files/myresult_skj_mid5_noLD_merged_norel_K*_rep*.out > /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_skj_output_files/admixture_cv_errors_skj_noLD_merged_norel.txt

qsub admixture_skj_mid5_noLD_merged_norel
#See Jaskiel_2025_2bRAD_Final_Workflow.Rmd for calculation of most probable K with Evanno's method (Evanno et al., 2005)

#-----------------------------------Yellowfin (neutral sites)
>ngsadmix_yft_mid5_noLD_merged_neutral
nano ngsadmix_yft_mid5_noLD_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsadmix_yft_mid5_noLD_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o ngsadmix_yft_mid5_noLD_merged_neutral.qlog # Name the file where to redirect standard output and error
for K in `seq 1 5` ;
do
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/NGSadmix -likes angsd_yft_mid5_noLD_merged_neutral.beagle.gz -K $K -P 1 -o /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_yft_output_files/ngsadmix_yft_mid5_noLD_merged_neutral_k${K};
done

qsub ngsadmix_yft_mid5_noLD_merged_neutral

# scp *qopt files to local machine, use admixturePlotting_v5.R to plot
bcftools convert -O v -o angsd_yft_mid5_noLD_merged_neutral.vcf angsd_yft_mid5_noLD_merged_neutral.bcf
plink --vcf angsd_yft_mid5_noLD_merged_neutral.vcf --make-bed --allow-extra-chr 0 --out angsd_yft_mid5_noLD_merged_neutral --const-fid 0

#Don't run this and remove once you're sure you don't need it
for K in `seq 1 5`; \
do admixture --cv angsd_yft_mid5_noLD_merged_neutral.bed $K | tee myresult_yft_mid5_noLD_merged_neutral_${K}.out; done
# use this to check K of least CV error:
grep -h CV myresult_yft_mid5_noLD_merged_neutral_*.out

#Trying multiple replicates per K and Evanno method for better confidence in K
>admixture_yft_mid5_noLD_merged_neutral
nano admixture_yft_mid5_noLD_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N admixture_yft_mid5_noLD_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o admixture_yft_mid5_noLD_merged_neutral.qlog # Name the file where to redirect standard output and error
module load admixture/1.3.0
for K in $(seq 1 5); do
  for REP in $(seq 1 10); do
    SEED=$((1000 + $RANDOM % 100000))  # generates a random 4–5 digit number
    admixture --cv=10 --seed=$SEED angsd_yft_mid5_noLD_merged_neutral.bed $K \
      | tee /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_yft_output_files/myresult_yft_mid5_noLD_merged_neutral_K${K}_rep${REP}.out
  done
done
grep "CV error" /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_yft_output_files/myresult_yft_mid5_noLD_merged_neutral_K*_rep*.out > /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_yft_output_files/admixture_cv_errors_yft_noLD_merged_neutral.txt

qsub admixture_yft_mid5_noLD_merged_neutral
#See Jaskiel_2025_2bRAD_Final_Workflow.Rmd for calculation of most probable K with Evanno's method (Evanno et al., 2005)

#----------------------------------Bigeye (Neutral sites)
>ngsadmix_bet_mid5_noLD_merged_neutral
nano ngsadmix_bet_mid5_noLD_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsadmix_bet_mid5_noLD_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o ngsadmix_bet_mid5_noLD_merged_neutral.qlog # Name the file where to redirect standard output and error
for K in `seq 1 5` ;
do
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/NGSadmix -likes angsd_bet_mid5_noLD_merged_neutral.beagle.gz -K $K -P 1 -o /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_bet_output_files/ngsadmix_bet_mid5_noLD_merged_neutral_k${K};
done

qsub ngsadmix_bet_mid5_noLD_merged_neutral

# scp *qopt files to local machine, use admixturePlotting_v5.R to plot
bcftools convert -O v -o angsd_bet_mid5_noLD_merged_neutral.vcf angsd_bet_mid5_noLD_merged_neutral.bcf
plink --vcf angsd_bet_mid5_noLD_merged_neutral.vcf --make-bed --allow-extra-chr 0 --out angsd_bet_mid5_noLD_merged_neutral --const-fid 0

#Don't run this and remove once you're sure you don't need it
for K in `seq 1 5`; \
do admixture --cv angsd_bet_mid5_noLD_merged_neutral.bed $K | tee myresult_bet_mid5_noLD_merged_neutral_${K}.out; done
# use this to check K of least CV error:
grep -h CV myresult_bet_mid5_noLD_merged_neutral_*.out

#Trying multiple replicates per K and Evanno method for better confidence in K
>admixture_bet_mid5_noLD_merged_neutral
nano admixture_bet_mid5_noLD_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N admixture_bet_mid5_noLD_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o admixture_bet_mid5_noLD_merged_neutral.qlog # Name the file where to redirect standard output and error
module load admixture/1.3.0
for K in $(seq 1 5); do
  for REP in $(seq 1 10); do
    SEED=$((1000 + $RANDOM % 100000))  # generates a random 4–5 digit number
    admixture --cv=10 --seed=$SEED angsd_bet_mid5_noLD_merged_neutral.bed $K \
      | tee /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_bet_output_files/myresult_bet_mid5_noLD_merged_neutral_K${K}_rep${REP}.out
  done
done
grep "CV error" /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_bet_output_files/myresult_bet_mid5_noLD_merged_neutral_K*_rep*.out > /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/ngsadmix_bet_output_files/admixture_cv_errors_bet_noLD_merged_neutral.txt

qsub admixture_bet_mid5_noLD_merged_neutral
#See Jaskiel_2025_2bRAD_Final_Workflow.Rmd for calculation of most probable K with Evanno's method (Evanno et al., 2005)

# ======================================= Basic Summary Stats ======================================
# - - - - - - - - - - - - - - - - - - Heterozygosity (neutral sites) - - - - - - - - - - - - - - - -
# Skipjack (195 individuals, 80% is 156), LD thinned sites, merged bams, no relatives
>angsd_het_skj_noLD_merged_norel
nano angsd_het_skj_noLD_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_het_skj_noLD_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_het_skj_noLD_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 156"
TODO="-makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doPost 1 -doGlf 2 -doCounts 1 -doMajorMinor 1 -dosnpstat 1 -doMaf 1"
angsd -bam bams_skj_merged_norel -sites finalsites_skj_unlinked -GL 1 -P 1 $TODO $FILTERS -out angsd_het_skj_noLD_merged_norel

qsub angsd_het_skj_noLD_merged_norel

# scp .beagle.gz and use R script to calculate heterozygosity and plot

# Yellowfin (96 individuals, 80% is 77), LD thinned sites, merged bams
>angsd_het_yft_noLD_merged_neutral
nano angsd_het_yft_noLD_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_het_yft_noLD_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_het_yft_noLD_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 77"
TODO="-makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doPost 1 -doGlf 2 -doCounts 1 -doMajorMinor 1 -dosnpstat 1 -doMaf 1"
angsd -bam bams_yft_merged -sites neutral_sites_yft_unlinked -GL 1 -P 1 $TODO $FILTERS -out angsd_het_yft_noLD_merged_neutral

qsub angsd_het_yft_noLD_merged_neutral

#Yellowfin Adults 
>angsd_het_yft_noLD_merged_neutral_adults
nano angsd_het_yft_noLD_merged_neutral_adults
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_het_yft_noLD_merged_neutral_adults # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_het_yft_noLD_merged_neutral_adults.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 23"
TODO="-makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doPost 1 -doGlf 2 -doCounts 1 -doMajorMinor 1 -dosnpstat 1 -doMaf 1"
angsd -bam bams_yft_merged_2020 -sites neutral_sites_yft_unlinked -GL 1 -P 1 $TODO $FILTERS -out angsd_het_yft_noLD_merged_neutral_adults

qsub angsd_het_yft_noLD_merged_neutral_adults

#Yellowfin Larvae
>angsd_het_yft_noLD_merged_neutral_larvae
nano angsd_het_yft_noLD_merged_neutral_larvae
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_het_yft_noLD_merged_neutral_larvae # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_het_yft_noLD_merged_neutral_larvae.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 54"
TODO="-makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doPost 1 -doGlf 2 -doCounts 1 -doMajorMinor 1 -dosnpstat 1 -doMaf 1"
angsd -bam bams_yft_merged_larvae -sites neutral_sites_yft_unlinked -GL 1 -P 1 $TODO $FILTERS -out angsd_het_yft_noLD_merged_neutral_larvae

qsub angsd_het_yft_noLD_merged_neutral_larvae

# Bigeye (57 individuals, 80% is 46) LD thinned, merged bams
>angsd_het_bet_noLD_merged_neutral
nano angsd_het_bet_noLD_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_het_bet_noLD_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_het_bet_noLD_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-maxHetFreq 0.5 -uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 46"
TODO="-makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doPost 1 -doGlf 2 -doCounts 1 -doMajorMinor 1 -dosnpstat 1 -doMaf 1"
angsd -bam bams_bet_merged -sites neutral_sites_bet_unlinked -GL 1 -P 1 $TODO $FILTERS -out angsd_het_bet_noLD_merged_neutral

qsub angsd_het_bet_noLD_merged_neutral

#Bigeye Adults
>angsd_het_bet_noLD_merged_neutral_adults
nano angsd_het_bet_noLD_merged_neutral_adults
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_het_bet_noLD_merged_neutral_adults # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_het_bet_noLD_merged_neutral_adults.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-maxHetFreq 0.5 -uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 19"
TODO="-makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doPost 1 -doGlf 2 -doCounts 1 -doMajorMinor 1 -dosnpstat 1 -doMaf 1"
angsd -bam bams_bet_merged_2020 -sites neutral_sites_bet_unlinked -GL 1 -P 1 $TODO $FILTERS -out angsd_het_bet_noLD_merged_neutral_adults

qsub angsd_het_bet_noLD_merged_neutral_adults

#Bigeye larvae
>angsd_het_bet_noLD_merged_neutral_larvae
nano angsd_het_bet_noLD_merged_neutral_larvae
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_het_bet_noLD_merged_neutral_larvae # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_het_bet_noLD_merged_neutral_larvae.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-maxHetFreq 0.5 -uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 26"
TODO="-makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doPost 1 -doGlf 2 -doCounts 1 -doMajorMinor 1 -dosnpstat 1 -doMaf 1"
angsd -bam bams_bet_merged_larvae -sites neutral_sites_bet_unlinked -GL 1 -P 1 $TODO $FILTERS -out angsd_het_bet_noLD_merged_neutral_larvae

qsub angsd_het_bet_noLD_merged_neutral_larvae

#Yellowfin Main pop
>angsd_het_yft_noLD_merged_neutral_main
nano angsd_het_yft_noLD_merged_neutral_main
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_het_yft_noLD_merged_neutral_main # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_het_yft_noLD_merged_neutral_main.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 63"
TODO="-makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doPost 1 -doGlf 2 -doCounts 1 -doMajorMinor 1 -dosnpstat 1 -doMaf 1"
angsd -bam bams_yft_merged_main -sites neutral_sites_yft_unlinked -GL 1 -P 1 $TODO $FILTERS -out angsd_het_yft_noLD_merged_neutral_main

qsub angsd_het_yft_noLD_merged_neutral_main

#Yellowfin Subpop
>angsd_het_yft_noLD_merged_neutral_subpop
nano angsd_het_yft_noLD_merged_neutral_subpop
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_het_yft_noLD_merged_neutral_subpop # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_het_yft_noLD_merged_neutral_subpop.qlog # Name the file where to redirect standard output and error
module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 14"
TODO="-makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doPost 1 -doGlf 2 -doCounts 1 -doMajorMinor 1 -dosnpstat 1 -doMaf 1"
angsd -bam bams_yft_merged_subpop -sites neutral_sites_yft_unlinked -GL 1 -P 1 $TODO $FILTERS -out angsd_het_yft_noLD_merged_neutral_subpop

qsub angsd_het_yft_noLD_merged_neutral_subpop

#------------------------------- Pairwise Fst between years of sampling within species (neutral sites) --------------------------------------
#Skipjack (ref)
#Doing SAF/SFS work
#2015 (17 individuals, so 80% is 14)
>sfs_skj_2015_merged_norel
nano sfs_skj_2015_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2015_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2015_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_skj_unlinked -b bams_skj_merged_norel_2015 -GL 1 -P 1 -minInd 14 $TODO -out skj_2015_merged_norel

qsub sfs_skj_2015_merged_norel

#2016 (66 individuals so 80% is 53)
>sfs_skj_2016_merged_norel
nano sfs_skj_2016_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2016_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2016_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_skj_unlinked -b bams_skj_merged_norel_2016 -GL 1 -P 1 -minInd 53 $TODO -out skj_2016_merged_norel

qsub sfs_skj_2016_merged_norel

#2017 (25 individuals so 80% is 20)
>sfs_skj_2017_merged_norel
nano sfs_skj_2017_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2017_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2017_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_skj_unlinked -b bams_skj_merged_norel_2017 -GL 1 -P 1 -minInd 20 $TODO -out skj_2017_merged_norel

qsub sfs_skj_2017_merged_norel

#2018 (71 individuals so 80% is 57)
>sfs_skj_2018_merged_norel
nano sfs_skj_2018_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2018_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2018_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_skj_unlinked -b bams_skj_merged_norel_2018 -GL 1 -P 1 -minInd 57 $TODO -out skj_2018_merged_norel

qsub sfs_skj_2018_merged_norel

#2019 (15 individuals so 80% is 12)
>sfs_skj_2019_merged_norel
nano sfs_skj_2019_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2019_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2019_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_skj_unlinked -b bams_skj_merged_norel_2019 -GL 1 -P 1 -minInd 12 $TODO -out skj_2019_merged_norel

qsub sfs_skj_2019_merged_norel

#2022 (only 1 individual so no minind)
>sfs_skj_2022_merged_norel
nano sfs_skj_2022_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2022_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2022_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_skj_unlinked -b bams_skj_2022 -GL 1 -P 1 -minInd 0 $TODO -out skj_2022_merged_norel

qsub sfs_skj_2022_merged_norel

# Writing down 2d-SFS priors - 2015 & 2016
>sfs_skj_2015_2016_merged_norel
nano sfs_skj_2015_2016_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2015_2016_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2015_2016_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS skj_2015_merged_norel.saf.idx skj_2016_merged_norel.saf.idx -fold 1 -P 1 > skj_2015_2016_merged_norel.sfs ; realSFS fst index skj_2015_merged_norel.saf.idx skj_2016_merged_norel.saf.idx -sfs skj_2015_2016_merged_norel.sfs -fstout skj_2015_2016_merged_norel

qsub sfs_skj_2015_2016_merged_norel

# Global Fst between populations 
realSFS fst stats skj_2015_2016_merged_norel.fst.idx
# Output:
FST.Unweight[nObs:26113]:0.006550 Fst.Weight:0.013112

# The difference between unweighted and weighted values is averaging versus ratio of sums method. 
# The weighted value is derived from the ratio of the separate per-locus sums of numerator and denominator values, 
# while the unweighted value is the average of per-locus values. [If that is unclear: weighted is sum(a)/sum(a+b), while unweighted is average(a/(a+b))].

# Writing down 2d-SFS priors - 2015 & 2017
>sfs_skj_2015_2017_merged_norel
nano sfs_skj_2015_2017_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2015_2017_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2015_2017_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS skj_2015_merged_norel.saf.idx skj_2017_merged_norel.saf.idx -fold 1 -P 1 > skj_2015_2017_merged_norel.sfs ; realSFS fst index skj_2015_merged_norel.saf.idx skj_2017_merged_norel.saf.idx -sfs skj_2015_2017_merged_norel.sfs -fstout skj_2015_2017_merged_norel

qsub sfs_skj_2015_2017_merged_norel

# Global Fst between populations 
realSFS fst stats skj_2015_2017_merged_norel.fst.idx
# Output:
FST.Unweight[nObs:26534]:0.010661 Fst.Weight:0.011377

# Writing down 2d-SFS priors - 2015 & 2018
>sfs_skj_2015_2018_merged_norel
nano sfs_skj_2015_2018_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2015_2018_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2015_2018_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS skj_2015_merged_norel.saf.idx skj_2018_merged_norel.saf.idx -fold 1 -P 1 > skj_2015_2018_merged_norel.sfs ; realSFS fst index skj_2015_merged_norel.saf.idx skj_2018_merged_norel.saf.idx -sfs skj_2015_2018_merged_norel.sfs -fstout skj_2015_2018_merged_norel

qsub sfs_skj_2015_2018_merged_norel

# Global Fst between populations 
realSFS fst stats skj_2015_2018_merged_norel.fst.idx
# Output:
FST.Unweight[nObs:26138]:0.006809 Fst.Weight:0.010103

# Writing down 2d-SFS priors - 2015 & 2019
>sfs_skj_2015_2019_merged_norel
nano sfs_skj_2015_2019_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2015_2019_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2015_2019_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS skj_2015_merged_norel.saf.idx skj_2019_merged_norel.saf.idx -fold 1 -P 1 > skj_2015_2019_merged_norel.sfs ; realSFS fst index skj_2015_merged_norel.saf.idx skj_2019_merged_norel.saf.idx -sfs skj_2015_2019_merged_norel.sfs -fstout skj_2015_2019_merged_norel

qsub sfs_skj_2015_2019_merged_norel

# Global Fst between populations 
realSFS fst stats skj_2015_2019_merged_norel.fst.idx
# Output:
FST.Unweight[nObs:21585]:0.020181 Fst.Weight:0.028928

# Writing down 2d-SFS priors - 2015 & 2022
>sfs_skj_2015_2022_merged_norel
nano sfs_skj_2015_2022_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2015_2022_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2015_2022_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS skj_2015_merged_norel.saf.idx skj_2022_merged_norel.saf.idx -fold 1 -P 1 > skj_2015_2022_merged_norel.sfs ; realSFS fst index skj_2015_merged_norel.saf.idx skj_2022_merged_norel.saf.idx -sfs skj_2015_2022_merged_norel.sfs -fstout skj_2015_2022_merged_norel

qsub sfs_skj_2015_2022_merged_norel

# Global Fst between populations 
realSFS fst stats skj_2015_2022_merged_norel.fst.idx
# Output:
FST.Unweight[nObs:26695]:0.052561 Fst.Weight:0.060200

# Writing down 2d-SFS priors - 2016 & 2017
>sfs_skj_2016_2017_merged_norel
nano sfs_skj_2016_2017_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2016_2017_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2016_2017_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS skj_2016_merged_norel.saf.idx skj_2017_merged_norel.saf.idx -fold 1 -P 1 > skj_2016_2017_merged_norel.sfs ; realSFS fst index skj_2016_merged_norel.saf.idx skj_2017_merged_norel.saf.idx -sfs skj_2016_2017_merged_norel.sfs -fstout skj_2016_2017_merged_norel

qsub sfs_skj_2016_2017_merged_norel

# Global Fst between populations 
realSFS fst stats skj_2016_2017_merged_norel.fst.idx
# Output:
FST.Unweight[nObs:26142]:0.007239 Fst.Weight:0.008733

# Writing down 2d-SFS priors - 2016 & 2018
>sfs_skj_2016_2018_merged_norel
nano sfs_skj_2016_2018_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2016_2018_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2016_2018_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS skj_2016_merged_norel.saf.idx skj_2018_merged_norel.saf.idx -fold 1 -P 1 > skj_2016_2018_merged_norel.sfs ; realSFS fst index skj_2016_merged_norel.saf.idx skj_2018_merged_norel.saf.idx -sfs skj_2016_2018_merged_norel.sfs -fstout skj_2016_2018_merged_norel

qsub sfs_skj_2016_2018_merged_norel

# Global Fst between populations 
realSFS fst stats skj_2016_2018_merged_norel.fst.idx
# Output:
FST.Unweight[nObs:26249]:0.003867 Fst.Weight:0.004239

# Writing down 2d-SFS priors - 2016 & 2019
>sfs_skj_2016_2019_merged_norel
nano sfs_skj_2016_2019_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2016_2019_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2016_2019_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS skj_2016_merged_norel.saf.idx skj_2019_merged_norel.saf.idx -fold 1 -P 1 > skj_2016_2019_merged_norel.sfs ; realSFS fst index skj_2016_merged_norel.saf.idx skj_2019_merged_norel.saf.idx -sfs skj_2016_2019_merged_norel.sfs -fstout skj_2016_2019_merged_norel

qsub sfs_skj_2016_2019_merged_norel

# Global Fst between populations 
realSFS fst stats skj_2016_2019_merged_norel.fst.idx
# Output:
FST.Unweight[nObs:21043]:0.025563 Fst.Weight:0.013789

# Writing down 2d-SFS priors - 2016 & 2022
>sfs_skj_2016_2022_merged_norel
nano sfs_skj_2016_2022_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2016_2022_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2016_2022_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS skj_2016_merged_norel.saf.idx skj_2022_merged_norel.saf.idx -fold 1 -P 1 > skj_2016_2022_merged_norel.sfs ; realSFS fst index skj_2016_merged_norel.saf.idx skj_2022_merged_norel.saf.idx -sfs skj_2016_2022_merged_norel.sfs -fstout skj_2016_2022_merged_norel

qsub sfs_skj_2016_2022_merged_norel

# Global Fst between populations 
realSFS fst stats skj_2016_2022_merged_norel.fst.idx
# Output:
FST.Unweight[nObs:25976]:-0.048256 Fst.Weight:0.094264

# Writing down 2d-SFS priors - 2017 & 2018
>sfs_skj_2017_2018_merged_norel
nano sfs_skj_2017_2018_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2017_2018_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2017_2018_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS skj_2017_merged_norel.saf.idx skj_2018_merged_norel.saf.idx -fold 1 -P 1 > skj_2017_2018_merged_norel.sfs ; realSFS fst index skj_2017_merged_norel.saf.idx skj_2018_merged_norel.saf.idx -sfs skj_2017_2018_merged_norel.sfs -fstout skj_2017_2018_merged_norel

qsub sfs_skj_2017_2018_merged_norel

# Global Fst between populations 
realSFS fst stats skj_2017_2018_merged_norel.fst.idx
# Output:
FST.Unweight[nObs:26218]:0.007537 Fst.Weight:0.008077

# Writing down 2d-SFS priors - 2017 & 2019
>sfs_skj_2017_2019_merged_norel
nano sfs_skj_2017_2019_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2017_2019_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2017_2019_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS skj_2017_merged_norel.saf.idx skj_2019_merged_norel.saf.idx -fold 1 -P 1 > skj_2017_2019_merged_norel.sfs ; realSFS fst index skj_2017_merged_norel.saf.idx skj_2019_merged_norel.saf.idx -sfs skj_2017_2019_merged_norel.sfs -fstout skj_2017_2019_merged_norel

qsub sfs_skj_2017_2019_merged_norel

# Global Fst between populations 
realSFS fst stats skj_2017_2019_merged_norel.fst.idx
# Output:
FST.Unweight[nObs:21335]:0.021885 Fst.Weight:0.035226

# Writing down 2d-SFS priors - 2017 & 2022
>sfs_skj_2017_2022_merged_norel
nano sfs_skj_2017_2022_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2017_2022_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2017_2022_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS skj_2017_merged_norel.saf.idx skj_2022_merged_norel.saf.idx -fold 1 -P 1 > skj_2017_2022_merged_norel.sfs ; realSFS fst index skj_2017_merged_norel.saf.idx skj_2022_merged_norel.saf.idx -sfs skj_2017_2022_merged_norel.sfs -fstout skj_2017_2022_merged_norel

qsub sfs_skj_2017_2022_merged_norel

# Global Fst between populations 
realSFS fst stats skj_2017_2022_merged_norel.fst.idx
# Output:
FST.Unweight[nObs:26412]:-0.029599 Fst.Weight:0.077059

# Writing down 2d-SFS priors - 2018 & 2019
>sfs_skj_2018_2019_merged_norel
nano sfs_skj_2018_2019_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2018_2019_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2018_2019_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS skj_2018_merged_norel.saf.idx skj_2019_merged_norel.saf.idx -fold 1 -P 1 > skj_2018_2019_merged_norel.sfs ; realSFS fst index skj_2018_merged_norel.saf.idx skj_2019_merged_norel.saf.idx -sfs skj_2018_2019_merged_norel.sfs -fstout skj_2018_2019_merged_norel

qsub sfs_skj_2018_2019_merged_norel

# Global Fst between populations 
realSFS fst stats skj_2018_2019_merged_norel.fst.idx
# Output:
FST.Unweight[nObs:21057]:0.026764 Fst.Weight:0.017780

# Writing down 2d-SFS priors - 2018 & 2022
>sfs_skj_2018_2022_merged_norel
nano sfs_skj_2018_2022_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2018_2022_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2018_2022_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS skj_2018_merged_norel.saf.idx skj_2022_merged_norel.saf.idx -fold 1 -P 1 > skj_2018_2022_merged_norel.sfs ; realSFS fst index skj_2018_merged_norel.saf.idx skj_2022_merged_norel.saf.idx -sfs skj_2018_2022_merged_norel.sfs -fstout skj_2018_2022_merged_norel

qsub sfs_skj_2018_2022_merged_norel

# Global Fst between populations 
realSFS fst stats skj_2018_2022_merged_norel.fst.idx
# Output:
FST.Unweight[nObs:26068]:-0.049228 Fst.Weight:0.098434

# Writing down 2d-SFS priors - 2019 & 2022
>sfs_skj_2019_2022_merged_norel
nano sfs_skj_2019_2022_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2019_2022_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2019_2022_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS skj_2019_merged_norel.saf.idx skj_2022_merged_norel.saf.idx -fold 1 -P 1 > skj_2019_2022_merged_norel.sfs ; realSFS fst index skj_2019_merged_norel.saf.idx skj_2022_merged_norel.saf.idx -sfs skj_2019_2022_merged_norel.sfs -fstout skj_2019_2022_merged_norel

qsub sfs_skj_2019_2022_merged_norel

# Global Fst between populations 
realSFS fst stats skj_2019_2022_merged_norel.fst.idx
# Output:
FST.Unweight[nObs:22051]:-0.072192 Fst.Weight:0.045822

#----------------------------------------------------------------------------------------
#Yellowfin
cd /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis
#FLAG: also try subpop vs main

#Doing SAF/SFS work
#2015 (8 individuals, so 80% is 6)
>sfs_yft_2015_merged_neutral
nano sfs_yft_2015_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2015_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2015_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_yft_unlinked -b bams_yft_merged_2015 -GL 1 -P 1 -minInd 6 $TODO -out yft_2015_merged_neutral

qsub sfs_yft_2015_merged_neutral

#2016 (4 individuals, so 80% is 3)
>sfs_yft_2016_merged_neutral
nano sfs_yft_2016_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2016_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2016_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_yft_unlinked -b bams_yft_merged_2016 -GL 1 -P 1 -minInd 3 $TODO -out yft_2016_merged_neutral

qsub sfs_yft_2016_merged_neutral

#2017 (5 individuals, so 80% is 4)
>sfs_yft_2017_merged_neutral
nano sfs_yft_2017_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2017_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2017_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_yft_unlinked -b bams_yft_merged_2017 -GL 1 -P 1 -minInd 4 $TODO -out yft_2017_merged_neutral

qsub sfs_yft_2017_merged_neutral

#2018 (43 individuals, so 80% is 34)
>sfs_yft_2018_merged_neutral
nano sfs_yft_2018_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2018_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2018_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_yft_unlinked -b bams_yft_merged_2018 -GL 1 -P 1 -minInd 34 $TODO -out yft_2018_merged_neutral

qsub sfs_yft_2018_merged_neutral

#2019 (1 individual, so no minind)
>sfs_yft_2019_merged_neutral
nano sfs_yft_2019_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2019_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2019_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_yft_unlinked -b bams_yft_2019 -GL 1 -P 1 -minInd 0 $TODO -out yft_2019_merged_neutral

qsub sfs_yft_2019_merged_neutral

#2020 (29 individuals, so 80% is 23)
>sfs_yft_2020_merged_neutral
nano sfs_yft_2020_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2020_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2020_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_yft_unlinked -b bams_yft_merged_2020 -GL 1 -P 1 -minInd 23 $TODO -out yft_2020_merged_neutral

qsub sfs_yft_2020_merged_neutral

#2022 (6 individuals, so 80% is 5)
>sfs_yft_2022_merged_neutral
nano sfs_yft_2022_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2022_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2022_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_yft_unlinked -b bams_yft_merged_2022 -GL 1 -P 1 -minInd 5 $TODO -out yft_2022_merged_neutral

qsub sfs_yft_2022_merged_neutral

# Writing down 2d-SFS priors - 2015 & 2016
>sfs_yft_2015_2016_merged_neutral
nano sfs_yft_2015_2016_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2015_2016_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2015_2016_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2015_merged_neutral.saf.idx yft_2016_merged_neutral.saf.idx -fold 1 -P 1 > yft_2015_2016_merged_neutral.sfs ; realSFS fst index yft_2015_merged_neutral.saf.idx yft_2016_merged_neutral.saf.idx -sfs yft_2015_2016_merged_neutral.sfs -fstout yft_2015_2016_merged_neutral

qsub sfs_yft_2015_2016_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2015_2016_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:101429]:0.018521 Fst.Weight:0.052742

# The difference between unweighted and weighted values is averaging versus ratio of sums method. 
# The weighted value is derived from the ratio of the separate per-locus sums of numerator and denominator values, 
# while the unweighted value is the average of per-locus values. [If that is unclear: weighted is sum(a)/sum(a+b), while unweighted is average(a/(a+b))].

# Writing down 2d-SFS priors - 2015 & 2017
>sfs_yft_2015_2017_merged_neutral
nano sfs_yft_2015_2017_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2015_2017_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2015_2017_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2015_merged_neutral.saf.idx yft_2017_merged_neutral.saf.idx -fold 1 -P 1 > yft_2015_2017_merged_neutral.sfs ; realSFS fst index yft_2015_merged_neutral.saf.idx yft_2017_merged_neutral.saf.idx -sfs yft_2015_2017_merged_neutral.sfs -fstout yft_2015_2017_merged_neutral

qsub sfs_yft_2015_2017_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2015_2017_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:97821]:0.055376 Fst.Weight:0.047808

# Writing down 2d-SFS priors - 2015 & 2018
>sfs_yft_2015_2018_merged_neutral
nano sfs_yft_2015_2018_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2015_2018_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2015_2018_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2015_merged_neutral.saf.idx yft_2018_merged_neutral.saf.idx -fold 1 -P 1 > yft_2015_2018_merged_neutral.sfs ; realSFS fst index yft_2015_merged_neutral.saf.idx yft_2018_merged_neutral.saf.idx -sfs yft_2015_2018_merged_neutral.sfs -fstout yft_2015_2018_merged_neutral

qsub sfs_yft_2015_2018_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2015_2018_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:96827]:0.031785 Fst.Weight:0.069079

# Writing down 2d-SFS priors - 2015 & 2019
>sfs_yft_2015_2019_merged_neutral
nano sfs_yft_2015_2019_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2015_2019_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2015_2019_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2015_merged_neutral.saf.idx yft_2019_merged_neutral.saf.idx -fold 1 -P 1 > yft_2015_2019_merged_neutral.sfs ; realSFS fst index yft_2015_merged_neutral.saf.idx yft_2019_merged_neutral.saf.idx -sfs yft_2015_2019_merged_neutral.sfs -fstout yft_2015_2019_merged_neutral

qsub sfs_yft_2015_2019_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2015_2019_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:97392]:0.081482 Fst.Weight:0.138098

# Writing down 2d-SFS priors - 2015 & 2020
>sfs_yft_2015_2020_merged_neutral
nano sfs_yft_2015_2020_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2015_2020_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2015_2020_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2015_merged_neutral.saf.idx yft_2020_merged_neutral.saf.idx -fold 1 -P 1 > yft_2015_2020_merged_neutral.sfs ; realSFS fst index yft_2015_merged_neutral.saf.idx yft_2020_merged_neutral.saf.idx -sfs yft_2015_2020_merged_neutral.sfs -fstout yft_2015_2020_merged_neutral

qsub sfs_yft_2015_2020_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2015_2020_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:100334]:0.035073 Fst.Weight:0.025052

# Writing down 2d-SFS priors - 2015 & 2022
>sfs_yft_2015_2022_merged_neutral
nano sfs_yft_2015_2022_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2015_2022_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2015_2022_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2015_merged_neutral.saf.idx yft_2022_merged_neutral.saf.idx -fold 1 -P 1 > yft_2015_2022_merged_neutral.sfs ; realSFS fst index yft_2015_merged_neutral.saf.idx yft_2022_merged_neutral.saf.idx -sfs yft_2015_2022_merged_neutral.sfs -fstout yft_2015_2022_merged_neutral

qsub sfs_yft_2015_2022_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2015_2022_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:101305]:0.024941 Fst.Weight:0.051843

# Writing down 2d-SFS priors - 2016 & 2017
>sfs_yft_2016_2017_merged_neutral
nano sfs_yft_2016_2017_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2016_2017_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2016_2017_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2016_merged_neutral.saf.idx yft_2017_merged_neutral.saf.idx -fold 1 -P 1 > yft_2016_2017_merged_neutral.sfs ; realSFS fst index yft_2016_merged_neutral.saf.idx yft_2017_merged_neutral.saf.idx -sfs yft_2016_2017_merged_neutral.sfs -fstout yft_2016_2017_merged_neutral

qsub sfs_yft_2016_2017_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2016_2017_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:99978]:0.044605 Fst.Weight:0.070586

# Writing down 2d-SFS priors - 2016 & 2018
>sfs_yft_2016_2018_merged_neutral
nano sfs_yft_2016_2018_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2016_2018_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2016_2018_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2016_merged_neutral.saf.idx yft_2018_merged_neutral.saf.idx -fold 1 -P 1 > yft_2016_2018_merged_neutral.sfs ; realSFS fst index yft_2016_merged_neutral.saf.idx yft_2018_merged_neutral.saf.idx -sfs yft_2016_2018_merged_neutral.sfs -fstout yft_2016_2018_merged_neutral

qsub sfs_yft_2016_2018_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2016_2018_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:99747]:0.003960 Fst.Weight:0.073581

# Writing down 2d-SFS priors - 2016 & 2019
>sfs_yft_2016_2019_merged_neutral
nano sfs_yft_2016_2019_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2016_2019_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2016_2019_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2016_merged_neutral.saf.idx yft_2019_merged_neutral.saf.idx -fold 1 -P 1 > yft_2016_2019_merged_neutral.sfs ; realSFS fst index yft_2016_merged_neutral.saf.idx yft_2019_merged_neutral.saf.idx -sfs yft_2016_2019_merged_neutral.sfs -fstout yft_2016_2019_merged_neutral

qsub sfs_yft_2016_2019_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2016_2019_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:99734]:0.383920 Fst.Weight:0.190273

# Writing down 2d-SFS priors - 2016 & 2020
>sfs_yft_2016_2020_merged_neutral
nano sfs_yft_2016_2020_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2016_2020_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2016_2020_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2016_merged_neutral.saf.idx yft_2020_merged_neutral.saf.idx -fold 1 -P 1 > yft_2016_2020_merged_neutral.sfs ; realSFS fst index yft_2016_merged_neutral.saf.idx yft_2020_merged_neutral.saf.idx -sfs yft_2016_2020_merged_neutral.sfs -fstout yft_2016_2020_merged_neutral

qsub sfs_yft_2016_2020_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2016_2020_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:103240]:0.007296 Fst.Weight:0.040054

# Writing down 2d-SFS priors - 2016 & 2022
>sfs_yft_2016_2022_merged_neutral
nano sfs_yft_2016_2022_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2016_2022_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2016_2022_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2016_merged_neutral.saf.idx yft_2022_merged_neutral.saf.idx -fold 1 -P 1 > yft_2016_2022_merged_neutral.sfs ; realSFS fst index yft_2016_merged_neutral.saf.idx yft_2022_merged_neutral.saf.idx -sfs yft_2016_2022_merged_neutral.sfs -fstout yft_2016_2022_merged_neutral

qsub sfs_yft_2016_2022_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2016_2022_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:104422]:0.087208 Fst.Weight:0.055431

# Writing down 2d-SFS priors - 2017 & 2018
>sfs_yft_2017_2018_merged_neutral
nano sfs_yft_2017_2018_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2017_2018_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2017_2018_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2017_merged_neutral.saf.idx yft_2018_merged_neutral.saf.idx -fold 1 -P 1 > yft_2017_2018_merged_neutral.sfs ; realSFS fst index yft_2017_merged_neutral.saf.idx yft_2018_merged_neutral.saf.idx -sfs yft_2017_2018_merged_neutral.sfs -fstout yft_2017_2018_merged_neutral

qsub sfs_yft_2017_2018_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2017_2018_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:95663]:0.064958 Fst.Weight:0.043955

# Writing down 2d-SFS priors - 2017 & 2019
>sfs_yft_2017_2019_merged_neutral
nano sfs_yft_2017_2019_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2017_2019_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2017_2019_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2017_merged_neutral.saf.idx yft_2019_merged_neutral.saf.idx -fold 1 -P 1 > yft_2017_2019_merged_neutral.sfs ; realSFS fst index yft_2017_merged_neutral.saf.idx yft_2019_merged_neutral.saf.idx -sfs yft_2017_2019_merged_neutral.sfs -fstout yft_2017_2019_merged_neutral

qsub sfs_yft_2017_2019_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2017_2019_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:96322]:0.038739 Fst.Weight:0.122689

# Writing down 2d-SFS priors - 2017 & 2020
>sfs_yft_2017_2020_merged_neutral
nano sfs_yft_2017_2020_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2017_2020_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2017_2020_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2017_merged_neutral.saf.idx yft_2020_merged_neutral.saf.idx -fold 1 -P 1 > yft_2017_2020_merged_neutral.sfs ; realSFS fst index yft_2017_merged_neutral.saf.idx yft_2020_merged_neutral.saf.idx -sfs yft_2017_2020_merged_neutral.sfs -fstout yft_2017_2020_merged_neutral

qsub sfs_yft_2017_2020_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2017_2020_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:98791]:0.104127 Fst.Weight:0.044772

# Writing down 2d-SFS priors - 2017 & 2022
>sfs_yft_2017_2022_merged_neutral
nano sfs_yft_2017_2022_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2017_2022_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2017_2022_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2017_merged_neutral.saf.idx yft_2022_merged_neutral.saf.idx -fold 1 -P 1 > yft_2017_2022_merged_neutral.sfs ; realSFS fst index yft_2017_merged_neutral.saf.idx yft_2022_merged_neutral.saf.idx -sfs yft_2017_2022_merged_neutral.sfs -fstout yft_2017_2022_merged_neutral

qsub sfs_yft_2017_2022_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2017_2022_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:99902]:0.066053 Fst.Weight:0.073657

# Writing down 2d-SFS priors - 2018 & 2019
>sfs_yft_2018_2019_merged_neutral
nano sfs_yft_2018_2019_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2018_2019_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2018_2019_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2018_merged_neutral.saf.idx yft_2019_merged_neutral.saf.idx -fold 1 -P 1 > yft_2018_2019_merged_neutral.sfs ; realSFS fst index yft_2018_merged_neutral.saf.idx yft_2019_merged_neutral.saf.idx -sfs yft_2018_2019_merged_neutral.sfs -fstout yft_2018_2019_merged_neutral

qsub sfs_yft_2018_2019_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2018_2019_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:95142]:0.038520 Fst.Weight:0.130300

# Writing down 2d-SFS priors - 2018 & 2020
>sfs_yft_2018_2020_merged_neutral
nano sfs_yft_2018_2020_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2018_2020_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2018_2020_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2018_merged_neutral.saf.idx yft_2020_merged_neutral.saf.idx -fold 1 -P 1 > yft_2018_2020_merged_neutral.sfs ; realSFS fst index yft_2018_merged_neutral.saf.idx yft_2020_merged_neutral.saf.idx -sfs yft_2018_2020_merged_neutral.sfs -fstout yft_2018_2020_merged_neutral

qsub sfs_yft_2018_2020_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2018_2020_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:98969]:0.014479 Fst.Weight:0.066756

# Writing down 2d-SFS priors - 2018 & 2022
>sfs_yft_2018_2022_merged_neutral
nano sfs_yft_2018_2022_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2018_2022_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2018_2022_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2018_merged_neutral.saf.idx yft_2022_merged_neutral.saf.idx -fold 1 -P 1 > yft_2018_2022_merged_neutral.sfs ; realSFS fst index yft_2018_merged_neutral.saf.idx yft_2022_merged_neutral.saf.idx -sfs yft_2018_2022_merged_neutral.sfs -fstout yft_2018_2022_merged_neutral

qsub sfs_yft_2018_2022_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2018_2022_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:99755]:0.000415 Fst.Weight:0.072277

# Writing down 2d-SFS priors - 2019 & 2020
>sfs_yft_2019_2020_merged_neutral
nano sfs_yft_2019_2020_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2019_2020_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2019_2020_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2019_merged_neutral.saf.idx yft_2020_merged_neutral.saf.idx -fold 1 -P 1 > yft_2019_2020_merged_neutral.sfs ; realSFS fst index yft_2019_merged_neutral.saf.idx yft_2020_merged_neutral.saf.idx -sfs yft_2019_2020_merged_neutral.sfs -fstout yft_2019_2020_merged_neutral

qsub sfs_yft_2019_2020_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2019_2020_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:98493]:0.109358 Fst.Weight:0.130806

# Writing down 2d-SFS priors - 2019 & 2022
>sfs_yft_2019_2022_merged_neutral
nano sfs_yft_2019_2022_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2019_2022_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2019_2022_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2019_merged_neutral.saf.idx yft_2022_merged_neutral.saf.idx -fold 1 -P 1 > yft_2019_2022_merged_neutral.sfs ; realSFS fst index yft_2019_merged_neutral.saf.idx yft_2022_merged_neutral.saf.idx -sfs yft_2019_2022_merged_neutral.sfs -fstout yft_2019_2022_merged_neutral

qsub sfs_yft_2019_2022_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2019_2022_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:99574]:0.606946 Fst.Weight:0.189460

# Writing down 2d-SFS priors - 2020 & 2022
>sfs_yft_2020_2022_merged_neutral
nano sfs_yft_2020_2022_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2020_2022_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2020_2022_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_2020_merged_neutral.saf.idx yft_2022_merged_neutral.saf.idx -fold 1 -P 1 > yft_2020_2022_merged_neutral.sfs ; realSFS fst index yft_2020_merged_neutral.saf.idx yft_2022_merged_neutral.saf.idx -sfs yft_2020_2022_merged_neutral.sfs -fstout yft_2020_2022_merged_neutral

qsub sfs_yft_2020_2022_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_2020_2022_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:103246]:-0.001726 Fst.Weight:0.034666

#------------------------ Bigeye reference-based analysis of genetic divergence (FST) across years
cd /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis

#Doing SAF/SFS work
#2016 (11 individuals, so 80% is 9)
>sfs_bet_2016_merged_neutral
nano sfs_bet_2016_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2016_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2016_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_bet_unlinked -b bams_bet_merged_2016 -GL 1 -P 1 -minInd 9 $TODO -out bet_2016_merged_neutral

qsub sfs_bet_2016_merged_neutral

#2017 (2 individuals, so set minind to 1)
>sfs_bet_2017_merged_neutral
nano sfs_bet_2017_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2017_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2017_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_bet_unlinked -b bams_bet_merged_2017 -GL 1 -P 1 -minInd 1 $TODO -out bet_2017_merged_neutral

qsub sfs_bet_2017_merged_neutral

#2018 (11 individuals, so set minind to 9)
>sfs_bet_2018_merged_neutral
nano sfs_bet_2018_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2018_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2018_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_bet_unlinked -b bams_bet_merged_2018 -GL 1 -P 1 -minInd 9 $TODO -out bet_2018_merged_neutral

qsub sfs_bet_2018_merged_neutral

#2019 (1 individual, so set minind to 0)
>sfs_bet_2019_merged_neutral
nano sfs_bet_2019_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2019_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2019_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_bet_unlinked -b bams_bet_merged_2019 -GL 1 -P 1 -minInd 0 $TODO -out bet_2019_merged_neutral

qsub sfs_bet_2019_merged_neutral

#2020 (24 individuals, so set minind to 19)
>sfs_bet_2020_merged_neutral
nano sfs_bet_2020_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2020_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2020_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_bet_unlinked -b bams_bet_merged_2020 -GL 1 -P 1 -minInd 19 $TODO -out bet_2020_merged_neutral

qsub sfs_bet_2020_merged_neutral

#2022 (8 individuals, so set minind to 6)
>sfs_bet_2022_merged_neutral
nano sfs_bet_2022_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2022_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2022_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_bet_unlinked -b bams_bet_merged_2022 -GL 1 -P 1 -minInd 6 $TODO -out bet_2022_merged_neutral

qsub sfs_bet_2022_merged_neutral

# Writing down 2d-SFS priors - 2016 & 2017
>sfs_bet_2016_2017_merged_neutral
nano sfs_bet_2016_2017_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2016_2017_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2016_2017_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS bet_2016_merged_neutral.saf.idx bet_2017_merged_neutral.saf.idx -fold 1 -P 1 > bet_2016_2017_merged_neutral.sfs ; realSFS fst index bet_2016_merged_neutral.saf.idx bet_2017_merged_neutral.saf.idx -sfs bet_2016_2017_merged_neutral.sfs -fstout bet_2016_2017_merged_neutral

qsub sfs_bet_2016_2017_merged_neutral

# Global Fst between populations 
realSFS fst stats bet_2016_2017_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:98797]:-0.010424 Fst.Weight:0.061151

# The difference between unweighted and weighted values is averaging versus ratio of sums method. 
# The weighted value is derived from the ratio of the separate per-locus sums of numerator and denominator values, 
# while the unweighted value is the average of per-locus values. [If that is unclear: weighted is sum(a)/sum(a+b), while unweighted is average(a/(a+b))].

# Writing down 2d-SFS priors - 2016 & 2018
>sfs_bet_2016_2018_merged_neutral
nano sfs_bet_2016_2018_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2016_2018_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2016_2018_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS bet_2016_merged_neutral.saf.idx bet_2018_merged_neutral.saf.idx -fold 1 -P 1 > bet_2016_2018_merged_neutral.sfs ; realSFS fst index bet_2016_merged_neutral.saf.idx bet_2018_merged_neutral.saf.idx -sfs bet_2016_2018_merged_neutral.sfs -fstout bet_2016_2018_merged_neutral

qsub sfs_bet_2016_2018_merged_neutral

# Global Fst between populations 
realSFS fst stats bet_2016_2018_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:98057]:0.023600 Fst.Weight:0.022591

# Writing down 2d-SFS priors - 2016 & 2019
>sfs_bet_2016_2019_merged_neutral
nano sfs_bet_2016_2019_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2016_2019_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2016_2019_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS bet_2016_merged_neutral.saf.idx bet_2019_merged_neutral.saf.idx -fold 1 -P 1 > bet_2016_2019_merged_neutral.sfs ; realSFS fst index bet_2016_merged_neutral.saf.idx bet_2019_merged_neutral.saf.idx -sfs bet_2016_2019_merged_neutral.sfs -fstout bet_2016_2019_merged_neutral

qsub sfs_bet_2016_2019_merged_neutral

# Global Fst between populations 
realSFS fst stats bet_2016_2019_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:93684]:0.126123 Fst.Weight:0.157868

# Writing down 2d-SFS priors - 2016 & 2020
>sfs_bet_2016_2020_merged_neutral
nano sfs_bet_2016_2020_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2016_2020_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2016_2020_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS bet_2016_merged_neutral.saf.idx bet_2020_merged_neutral.saf.idx -fold 1 -P 1 > bet_2016_2020_merged_neutral.sfs ; realSFS fst index bet_2016_merged_neutral.saf.idx bet_2020_merged_neutral.saf.idx -sfs bet_2016_2020_merged_neutral.sfs -fstout bet_2016_2020_merged_neutral

qsub sfs_bet_2016_2020_merged_neutral

# Global Fst between populations 
realSFS fst stats bet_2016_2020_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:96858]:0.021080 Fst.Weight:0.019450

# Writing down 2d-SFS priors - 2016 & 2022
>sfs_bet_2016_2022_merged_neutral
nano sfs_bet_2016_2022_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2016_2022_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2016_2022_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS bet_2016_merged_neutral.saf.idx bet_2022_merged_neutral.saf.idx -fold 1 -P 1 > bet_2016_2022_merged_neutral.sfs ; realSFS fst index bet_2016_merged_neutral.saf.idx bet_2022_merged_neutral.saf.idx -sfs bet_2016_2022_merged_neutral.sfs -fstout bet_2016_2022_merged_neutral

qsub sfs_bet_2016_2022_merged_neutral

# Global Fst between populations 
realSFS fst stats bet_2016_2022_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:98192]:0.017628 Fst.Weight:0.038353

# Writing down 2d-SFS priors - 2017 & 2018
>sfs_bet_2017_2018_merged_neutral
nano sfs_bet_2017_2018_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2017_2018_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2017_2018_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS bet_2017_merged_neutral.saf.idx bet_2018_merged_neutral.saf.idx -fold 1 -P 1 > bet_2017_2018_merged_neutral.sfs ; realSFS fst index bet_2017_merged_neutral.saf.idx bet_2018_merged_neutral.saf.idx -sfs bet_2017_2018_merged_neutral.sfs -fstout bet_2017_2018_merged_neutral

qsub sfs_bet_2017_2018_merged_neutral

# Global Fst between populations 
realSFS fst stats bet_2017_2018_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:104686]:0.045416 Fst.Weight:0.078828

# Writing down 2d-SFS priors - 2017 & 2019
>sfs_bet_2017_2019_merged_neutral
nano sfs_bet_2017_2019_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2017_2019_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2017_2019_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS bet_2017_merged_neutral.saf.idx bet_2019_merged_neutral.saf.idx -fold 1 -P 1 > bet_2017_2019_merged_neutral.sfs ; realSFS fst index bet_2017_merged_neutral.saf.idx bet_2019_merged_neutral.saf.idx -sfs bet_2017_2019_merged_neutral.sfs -fstout bet_2017_2019_merged_neutral

qsub sfs_bet_2017_2019_merged_neutral

# Global Fst between populations 
realSFS fst stats bet_2017_2019_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:98039]:0.386715 Fst.Weight:0.211814

# Writing down 2d-SFS priors - 2017 & 2020
>sfs_bet_2017_2020_merged_neutral
nano sfs_bet_2017_2020_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2017_2020_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2017_2020_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS bet_2017_merged_neutral.saf.idx bet_2020_merged_neutral.saf.idx -fold 1 -P 1 > bet_2017_2020_merged_neutral.sfs ; realSFS fst index bet_2017_merged_neutral.saf.idx bet_2020_merged_neutral.saf.idx -sfs bet_2017_2020_merged_neutral.sfs -fstout bet_2017_2020_merged_neutral

qsub sfs_bet_2017_2020_merged_neutral

# Global Fst between populations 
realSFS fst stats bet_2017_2020_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:102937]:-0.005643 Fst.Weight:0.073877

# Writing down 2d-SFS priors - 2017 & 2022
>sfs_bet_2017_2022_merged_neutral
nano sfs_bet_2017_2022_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2017_2022_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2017_2022_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS bet_2017_merged_neutral.saf.idx bet_2022_merged_neutral.saf.idx -fold 1 -P 1 > bet_2017_2022_merged_neutral.sfs ; realSFS fst index bet_2017_merged_neutral.saf.idx bet_2022_merged_neutral.saf.idx -sfs bet_2017_2022_merged_neutral.sfs -fstout bet_2017_2022_merged_neutral

qsub sfs_bet_2017_2022_merged_neutral

# Global Fst between populations 
realSFS fst stats bet_2017_2022_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:104976]:0.300577 Fst.Weight:0.105598

# Writing down 2d-SFS priors - 2018 & 2019
>sfs_bet_2018_2019_merged_neutral
nano sfs_bet_2018_2019_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2018_2019_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2018_2019_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS bet_2018_merged_neutral.saf.idx bet_2019_merged_neutral.saf.idx -fold 1 -P 1 > bet_2018_2019_merged_neutral.sfs ; realSFS fst index bet_2018_merged_neutral.saf.idx bet_2019_merged_neutral.saf.idx -sfs bet_2018_2019_merged_neutral.sfs -fstout bet_2018_2019_merged_neutral

qsub sfs_bet_2018_2019_merged_neutral

# Global Fst between populations 
realSFS fst stats bet_2018_2019_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:96928]:0.362844 Fst.Weight:0.196647

# Writing down 2d-SFS priors - 2018 & 2020
>sfs_bet_2018_2020_merged_neutral
nano sfs_bet_2018_2020_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2018_2020_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2018_2020_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS bet_2018_merged_neutral.saf.idx bet_2020_merged_neutral.saf.idx -fold 1 -P 1 > bet_2018_2020_merged_neutral.sfs ; realSFS fst index bet_2018_merged_neutral.saf.idx bet_2020_merged_neutral.saf.idx -sfs bet_2018_2020_merged_neutral.sfs -fstout bet_2018_2020_merged_neutral

qsub sfs_bet_2018_2020_merged_neutral

# Global Fst between populations 
realSFS fst stats bet_2018_2020_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:102463]:0.011152 Fst.Weight:0.021051

# Writing down 2d-SFS priors - 2018 & 2022
>sfs_bet_2018_2022_merged_neutral
nano sfs_bet_2018_2022_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2018_2022_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2018_2022_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS bet_2018_merged_neutral.saf.idx bet_2022_merged_neutral.saf.idx -fold 1 -P 1 > bet_2018_2022_merged_neutral.sfs ; realSFS fst index bet_2018_merged_neutral.saf.idx bet_2022_merged_neutral.saf.idx -sfs bet_2018_2022_merged_neutral.sfs -fstout bet_2018_2022_merged_neutral

qsub sfs_bet_2018_2022_merged_neutral

# Global Fst between populations 
realSFS fst stats bet_2018_2022_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:104295]:0.017616 Fst.Weight:0.030701

# Writing down 2d-SFS priors - 2019 & 2020
>sfs_bet_2019_2020_merged_neutral
nano sfs_bet_2019_2020_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2019_2020_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2019_2020_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS bet_2019_merged_neutral.saf.idx bet_2020_merged_neutral.saf.idx -fold 1 -P 1 > bet_2019_2020_merged_neutral.sfs ; realSFS fst index bet_2019_merged_neutral.saf.idx bet_2020_merged_neutral.saf.idx -sfs bet_2019_2020_merged_neutral.sfs -fstout bet_2019_2020_merged_neutral

qsub sfs_bet_2019_2020_merged_neutral

# Global Fst between populations 
realSFS fst stats bet_2019_2020_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:95607]:0.149189 Fst.Weight:0.180572

# Writing down 2d-SFS priors - 2019 & 2022
>sfs_bet_2019_2022_merged_neutral
nano sfs_bet_2019_2022_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2019_2022_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2019_2022_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS bet_2019_merged_neutral.saf.idx bet_2022_merged_neutral.saf.idx -fold 1 -P 1 > bet_2019_2022_merged_neutral.sfs ; realSFS fst index bet_2019_merged_neutral.saf.idx bet_2022_merged_neutral.saf.idx -sfs bet_2019_2022_merged_neutral.sfs -fstout bet_2019_2022_merged_neutral

qsub sfs_bet_2019_2022_merged_neutral

# Global Fst between populations 
realSFS fst stats bet_2019_2022_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:97052]:0.652963 Fst.Weight:0.224265

# Writing down 2d-SFS priors - 2020 & 2022
>sfs_bet_2020_2022_merged_neutral
nano sfs_bet_2020_2022_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2020_2022_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2020_2022_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS bet_2020_merged_neutral.saf.idx bet_2022_merged_neutral.saf.idx -fold 1 -P 1 > bet_2020_2022_merged_neutral.sfs ; realSFS fst index bet_2020_merged_neutral.saf.idx bet_2022_merged_neutral.saf.idx -sfs bet_2020_2022_merged_neutral.sfs -fstout bet_2020_2022_merged_neutral

qsub sfs_bet_2020_2022_merged_neutral

# Global Fst between populations 
realSFS fst stats bet_2020_2022_merged_neutral.fst.idx
# Output:
FST.Unweight[nObs:102790]:0.002733 Fst.Weight:0.030966

#------------------------------ PW Fst between adults and larvae ---------------------------
#Yellowfin
#Larvae (67 so 80% is 54)
>sfs_yft_larvae_merged_neutral
nano sfs_yft_larvae_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_larvae_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_larvae_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_yft_unlinked -b bams_yft_merged_larvae -GL 1 -P 1 -minInd 54 $TODO -out yft_larvae_merged_neutral

qsub sfs_yft_larvae_merged_neutral

#For adults, you can just use yft_2020_merged_neutral.saf.idx because those are the only adults and you produced this already

#Fst adults vs larvae
>sfs_yft_adults_v_larvae_merged_neutral
nano sfs_yft_adults_v_larvae_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_adults_v_larvae_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_adults_v_larvae_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_larvae_merged_neutral.saf.idx yft_2020_merged_neutral.saf.idx -fold 1 -P 1 > yft_adults_v_larvae_merged_neutral.sfs ; realSFS fst index yft_larvae_merged_neutral.saf.idx yft_2020_merged_neutral.saf.idx -sfs yft_adults_v_larvae_merged_neutral.sfs -fstout yft_adults_v_larvae_merged_neutral

qsub sfs_yft_adults_v_larvae_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_adults_v_larvae_merged_neutral.fst.idx

output:
FST.Unweight[nObs:99418]:0.009726 Fst.Weight:0.033441

#Bigeye
#Larvae (33 so 80% is 26)
>sfs_bet_larvae_merged_neutral
nano sfs_bet_larvae_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_larvae_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_larvae_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_bet_unlinked -b bams_bet_merged_larvae -GL 1 -P 1 -minInd 26 $TODO -out bet_larvae_merged_neutral

qsub sfs_bet_larvae_merged_neutral

#For adults, you can just use bet_2020_merged_neutral.saf.idx because those are the only adults and you produced this already

#Fst adults vs larvae
>sfs_bet_adults_v_larvae_merged_neutral
nano sfs_bet_adults_v_larvae_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_adults_v_larvae_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_adults_v_larvae_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS bet_larvae_merged_neutral.saf.idx bet_2020_merged_neutral.saf.idx -fold 1 -P 1 > bet_adults_v_larvae_merged_neutral.sfs ; realSFS fst index bet_larvae_merged_neutral.saf.idx bet_2020_merged_neutral.saf.idx -sfs bet_adults_v_larvae_merged_neutral.sfs -fstout bet_adults_v_larvae_merged_neutral

qsub sfs_bet_adults_v_larvae_merged_neutral

# Global Fst between populations 
realSFS fst stats bet_adults_v_larvae_merged_neutral.fst.idx

output:
FST.Unweight[nObs:102340]:0.010534 Fst.Weight:0.012001

#-------------------------- PW Fst between main yft pop and subpop -----------------------

#Main pop (79 so 80% is 63)
>sfs_yft_main_merged_neutral
nano sfs_yft_main_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_main_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_main_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_yft_unlinked -b bams_yft_merged_main -GL 1 -P 1 -minInd 63 $TODO -out yft_main_merged_neutral

qsub sfs_yft_main_merged_neutral

#Subpop (17 so 80% is 14)
>sfs_yft_subpop_merged_neutral
nano sfs_yft_subpop_merged_neutral
#!/bin/bash -l
#$ -V # inherit the subpopmission environment
#$ -cwd # start job in subpopmission directory
#$ -N sfs_yft_subpop_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_subpop_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_yft_unlinked -b bams_yft_merged_subpop -GL 1 -P 1 -minInd 14 $TODO -out yft_subpop_merged_neutral

qsub sfs_yft_subpop_merged_neutral

#Fst main pop vs subpop
>sfs_yft_main_v_subpop_merged_neutral
nano sfs_yft_main_v_subpop_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_main_v_subpop_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_main_v_subpop_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
realSFS yft_main_merged_neutral.saf.idx yft_subpop_merged_neutral.saf.idx -fold 1 -P 1 > yft_main_v_subpop_merged_neutral.sfs ; realSFS fst index yft_main_merged_neutral.saf.idx yft_subpop_merged_neutral.saf.idx -sfs yft_main_v_subpop_merged_neutral.sfs -fstout yft_main_v_subpop_merged_neutral

qsub sfs_yft_main_v_subpop_merged_neutral

# Global Fst between populations 
realSFS fst stats yft_main_v_subpop_merged_neutral.fst.idx

output:
FST.Unweight[nObs:91984]:0.057886 Fst.Weight:0.390871

#------------------------------------------ Tajima's D ------------------------------------------
##ANGSD EXAMPLE
#from https://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests
module load angsd
module list
#angsd/0.935

#Step 1: Finding a 'global estimate' of the SFS
#First estimate the site allele frequency likelihood
./angsd -bam bam.filelist -doSaf 1 -anc chimpHg19.fa -GL 1 -P 24 -out out
#Obtain the maximum likelihood estimate of the SFS using the realSFS program found in the misc subfolder
./misc/realSFS out.saf.idx -P 24 > out.sfs
#Step 2: Calculate the thetas for each site
realSFS saf2theta out.saf.idx -sfs out.sfs -outname out
#The output from the above command are two files out.thetas.gz and out.thetas.idx. A formal description of these files can be found in the doc/formats.pdf in the angsd package. It is possible to extract the logscale persite thetas using the ./thetaStat print program.
thetaStat print out.thetas.idx 2>/dev/null |head
#Step 3a: Estimate Tajimas D and other statistics
#calculate Tajimas D
./misc/thetaStat do_stat out.thetas.idx

#After filtering for LD
#Skipjack (195 individuals, so 80% is 156)
>sfs_skj_dem_merged_norel
nano sfs_skj_dem_merged_norel
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_dem_merged_norel # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_dem_merged_norel.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_skj_unlinked -b bams_skj_merged_norel -GL 1 -P 1 -minInd 156 $TODO -out skj_dem_merged_norel

qsub sfs_skj_dem_merged_norel

#Yellowfin (neutral)
>sfs_yft_dem_merged_neutral
nano sfs_yft_dem_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_dem_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_dem_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_yft_unlinked -b bams_yft_merged -GL 1 -P 1 -minInd 77 $TODO -out yft_dem_merged_neutral

qsub sfs_yft_dem_merged_neutral

#Bigeye
>sfs_bet_dem_merged_neutral
nano sfs_bet_dem_merged_neutral
#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_dem_merged_neutral # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_dem_merged_neutral.qlog # Name the file where to redirect standard output and error
module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites neutral_sites_bet_unlinked -b bams_bet_merged -GL 1 -P 1 -minInd 46 $TODO -out bet_dem_merged_neutral

qsub sfs_bet_dem_merged_neutral

##Obtain the maximum likelihood estimate of the SFS using the realSFS program found in the misc subfolder
#Ran this code to make sfs files already, no need to rerun 
realSFS skj_dem_merged_norel.saf.idx -fold 1 >skj_dem_merged_norel_folded.sfs #28449
realSFS yft_dem_merged_neutral.saf.idx -fold 1 >yft_dem_merged_folded_neutral.sfs #
realSFS bet_dem_merged_neutral.saf.idx -fold 1 >bet_dem_merged_folded_neutral.sfs #

#Step 2: Calculate the thetas for each site
realSFS saf2theta skj_dem_merged_norel.saf.idx -sfs skj_dem_merged_norel_folded.sfs -outname skj_theta_noLD_fold_merged_norel -fold 1
realSFS saf2theta yft_dem_merged_neutral.saf.idx -sfs yft_dem_merged_folded_neutral.sfs -outname yft_theta_noLD_fold_merged_neutral -fold 1
realSFS saf2theta bet_dem_merged_neutral.saf.idx -sfs bet_dem_merged_folded_neutral.sfs -outname bet_theta_noLD_fold_merged_neutral -fold 1

#The output from the above command are two files out.thetas.gz and out.thetas.idx. A formal description of these files can be found in the doc/formats.pdf in the angsd package. It is possible to extract the logscale persite thetas using the ./thetaStat print program.
thetaStat print skj_theta_noLD_fold_merged_norel.thetas.idx 2>/dev/null |head
thetaStat print yft_theta_noLD_fold_merged_neutral.thetas.idx 2>/dev/null |head 
thetaStat print bet_theta_noLD_fold_merged_neutral.thetas.idx 2>/dev/null |head

#Step 3a: Estimate Tajimas D and other statistics
#calculate Tajimas D

#skj 
thetaStat do_stat skj_theta_noLD_fold_merged_norel.thetas.idx
cat skj_theta_noLD_fold_merged_norel.thetas.idx.pestPG

#yft
thetaStat do_stat yft_theta_noLD_fold_merged_neutral.thetas.idx
cat yft_theta_noLD_fold_merged_neutral.thetas.idx.pestPG

#bet
thetaStat do_stat bet_theta_noLD_fold_merged_neutral.thetas.idx
cat bet_theta_noLD_fold_merged_neutral.thetas.idx.pestPG


#---Adults vs Larvae 

#yft larvae and adults
realSFS yft_larvae_merged_neutral.saf.idx -fold 1 >yft_larvae_merged_folded_neutral.sfs
realSFS yft_2020_merged_neutral.saf.idx -fold 1 >yft_2020_merged_folded_neutral.sfs

#bet larvae and adults
realSFS bet_larvae_merged_neutral.saf.idx -fold 1 >bet_larvae_merged_folded_neutral.sfs
realSFS bet_2020_merged_neutral.saf.idx -fold 1 >bet_2020_merged_folded_neutral.sfs

#Calculate thetas
realSFS saf2theta yft_larvae_merged_neutral.saf.idx -sfs yft_larvae_merged_folded_neutral.sfs -outname yft_larvae_merged_theta_fold_neutral -fold 1
realSFS saf2theta yft_2020_merged_neutral.saf.idx -sfs yft_2020_merged_folded_neutral.sfs -outname yft_2020_merged_theta_fold_neutral -fold 1

realSFS saf2theta bet_larvae_merged_neutral.saf.idx -sfs bet_larvae_merged_folded_neutral.sfs -outname bet_larvae_merged_theta_fold_neutral -fold 1
realSFS saf2theta bet_2020_merged_neutral.saf.idx -sfs bet_2020_merged_folded_neutral.sfs -outname bet_2020_merged_theta_fold_neutral -fold 1

#Calculate stats
thetaStat do_stat yft_larvae_merged_theta_fold_neutral.thetas.idx
cat yft_larvae_merged_theta_fold_neutral.thetas.idx.pestPG

thetaStat do_stat yft_2020_merged_theta_fold_neutral.thetas.idx
cat yft_2020_merged_theta_fold_neutral.thetas.idx.pestPG

thetaStat do_stat bet_larvae_merged_theta_fold_neutral.thetas.idx
cat bet_larvae_merged_theta_fold_neutral.thetas.idx.pestPG

thetaStat do_stat bet_2020_merged_theta_fold_neutral.thetas.idx
cat bet_2020_merged_theta_fold_neutral.thetas.idx.pestPG

#---Main pop vs subpop
realSFS yft_main_merged_neutral.saf.idx -fold 1 > yft_main_merged_folded_neutral.sfs
realSFS yft_subpop_merged_neutral.saf.idx -fold 1 > yft_subpop_merged_folded_neutral.sfs

#Calculate thetas
realSFS saf2theta yft_main_merged_neutral.saf.idx -sfs yft_main_merged_folded_neutral.sfs -outname yft_main_merged_theta_fold_neutral -fold 1
realSFS saf2theta yft_subpop_merged_neutral.saf.idx -sfs yft_subpop_merged_folded_neutral.sfs -outname yft_subpop_merged_theta_fold_neutral -fold 1

#Stats
thetaStat do_stat yft_main_merged_theta_fold_neutral.thetas.idx
thetaStat do_stat yft_subpop_merged_theta_fold_neutral.thetas.idx

#- Output in the ./thetaStat print thetas.idx are the log scaled per site estimates of the thetas
#- Output in the pestPG file are the sum of the per site estimates for a region
#The .pestPG file is a 14 column file (tab seperated). The first column contains information about the region. The second and third column is the reference name and the center of the window.
#We then have 5 different estimators of theta, these are: Watterson, pairwise, FuLi, fayH, L. And we have 5 different neutrality test statistics: Tajima's D, Fu&Li F's, Fu&Li's D, Fay's H, Zeng's E. The final column is the effetive number of sites with data in the window.
#Format is:
#(indexStart,indexStop)(posStart,posStop)(regStat,regStop) chrname wincenter tW tP tF tH tL tajD fulif fuliD fayH zengsE numSites
#Most likely you are just interest in the wincenter (column 3) and the column 9 which is the Tajima's D statistic.
#The first 3 columns relates to the region. The next 5 columns are 5 different estimators of theta, and the next 5 columns are neutrality test statistics. The final column is the number of sites with data in the region.

#See Jaskiel_2026_2bRAD_Final_Workflow.Rmd for plotting