# HOWTO: Transcriptome assembly (Illumina reads)

This is a rough howto on producing a transcriptome assembly from Illumina RNASeq data, with scripts specific to the Cedar cluster at Compute Canada, which uses the SLURM-scheduler.

See the [general Compute Canada wiki](https://docs.computecanada.ca/wiki/Compute_Canada_Documentation) and the [Cedar wiki](https://docs.computecanada.ca/wiki/Cedar) for additional info on job scheduling, ressource allocation, and [available software](https://docs.computecanada.ca/wiki/Available_software).


### Step 1: FASTQC
[FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) runs some useful metrics on raw reads, including overrepresented sequence motifs (e.g. adapters), and quality of reads.

    #!/bin/bash
    #SBATCH --time=00:20:00
    #SBATCH --account=def-PI
    #SBATCH --mem-per-cpu=2000M
    module load fastqc

    for file in *fastq.gz ; do fastqc fastq ${file}; done


### Step 2: Rcorrector
[Rcorrector](https://github.com/mourisl/Rcorrector) (Song and Florea 2015) provides error correction in Illumina RNA-seq data. Download and install locally in your home directory.

    #!/bin/bash
    #SBATCH --time=12:00:00
    #SBATCH --account=def-PI
    #SBATCH --mem-per-cpu=16000M
    #SBATCH --cpus-per-task=8
    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

    perl ~/progs/rcorrector/run_rcorrector.pl -p R1.fastq.gz R2.fastq.gz -t 8


### Step 3: trimmomatic
[Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) (Bolger et al. 2014) trims sequencing adapters from the corrected reads, and does some quality-trimming as well. If your RNAseq-data was generated using the SmartSeq2 protocol (Picelli et al. 2014), include the SmartSeq2 adapter in your adapter-file.

###### adapters_nex_smarter.fasta:
    >SMARTER-sequence
    AAGCAGTGGTATCAACGCAGAGTAC
    >Transposase/1
    CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
    >Transposase_rc/2
    CTGTCTCTTATACACATCTGACGCTGCCGACGA

###### Running trimmomatic:
    #!/bin/bash
    #SBATCH --time=16:00:00
    #SBATCH --account=def-PI
    #SBATCH --mem-per-cpu=4000M
    #SBATCH --cpus-per-task=8
    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
    module load java/13.0.1

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 R1.cor.fq.gz R2.cor.fq.gz R1_cor.trim_pair.fastq.gz R1_cor.trim_unpair.fastq.gz R2_cor.trim_pair.fastq.gz R2_cor.trim_unpair.fastq.gz ILLUMINACLIP:adapters_nex_smarter.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

You can run FASTQC again if you'd like.


### Step 4: Combine unpaired reads (trimmed and corrected)
    gunzip R*_cor.trim_unpair.fastq.gz
    cat R*_cor.trim_unpair.fastq > SAMPLE_cor.trim_unpair.fastq

### Step 5: rnaSPAdes
Assembly using [rnaSPAdes](http://cab.spbu.ru/software/rnaspades/) (Bushmanova et al. 2018) – this seems to perform better than Trinity, at least for single-cell assemblies and when using rcorrector.

    #!/bin/bash
    #SBATCH --time=48:00:00
    #SBATCH --account=def-PI
    #SBATCH --mem-per-cpu=250000M
    #SBATCH --cpus-per-task=12
    #SBATCH --mail-user=YOUREMAIL@BLA.COM
    #SBATCH --mail-type=ALL
    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
    module load gcc/7.3.0
    module load spades/3.13.0

    spades.py --rna -1 R1_cor.trim_pair.fastq.gz -2 R2_cor.trim_pair.fastq.gz -s R1-R2_cor.trim_unpair.fastq.gz -o SAMPLE_rna-spades -t 12
