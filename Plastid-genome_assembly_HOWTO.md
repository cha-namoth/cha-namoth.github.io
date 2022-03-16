# HOWTO: Plastid-genome assembly (one of many different ways)

## General bits and bobs

You can run all the commands with `screen`, or with `nohup`.

##### screen
To run in `screen`, just type `screen -S [desired name of your session]`, and it will open a new screen session for you. Run your command, and to exit that screen and to keep the command running in the background:

1) press CTRL+A, release the keys
2) press D to detach.
3) to reconnect to that session/screen, type `screen -r [session name]`
4) if you want to close a session (it also terminates any command you've been running in there), type `exit` when you're in that screen.

##### nohup
For `nohup`, you'll have to first make a quick bash-script with the command you want to execute, with your favorite editor (e.g. nano). You can also do this on your local computer and upload the script to the server, just make sure it's in plain text – i.e. don't use Word (TextEdit kind of works, but I recommend either SublimeText or Atom):

    nano [script].sh
    add: `#!/bin/bash` in first line
    copy desired command below it
    save and exit script
    chmod +x [script].sh (makes it executable)
    nohup ./[script].sh &

Your script is now running the command in the background and you can do other stuff!



## Trimming your reads
Trim off the sequencing adapters and low-quality parts of the reads with trimmomatic (http://www.usadellab.org/cms/index.php?page=trimmomatic).

    #!/bin/bash
    java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -phred33 \
    NeomP_S1_S1_R1_001.fastq.gz \
    NeomP_S1_S1_R2_001.fastq.gz \
    NeomP_S1_R1.trim_pair.fastq.gz \
    NeomP_S1_R1.trim_unpair.fastq.gz \
    NeomP_S1_R2.trim_pair.fastq.gz \
    NeomP_S1_R2.trim_unpair.fastq.gz \
    ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:25


PE = use paired-end mode\
-threads 4 = use 4 threads\
-phred33 = use phred33 to interpret quality scores\
This is followed by 2 lines of input files (R1 and R2), and then 4 lines of output files (2 paired, 2 unpaired)\
ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 = adapter file, plus instructions how stringent to be in trimming\
SLIDINGWINDOW:4:20 = Scan reads with a 4-base sliding window, trimming when average quality per base drops below 20\
MINLEN:25 = discard any reads shorter than 25bp after trimming

**IMPORTANT:** The `NexteraPE-PE.fa` file containing adapter sequences to be trimmed needs to be in the same directory you're running this script from.


## Genome Assembly with SPAdes

### A bit of background on assemblers
As hinted at above, there are several genome assemblers available, and it's an evolving field. We will use SPAdes (https://cab.spbu.ru/files/release3.15.0/manual.html) since it appears to generally do a good job of assemblies. It also depends what you are assembling – eukaryotic genomes can be very hard to assemble due to repeating and non-coding regions, something that most bacterial genomes don't suffer from. To complicate things even more, it also depends on what your starting material looks like: We have both short (Illumina) and long reads (MinION/Nanopore), and SPAdes can utilize both of them. Many assemblers can only use one or the other (e.g. Velvet for Illumina).

MinION/Nanopore reads are also kind of special, since they're super long, but also much more full of sequencing errors. Special assemblers that ONLY use these long reads exist and are getting better and better (e.g. canu, flye). They usually still require 'polishing', which basically just means correcting sequencing errors on the final assemblies using a variety of different methods. But we'll be only using SPAdes for now.

Additional reading:
https://bioinformatics.uconn.edu/bacterial-genome-assembly-tutorial/


### Assembly with SPAdes

    #!/bin/bash
    spades.py -o [OUTPUT-FOLDER] \
    --pe1-1 Sample.R1.trim_pair.fastq.gz \
    --pe1-2 Sample.R2.trim_pair.fastq.gz \
    --pe1-s Sample.trim_unpair.fastq.gz \
    --nanopore Sample.nanopore.fastq.gz \
    -t 16 -m 250 --sc

-o = desired output folder\
--pe1-1 = fastq-file (paired end reads, R1)\
--pe1-2 = fastq-file (paired end reads, R2)\
--pe1-s = fastq-file (unpaired reads)\
--nanopore = nanopore reads\
-t = number of threads to use\
-m = amount of memory in GB (250 will be probably enough and is the default)\
--sc = single-cell mode (accounts for uneven coverage from whole-genome amplification)

The assembly might take a few hours to several days, depending on how many reads will be assembled. The final assembly will be `contigs.fasta` in your desired output folder.

The fasta-headers in a genome assembly vary with the exact assembly program used, but for SPAdes they have this format:

    >NODE_1_length_39736_cov_60.147098
    AATAAGGTATTCGATTGCTGTGAAAGCTCTTTTTGAGGATGGGAGAATAATGGAAATTTT
    CGCACAAGGGATTCAAGTAGAAAATTCTAAAAATATGACAAACTTACAACCATTCTCAAT
    >NODE_2_length_19886_cov_92.469417
    AGCAATATTCAGAGGATTGGAACATGTTGATAAAGAAGTTTTCACAGCAATCGAATACCT
    TATTAGAGAAGTAGAGTATAAAACCATAAATATTTCTGAACTTTTCTACTTGTTGATGTT

'NODE_1' = the contig number \
'length_39736' = its length \
'cov_60.147098' = its relative coverage

The lower the NODE-number, the longer it is – NODE_1 is always the longest in a SPAdes-assembly. Generally the fewer, but the longer the contigs the better, as a rough rule for genome assemblies.



## Assembly evaluation (QUAST and BUSCO)

Once the assembly has finished, we can evaluate it to give us a better idea of, for example, how many contigs have been produced, what is the biggest contig, etc. QUAST is a good tool for this.


##### General assembly metrics: QUAST

    conda activate QUAST
    quast -t 2 -o [OUTPUT_QUAST] contigs.fasta
    conda deactivate

The report files in the output folder have all the relevant results. For example, maybe we can see that the largest contig is >100.000bp – that might well be a part of the plastid genome. Likely we'll have some larger contigs, and tons of smaller ones.
´

##### Genome completeness: BUSCO

    conda activate busco_v5
    busco -c 8 -m geno -i contigs.fasta -l [LINEAGE] -o [OUTPUT_BUSCO]
    conda deactivate

This is another tool we often use to see how 'complete' a genome or transcriptome is. BUSCO looks for core conserved genes, and uses that presence/absence as a proxy to how complete the genome is (the more BUSCOs recovered, the more complete). This is somewhat of an imperfect proxy, since many genes are just not that conserved across all eukaryotes, so BUSCO has several lineage-specific datasets to choose from (the list also the number of genes): https://busco.ezlab.org/list_of_lineages.html

Here we'll run `eukaryota_odb10` and `chlorophyta_odb10` (run separately).




## Assembly exploration with BlobToolKit

First we want to identify all the contigs (assembled sequencing reads, making longer sequences) that are bacterial. To do this, we routinely use a tool called 'BlobTools' (https://blobtoolkit.genomehubs.org). It's interactive, and you have a visual representation of what each contig was identified as (colour), its length (size of the circle), GC-content (x-axis), and relative coverage (y-axis; how many reads mapped back to that contig, the more the higher the coverage).
In most blobplots it's very obvious that there are many big bacterial contigs (EXAMPLE). Blobtools actually lets you select any contigs identified as bacteria (or anything else for that matter), and extract their headers.

I also wrote something on how to actually use blobtoolkit: https://github.com/cha-namoth/cha-namoth.github.io/blob/master/BlobToolKitv2_conda_HOWTO.md

You can start with 'Prep your data' in the guide linked above, since everything else is already installed.
At Step 3, use only the minimap2-approach, since we both have long and short read data, and bowtie2 can't handle long reads.

At the very last step, be sure to use Firefox when accessing http://localhost:8080, nothing else will work properly.


## Plastid genome assembly with NOVOplasty




## Plastid genome annotation with MFannot

There are several automatic plastid genome annotators out there, including GeSeq: https://chlorobox.mpimp-golm.mpg.de/geseq.html
MFannot seems to perform reasonably well with extreme genomes like those from apicomplexans, and they even have a web interface if that is easier: https://megasun.bch.umontreal.ca/cgi-bin/mfannot/mfannotInterface.pl

To run MFannot on Jezero, we need to first setup Docker and start a Docker-container for MFannot. For this, regular users (ie. non-sudo) need to run this rootless setup first:

`/usr/bin/dockerd-rootless-setuptool.sh install`

Then start the docker-service: `systemctl --user start docker`


Go to your directory where you want to run MFannot from:

    cd /path/to/directory/with/files
    docker run -it --mount type=bind,src=$(pwd),dst=$(pwd) -w $(pwd) nbeck/mfannot

This brings you to a new prompt (starts with `root`). You are now in the Docker container, and you only have access to the files in that particular directory you started the Docker command from (see above). Run this to annotate your fasta:

    mfannot --sqn --tbl input.fasta

`--sqn` = produce a `.sqn` file (can be converted to a flat GenBank file)\
`--tbl` = produce a `.tbl` file (table, easy to read)\
default is a `.new` file, which can be converted to GFF (see below)

Type `exit` to quit the Docker container.

The Docker container is still running – stop and remove it, but first show the ID and names for active containers:

    docker ps -a
    docker stop [container-ID/name]
    docker rm [container-ID/name]

### Converting MFannot output
Convert `.new` to GFF (can import into Geneious):

    agat_convert_mfannot2gff.pl -m [mfannot.new] -o [mfannot.gff]

Convert `.sqn` to flat-file GenBank (can import into Geneious):

    asn2gb -i [mfannot.sqn] -o [mfannot.gb]
