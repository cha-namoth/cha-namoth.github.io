# BlobToolKit v2 README


Also see tutorials here: https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/

## Installation
###### Install Conda:
`curl https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh > Miniconda3.sh`

Make the script executable:
`chmod +x Miniconda3.sh`

Logout, log back in.


###### Install BlobToolKit:
Create new conda environment and install BlobToolKit:

	conda create -n btk_env -c conda-forge -y python=3.6 docopt psutil pyyaml ujson tqdm nodejs=10 yq
	conda activate btk_env
	conda install -c bioconda -y pysam seqtk
	conda install -c conda-forge -y geckodriver selenium pyvirtualdisplay
	pip install fastjsonschema

	mkdir -p ~/progs/blobtoolkit
	cd ~/progs/blobtoolkit
	git clone https://github.com/blobtoolkit/blobtools2
	git clone https://github.com/blobtoolkit/viewer
	git clone https://github.com/blobtoolkit/specification
	git clone https://github.com/blobtoolkit/insdc-pipeline

	cd viewer
	npm install
	cd ..

(there might be a few errors during `npm install` – it's probably fine)


###### Fetch NCBI taxdump:
	mkdir -p taxdump
	cd taxdump
	curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -

## Prep your data
#### Step 1: Run megablast against nt

Run first, then blastx overnight. Megablast vs nt: Megablast is optimized for highly similar sequences

	blastn \
	-task megablast \
	-query Trinity.fasta \
	-db /scratch/NCBI_NT/nt \
	-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
	-culling_limit 5 \
	-num_threads 12 \
	-evalue 1e-25 \
	-max_target_seqs 5 > Trinity_vs_nt.blastn

#### Step 2: Run diamond blastx against uniprot reference proteomes
uniprot db was converted for diamond by Filip, the taxids file was provided by Filip as well and can be found under `/scratch/uniprot/`.
Note that diamond blastx is VERY SLOW – run overnight with 24 threads.

Diamond blastx:

	/opt/bin/diamond blastx \
	--query Trinity.fasta \
	--db /scratch/uniprot/uniprot_ref_proteomes.diamond.dmnd \
	--outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
	--sensitive \
	--max-target-seqs 1 \
	--evalue 1e-25 \
	--threads 8 \
	--out Trinity.fasta.vs.uniprot_ref.mts1.1e25.out

Taxify diamond results (this currently doesn't work in blobtoolsv2 – working on it. In the meantime, taxify the uniprot-output with blobtoolsv1):

	/opt/blobtools/blobtools taxify -f Trinity.fasta.vs.uniprot_ref.mts1.1e25.out -m /scratch/uniprot/uniprot_ref_proteomes.taxids -s 0 -t 2

*TODO: Taxify UniProt database with NCBI taxdump: https://blobtoolkit.genomehubs.org/install/#databases*

#### Step 3: Map trimmed reads to the assembly
Mapping is used to assess the coverage of reads on contigs.

Index the contigs for Bowtie2:

	bowtie2-build Trinity.fasta Trinity.fasta

bowtie2 mapping of merged reads (slow, run on as many threads as possible with -p):

	bowtie2 -x Trinity.fasta -U ../genus_species.paired.trim.fastq -S genus_species.sam -p 10

*TODO: Add mapping via minimap2 (much faster)*


## Using Blobtools v2 on your data
To use, ssh into Soyouz with ports forwarded  (change user@server.ca):
`ssh -L 8080:127.0.0.1:8080 -L 8000:127.0.0.1:8000 user@server.ca`

Activate newly created Conda blobtoolkit environment:
`conda activate btk_env`


Create BlobDir:

	~/progs/blobtoolkit/blobtools2/blobtools create \
	--fasta ~/path/to/assembly/assembly_transcripts.fasta \
	--taxid 1911741 \
	--taxdump ~/progs/blobtoolkit/taxdump \
	~/path/to/blobdir/

The taxid is from NCBI, add it if your organism has one. The last line is the path to where your BlobDir directory and all relevant data will be created (does not need to exist). This will create a couple of .json files in the BlobDir, and we will add to this in the next few commands.


Add taxonomic hits:

	~/progs/blobtoolkit/blobtools2/blobtools add \
	--hits ~/path/to/assembly/transcripts_vs_nt.blastn \
	--hits ~/path/to/assembly/transcripts_vs_uniprot_ref.mts1.1e25.taxified.out \
	--taxrule bestsumorder \
	--taxdump ~/progs/blobtoolkit/taxdump \
	~/path/to/blobdir/


Add coverage:

	~/progs/blobtoolkit/blobtools2/blobtools add \
	--cov ~/path/to/assembly/SPO2_transcripts_mapped.sorted.bam \
	~/path/to/blobdir/


Add BUSCO hits:

	~/progs/blobtoolkit/blobtools2/blobtools add \
	--busco ~/path/to/assembly/SAMPLE_BUSCO_full_table.tsv \
	~/path/to/blobdir/


Start viewer:

	~/progs/blobtoolkit/blobtools2/blobtools host --port 8080 \
	--api-port 8000 \
	--hostname localhost \
	--viewer ~/progs/blobtoolkit/viewer \
	~/path/to/blobdir/


Then on own computer, in Firefox, open (might take a while):
`http://localhost:8080`
