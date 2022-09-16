# BlobToolKit v2 README


Also see tutorials here: https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/


## Prep your data

**As of September 2022 this is broken:**
1) The blast formats need to be updated (see https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/adding-data-to-a-dataset/adding-hits/#hits), basically the number of columns expected from blobtoolkit differs from the columns provided
2) There's something wrong with the host command. Not sure what.
3) You can try running a local installation on your computer, as far as I can tell this is only an option on OS X/Linux machines (https://github.com/blobtoolkit/blobtoolkit)
4) What semi-works is the installation via `pip install blobtoolkit` on Jezero, at least until the last step (`blobtools host --port 8080 --api-port 8000 --hostname localhost ~/blobdir/`)

#### Step 1: Run megablast against nt

Run first, then blastx overnight. Megablast vs nt: Megablast is optimized for highly similar sequences

	blastn \
	-task megablast \
	-query contigs.fasta \
	-db /Data/databases/NT_blast/nt \
	-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
	-culling_limit 5 \
	-num_threads 12 \
	-evalue 1e-25 \
	-max_target_seqs 5 > contigs_vs_nt.blastn

#### Step 2: Run diamond blastx against uniprot reference proteomes
Note that diamond blastx can be VERY SLOW – run overnight with 12-14 threads.

Diamond blastx:

	diamond blastx \
	--query contigs.fasta \
	--max-target-seqs 1 \
	--sensitive \
	--threads 12 -M 60 \
	--db /Data/databases/uniprot_ref_diamond/uniprot_ref_proteomes.dmnd \
	--evalue 1e-25 \
	--outfmt 6 \
	--out contigs_vs_uniprot_ref.mts1.1e25.out


###### If Uniprot results are not taxified (example for Jezero):
You'll have to use blobtoolsv1 to do this.

	conda activate blobtools_V1
	/opt/blobtools/blobtools taxify -f contigs_vs_uniprot_ref.mts1.1e25.out -m /Data/databases/uniprot_ref_diamond/uniprot_ref_proteomes.taxids -s 0 -t 2


#### Step 3: Map trimmed reads to the assembly
Mapping is used to assess the coverage of reads on contigs. You can do this either with minimap2 or bowtie2 (If you have Nanopore or PacBio reads you'll need to use minimap2). Bowtie2 instructions are at the very bottom.

##### minimap2 mapping of reads (short and long):

	minimap2 -t 12 \
	-a assembly.fasta \
	genus_species_R1.trim_pair.fastq.gz \
	genus_species_R2.trim_pair.fastq.gz \
	genus_species.trim_unpair.fastq.gz \
	nanopore_basecalled.fastq.gz \
	> genus_species.minimap2.sam


##### Convert SAM to BAM and sort (applies to both minimap2 and bowtie2):

	samtools view -@ 10 -bS -o genus_species.bam genus_species.sam
	samtools sort -@ 10 -o genus_species.sorted.bam genus_species.bam


## Using Blobtools v2 on your data
To use, ssh into Soyouz with ports forwarded  (change user@server.ca):
`ssh -L 8080:127.0.0.1:8080 -L 8000:127.0.0.1:8000 user@server.ca`

Activate newly created Conda blobtoolkit environment:
`conda activate btk_env` (Soyouz)
`conda activate blobtoolkit` (Jezero)

Create BlobDir:

	blobtools create \
	--fasta assembly.fasta \
	--taxid 1911741 \
	--taxdump /opt/blobtoolkit/taxdump/ \
	~/path/to/blobdir/

The taxid is from NCBI, you can add it if your organism has one. The last line is the path to where your BlobDir directory and all relevant data will be created (does not need to exist). This will create a couple of .json files in the BlobDir, and we will add to this in the next few commands.


Add taxonomic hits:

	blobtools add \
	--hits transcripts_vs_nt.blastn \
	--hits transcripts_vs_uniprot_ref.mts1.1e25.taxified.out \
	--taxrule bestsumorder \
	--taxdump /opt/blobtoolkit/taxdump/ \
	~/path/to/blobdir/


Add coverage:

	blobtools add \
	--cov ~/path/to/assembly/genus_species.sorted.bam \
	~/path/to/blobdir/


Add BUSCO hits (optional):

	blobtools add \
	--busco full_table_SAMPLE_BUSCO.tsv \
	~/path/to/blobdir/


Start viewer:

	conda deactivate
	conda activate btk_env
	/opt/blobtools2/blobtools2/blobtools host --port 8080 \
	--api-port 8000 \
	--hostname localhost \
	--viewer /opt/blobtoolkit/viewer/ \
	~/path/to/blobdir/

Then on own computer, in Firefox or Chrome, open `http://localhost:8080` (might take a while)


## Installation
###### Install Miniconda (if not already installed):
`curl https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh > Miniconda3.sh`

Make the script executable: `chmod +x Miniconda3.sh`
Execute script: `./Miniconda3.sh`

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


### ALTERNATIVE: map short reads only with bowtie2

Index the contigs for Bowtie2:

	bowtie2-build assembly.fasta assembly.fasta

bowtie2 mapping of merged reads (slow, run on as many threads as possible with -p):

	bowtie2 -x assembly.fasta -U ../genus_species.paired.trim.fastq -S genus_species.sam -p 10


### Notes to incorporate
Install instructions:

	sudo apt update && sudo apt-get -y install firefox xvfb
	conda create -n btk_env -c conda-forge -y python=3.6 docopt psutil pyyaml ujson tqdm nodejs=10 yq
	conda activate btk_env
	conda install -c bioconda -y pysam seqtk
	conda install -c conda-forge -y geckodriver selenium pyvirtualdisplay
	pip install fastjsonschema

	sudo mkdir /opt/blobtools2
	cd /opt/blobtools2;
	sudo git clone https://github.com/blobtoolkit/blobtools2;
	sudo git clone https://github.com/blobtoolkit/viewer;
	sudo git clone https://github.com/blobtoolkit/specification;
	sudo git clone https://github.com/blobtoolkit/insdc-pipeline;

	cd viewer
	sudo su -    #switch to root to get npm to work
	/opt/Anaconda3/envs/btk_env/bin/npm install    #need full path to npm when in sudo

	Testing:
	ssh -L 8080:127.0.0.1:8080 -L 8000:127.0.0.1:8000 username@server.ca
	conda activate btk_env
	/opt/blobtools2/blobtools2/blobtools host --port 8080 --api-port 8000 --hostname localhost --viewer /opt/blobtools2/viewer test_BlobDir
