# BlobToolKit v2 README


Also see tutorials here: https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/

## Changes as of November 2022

1) The blastn and diamondblast commands are reformatted, and there's a new uniprot database version, which removes the previously required 'taxify' step.

2) There is no need for conda anymore. Follow the steps below with `pipenv`.

3) If you have a Mac or Unix machine, you can still always run a local blobtoolkit installation on your computer if the final viewing step fails: https://github.com/blobtoolkit/blobtoolkit
After install, download your blobdirs and host them via: `blobtools view --local path/to/blobdir`


## Prep your data

#### Step 1: Run megablast against nt

Run first, then blastx overnight. Megablast vs nt: Megablast is optimized for highly similar sequences

	blastn \
	-task megablast \
	-query contigs.fasta \
	-db /Data/databases/NT_blast/nt \
	-outfmt '6 qseqid staxids bitscore std' \
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
	--db /Data/databases/uniprot_ref_diamond_2022/reference_proteomes.dmnd \
	--evalue 1e-25 \
	--outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
	--out contigs_vs_uniprot_ref.mts1.1e25.out



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
To use, ssh into Jezero with ports forwarded  (change user@server.ca):
`ssh -L 8080:127.0.0.1:8080 -L 8000:127.0.0.1:8000 user@server.ca`


### Activate conda blobtoolkit environment

	conda activate blobtoolkit-pip


Create BlobDir:

	/opt/Anaconda3/envs/blobtoolkit-pip/bin/blobtools create \
	--fasta assembly.fasta \
	--taxid 1911741 \
	--taxdump /opt/blobtoolkit/taxdump_new/ \
	~/path/to/blobdir/

The taxid is from NCBI, you can add it if your organism has one. The last line is the path to where your BlobDir directory and all relevant data will be created (does not need to exist). This will create a couple of .json files in the BlobDir, and we will add to this in the next few commands.


Add taxonomic hits:

	/opt/Anaconda3/envs/blobtoolkit-pip/bin/blobtools add \
	--hits transcripts_vs_nt.blastn \
	--hits transcripts_vs_uniprot_ref.mts1.1e25.out \
	--taxrule bestsumorder \
	--taxdump /opt/blobtoolkit/taxdump_new/ \
	~/path/to/blobdir/


Add coverage:

	/opt/Anaconda3/envs/blobtoolkit-pip/bin/blobtools add \
	--cov ~/path/to/assembly/genus_species.sorted.bam \
	~/path/to/blobdir/


Add BUSCO hits (optional):

	/opt/Anaconda3/envs/blobtoolkit-pip/bin/blobtools add \
	--busco full_table_SAMPLE_BUSCO.tsv \
	~/path/to/blobdir/


Start viewer:

	/opt/Anaconda3/envs/blobtoolkit-pip/bin/blobtools host --port 8080 \
	--api-port 8000 \
	--hostname localhost \
	~/path/to/blobdir/

Then on own computer, in Firefox or Chrome, open `http://localhost:8080` (might take a while)


### ALTERNATIVE way to start viewer if above gives you an 'can't access API' error

	conda deactivate
	conda deactivate
	conda activate blobtoolkit-pip
	/opt/Anaconda3/envs/blobtoolkit-pip/bin/blobtools view --remote path_to_your_blobtools_directory

It will say for remote access use: `ssh -L 8003:127.0.0.1:8003 -L 8002:127.0.0.1:8002 username@remote_host`, this can be different every time so copy what it says when you run it.
It will also say view dataset at ... and provide a link, copy that link and open in Chrome or Firefox.
Open a new terminal and log in using the proxy prompt from the previous terminal (make sure to change username@remote_host).
Once the new terminal is logged in with the appropriate proxy the viewer in browser should work.


## Installing blobtoolkit locally (OLD, DEPRECATED AS OF OCT 2023)

Install the most recent version of blobtoolkit locally (i.e. only for yourself), in a virtual pip environment. The installation step will take a few minutes:

	conda deactivate
	curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
	pipenv install blobtoolkit==4.1.5 numpy==1.24.4 pandas==2.0.3 scipy=1.10.1
	pipenv shell

The last command will enter the virtual environment. From here, proceed with the commands below as per usual.


## Installation (OLD, DEPRECATED)
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


### If Uniprot results are not taxified (OLD – YOU DO NOT USUALLY NEED TO DO THIS):
You'll have to use blobtoolsv1 to do this.

	conda activate blobtools_V1
	/opt/blobtools/blobtools taxify -f contigs_vs_uniprot_ref.mts1.1e25.out -m /Data/databases/uniprot_ref_diamond/uniprot_ref_proteomes.taxids -s 0 -t 2

### OLD instructions for starting viewer

	conda deactivate
	conda activate btk_env
	/opt/blobtools2/blobtools2/blobtools host --port 8080 \
	--api-port 8000 \
	--hostname localhost \
	--viewer /opt/blobtoolkit/viewer/ \
	~/path/to/blobdir/

Then on own computer, in Firefox or Chrome, open `http://localhost:8080` (might take a while)


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
