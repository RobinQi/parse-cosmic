parse-cosmic
============

This is a script to carefully parse through and standardize somatic variants downloadable from Sanger's [COSMIC database](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/). As inputs, you will need 2 tab-delimited files from COSMIC, reference fasta files for Build36 (hg18) and Build37 (hg19), and a liftOver chain file to convert Build36 variant loci to Build37 where necessary. The script also requires liftOver, samtools, and the bedtools' binaries in your system PATH. For more details, keep reading.

Dependencies
------------

Make sure that **liftOver** is accessible to this script via command-line. Source code or binaries for OSX/Linux can be [downloaded here](http://hgdownload.cse.ucsc.edu/admin/exe/).

You will also need **samtools** and **bedtools**. Debian or Ubuntu users can install them as follows:

    sudo apt-get install samtools bedtools

For Fedora, CentOS6 or RHEL6 users:

    sudo yum install samtools BEDTools

All other users can compile samtools and bedtools binaries from source code:

    git clone https://github.com/arq5x/bedtools.git
    cd bedtools
    make
    cd ..
    git clone https://github.com/samtools/samtools.git
    cd samtools
    make
    cd ..

After compilation, the binaries can be moved somewhere that your system PATH can find it. For example:

    sudo mv bedtools/bin/* /usr/local/bin/
    sudo mv samtools/samtools /usr/local/bin/

Input data
----------

Get the latest tab-delimited file from COSMIC that lists somatic mutations, but excludes fusions:

    wget ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicCompleteExport_v*.tsv.gz
    gunzip CosmicCompleteExport_v*.tsv.gz

Get the file that lists sequences for each insertion (COSMIC stores this separately, dunno why):

    wget ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicInsMutExport_v*.tsv.gz
    gunzip CosmicInsMutExport_v*.tsv.gz

Get the liftOver chain to convert Build36 to Build37:

    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
    gunzip hg18ToHg19.over.chain.gz

Get Build37 and Build36 reference fasta files from Ensembl, and generate indexes for them:

    wget ftp://ftp.ensembl.org/pub/release-72/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.72.dna.primary_assembly.fa.gz
    gunzip Homo_sapiens.GRCh37.72.dna.primary_assembly.fa.gz
    samtools faidx Homo_sapiens.GRCh37.72.dna.primary_assembly.fa
    wget ftp://ftp.ensembl.org/pub/release-54/fasta/homo_sapiens/dna/Homo_sapiens.NCBI36.54.dna.toplevel.fa.gz
    gunzip Homo_sapiens.NCBI36.54.dna.toplevel.fa.gz
    samtools faidx Homo_sapiens.NCBI36.54.dna.toplevel.fa

Notes
-----

Type `perl parse_cosmic.pl --help` for usage, and `perl parse_cosmic.pl --man` for the full manual.
