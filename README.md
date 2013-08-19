parse-cosmic
============

A script to help carefully parse through and standardize somatic variant lists from COSMIC

Download the massive tab-delim file that contains curated mutation data from COSMIC, except for fusions:

    wget ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicCompleteExport_v66_250713.tsv.gz
    gunzip CosmicCompleteExport_v66_250713.tsv.gz

Download the tab-delim file that contains the inserted sequence for each insertion (COSMIC stores this separately for whatever reason):

    wget ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicInsMutExport_v66_250713.tsv.gz
    gunzip CosmicInsMutExport_v66_250713.tsv.gz

Run the script to parse out mutations of genome-wide screened cases, and prep per-sample Build37 variant lists:

    mkdir var_files
    perl bin/parse_cosmic.pl CosmicCompleteExport_v66_250713.tsv CosmicInsMutExport_v66_250713.tsv var_files
