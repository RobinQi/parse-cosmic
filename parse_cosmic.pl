#!/usr/bin/perl -w

use strict;
use warnings;
use IO::File;
use File::Temp;
use Getopt::Long qw( GetOptions );
use Pod::Usage qw( pod2usage );

# Default file paths and constants
my $max_indel_length = 100;
my $max_muts_per_sample = 1500;
my %valid_chrom = map {($_,1)} ( 1..22, qw( X Y MT ));
my %complement = qw( A T T A C G G C N N );
my $ref_seq_b36 = "/ifs/e63data/leew1/ref/resource/tcga-human-build36-wugsc-1.0.fasta";
my $ref_seq_b37 = "/ifs/e63data/leew1/ref/resource/GRCh37-lite.fa";
my $hg18_to_hg19_chain = "/ifs/e63data/sander-lab/kandoth/srv/hg18ToHg19.over.chain";

# Check for missing or crappy arguments
unless( @ARGV and $ARGV[0]=~m/^-/ ) {
    pod2usage( -verbose => 0, -message => "$0: Missing or invalid arguments!\n", -exitval => 2 );
}

# Parse options and print usage if there is a syntax error, or if usage was explicitly requested
my ( $man, $help ) = ( 0, 0 );
my ( $cosmic_full_file, $cosmic_ins_file, $output_dir );
GetOptions(
    'help!' => \$help,
    'man!' => \$man,
    'complete-cosmic-file=s' => \$cosmic_full_file,
    'inserted-sequence-file=s' => \$cosmic_ins_file,
    'output-dir=s' => \$output_dir,
    'b36-fasta=s' => \$ref_seq_b36,
    'b37-fasta=s' => \$ref_seq_b37,
    'liftover-chain-hg18-to-hg19=s' => \$hg18_to_hg19_chain
) or pod2usage( -verbose => 1, -input => \*DATA, -exitval => 2 );

pod2usage( -verbose => 1, -input => \*DATA, -exitval => 0 ) if( $help );
pod2usage( -verbose => 2, -input => \*DATA, -exitval => 0 ) if( $man );

# Make sure all our input files exist and are non-zero
( -s $ref_seq_b36 ) or die "Build36 reference fasta missing!\nPath: $ref_seq_b36\n";
( -s $ref_seq_b37 ) or die "Build37 reference fasta missing!\nPath: $ref_seq_b37\n";
( -s $hg18_to_hg19_chain ) or die "liftOver chain missing!\nPath: $hg18_to_hg19_chain\n";

# Load up the inserted sequences of insertions, which COSMIC stores separately, for whatever reason
my %ins_nucs = map{chomp; split(/\t/)} grep {m/^\w+\t[ACGTacgt]+$/} `cut -f 9,11 $cosmic_ins_file`;

# From the main COSMIC file, parse out SNVs and small indels, into a minimal 5-column format
my %vars = ();
my %tissue_site_counts = ();
my $cosmicFh = IO::File->new( $cosmic_full_file );
warn "Parsing COSMIC DB for usable SNVs and small-indels in GWS cases...\n";
while( my $line = $cosmicFh->getline ) {
    next if( $line =~ m/^Gene name/ );
    chomp( $line );
    my @col = split( /\t/, $line );
    my ( $samp_id, $samp_site, $samp_hist, $samp_subtype, $gws, $origin ) = @col[4,6,8,9,10,23];
    my ( $mut_id, $mut_nuc, $mut_aa, $mut_type, $zygos, $hg18_locus, $hg18_strand, $hg19_locus,
        $hg19_strand, $status ) = @col[11..20];

    # Skip variants of samples that were not "Genome-wide screened"
    next unless( $gws eq 'y' );

    # Skip variants belonging to recurrent tumors, metastases, adjacent, etc.
    # We don't want to report the same variant from the same patient more than once
    next unless( $origin =~ m/^(primary|secondary|surgery)/ );

    # Skip known germline variants or variants from samples annotated as normals
    next if( $status =~ m/germline/ or $samp_hist =~ m/NORMAL/ or $samp_subtype =~ m/normal/);

    # Skip mutation types that are too hard to annotate
    next if( $mut_type =~ m/^Whole gene deletion/ );

    # Skip variants that don't have both nucleotide change and loci available
    next unless( $mut_nuc and ( $hg18_locus or $hg19_locus ));

    # Pull b36 or b37 loci of the variant giving preference to b36, which we will liftOver later
    my $build = ( $hg18_locus ? "b36" : "b37" );
    my $locus = ( $hg18_locus ? $hg18_locus : $hg19_locus );
    my ( $chr, $start, $stop ) = $locus =~ m/^(\w+):(\d+)-(\d+)$/;
    die "Cannot identify locus for:\n$line\n" unless( $chr and $start and $stop );
    $chr =~ s/^(23|x)$/X/; $chr =~ s/^(24|y)$/Y/; $chr =~ s/^(25|M|m|mt)$/MT/;
    die "Invalid chrom name in:\n$line\n" unless( $valid_chrom{$chr} );

    # Fetch the reference sequence at this locus to do some QC
    my $ref_seq = ( $hg18_locus ? $ref_seq_b36 : $ref_seq_b37 );
    my $fetched_ref = `samtools faidx $ref_seq $chr:$start-$stop | grep -v ^\\>`;
    chomp( $fetched_ref );

    # Try to find out what kind of variant this is, and convert it to a minimal 5-column format
    my ( $ref, $var, @tmp );
    if( @tmp = $mut_nuc =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*([ACGTacgt]+)>([ACGTacgt]+)$/ ) { # SNV
        ( $ref, $var ) = ( uc( $tmp[0] ), uc( $tmp[1] ));
    }
    elsif( @tmp = $mut_nuc =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*ins([ACGTacgt]+)$/ ) { # insertion
        ( $ref, $var ) = ( "-", uc( $tmp[0] ));
        my $ins_nuc = uc( $ins_nucs{$mut_id} );
        if( length( $var ) > $max_indel_length ) {
            warn "Skipped: Long insertion >$max_indel_length bps in $mut_nuc at $build locus: $locus\n";
            next;
        }
        elsif( $ins_nuc and $var ne $ins_nuc ) {
            warn "Skipped: Inserted sequence per COSMIC ($ins_nuc) differs from $mut_nuc at $build locus: $locus\n";
            next;
        }
    }
    elsif( @tmp = $mut_nuc =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*del([ACGTacgt]+)$/ ) { # deletion
        ( $ref, $var ) = ( uc( $tmp[0] ), "-" );
        if( length( $ref ) > $max_indel_length ) {
            warn "Skipped: Long deletion >$max_indel_length bps in $mut_nuc at $build locus: $locus\n";
            next;
        }
        elsif( length( $ref ) != ( $stop - $start + 1 )) {
            warn "Skipped: Length of deleted sequence in $mut_nuc doesn't match $build locus: $locus\n";
            next;
        }
    }
    elsif( @tmp = $mut_nuc =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*del(\d+)$/ ) { # deletion without sequence specified
        my $del_length = $tmp[0];
        if( $del_length > $max_indel_length ) {
            warn "Skipped: Long deletion >$max_indel_length bps in $mut_nuc at $build locus: $locus\n";
            next;
        }

        if( length( $fetched_ref ) != $del_length ) {
            warn "Skipped: Length of deleted sequence in $mut_nuc doesn't match $build locus: $locus\n";
            next;
        }
        ( $ref, $var ) = ( uc( $fetched_ref ), "-" );
    }

    unless( $ref and $var ) {
        warn "Skipped: Cannot identify ref/var from $mut_nuc at $build locus: $locus\n";
        next;
    }

    unless(( length($ref)==1 && length($var)==1 ) || ( length($ref)==2 && length($var)==2 ) ||
           ( $ref eq "-" && $var=~m/[ACGT]+/ ) || ( $ref=~m/[ACGT]+/ && $var eq "-" )) {
        warn "Skipped: Cannot annotate $ref/$var at $build locus: $locus\n";
        next;
    }

    # Reverse complement the ref/var if the variant is from the negative strand
    my $locus_strand = ( $hg18_strand ? $hg18_strand : $hg19_strand );
    my $rc_ref = reverse map {$complement{$_}} split( //, $ref );
    if( $locus_strand eq "-" and $rc_ref eq $fetched_ref ) {
        if( $ref eq "-" && $var =~ m/[ACGT]+/ ) {
            $var = reverse map {$complement{$_}} split( //, $var );
        }
        elsif( $ref =~ m/[ACGT]+/ && $var eq "-" ) {
            $ref = $rc_ref;
        }
        else {
            $ref = $rc_ref;
            $var = reverse map {$complement{$_}} split( //, $var );
        }
    }

    # Skip SNVs and deletions with reference alleles that don't match the reference fasta
    if( $ref ne $fetched_ref and $ref ne "-" ) {
        warn "Skipped: Variant's ref allele $ref is $fetched_ref in the $build fasta at locus: $locus\n";
        next;
    }

    # If we got to this point, then the variant can be annotated. Store info about it into a hash
    $samp_site = $samp_hist if( $samp_site eq "NS" );
    ++$tissue_site_counts{$samp_site} unless( defined $vars{$samp_id} );
    $vars{$samp_id}{$build}{"$chr\t$start\t$stop\t$ref\t$var\n"} = 1;
}
$cosmicFh->close;

# Create a variant file for each sample, in the 5-column format
my ( $samp_count, $tot_mut_count ) = ( 0, 0 );
warn "\nRunning liftOver and prepping Build37 files for transcript annotation...\n";
foreach my $samp_id ( keys %vars ) {

    # Make temporary files to use with liftOver
    my $build36_fh = File::Temp->new;
    my $build37_fh = File::Temp->new;
    my $unmapped_fh = File::Temp->new;
    ( $build36_fh and $build37_fh and $unmapped_fh ) or die "Couldn't create a temp file. $!";

    # Write build36 vars to a file in 5-column format, but use 0-based start since liftOver expects BED format
    my $mut_count = 0;
    if( defined $vars{$samp_id}{b36} && scalar( keys %{$vars{$samp_id}{b36}} ) > 0 ) {
        $mut_count += scalar( keys %{$vars{$samp_id}{b36}} );
        foreach my $line ( keys %{$vars{$samp_id}{b36}} ) {
            my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
            $build36_fh->print( join( "\t", $chr, $start-1, $stop, $ref, $var ));
        }

        # Run liftOver to convert the Build36 variants to Build37
        my $lift_cmd = join( " ", "liftOver", $build36_fh, $hg18_to_hg19_chain, $build37_fh, $unmapped_fh );
        print STDERR "Running liftOver as follows:\n$lift_cmd\n";
        system( $lift_cmd ) or die "Failed to run 'liftOver' for variants in $samp_id\n";
    }

    # Append additional b37 loci to the liftOver'ed file, and remove any duplicates
    if( defined $vars{$samp_id}{b37} && scalar( keys %{$vars{$samp_id}{b37}} ) > 0 ) {
        $mut_count += scalar( keys %{$vars{$samp_id}{b37}} );
        foreach my $line ( keys %{$vars{$samp_id}{b37}} ) {
            my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
            $build37_fh->print( join( "\t", $chr, $start-1, $stop, $ref, $var ));
        }
    }

    # Remove duplicate variants that are identifiable only after liftOver
    my $build37_file = $build37_fh->filename;
    print `sort -u $build37_file -o $build37_file`;

    # Sort by loci, and write to output dir, unless this sample is hypermutated
    if( -s $build37_file && $mut_count <= $max_muts_per_sample ) {
        ++$samp_count;
        $tot_mut_count += $mut_count;
        print `joinx sort -s -i $build37_file -o $output_dir/$samp_id.var`;
    }
    else {
        warn "Skipped: Hypermutated sample $samp_id with $mut_count variants\n";
    }
    exit;
}

print "\nFetched $tot_mut_count SNVs and small indels across $samp_count samples from the following tissue types:\n";
foreach my $key ( sort {$tissue_site_counts{$b} <=> $tissue_site_counts{$a}} keys %tissue_site_counts ) {
    print $tissue_site_counts{$key} . "\t$key\n";
}

__DATA__

=head1 NAME

 parse_cosmic.pl - Parse the giant tab-delimited lists of somatic variants downloadable from COSMIC

=head1 SYNOPSIS

 perl parse_cosmic.pl --complete-cosmic-file CosmicCompleteExport.tsv --inserted-sequence-file CosmicInsMutExport.tsv --output-dir var_files

=head1 OPTIONS

=over 8

=item B<--complete-cosmic-file>

 The CosmicCompleteExport file from ftp.sanger.ac.uk/pub/CGP/cosmic/data_export

=item B<--inserted-sequence-file>

 The CosmicInsMutExport file from ftp.sanger.ac.uk/pub/CGP/cosmic/data_export

=item B<--output-dir>

 An output directory for per-sample variant lists

=item B<--b36-fasta>

 A Build36 (HG18) reference fasta sequence with an accompanying fai index (Default: /ifs/e63data/leew1/ref/resource/tcga-human-build36-wugsc-1.0.fasta)

=item B<--b37-fasta>

 A Build37 (HG19) reference fasta sequence with an accompanying fai index (Default: /ifs/e63data/leew1/ref/resource/GRCh37-lite.fa)

=item B<--liftover-chain-hg18-to-hg19>

 A liftOver chain file to convert hg18 loci to hg19 (Default: /ifs/e63data/sander-lab/kandoth/srv/hg18ToHg19.over.chain)

=item B<--help>

 Prints a brief help message and quits

=item B<--man>

 Prints the detailed manual

=back

=head1 DESCRIPTION

 Takes the giant tab-delimited lists of somatic variants downloadable from COSMIC, fix their various issues, and generate a standardized Build37 variant list in 5-column format

 Here is a summary of what this script does, when parsing out COSMIC variants:

  1. Use only variants with a value of "y" under column Genome-wide screen. This indicates that whole-genome or exome-seq was performed
  2. Use only samples annotated as "primary" or "secondary" under column "Tumour origin"
  3. Skip germline variants or variants from samples annotated as normals
  4. Skip variants that don't have both nucleotide change and loci provided
  5. Skip insertions and deletions longer than 100 base-pairs
  6. Skip deletions where the length of deleted base-pairs doesn't fit the provided genomic loci
  7. Skip complex indels which are not supported in MAF format
  8. If variant loci is only reported on Build36, liftOver to Build37
  9. If reference/variant alleles are reported on the reverse strand, standardize it to the forward strand instead
 10. Skip variants where the reported reference allele doesn't match the reference fasta sequence at that locus
 11. After all the above steps, skip hypermutated samples that report more than 1500 variants

=head1 AUTHOR

 Cyriac Kandoth (ckandoth@gmail.com)

=head1 LICENSE

 LGPLv3 (The Genome Institute, Washington University, St Louis, MO 63108, USA)

=cut
