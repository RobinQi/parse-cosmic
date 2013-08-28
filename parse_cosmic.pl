#!/usr/bin/perl -w

use strict;
use warnings;
use IO::File;
use File::Temp qw( tempdir );
use Getopt::Long qw( GetOptions );
use Pod::Usage qw( pod2usage );

# Default file paths and constants
my $max_indel_length = 100;
my $max_muts_per_sample = 1500;
my %valid_chrom = map{($_,1)} ( 1..22, qw( X Y MT ));
my %complement = qw( A T T A C G G C N N - - );
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
    if( @tmp = $mut_nuc =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*([ACGTacgt]+)>([ACGTacgt]+)$/ ) { # SNP/DNP/ONP
        ( $ref, $var ) = ( uc( $tmp[0] ), uc( $tmp[1] ));
    }
    elsif( @tmp = $mut_nuc =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*ins([ACGTacgt]+)$/ ) { # Insertion
        ( $ref, $var ) = ( "-", uc( $tmp[0] ));
        unless( defined $ins_nucs{$mut_id} ) {
            warn "Skipped: Sequence unavailable for insertion with mutation ID# $mut_id at $build locus: $locus\n";
            next;
        }
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
    elsif( @tmp = $mut_nuc =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*del([ACGTacgt]+)$/ ) { # Deletion
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
    elsif( @tmp = $mut_nuc =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*del(\d+)$/ ) { # Deletion without sequence specified
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

    unless(( length($ref) == length($var) && length($ref) <= 3 ) ||
           ( $ref eq "-" && $var=~m/[ACGT]+/ ) || ( $ref=~m/[ACGT]+/ && $var eq "-" )) {
        warn "Skipped: Cannot annotate $ref/$var at $build locus: $locus\n";
        next;
    }

    # Reverse complement the ref/var if the variant is from the negative strand
    my $locus_strand = ( $hg18_strand ? $hg18_strand : $hg19_strand );
    my $rc_ref = reverse map{$complement{$_}} split( //, $ref );
    # When possible, trust the fetched_ref to determine strand. COSMIC's strand info is sometimes wrong
    if( $rc_ref eq $fetched_ref or ( $ref eq "-" and $locus_strand eq "-" )) {
        if( $ref eq "-" && $var =~ m/[ACGT]+/ ) {
            $var = reverse map{$complement{$_}} split( //, $var );
        }
        elsif( $ref =~ m/[ACGT]+/ && $var eq "-" ) {
            $ref = $rc_ref;
        }
        else {
            $ref = $rc_ref;
            $var = reverse map{$complement{$_}} split( //, $var );
        }
    }

    # Skip SNVs and deletions with reference alleles that don't match the reference fasta
    if( $ref ne $fetched_ref and $ref ne "-" ) {
        warn "Skipped: Variant's ref allele $ref is $fetched_ref in the $build fasta at locus: $locus\n";
        next;
    }

    # Save information on the tumor type/histology where available
    $samp_site = $samp_hist if( $samp_site eq "NS" );
    ++$tissue_site_counts{$samp_site} unless( defined $vars{$samp_id} );

    # If we got to this point, then the variant can be annotated. Store info about it into a hash
    # Also add a chr-prefix and use 0-based start loci, to simplify liftOver and sortBed done later
    $start--;
    $vars{$samp_id}{$build}{"chr$chr\t$start\t$stop\t$ref\t$var\n"} = 1;
}
$cosmicFh->close;

# Create a variant file for each sample, in the 5-column format
my ( $samp_count, $tot_mut_count ) = ( 0, 0 );
warn "\nRunning liftOver and prepping Build37 files for transcript annotation...\n";
foreach my $samp_id ( keys %vars ) {

    # Make temporary files to use with liftOver
    my $tmp_dir = tempdir;
    ( -e $tmp_dir ) or die "Couldn't create a temporary directory!";
    my $build36_file = "$tmp_dir/build36";
    my $build37_file = "$tmp_dir/build37";
    my $unmapped_file = "$tmp_dir/unmapped";

    # Write build36 vars to a file in 5-column format, but use 0-based start since liftOver expects BED format
    if( defined $vars{$samp_id}{b36} && scalar( keys %{$vars{$samp_id}{b36}} ) > 0 ) {
        my $build36_fh = IO::File->new( $build36_file, ">" );
        foreach my $line ( keys %{$vars{$samp_id}{b36}} ) {
            $build36_fh->print( $line );
        }
        $build36_fh->close;

        # Run liftOver to map the Build36 variants to Build37. Discard those that cannot be mapped
        print `liftOver $build36_file $hg18_to_hg19_chain $build37_file $unmapped_file >/dev/null 2>&1`;
    }

    # Append additional b37 loci to the liftOver'ed file, and remove any duplicates
    if( defined $vars{$samp_id}{b37} && scalar( keys %{$vars{$samp_id}{b37}} ) > 0 ) {
        my $build37_fh = IO::File->new( $build37_file, ">>" );
        foreach my $line ( keys %{$vars{$samp_id}{b37}} ) {
            $build37_fh->print( $line );
        }
        $build37_fh->close;
    }

    # Remove duplicate variants that are identifiable only after liftOver, and count those remaining
    print `sort -u $build37_file -o $build37_file`;
    my ( $mut_count ) = map{split(/\s+/)}`wc -l $build37_file`;

    # Sort by loci, restore 1-based start loci, and write to output dir, unless this sample is hypermutated
    if( -s $build37_file && $mut_count <= $max_muts_per_sample ) {
        ++$samp_count;
        $tot_mut_count += $mut_count;
        unlink( "$output_dir/$samp_id.var" ) if( -e "$output_dir/$samp_id.var" );
        print `sed 's/^chr//' $build37_file | sortBed | awk '{OFS="\\t"; ++\$2; print}' > $output_dir/$samp_id.var`;
    }
    elsif( ! -s $build37_file ) {
        die "Unexpected empty file $build37_file!";
    }
    else {
        warn "Skipped: Hypermutated sample $samp_id with $mut_count variants\n";
    }
}

print "\nFetched $tot_mut_count SNVs and small indels across $samp_count samples from the following tissue types:\n";
foreach my $key ( sort {$tissue_site_counts{$b} <=> $tissue_site_counts{$a}} keys %tissue_site_counts ) {
    print $tissue_site_counts{$key} . "\t$key\n";
}

__DATA__

=head1 NAME

 parse_cosmic.pl - Parse the giant tab-delimited lists of somatic variants downloadable from COSMIC

=head1 SYNOPSIS

 perl parse_cosmic.pl --complete-cosmic-file CosmicCompleteExport.tsv --inserted-sequence-file CosmicInsMutExport.tsv --output-dir var_files --b36-fasta Homo_sapiens.NCBI36.54.dna.toplevel.fa --b37-fasta Homo_sapiens.GRCh37.72.dna.primary_assembly.fa --liftover-chain-hg18-to-hg19 hg18ToHg19.over.chain

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
