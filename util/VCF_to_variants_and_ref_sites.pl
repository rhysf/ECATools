#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use read_FASTA;
use genome_array;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: perl $0 -f <reference FASTA> -v <VCF>\n
Optional: -s\tSettings [1]
            \t\t1 = homozygous only
            \t\t2 = homozygous or heterozygous & ignore indels
            \t\t3 = all variants including indels
          -m\tMin read depth [4]
          -e\tExclude variants on contig/chromosome (y/n) [n]
          -z\tExclude variants not on contig/chromosome (y/n) [n]
          -o\tOutput directory [ECATools_output]\n";
our($opt_e, $opt_f, $opt_m, $opt_o, $opt_s, $opt_v, $opt_z);
getopt('efmosvz');
die $usage unless ($opt_f && $opt_v);
if(!defined $opt_s) { $opt_s = 1; }
if(!defined $opt_m) { $opt_m = 4; }
if(!defined $opt_e) { $opt_e = 'N'; }
if(!defined $opt_z) { $opt_z = 'N'; }
if(!defined $opt_o) { $opt_o = 'ECATools_output'; }

# outfile
if(! -d $opt_o) { `mkdir $opt_o`; }
my($filename, $dirs, $suffix) = fileparse($opt_v);
my $settings = "m-$opt_m-s-$opt_s-e-$opt_e-z-$opt_z";
my $outfile1 = "$opt_o/$filename-$settings-reference-bases.tab";
my $outfile2 = "$opt_o/$filename-$settings-variant-bases.tab";

if((-e $outfile1) || (-e $outfile2)) {
        warn "Skipping as $outfile1 and/or $outfile2 already made\n";
} else {
        # Save sequences
        my $fasta = fastafile::fasta_to_struct($opt_f);

        # Save genome to array
        my $genome_array = genomearray::make_genome_hash_array_from_seq_struct($fasta);

        # Convert ref bases to 1 and variants to 2
        $genome_array = genomearray::fill_genome_hash_array_from_vcf($genome_array, $opt_v, $opt_s, $opt_e, $opt_z, $opt_m);

        # Print tab files for locations of 1s (reference)
        genomearray::print_tab_file_from_regions_in_genome_hash($fasta, $genome_array, 1, $outfile1);

        # Print tab files for locations of 1s (variant)
        genomearray::print_tab_file_from_regions_in_genome_hash($fasta, $genome_array, 2, $outfile2);
}
