#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/modules";
use read_FASTA;
use read_Tab;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: perl $0 -f <Genome.FASTA> -a <ECA.tab> > output\n
Note: Expects output=opt_f-EVCA-positions.fasta\n";
our($opt_f, $opt_a);
getopt('fa');
die $usage unless (($opt_f) && ($opt_a));

# Input FASTA
my $fasta = fastafile::fasta_to_struct($opt_f);

# Save hash of arrays
warn "Save hash of arrays (contig -> nt array)...\n";
my %hash_of_chr_arrays;
foreach my $contig(keys %{$$fasta{'seq'}}) {
	my $seq = $$fasta{'seq'}{$contig};
	my @nt = split //, $seq;
	$hash_of_chr_arrays{$contig} = \@nt;
}

# We want to keep order of file (so cant save the file to memory and sort)
# Also cant assume that all of a single contig will be together in the EVCA file
warn "Going through each position in $opt_a (takes a long time)...\n";
my $extracted_seq = '';
open my $fh, '<', $opt_a or die "Cannot open $opt_a : $!\n";
TAB: while(my $line=<$fh>) {
	chomp $line;
	my @bits = split /\t/, $line;
	my ($contig, $pos) = @bits;
	warn "$contig $pos...\n";
	die "No $contig found in fasta file\n" if(!defined $hash_of_chr_arrays{$contig});
	my $nt_of_interest = $hash_of_chr_arrays{$contig}[($pos - 1)];
	$extracted_seq .= $nt_of_interest;
}
print ">$opt_f\n$extracted_seq\n";
