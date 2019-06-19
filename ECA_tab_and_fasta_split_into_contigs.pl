#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/modules";
use read_FASTA;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: perl $0 -f <ECA.FASTA> -e <ECA.tab (contig tab position)>\n
Notes: produces the following outfiles:
          [opt_f]-[contigs].fasta\n";
our($opt_e, $opt_f);
getopt('ef');
die $usage unless ($opt_e && $opt_f);

# Save sequences
my $fasta = fastafile::fasta_to_struct($opt_f);

# Save tabs (order is important)
#my $eca_positions = tabfile::save_columns_to_one_hash($opt_e, 0, 1);
my $count_of_sites = 0;
my $length_of_contig = 0;
my $new_contig;
open my $fh, '<', $opt_e or die "Cannot open $opt_e : $!\n";
while(my $line=<$fh>) {
	chomp $line;
	my @bits = split /\t/, $line;
	my ($contig, $position) = @bits;

	# New contig?
	if((!defined $new_contig) || ($new_contig ne $contig)) {
		
		# Init?
		if($length_of_contig ne 0) { 
			# extract
			my $start_of_contig = ($count_of_sites - $length_of_contig);
			warn "Contig $new_contig $start_of_contig - $count_of_sites (length = $length_of_contig)\n";	
			my $output = ($opt_f . "-" . $new_contig . ".fasta");
			open my $ofh, '>', $output or die "Cannot open $output : $!\n";
			foreach my $isolates(keys %{$$fasta{'seq'}}) {
				my $seq = substr $$fasta{'seq'}{$isolates}, $start_of_contig, ($length_of_contig);
				print $ofh ">$isolates\n$seq\n";
			}
			close $ofh;
		}
		# Init
		$length_of_contig = 0;
		$new_contig = $contig;
	}
	$count_of_sites++;
	$length_of_contig++;
}

# extract (last contig)
my $start_of_contig = ($count_of_sites - $length_of_contig);
warn "Contig $new_contig $start_of_contig - $count_of_sites (length = $length_of_contig)\n";	
my $output = ($opt_f . "-" . $new_contig . ".fasta");
open my $ofh, '>', $output or die "Cannot open $output : $!\n";
foreach my $isolates(keys %{$$fasta{'seq'}}) {
	my $seq = substr $$fasta{'seq'}{$isolates}, $start_of_contig, ($length_of_contig);
	print $ofh ">$isolates\n$seq\n";
}
close $ofh;
