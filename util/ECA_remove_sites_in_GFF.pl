#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib "$Bin/modules";
use read_Tab;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: perl $0 <EVCA.tab (supercontig tab position)> <Combined Gubbins recombination prediction.GFF (or any!))> > output\n";
die $usage if (@ARGV !=2);
my ($eca, $gff) = @ARGV;

# Save ECA list
my $eca_positions = tabfile::save_columns_to_one_hash($eca, 0, 1);

# go through GFF and delete positions from eca_positions
open my $fh, '<', $gff or die "Cannot open $gff : $!\n";
while(my $line=<$fh>) {
	chomp $line;
	my @bits = split /\t/, $line;
	my ($contig, $start, $stop, $info) = ($bits[0], $bits[3], $bits[4], $bits[8]);
	for(my $i=$start; $i<$stop; $i++) {
		if(defined $$eca_positions{$contig}{$i}) { delete $$eca_positions{$contig}{$i}; }
	}
}
close $fh;

# print
foreach my $contig(keys %{$eca_positions}) {
	foreach my $pos(keys %{$$eca_positions{$contig}}) {
		print "$contig\t$pos\n";
	}
}
