#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/modules";
use read_VCF_lines;
use read_Tab;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: perl $0 -v <Name tab Type(VCF/PILEUP) tab Location> -a <ECA file>\n
Optional: -s\tSettings (1 = homozygous or heterozygous & ignore indels, 2 = homozygous only) [2]
          -i\tInclude headers (y/n) [y]
          -g\tGenome length for warnings, currently set at gscidC R265 post-Pilon length exc. Ns [17260762] 
	  -e\tExclude variants on contig/chromosome []
          -z\tExclude variants not on contig/chromosome []\n
Outputs:  Produces the following outfiles:
          ECVA-positions-from-[settings].tab
          ECVA-positions-from-[settings].fasta
          ECVA-positions-from-[settings].pairwise\n";
our($opt_v, $opt_a, $opt_s, $opt_i, $opt_n, $opt_c, $opt_g, $opt_e, $opt_z);
getopt('vasincmigez');
die $usage unless (($opt_v) && ($opt_a));
if(!defined $opt_s) { $opt_s = 2; }
die "opt_s must be 1 or 2" if (($opt_s ne 1) && ($opt_s ne 2));
if(!defined $opt_n) { $opt_n = 3; }
if(!defined $opt_i) { $opt_i = 'y'; }

# Output file suffix (onto VCFs)
my ($filename, $dirs, $suffix) = fileparse($opt_v);
my $settings = ("v-" . $filename . "-s-" . $opt_s);
if($opt_e) { $settings .= ("-e-" . $opt_e); }
if($opt_z) { $settings .= ("-z-" . $opt_z); }
$settings .= ".vcf";

# Input files
my $Name_Type_Location_file = tabfile::save_columns_to_column_hash($opt_v, 0, 1, 2);
die "No input specified with $opt_v: $!\n" if(!defined $Name_Type_Location_file);

# Find verified positions 
warn "Find ECA positions...\n";
my ($verified_eca_positions) = &save_positions_from_ECA_file($opt_a, $opt_e, $opt_z);

# Print positions from VCFs 
my (%found, %NR_count);
foreach my $isolate(keys %{$Name_Type_Location_file}) {
	my $VCF = $$Name_Type_Location_file{$isolate}{'VCF'};
	my $outfile = ($VCF . $settings);
	open my $fh, '<', $VCF or die "Cannot open $VCF: $!\n";
	open my $ofh, '>', $outfile or die "Cannot open $outfile: $!\n";
	warn "Reading $VCF, printing ECA positions to $outfile...\n";
	VCF: while(my $line=<$fh>) {
		chomp $line;
		my $VCF_line = vcflines::read_VCF_lines($line);
		if(($line =~ m/^\#/) && ($opt_i eq 'y')) { print $ofh "$line\n"; }
		next VCF if($$VCF_line{'next'} eq 1);
		next VCF if(!defined $$verified_eca_positions{$$VCF_line{'supercontig'}}{$$VCF_line{'position'}});
		next VCF if(($opt_s eq 2) && ($$VCF_line{'base_type0'} eq 'heterozygous'));
		next VCF if($$VCF_line{'base_type0'} !~ m/heterozygous|snp|reference/);

		# Print VCF line
		print $ofh "$line\n";
	}
	close $fh;
	close $ofh;
}

sub save_positions_from_ECA_file {
	my ($file, $exclude, $include) = @_;
	
	my (%verified_eca_positions);
	my $count = 0;
	warn "Saving positions found in $file...\n";
	open my $fh, '<', $file or die "Cannot open $file: $!\n";
	ECA: while(my $line = <$fh>) {
	   	chomp $line;
	       	my @bits = split /\t/, $line;
		my ($contig, $pos) = @bits;
		next ECA if(($exclude) && ($contig ne $exclude));
		next ECA if(($include) && ($contig eq $include));
		$verified_eca_positions{$contig}{$pos} = 1;
		$count++;
	}
	close $fh;
	warn "$count entirely covered in all from file\n";
	if($exclude) { warn "(over all but $exclude)\n"; }
	if($include) { warn "(over only $include)\n"; }
	return \%verified_eca_positions;
}
