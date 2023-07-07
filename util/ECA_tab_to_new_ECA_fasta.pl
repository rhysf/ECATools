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
my $usage = "Usage: perl $0 -v <Name tab Type(VCF/PILEUP) tab Location> -a <ECA.tab (contig tab position)>\n
Optional: -s\tSettings (1 = homozygous or heterozygous & ignore indels, 2 = homozygous only) [2]
          -i\tInclude invariant sites (y/n) [y]
          -m\tMin read depth [4]
          -g\tGenome length for warnings, currently set at gscidC R265 post-Pilon length exc. Ns [17260762] 
	  -e\tExclude variants on contig/chromosome []
          -z\tExclude variants not on contig/chromosome []\n
Output:   -o\tOutput [./ECVA-positions-from-EVCA-file-opt_a-[settings].fasta]\n
Warning: Currently makes the order of nucleotides arbitrary - so - this could be a problem if you want to find out genomic loci specific info (e.g. recombination)\n";
our($opt_a, $opt_c, $opt_e, $opt_g, $opt_i, $opt_m, $opt_n, $opt_o, $opt_s, $opt_v, $opt_z);
getopt('acegimnosvz');
die $usage unless (($opt_v) && ($opt_a));
if(!defined $opt_s) { $opt_s = 2; }
die "opt_s must be 1 or 2" if (($opt_s ne 1) && ($opt_s ne 2));
if(!defined $opt_n) { $opt_n = 3; }
if(!defined $opt_m) { $opt_m = 4; }
if(!defined $opt_g) { $opt_g = 17260762; }
if(!defined $opt_i) { $opt_i = 'y'; }

# Output files
my ($filename, $dirs, $suffix) = fileparse($opt_v);
my $settings = ("v-" . $filename . "-m-" . $opt_m . "-s-" . $opt_s . "-i-" . $opt_i);
if($opt_e) { $settings .= ("-e-" . $opt_e); }
if($opt_z) { $settings .= ("-z-" . $opt_z); }
if(!defined $opt_o) { $opt_o = 'EVCA-positions-from-EVCA-file-opt_a-' . $settings . '.fasta'; }

# Input file 1 (Name tab Type tab Location)
my $Name_Type_Location_file = tabfile::save_columns_to_column_hash($opt_v, 0, 1, 2);
die "No input specified with $opt_v: $!\n" if(!defined $Name_Type_Location_file);

# Input file 2 (Save ECA positions: contig tab position = 1)
my $verified_eca_positions = &Tab_to_ECA_positions($opt_a, $opt_e, $opt_z);

# Read each VCF and save all positions in ECA positions file (ignore ambiguous and indels) 
my ($found, $NR_count) = &save_supercontig_tab_positions_from_VCF($Name_Type_Location_file, $opt_s, $opt_e, $opt_z, $verified_eca_positions);

# Get all ECA homozygous (and possibly heterozygous positions with ambiguity characters) for FASTA and save these for pairwise comparisons
if($opt_s eq 1) { warn "Getting all EVCA homozygous and heterzygous positions with ambiguity characters...\n"; }
elsif($opt_s eq 2) { warn "Getting all EVCA homozygous snps...\n"; }
my ($verified_eca_positions_per_isolate, $comparisons_snps, $comparisons_hets) = &save_ECA_bases($Name_Type_Location_file, $verified_eca_positions);

# Print FASTA file for each file
warn "Printing FASTA for the ECA positions in each of the files...\n";
open my $ofh, '>', $opt_o or die "Cannot open $opt_o: $!\n";
foreach my $isolate(keys %{$Name_Type_Location_file}) {
	my $VCF = $$Name_Type_Location_file{$isolate}{'VCF'};
	my $N_count = 0;
	print $ofh ">$isolate\n";
	foreach my $contig(keys %{$verified_eca_positions_per_isolate}) {
		foreach my $position(keys %{$$verified_eca_positions_per_isolate{$contig}}) {
			die "No VCF $VCF for $contig $position in isolate $isolate\n" if(!defined $VCF);	
			if(defined $$verified_eca_positions_per_isolate{$contig}{$position}{$VCF}) {
				print $ofh $$verified_eca_positions_per_isolate{$contig}{$position}{$VCF};
			} else { 
				#warn "VCF $VCF defined, but no position saved at $contig $position - adding N\n"; 
				$N_count++;
				print $ofh "N";
			}
		}
	}
	print $ofh "\n";
	warn "Number of N's inserted into $isolate ECA sequence = $N_count\n";
}
close $ofh;

sub save_supercontig_tab_positions_from_VCF {
	my ($Name_Type_Location_file, $settings, $exclude, $include, $verified_positions_already) = @_;
	my (%found, %NR_count);
	warn "save_supercontig_tab_positions_from_VCF: $settings\n";
	if($exclude) { warn "	-e) $exclude\n"; }
	if($include) { warn "	-z) $include\n"; }
	foreach my $isolate(keys %{$Name_Type_Location_file}) {
		my $VCF = $$Name_Type_Location_file{$isolate}{'VCF'};
		warn "Opening and saving from $VCF...\n";
		open my $fh, '<', $VCF or die "Cannot open $VCF: $!\n";
		VCF: while(my $line=<$fh>) {
			chomp $line;
			my $VCF_line = vcflines::read_VCF_lines($line);
			next VCF if($$VCF_line{'next'} eq 1);
			next VCF if(($exclude) && ($$VCF_line{'supercontig'} ne $exclude));
			next VCF if(($include) && ($$VCF_line{'supercontig'} eq $include));
			next VCF if(!defined $$verified_positions_already{$$VCF_line{'supercontig'}}{$$VCF_line{'position'}});

			my $POSITION = ($$VCF_line{'supercontig'} . "\t" . $$VCF_line{'position'});
			$found{$POSITION}++;
			$NR_count{$$VCF_line{'base_type0'}}{$POSITION}++;
			$NR_count{'all'}++;
		}
	}
	return (\%found, \%NR_count);
}

sub Tab_to_ECA_positions {
	my ($file, $exclude, $include) = @_;
	my %verified_eca_positions;
	my $count = 0;
	warn "Saving positions found in all files from $opt_a...\n";
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
	warn "$count entirely covered and verified in all\n";
	if($exclude) { warn "(over all but $exclude)\n"; }
	if($include) { warn "(over only $include)\n"; }
	return \%verified_eca_positions;
}

sub save_ECA_bases {
	my ($Name_Type_Location_file, $verified_eca_positions) = @_;
	my (%comparisons_snps, %comparisons_hets, %verified_eca_positions_per_isolate);
	foreach my $isolate(keys %{$Name_Type_Location_file}) {
		my $VCF = $$Name_Type_Location_file{$isolate}{'VCF'};
		warn "Saving nucleotides from $VCF...\n";
		open my $fh, '<', $VCF or die "Cannot open $VCF: $!\n";
		VCF: while(my $line=<$fh>) {
			chomp $line;
			my $VCF_line = vcflines::read_VCF_lines($line);
			my ($supercontig, $position, $consensus, $base_type, $amb_char) = ($$VCF_line{'supercontig'}, $$VCF_line{'position'}, $$VCF_line{'0base1'}, $$VCF_line{'base_type0'}, $$VCF_line{'amb_char0'});
			next VCF if($$VCF_line{'next'} eq 1);

			# ignore if not an ECA position
			next VCF if(!defined $$verified_eca_positions{$supercontig}{$position});

			# Save nucleotides
			if(($opt_s eq 1) && ($base_type eq 'heterozygous')) { 
				$verified_eca_positions_per_isolate{$supercontig}{$position}{$VCF} = $amb_char; 
				$comparisons_hets{$VCF}{$supercontig}{$position} = $amb_char;
			}
			if($base_type =~ m/snp|reference/) { 	
				die "What is this $consensus in $VCF $supercontig $position" if($consensus !~ m/A|T|C|G/i);
				$verified_eca_positions_per_isolate{$supercontig}{$position}{$VCF} = $consensus; 
				$comparisons_snps{$VCF}{$supercontig}{$position} = $consensus;
			}
		}
		close $fh;
	}
	return (\%verified_eca_positions_per_isolate, \%comparisons_snps, \%comparisons_hets);
}
