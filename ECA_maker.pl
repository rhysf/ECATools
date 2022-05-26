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
my $usage = "Usage: perl $0 -v <Name tab Type(VCF/PILEUP) tab Location>\n
Optional: -s\tSettings [2]
             \t\t1 = homozygous or heterozygous & ignore indels
	     \t\t2 = homozygous only
	     \t\t3 = all variants including indels (partly experimental)
          -m\tMin read depth [4]
          -i\tInclude invariant sites (y/n) [n]
          -g\tGenome length for warnings [] 
	  -e\tExclude variants on contig/chromosome [N]
          -z\tExclude variants not on contig/chromosome [N]
	  -a\tExclude ambigious sites (y/n) [y]\n
Notes: produces the following outfiles:
          ECVA-positions-from-[settings].tab
          ECVA-positions-from-[settings].fasta
          ECVA-positions-from-[settings].pairwise\n";
our($opt_v, $opt_s, $opt_c, $opt_m, $opt_i, $opt_g, $opt_e, $opt_z, $opt_a);
getopt('vscmigeza');
die $usage unless ($opt_v);
if(!defined $opt_s) { $opt_s = 2; }
if(!defined $opt_m) { $opt_m = 4; }
if(!defined $opt_i) { $opt_i = 'n'; }
if(!defined $opt_a) { $opt_a = 'y'; }
if(!defined $opt_e) { $opt_e = 'N'; }
if(!defined $opt_z) { $opt_z = 'N'; }
die "opt_s must be 1, 2 or 3" if(($opt_s ne 1) && ($opt_s ne 2) && ($opt_s ne 3));
die "opt_i must equal n or y: $opt_i\n" if($opt_i !~ m/[ny]/i);

# Output files
my ($filename, $dirs, $suffix) = fileparse($opt_v);
my $settings = ("v-" . $filename . "-m-" . $opt_m . "-s-" . $opt_s . "-i-" . $opt_i);
if($opt_e) { $settings .= ("-e-" . $opt_e); }
if($opt_z) { $settings .= ("-z-" . $opt_z); }
my $Outfile_ECA_positions = 'EVCA-positions-from-' . $settings . '.tab';
my $Outfile_ECA_FASTA = 'EVCA-positions-from-' . $settings . '.fasta';
my $Outfile_ECA_pairwise = 'EVCA-positions-from-' . $settings . '.pairwise';

# Check input files
my $Name_Type_Location_file = tabfile::save_columns_to_column_hash($opt_v, 0, 1, 2);
die "No input specified with $opt_v: $!\n" if(!defined $Name_Type_Location_file);
my $check1 = &Name_Type_Location_files_check_exist($Name_Type_Location_file);

# Read each VCF and save almost all variants (ignore ambiguous and indels, and het when unwanted) and reference if invariants wanted
warn "Get polymorphisms from VCFs...\n";
my ($contig_tab_position_to_tally, $NR_count);
($contig_tab_position_to_tally, $NR_count) = &ECA_save_supercontig_tab_positions_from_VCF($Name_Type_Location_file, $contig_tab_position_to_tally, $NR_count, $opt_s, $opt_e, $opt_z, $opt_m, '1st_pass');

# Get reference bases from sites (that have at least 1 variant)
warn "Get reference bases from sites that have at least 1 variant...\n"; 
($contig_tab_position_to_tally, $NR_count) = &ECA_save_supercontig_tab_positions_from_VCF($Name_Type_Location_file, $contig_tab_position_to_tally, $NR_count, $opt_s, $opt_e, $opt_z, $opt_m, '2nd_pass');

# Count unique loci
if(defined $opt_g) { my $count_at_unique_loci = &tallies_from_VCF($contig_tab_position_to_tally, $NR_count, $opt_g); }

# Find verified positions 
warn "Find ECA positions...\n";
my ($verified_eca_positions) = &find_ECA_positions_from_contig_tab_position_to_tally($contig_tab_position_to_tally, scalar(keys(%{$Name_Type_Location_file})), $opt_e, $opt_z, $opt_i);

# Get all ECA homozygous positions (and possibly heterozygous positions with ambiguity characters) for FASTA and save these for pairwise comparisons
if($opt_s eq 1) { warn "Saving ECA homozygous and heterzygous bases (with ambiguity characters)...\n"; }
elsif($opt_s eq 2) { warn "Saving ECA homozygous bases...\n"; }
my ($verified_eca_positions_per_isolate, $comparisons_snps, $comparisons_hets, $indel_lengths) = &save_ECA_bases($Name_Type_Location_file, $verified_eca_positions);

# Print FASTA file for each file and position in Tab in same order
warn "Printing FASTA and Tab files...\n";
my $isolate_number = 0;
open my $ofh1, '>', $Outfile_ECA_positions or die "Cannot open $Outfile_ECA_positions: $!\n";
open my $ofh2, '>', $Outfile_ECA_FASTA or die "Cannot open $Outfile_ECA_FASTA: $!\n";
foreach my $isolate(keys %{$Name_Type_Location_file}) {
	my $VCF = $$Name_Type_Location_file{$isolate}{'VCF'};
	print $ofh2 ">$isolate\n";
	foreach my $contig(sort keys %{$verified_eca_positions_per_isolate}) {
		POSITIONS: foreach my $position(sort { $a <=> $b } keys %{$$verified_eca_positions_per_isolate{$contig}}) {
			die "VCF $VCF defined, but no position saved at contig) $contig position) $position\n" if(!defined $$verified_eca_positions_per_isolate{$contig}{$position}{$VCF});

			#warn "$contig $position $base\n";

			# Max length for indels
			my $max_length_insertion = $$indel_lengths{$contig}{$position}{'INS'};
			my $max_length_deletion  = $$indel_lengths{$contig}{$position}{'DEL'};

			# Invariant?
			if(($opt_i =~ m/n/) && ($max_length_deletion eq 1)) {
				my %bases_found;
				foreach my $vcfs(keys %{$$verified_eca_positions_per_isolate{$contig}{$position}}) {
					if(!defined $$verified_eca_positions_per_isolate{$contig}{$position}{$vcfs}) {
						warn "ERROR: No invariant base saved for $contig $position $vcfs\n";
					}
					my $base_found = $$verified_eca_positions_per_isolate{$contig}{$position}{$vcfs};
					$bases_found{$base_found}++;
				}
				next POSITIONS if(scalar(keys(%bases_found)) < 2);
			}
			# End of new invariant check

			# Base
			my $base = $$verified_eca_positions_per_isolate{$contig}{$position}{$VCF};

			# Add dashes to base if an insertion found
			if(length($base) < $max_length_insertion) {
				for(my $i = length($base); $i < $max_length_insertion; $i++) {
					$base .= '-';
				}
			}

			# Add deleted bases to base if a deletion found
			if(length($base) < $max_length_deletion) {
				my $extra_base = $$indel_lengths{$contig}{$position}{'DEL_BASES'};
				$base .= $extra_base;
			}

			# Print base location
			if($isolate_number eq 0) { print $ofh1 "$contig\t$position\n"; }
			
			# Print base in FASTA
			print $ofh2 $base;
		}
	}
	print $ofh2 "\n";
	$isolate_number = 1;
}
close $ofh1;
close $ofh2;

# Calculate pairwise comparisons 
warn "Calculate pairwise comparisons...\n";
open my $ofh3, '>', $Outfile_ECA_pairwise or die "Cannot open $Outfile_ECA_pairwise: $!\n";
VCF1: foreach my $isolate(keys %{$Name_Type_Location_file}) {
	my $VCF1 = $$Name_Type_Location_file{$isolate}{'VCF'};
	VCF2: foreach my $isolate(keys %{$Name_Type_Location_file}) {
		my $VCF2 = $$Name_Type_Location_file{$isolate}{'VCF'};
		next VCF2 if ($VCF2 eq $VCF1);
		
		my ($overlap_snps, $overlap_hets) = (0, 0);
		my ($total_snps_vcf1, $total_snps_vcf2) = (0, 0);
		my ($total_hets_vcf1, $total_hets_vcf2) = (0, 0);

		# Overlap and total for snps
		foreach my $contig(keys %{$$comparisons_snps{$VCF1}}) {
			$total_snps_vcf1 += scalar(keys(%{$$comparisons_snps{$VCF1}{$contig}}));
			foreach my $position(keys %{$$comparisons_snps{$VCF1}{$contig}}) {
				my $snp = $$comparisons_snps{$VCF1}{$contig}{$position};
				if((defined $$comparisons_snps{$VCF2}{$contig}{$position}) && ($$comparisons_snps{$VCF2}{$contig}{$position} eq $snp)) { $overlap_snps++; }
			}
		}
		foreach my $contig(keys %{$$comparisons_snps{$VCF2}}) { $total_snps_vcf2 += scalar(keys(%{$$comparisons_snps{$VCF2}{$contig}})); }
		my $percent_hom_shared = 0;
		if(($total_snps_vcf1 + $total_snps_vcf2) ne 0) {
			$percent_hom_shared = (($overlap_snps / (($total_snps_vcf1 + $total_snps_vcf2) / 2)) * 100);
		}
		
		# Overlap and total for hets
		if($opt_s eq 1) {
			foreach my $contig(keys %{$$comparisons_hets{$VCF1}}) {
				$total_hets_vcf1 += scalar(keys(%{$$comparisons_hets{$VCF1}{$contig}}));
				foreach my $position(keys %{$$comparisons_hets{$VCF1}{$contig}}) {
					my $het = $$comparisons_hets{$VCF1}{$contig}{$position};
					if((defined $$comparisons_hets{$VCF2}{$contig}{$position}) && ($$comparisons_hets{$VCF2}{$contig}{$position} eq $het)) { $overlap_hets++; }
				}
			}
			foreach my $contig(keys %{$$comparisons_hets{$VCF2}}) { $total_hets_vcf2 += scalar(keys(%{$$comparisons_hets{$VCF2}{$contig}})); }
		
			my $percent_het_shared = 0;
			if(($total_hets_vcf1 + $total_hets_vcf2) ne 0) {
				$percent_het_shared = (($overlap_hets / (($total_hets_vcf1 + $total_hets_vcf2) / 2)) * 100);
			}
			print $ofh3 "$VCF1 - $VCF2 = $overlap_snps shared snps ($percent_hom_shared\%) and $overlap_hets shared hets ($percent_het_shared\%)\n";
		}
		if($opt_s eq 2) { print $ofh3 "$VCF1 - $VCF2 = $overlap_snps shared snps ($percent_hom_shared\%)\n"; }	
	}
	print $ofh3 "\n";
}
close $ofh3;

# Convert FASTA to Nexus
my $fasta_to_nexus_CMD = "perl $Bin/FASTA-parser.pl -s $Outfile_ECA_FASTA -p nexus -d dna -g n > $Outfile_ECA_FASTA.nex";
system($fasta_to_nexus_CMD);

sub tallies_from_VCF {
	my ($all, $nr, $genome_length) = @_;
	warn "tallies_from_VCF...\n";
	my ($sum, $count_at_unique_loci, $count_snp_at_unique_loci, $count_ref_at_unique_loci, $count_het_at_unique_loci) = (0, 0, 0, 0, 0);
	$sum += $$nr{'all'};
	$count_at_unique_loci += scalar(keys(%{$all}));
	$count_snp_at_unique_loci += scalar(keys(%{$$nr{'snp'}}));
	$count_ref_at_unique_loci += scalar(keys(%{$$nr{'reference'}}));
	$count_het_at_unique_loci += scalar(keys(%{$$nr{'heterozygous'}}));
		
	my $snp_per_kb = (($count_snp_at_unique_loci / $genome_length) * 1000);
	my $ref_per_kb = (($count_ref_at_unique_loci / $genome_length) * 1000);
	my $het_per_kb = (($count_het_at_unique_loci / $genome_length) * 1000);
	
	warn "$sum total positions found in VCFs\n";
	warn "$count_at_unique_loci / $sum at unique loci (non-redundant)\n";
	warn "Unique loci (non-redundant; NR), and per Kb using genome length $genome_length nt:\n";
	warn "$count_snp_at_unique_loci NR SNPs = $snp_per_kb/Kb\n";
	warn "$count_ref_at_unique_loci NR REFs = $ref_per_kb/Kb\n";
	warn "$count_het_at_unique_loci NR HETs = $het_per_kb/Kb\n";
	return $count_at_unique_loci;
}

sub Name_Type_Location_files_check_exist {
	my $Name_Type_Location_file = $_[0];
	warn "Name_Type_Location_files_check_exist...\n";
	foreach my $isolate(keys %{$Name_Type_Location_file}) {
		foreach my $file_type(keys %{$$Name_Type_Location_file{$isolate}}) {
			my $file = $$Name_Type_Location_file{$isolate}{$file_type};
			warn "Checking $file...\n";
			open my $fh, '<', $file or die "Cannot open $file. Check Name_Type_Location file: $!\n";
			close $fh;
		}
	}
	warn "Name_Type_Location_files_check_exist: Yes\n";
	return 1;
}

sub ECA_save_supercontig_tab_positions_from_VCF {
	my ($Name_Type_Location_file, $contig_tab_position_to_tally, $NR_count, $settings, $exclude, $include, $min_depth, $pass) = @_;

	# Save positions
	warn "ECA_save_supercontig_tab_positions_from_VCF (setting=$settings, exclude=$exclude, include=$include, min depth=$min_depth, pass=$pass)...\n";
	foreach my $isolate(keys %{$Name_Type_Location_file}) {
		my $VCF = $$Name_Type_Location_file{$isolate}{'VCF'};
		my $depth_flag_count = 0;
		warn "Opening and saving from $VCF...\n";
		open my $fh, '<', $VCF or die "Cannot open $VCF: $!\n";
		VCF: while(my $line=<$fh>) {
			chomp $line;
			my $VCF_line = vcflines::read_VCF_lines($line);
			my $supercontig = $$VCF_line{'supercontig'};
			my $position = $$VCF_line{'position'};
			my $base_type = $$VCF_line{'base_type0'};

			# Ignore ambigious, or contigs of uninterest
			next VCF if($$VCF_line{'next'} eq 1);
			next VCF if(($exclude ne 'N') && ($supercontig eq $exclude));
			next VCF if(($include ne 'N') && ($supercontig ne $include));
			my $sc_and_position = ($supercontig . "\t" . $position);

			# Ignore heterozygous sites if setting 2 (only homozygous)
			# Only keep HET/SNP on 1st pass
			# Only keep REF on sites previously found on 2nd pass
			next VCF if(($settings eq 2) && ($base_type eq 'heterozygous'));
			next VCF if(($settings ne 3) && ($base_type =~ m/insertion|deletion/));
			if($pass eq '1st_pass') {
				next VCF if($base_type !~ m/heterozygous|snp|insertion|deletion/);
				#next VCF if($base_type !~ m/heterozygous|snp/);
			} else {
				next VCF if($base_type !~ m/reference/);
				# Ignore position not previously found
				next VCF if(!defined $$contig_tab_position_to_tally{$sc_and_position}); 
			}

			# Min depth filtering
			if((defined $$VCF_line{'DP0'}) && ($$VCF_line{'DP0'} ne '?'))  {
				next if($$VCF_line{'DP0'} < $min_depth);
			}
			# Alert this wont be done max once per VCF
			else {
				if($depth_flag_count eq 0) {
					$depth_flag_count = 1;
					warn "Depth flag not identified in VCF. Ignoring -m\n";
				}
			}

			# Save positions in hashes
			$$contig_tab_position_to_tally{$sc_and_position}++;
			$$NR_count{$base_type}{$sc_and_position}++;
			$$NR_count{'all'}++;
		}
		close $fh;
	}
	return ($contig_tab_position_to_tally, $NR_count);
}

sub find_ECA_positions_from_contig_tab_position_to_tally {
	my ($contig_tab_position_to_tally, $num_of_positions_found, $exclude, $include, $invariants) = @_;
	my $count = 0;
	my (%verified_eca_positions);
	foreach my $lines(sort keys %{$contig_tab_position_to_tally}) {
		my $count_found = $$contig_tab_position_to_tally{$lines};
		if($count_found eq $num_of_positions_found) {
			my @line_parts = split /\t/, $lines;
			my ($contig, $pos) = @line_parts;
			$verified_eca_positions{$contig}{$pos} = 1;
			$count++;
		}
	}
	warn "$count entirely covered in all (ECA) (inc. invariants = $invariants)\n";
	if($exclude ne 'N') { warn "(over all but $exclude)\n"; }
	if($include ne 'N') { warn "(over only $include)\n"; }
	return \%verified_eca_positions;
}

sub save_ECA_bases {
	my ($Name_Type_Location_file, $verified_eca_positions) = @_;
	my (%comparisons_snps, %comparisons_hets, %verified_eca_positions_per_isolate);
	my %indel_lengths;
	foreach my $isolate(keys %{$Name_Type_Location_file}) {
		my $VCF = $$Name_Type_Location_file{$isolate}{'VCF'};
		warn "Saving nucleotides from $VCF...\n";
		open my $fh, '<', $VCF or die "Cannot open $VCF: $!\n";
		VCF: while(my $line=<$fh>) {
			chomp $line;
			my $VCF_line = vcflines::read_VCF_lines($line);
			my $supercontig = $$VCF_line{'supercontig'};
			my $position = $$VCF_line{'position'};
			my $ref_base = $$VCF_line{'reference_VCF_format'};
			my $consensus = $$VCF_line{'0base1'};
			my $base_type = $$VCF_line{'base_type0'};
			my $amb_char = $$VCF_line{'amb_char0'};

			# ignore if header
			next VCF if($$VCF_line{'next'} eq 1);

			# ignore if not an ECA position
			next VCF if(!defined $$verified_eca_positions{$supercontig}{$position});
			
			# length for indels
			if(!defined $indel_lengths{$supercontig}{$position}{'INS'}) { $indel_lengths{$supercontig}{$position}{'INS'} = 1; }
			if(!defined $indel_lengths{$supercontig}{$position}{'DEL'}) { $indel_lengths{$supercontig}{$position}{'DEL'} = 1; }

			# Save nucleotides
			if(($opt_s eq 1) && ($base_type eq 'heterozygous')) { 
				$verified_eca_positions_per_isolate{$supercontig}{$position}{$VCF} = $amb_char; 
				$comparisons_hets{$VCF}{$supercontig}{$position} = $amb_char;

			}
			elsif($base_type =~ m/snp|reference/) { 
				$consensus = substr $consensus, 0, 1;
				die "ERROR: unexpected base found in $VCF $supercontig $position: $consensus\n" if($consensus !~ m/[ATCG]/i);
				$verified_eca_positions_per_isolate{$supercontig}{$position}{$VCF} = $consensus; 
				$comparisons_snps{$VCF}{$supercontig}{$position} = $consensus;
			}
			elsif($base_type =~ m/insertion/) {

				$verified_eca_positions_per_isolate{$supercontig}{$position}{$VCF} = $consensus; 
				$comparisons_snps{$VCF}{$supercontig}{$position} = $consensus;

				# Longest insertion
				my $length = length($consensus);
				if($length > $indel_lengths{$supercontig}{$position}{'INS'}) { $indel_lengths{$supercontig}{$position}{'INS'} = $length; }
			}
			elsif($base_type =~ m/deletion/) {
				
				# how many dashes to add?
				my $dashes = '';
				for(my $i=0; $i < (length($ref_base) - length($consensus)); $i++) { $dashes .= '-'; }

				# save
				$verified_eca_positions_per_isolate{$supercontig}{$position}{$VCF} = ($consensus . $dashes); 
				$comparisons_snps{$VCF}{$supercontig}{$position} = ($consensus . $dashes);

				# Longest deletion
				my $length = (length($ref_base) - length($consensus)) + 1;
				if($length > $indel_lengths{$supercontig}{$position}{'DEL'}) { 
					$indel_lengths{$supercontig}{$position}{'DEL'} = $length; 
					my $del_base = substr $ref_base, length($consensus), length($ref_base);
					$indel_lengths{$supercontig}{$position}{'DEL_BASES'} = $del_base; 
				}
			}
			else {
				warn "WARNING: base found unaccounted for in $VCF $supercontig $position: $base_type\n";
			}
		}
		close $fh;
	}
	return (\%verified_eca_positions_per_isolate, \%comparisons_snps, \%comparisons_hets, \%indel_lengths);
}
