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
my $usage = "Usage: perl $0 -f <reference FASTA> -v <name_VCF_location file>\n
Optional: -r\tSteps to run [123]
            \t\t1 = Print VCF regions, according to settings
            \t\t2 = Find and print ECA positions
            \t\t3 = Print ECA FASTA and pairwise comparisons
          -s\tSettings [1]
            \t\t1 = homozygous only
            \t\t2 = homozygous or heterozygous & ignore indels
            \t\t3 = all variants (including indels)
          -m\tMin read depth [4]
          -i\tInclude invariant sites (y/n) [n]
          -e\tExclude variants on contig/chromosome (y/n) [n]
          -z\tExclude variants not on contig/chromosome (y/n) [n]
          -a\tExclude ambigious sites (y/n) [y]\n
          -o\tOutput directory [ECATools_output]\n
Parallel: -g\tRun commands on the grid (y/n) [n]
          -p\tPlatform (UGER, LSF, GridEngine) [UGER]
          -q\tQueue name [short]\n
Notes:    If running in parallel, need to complete steps 1 separately\n";
our($opt_a, $opt_c, $opt_e, $opt_f, $opt_g, $opt_i, $opt_m, $opt_o, $opt_p, $opt_q, $opt_r, $opt_s, $opt_v, $opt_z);
getopt('acefgimopqrsvz');
die $usage unless ($opt_f && $opt_v);
if(!defined $opt_a) { $opt_a = 'y'; }
if(!defined $opt_e) { $opt_e = 'n'; }
if(!defined $opt_i) { $opt_i = 'n'; }
if(!defined $opt_m) { $opt_m = 4; }
if(!defined $opt_o) { $opt_o = 'ECATools_output'; }
if(!defined $opt_r) { $opt_r = '123'; }
if(!defined $opt_s) { $opt_s = 1; }
if(!defined $opt_z) { $opt_z = 'n'; }
if(!defined $opt_g) { $opt_g = 'n'; }
if(!defined $opt_p) { $opt_p = 'UGER'; }
if(!defined $opt_q) { $opt_q = 'short'; }
die "opt_s must be 1, 2 or 3" if($opt_s !~ m/[123]/);
die "opt_i must equal n or y: $opt_i\n" if($opt_i !~ m/[ny]/i);
warn "$0: -f $opt_f -v $opt_v -r $opt_r -s $opt_s -m $opt_m -i $opt_i -e $opt_e -z $opt_z -a $opt_a -o $opt_o -g $opt_g -p $opt_p -q $opt_q\n";

# Output files
my ($filename, $dirs, $suffix) = fileparse($opt_v);
my $outfile_settings = ("v-" . $filename . "-m-" . $opt_m . "-s-" . $opt_s . "-i-" . $opt_i);
if($opt_e) { $outfile_settings .= ("-e-" . $opt_e); }
if($opt_z) { $outfile_settings .= ("-z-" . $opt_z); }
if(! -d $opt_o) { `mkdir $opt_o`; }
my $outfile_ECA_positions = $opt_o . '/ECA-positions-from-' . $outfile_settings . '.tab';
my $outfile_ECA_FASTA = $opt_o . '/ECA-positions-from-' . $outfile_settings . '.fasta';
my $outfile_ECA_pairwise = $opt_o . '/ECA-positions-from-' . $outfile_settings . '.pairwise';

# Dependencies
my $vcf_to_variants_and_ref_sites = "$Bin/util/VCF_to_variants_and_ref_sites.pl";
my $Run_Commands_python = "$Bin/util/Run_cmds_on_grid.py";
foreach($vcf_to_variants_and_ref_sites, $Run_Commands_python) { die "Error: Dependency script $_ not found" if(! -e $_); }

# Check input files
my $Name_Type_Location_hash = tabfile::save_columns_to_column_hash($opt_v, 0, 1, 2);
die "No input specified with $opt_v: $!\n" if(!defined $Name_Type_Location_hash);
my $check1 = &check_files_exist_in_Name_Type_Location_files($Name_Type_Location_hash);
my $number_of_VCFs = scalar(keys(%{$Name_Type_Location_hash}));

# Print VCF regions, according to settings
if($opt_r =~ m/1/) {
	my $settings = "-f $opt_f -s $opt_s -m $opt_m -e $opt_e -z $opt_z -o $opt_o";
	my $outfile_parallel_commands = "$opt_o/vcf_to_variants_and_ref_sites.cmds";
	&print_supercontig_tab_positions_from_VCF($Name_Type_Location_hash, $settings, $outfile_parallel_commands, $opt_p, $opt_q);
	die "ECATools: Wait for step 1 to complete before running steps 2 and 3\n" if($opt_g eq 'y');
}

# Find and print ECA sites
if($opt_r =~ m/2/) {

	# Save all variant sites
	my $settings = "m-$opt_m-s-$opt_s-e-$opt_e-z-$opt_z";
	my $variants_struct = &save_contig_position_to_number_VCFs_found("$opt_o/*$settings-variant-bases.tab");

	# Save all reference sites
	my $reference_struct = &save_contig_position_to_number_VCFs_found("$opt_o/*$settings-reference-bases.tab");

	# Find ECA positions
	my $verified_eca_positions = &find_and_print_ECA_positions($variants_struct, $reference_struct, $number_of_VCFs, $opt_i, $outfile_ECA_positions);
}

# Print ECA sites
if($opt_r =~ m/3/) {

	# Save ECA positions for FASTA and pairwise comparisons
	my $verified_eca_positions = tabfile::save_columns_to_one_hash($outfile_ECA_positions, 0, 1);
	my ($verified_eca_positions_per_isolate, $comparisons_snps, $comparisons_hets, $indel_lengths) = &save_ECA_bases($Name_Type_Location_hash, $verified_eca_positions, $opt_s);

	# Print FASTA file for each file and position in Tab in same order
	&print_fasta_for_ECA_sites($Name_Type_Location_hash, $verified_eca_positions_per_isolate, $indel_lengths, $opt_i, $outfile_ECA_FASTA);

	# Calculate pairwise comparisons
	if($opt_i eq 'y') { 
		&calculate_pairwise_differences($Name_Type_Location_hash, $comparisons_snps, 'SNPs', $outfile_ECA_pairwise);
		
		# heterozygous
		if($opt_s =~ m/23/) {
			&calculate_pairwise_differences($Name_Type_Location_hash, $comparisons_hets, 'heterozygous_positions', $outfile_ECA_pairwise);
		}
	}
}

### sub routines

sub print_fasta_for_ECA_sites {
	my ($Name_Type_Location_hash, $verified_eca_positions_per_isolate, $indel_lengths, $include_invariant_sites, $outfile_ECA_FASTA) = @_;

	# Outfile
	warn "print_fasta_for_ECA_sites: $outfile_ECA_FASTA\n";
	open my $ofh, '>', $outfile_ECA_FASTA or die "Cannot open $outfile_ECA_FASTA: $!\n";

	# Go through each isolate to print (sort keys for contig and sort { $a <=> $b } keys must be consistent with find_and_print_ECA_positions)
	foreach my $isolate(keys %{$Name_Type_Location_hash}) {
		my $VCF = $$Name_Type_Location_hash{$isolate}{'VCF'};

		# header
		print $ofh ">$isolate\n";

		# save sequence
		my $sequence = '';
		foreach my $contig(sort keys %{$verified_eca_positions_per_isolate}) {
			POSITIONS: foreach my $position(sort { $a <=> $b } keys %{$$verified_eca_positions_per_isolate{$contig}}) {
				die "Error: No position saved for $contig pos. $position in $VCF\n" if(!defined $$verified_eca_positions_per_isolate{$contig}{$position}{$VCF});

				# Base
				my $base = $$verified_eca_positions_per_isolate{$contig}{$position}{$VCF};

				# Add dashes to base if an insertion found
				if(defined $$indel_lengths{$contig}{$position}{'INS'}) {
					my $max_length_insertion = $$indel_lengths{$contig}{$position}{'INS'};

					next if(length($base) >= $max_length_insertion);

					for(my $i = length($base); $i < $max_length_insertion; $i++) {
						$base .= '-';
					}
				}

				# Add deleted bases to base if a deletion found
				if(defined $$indel_lengths{$contig}{$position}{'DEL'}) {
					my $max_length_deletion  = $$indel_lengths{$contig}{$position}{'DEL'};
					if(length($base) < $max_length_deletion) {
						my $extra_base = $$indel_lengths{$contig}{$position}{'DEL_BASES'};
						$base .= $extra_base;
					}
				}
				
				# save sequence
				$sequence .= $base;
			}
		}

		# print sequence
		$sequence =~ s/(\S{60})/$1\n/g;
		print $ofh "$sequence\n";
	}
	close $ofh;
	return 1;
}

sub calculate_pairwise_differences {
	my ($Name_Type_Location_hash, $comparisons_variants, $comparison_type, $outfile) = @_;

	warn "calculate_pairwise_differences: $comparison_type to $outfile\n";

	# outfile
	open my $ofh, '>>', $outfile or die "Cannot open $outfile: $!\n";
	print $ofh "file1\tfile2\tshared $comparison_type\tshared $comparison_type (\%)\n";

	# go through pairwise VCF1 and VCF2 and find overlap
	VCF1: foreach my $isolate(keys %{$Name_Type_Location_hash}) {
		my $VCF1 = $$Name_Type_Location_hash{$isolate}{'VCF'};
		VCF2: foreach my $isolate(keys %{$Name_Type_Location_hash}) {
			my $VCF2 = $$Name_Type_Location_hash{$isolate}{'VCF'};
			next VCF2 if ($VCF2 eq $VCF1);
			
			my ($overlap_snps) = 0;
			my ($total_snps_vcf1, $total_snps_vcf2) = (0, 0);

			# Overlap and total
			foreach my $contig(keys %{$$comparisons_variants{$VCF1}}) {
				$total_snps_vcf1 += scalar(keys(%{$$comparisons_variants{$VCF1}{$contig}}));
				POS: foreach my $position(keys %{$$comparisons_variants{$VCF1}{$contig}}) {
					my $vcf1_snp = $$comparisons_variants{$VCF1}{$contig}{$position};
					next POS if(!defined $$comparisons_variants{$VCF2}{$contig}{$position});
					my $vcf2_snp = $$comparisons_variants{$VCF2}{$contig}{$position};

					if($vcf1_snp eq $vcf2_snp) { $overlap_snps++; }
				}
			}
			foreach my $contig(keys %{$$comparisons_variants{$VCF2}}) { $total_snps_vcf2 += scalar(keys(%{$$comparisons_variants{$VCF2}{$contig}})); }
			my $percent_shared = 0;
			if(($total_snps_vcf1 + $total_snps_vcf2) ne 0) {
				$percent_shared = (($overlap_snps / (($total_snps_vcf1 + $total_snps_vcf2) / 2)) * 100);
			}

			print $ofh "$VCF1\t$VCF2\t$overlap_snps\t$percent_shared\n";
		}
		print $ofh "\n";
	}
	close $ofh;
	return 1;
}

sub check_files_exist_in_Name_Type_Location_files {
	my $Name_Type_Location_file = $_[0];
	warn "check_files_exist_in_Name_Type_Location_files...\n";
	foreach my $isolate(keys %{$Name_Type_Location_file}) {
		foreach my $file_type(keys %{$$Name_Type_Location_file{$isolate}}) {
			my $file = $$Name_Type_Location_file{$isolate}{$file_type};
			warn "Checking $file...\n";
			if(! -e $file) { die "Cannot open $file. Check Name_Type_Location file: $!\n"; }
		}
	}
	warn "check_files_exist_in_Name_Type_Location_files: Yes\n";
	return 1;
}

sub save_ECA_bases {
	my ($Name_Type_Location_file, $verified_eca_positions, $settings) = @_;
	my (%comparisons_snps, %comparisons_hets, %verified_eca_positions_per_isolate);
	my %indel_lengths;

	# warn message
	if($settings eq 1) { warn "save_ECA_bases: Saving ECA homozygous bases...\n"; }
	if($settings eq 2) { warn "save_ECA_bases: Saving ECA homozygous and heterzygous bases (with ambiguity characters)...\n"; }
	if($settings eq 3) { warn "save_ECA_bases: Saving ECA homozygous, heterzygous bases (with ambiguity characters) and indels...\n"; }

	foreach my $isolate(keys %{$Name_Type_Location_file}) {
		my $VCF = $$Name_Type_Location_file{$isolate}{'VCF'};
		warn "save_ECA_bases: $VCF\n";

		# sometimes 2 vcf lines for the same position are given (Pilon). Only consider the first one
		my $last_contig_and_position;

		# Run through VCF's extracting sequences for EXA bases
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

			# Ignore if header
			next VCF if($$VCF_line{'next'} eq 1);

			# Ignore if not an ECA position
			next VCF if(!defined $$verified_eca_positions{$supercontig}{$position});

			# Ignore ambigious positions
			next VCF if($base_type eq 'ambiguous');

			# Ignore lines given twice
			if(!defined $last_contig_and_position) { $last_contig_and_position = "$supercontig\t$position"; }
			else {
				my $current_contig_and_position = "$supercontig\t$position";
				if($current_contig_and_position eq $last_contig_and_position) {
					$last_contig_and_position = $current_contig_and_position;
					next VCF;
				}
				$last_contig_and_position = $current_contig_and_position;
			}

			# Save consensus sequences

			# homozygous
			if(($base_type eq 'snp') || ($base_type eq 'reference')) {
				$consensus = substr $consensus, 0, 1;
				die "ERROR: unexpected base found in $VCF $supercontig $position: $consensus\n" if($consensus !~ m/[ATCG]/i);
				$verified_eca_positions_per_isolate{$supercontig}{$position}{$VCF} = $consensus; 
				$comparisons_snps{$VCF}{$supercontig}{$position} = $consensus;
			}

			# heterozygous
			elsif(($settings =~ m/[23]/) && ($base_type eq 'heterozygous')) { 
				$verified_eca_positions_per_isolate{$supercontig}{$position}{$VCF} = $amb_char; 
				$comparisons_hets{$VCF}{$supercontig}{$position} = $amb_char;
			}

			# insertion
			elsif($base_type =~ m/insertion/) {
				$verified_eca_positions_per_isolate{$supercontig}{$position}{$VCF} = $consensus; 
				$comparisons_snps{$VCF}{$supercontig}{$position} = $consensus;

				# Longest insertion
				my $length = length($consensus);
				if(!defined $indel_lengths{$supercontig}{$position}{'INS'}) { $indel_lengths{$supercontig}{$position}{'INS'} = 1; }
				if($length > $indel_lengths{$supercontig}{$position}{'INS'}) { $indel_lengths{$supercontig}{$position}{'INS'} = $length; }
			}

			# deletion
			elsif($base_type =~ m/deletion/) {

				# how many dashes to add?
				my $dashes = '';
				for(my $i=0; $i < (length($ref_base) - length($consensus)); $i++) { $dashes .= '-'; }

				# save
				$verified_eca_positions_per_isolate{$supercontig}{$position}{$VCF} = ($consensus . $dashes); 
				$comparisons_snps{$VCF}{$supercontig}{$position} = ($consensus . $dashes);

				# Longest deletion
				my $length = (length($ref_base) - length($consensus)) + 1;
				if(!defined $indel_lengths{$supercontig}{$position}{'DEL'}) { $indel_lengths{$supercontig}{$position}{'DEL'} = 1; }
				if($length > $indel_lengths{$supercontig}{$position}{'DEL'}) { 
					$indel_lengths{$supercontig}{$position}{'DEL'} = $length; 
					my $del_base = substr $ref_base, length($consensus), length($ref_base);
					$indel_lengths{$supercontig}{$position}{'DEL_BASES'} = $del_base; 
				}
			}
			else {
				die "ERROR: base unaccounted for: $VCF : $line ($base_type)\n";
			}
		}
		close $fh;
	}
	return (\%verified_eca_positions_per_isolate, \%comparisons_snps, \%comparisons_hets, \%indel_lengths);
}

sub find_and_print_ECA_positions {
	my ($contig_position_to_number_VCFs_variants, $contig_position_to_number_VCFs_reference, $number_of_vcfs, $include_invariants, $outfile) = @_;

	# outfile
	open my $ofh, '>', $outfile or die "Cannot open $outfile: $!";

	warn "find_and_print_ECA_positions: $number_of_vcfs, $outfile...\n";
	my $count = 0;
	my $exclude_invariants = 0;
	my (%verified_eca_positions);
	foreach my $contig(sort keys %{$contig_position_to_number_VCFs_variants}) {
		POSITION: foreach my $position(sort { $a <=> $b } keys %{$$contig_position_to_number_VCFs_variants{$contig}}) {

			# variants found (>= 1)
			my $variants_found = $$contig_position_to_number_VCFs_variants{$contig}{$position};

			# ref bases found
			my $ref_found = 0;
			if(defined $$contig_position_to_number_VCFs_reference{$contig}{$position}) {
				$ref_found += $$contig_position_to_number_VCFs_reference{$contig}{$position};
			} 

			# invariants
			else {
				if($include_invariants eq 'n') {
					$exclude_invariants++;
					next POSITION;
				}
			}

			# ECA?
			if(($variants_found + $ref_found) eq $number_of_vcfs) {
				$count++;

				# Print base location
				#print $ofh "$contig\t$position\t$ref_found\t$variants_found\n";
				print $ofh "$contig\t$position\n";
			}
		}
	}
	close $ofh;
	if($include_invariants eq 'n') { warn "find_and_print_ECA_positions: $exclude_invariants invariants excluded\n"; }
	warn "find_and_print_ECA_positions: $count entirely covered in all (ECA)\n";
	return \%verified_eca_positions;
}

sub save_contig_position_to_number_VCFs_found {
	my ($find_files) = $_[0];
	my @variant_files = <"$find_files">;
	
	my %contig_position_to_number_VCFs;
	foreach my $file(@variant_files) {
		warn "save_contig_position_to_number_VCFs_found: saving $file...\n";
		
		open my $fh, '<', $file or die "Cannot open $file : $!";
		my $contig;
		POSITION: while(my $line = <$fh>) {
			chomp $line;
			if($line =~ m/^##/) {
				$line =~ s/^##//;
				$contig = $line;
				next POSITION;
			}

			# positions
			my @bits = split /\t/, $line;
			my ($start, $stop) = @bits;
			for(my $i = $start; $i < $stop; $i++) {
				$contig_position_to_number_VCFs{$contig}{$i}++;
			}
		}
	}
	return \%contig_position_to_number_VCFs;
}

sub print_supercontig_tab_positions_from_VCF {
	my ($Name_Type_Location_hash, $settings, $cmds_outfile, $platform, $queue) = @_;

	my @cmds;
	foreach my $isolate(keys %{$Name_Type_Location_hash}) {
		my $VCF = $$Name_Type_Location_hash{$isolate}{'VCF'};
		my $cmd = "perl $vcf_to_variants_and_ref_sites -v $VCF $settings";

		# run on grid?
		if($opt_g eq 'n') { system($cmd); }
		else {
			push @cmds, $cmd;
		}
	}
	# Run paralel cmds:
	if($opt_g ne 'n') { &run_cmds_on_grid(\@cmds, $cmds_outfile, $platform, $queue); }
	return 1;
}

sub run_cmds_on_grid {
	my ($cmds, $cmds_file, $platform, $queue) = @_;
	open my $ofh, '>', $cmds_file or die "Cannot open $cmds_file : $!\n";
	print $ofh @{$cmds};
	close $ofh;
	warn "Running commands in parallel...\n";
	my $cmd_parallel = "$Run_Commands_python --platform $platform --queue $queue --mem 2 --throttle_nodes 95 --cmds_per_node 10 $cmds_file 2>&1 | tee $cmds_file.log";
	process_cmd($cmd_parallel);	
}
