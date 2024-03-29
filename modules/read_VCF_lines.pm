package vcflines;
use strict;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);
use Data::Dumper;

### rfarrer@broadinstitute.org

sub read_VCF_lines {
	my $VCF_line = $_[0];
	my @bits = split /\t/, $VCF_line;
	my %VCF_info;

	# Save header info including isolate names
	if($VCF_line =~ m/^\#/) { 
		my $VCF_struct = &VCF_header_to_struct($VCF_line, \%VCF_info); 
		return $VCF_struct; 
	} else {
		$VCF_info{'next'}=0;
		$VCF_info{'header'}='N';
	}

	# Initial quality check
	if(@bits < 9) {
		warn "$0: Bad VCF with < 9 columns: $VCF_line\n";
		$VCF_info{'next'}=1; 
		return \%VCF_info;
	}

	# Multi VCF
	if(@bits > 10) {
		for(my $i=9; $i < scalar(@bits); $i++) {
			my $sample_info_id = ('sample_info' . ($i - 9));
			#warn "$sample_info_id = $i = $bits[$i]\n";
			$VCF_info{$sample_info_id} = $bits[$i];
		}
	}

	# Parts continued
	$VCF_info{'supercontig'}          = $bits[0];
	$VCF_info{'position'}             = $bits[1];
	$VCF_info{'id'}                   = $bits[2];
	$VCF_info{'reference_VCF_format'} = $bits[3];
	$VCF_info{'consensus_VCF_format'} = $bits[4];
	$VCF_info{'cons_qual'}            = $bits[5];
	$VCF_info{'filter'}               = $bits[6];
	$VCF_info{'info'}                 = $bits[7];
	$VCF_info{'format'}               = $bits[8];

	# Split format parts, and save available names 
	my @format_parts = split /:/, $VCF_info{'format'};
	my %format_part_ids_available;
	foreach(@format_parts) { $format_part_ids_available{$_} = 1; }
	$VCF_info{'number_of_samples'} = (scalar(@bits) - 9);

	# Sample_info deals with Multi VCF as well
	for(my $i=9; $i < scalar(@bits); $i++) {

		# Save sample_info[isolate number]
		my $isolate_number = ($i - 9);
		my $sample_info_id = ('sample_info' . $isolate_number);
		$VCF_info{$sample_info_id} = $bits[$i];
		my @sample_info_parts = split /:/, $VCF_info{$sample_info_id};

		# Check format parts matches sample info parts
		# Die if less format parts than sample info parts
		die "$VCF_line: scalar(@sample_info_parts) < scalar(@format_parts) (sample info < format parts). Have not coded for this eventuality!\n" if(scalar(@format_parts) < scalar(@sample_info_parts));
		# Subset sample_info_parts by format_parts
		if(scalar(@format_parts) > scalar(@sample_info_parts)) {
			my @reduced_sample_info_parts;
			for(my $i=0; $i<scalar(@format_parts); $i++) {
				push @reduced_sample_info_parts, $sample_info_parts[$i];
			}
			@sample_info_parts = @reduced_sample_info_parts;
		}
		die "Format and Sample Info do not match: $VCF_line\n" if(scalar(@format_parts) ne scalar(@sample_info_parts));

		# Save genotype, depth etc.
		my ($GT_id, $DP_id, $base_type_ID, $amb_char_ID) = ("GT$isolate_number", "DP$isolate_number", "base_type$isolate_number", "amb_char$isolate_number");
		for(my $f=0; $f<scalar(@format_parts); $f++) {
			my $format_part = $format_parts[$f];
			my $sample_part = $sample_info_parts[$f];
			my $format_part_id = ($format_part . $isolate_number);
			die "$format_part_id overwrites known format id on $VCF_line\n" if(defined $format_part_ids_available{$format_part_id});
			$VCF_info{$format_part_id} = $sample_part;
		}
		die "Unable to find genotype ($GT_id): $VCF_line\n" if(!defined $VCF_info{$GT_id});
		if(!defined $VCF_info{$DP_id}) { $VCF_info{$DP_id} = '?'; }

		# Check alleles are different
		if($VCF_info{$GT_id} =~ m/[\/\|]/) {
			my @allele_parts = split /[\/\|]/, $VCF_info{$GT_id};
			die "Have not coded for multiple alleles for $VCF_line in check alleles are different\n" if(scalar(@allele_parts) ne 2);
			if($allele_parts[0] eq $allele_parts[1]) { 
				$VCF_info{$GT_id} = $allele_parts[0]; 
			}
		}

		# Determine base (base2 only if diploid)
		my ($base1, $base2, $base_type) = &VCF_struct_determine_bases_and_base_type(\%VCF_info, $GT_id);
		$VCF_info{$base_type_ID}= $base_type;
		$VCF_info{($isolate_number . 'base1')} = $base1;
		$VCF_info{($isolate_number . 'base2')} = $base2;
		
		# Ambiguity character
		if($VCF_info{$base_type_ID} eq 'heterozygous') { $VCF_info{$amb_char_ID} = &get_ambiguity_char($base1, $base2); }
	}	

	# Phased (by my scripts)
	if($VCF_info{'info'} =~ m/PHASE\=(.+)/) {
		$VCF_info{'phased'} = 1;
		my $phase_name = $1;
		$phase_name =~ s/\;//g;
		$VCF_info{'phase_group'} = $phase_name;
		#$VCF_info{'base_type'}='phased'; 
	}

	# Return
	return \%VCF_info;
}

############### Local subroutines

sub VCF_header_to_struct {
	my ($VCF_line, $VCF_struct) = @_;
	my @bits = split /\t/, $VCF_line;
	$$VCF_struct{'next'}=1; 
	$$VCF_struct{'header'}='Y';
	if($VCF_line =~ m/^\#CHROM\tPOS\tID\tREF/) {
		for(my $i=9; $i < scalar(@bits); $i++) {
			$$VCF_struct{'isolate_names'}{($i - 9)} = $bits[$i];
		}
	}
	return $VCF_struct;
}

sub get_ambiguity_char {
	my ($base1, $base2) = @_;
	my $ambiguity_char;
				
	# K
	if(($base1 eq 'T') && ($base2 eq 'G')) { $ambiguity_char = 'K'; }
	elsif(($base1 eq 'G') && ($base2 eq 'T')) { $ambiguity_char = 'K'; }
							
	# M
	elsif(($base1 eq 'A') && ($base2 eq 'C')) { $ambiguity_char = 'M'; }
	elsif(($base1 eq 'C') && ($base2 eq 'A')) { $ambiguity_char = 'M'; }
										
	# R
	elsif(($base1 eq 'A') && ($base2 eq 'G')) { $ambiguity_char = 'R'; }
	elsif(($base1 eq 'G') && ($base2 eq 'A')) { $ambiguity_char = 'R'; }
													
	# Y
	elsif(($base1 eq 'T') && ($base2 eq 'C')) { $ambiguity_char = 'Y'; }
	elsif(($base1 eq 'C') && ($base2 eq 'T')) { $ambiguity_char = 'Y'; }
																
	# S
	elsif(($base1 eq 'G') && ($base2 eq 'C')) { $ambiguity_char = 'S'; }
	elsif(($base1 eq 'C') && ($base2 eq 'G')) { $ambiguity_char = 'S'; }
																			
	# W
	elsif(($base1 eq 'A') && ($base2 eq 'T')) { $ambiguity_char = 'W'; }
	elsif(($base1 eq 'T') && ($base2 eq 'A')) { $ambiguity_char = 'W'; }
	return $ambiguity_char;
}

sub VCF_struct_determine_bases_and_base_type {
	my($VCF_struct, $GT_id) = @_;
	my ($base1, $base2, $base_type);
	$base2 = 'None';

	# Variables from struct that are needed
	foreach my $needed_vars($GT_id, 'reference_VCF_format', 'consensus_VCF_format') {
		die "Error: $needed_vars not defined in VCF_struct\n" if(!defined $$VCF_struct{$needed_vars});
	}
	my $GT = $$VCF_struct{$GT_id};
	my $ref_base = $$VCF_struct{'reference_VCF_format'};
	my $cons_base = $$VCF_struct{'consensus_VCF_format'};

	# Save bases and GT parts
	my @consensus_bases = split /,/, $cons_base;
	my @all_bases = ($ref_base);
	push(@all_bases, @consensus_bases);
	my @GT_parts = split /[\/\|]/, $GT;
	my @GT_bases;
	foreach my $GT_part(@GT_parts) {

		# Ambiguous GT (all things non numerical such as . or *)
		if($GT_part !~ m/^\d+$/) {
			if($GT_part ne '.') { warn "something ambiguous: $GT_part\n"; }
			$base1 = 'N';
			$base_type = 'ambiguous';
			return ($base1, $base2, $base_type);
		}
		die "VCF_struct_determine_bases_and_base_type: Error. No $GT_part defined from ref $ref_base and cons $cons_base\n" if(!defined $all_bases[$GT_part]);

		# Ambiguous bases
		if($all_bases[$GT_part] =~ m/[N\.\*]/) {
			$base1 = 'N';
			$base_type = 'ambiguous';
			return ($base1, $base2, $base_type);
		}

		# Save bases from GT
		push @GT_bases, $all_bases[$GT_part];
	}

	# polyploid?
	die "VCF_struct_determine_bases_and_base_type: polyploid not implemented yet: $GT\n" if(scalar(@GT_parts) > 2);

	# Homozygous ref-calls
	if($GT eq 0) {
		$base1 = $ref_base;
		$base_type = 'reference';
		return ($base1, $base2, $base_type);
	}

	# Homozygous non-reference
	elsif(($GT ne 0) && (scalar(@GT_bases) eq 1)) {
		my $consensus = $GT_bases[0];

		# Homozygous SNP
		if(length($ref_base) eq length($consensus)) { 

			# SNP
			if((length($ref_base) eq 1) && (length($consensus) eq 1)) { 
				$base1 = $consensus;
				$base_type = 'snp';
				return ($base1, $base2, $base_type);
			}

			# SNP(s) disguised as an indel
			if(length($ref_base) eq length($consensus)) {
				my @bases_reference = split //, $ref_base;
				my @bases_consensus = split //, $consensus;
				my $snp_count = 0;
				my ($ref_base_saved, $cons_base_saved);
				for(my $i=0; $i<scalar(@bases_reference); $i++) {
					my $ref_base = $bases_reference[$i];
					my $cons_base = $bases_consensus[$i];
					if($ref_base ne $cons_base) { 
						$ref_base_saved = $ref_base;
						$cons_base_saved = $cons_base;
						$snp_count++; 
					}
				}
				if($snp_count eq 0) {
					$base1 = $ref_base;
					$base_type = 'reference';
					return ($base1, $base2, $base_type);
				}
				elsif($snp_count eq 1) {
					$base1 = $ref_base_saved;
					$base2 = $cons_base_saved;
					$base_type = 'snp';
					return ($base1, $base2, $base_type);
				}
				if($snp_count > 1) {
					$base1 = $consensus;
					$base_type = ('snp_multi' . $snp_count);
					return ($base1, $base2, $base_type);
				}
			}
	
			# Ambiguous?
			warn "Nothing found for this apparant homozygous snp:\n";
			warn Dumper($VCF_struct);
			$base1 = 'N';
			$base_type = 'ambiguous';
			return ($base1, $base2, $base_type);
		}

		# Homozygous indel
		if(length($ref_base) ne length($consensus)) {

			# Deletion (maybe with snps in there too!)
			if(length($ref_base) > length($consensus)) { 
				$base1 = $consensus;
				$base_type = 'deletion';
				return ($base1, $base2, $base_type);
			}
			if(length($ref_base) eq length($consensus)) { 
				$base1 = $consensus;
				$base_type = 'deletion';
				return ($base1, $base2, $base_type);
			}

			# Insertion (maybe with snps in there too!)
			if(length($ref_base) < length($consensus)) { 
				$base1 = $consensus;
				$base_type = 'insertion';
				return ($base1, $base2, $base_type);
			}

			# Ambiguous?
			warn "Nothing found for this apparent homozygous indel:\n";
			warn Dumper($VCF_struct);
			$base1 = 'N';
			$base_type = 'ambiguous';
			return ($base1, $base2, $base_type);
		}
	}

	# Bi-allelic heterozygous positions & indels
	elsif(scalar(@GT_bases) eq 2) {
		$base1 = $GT_bases[0];
		$base2 = $GT_bases[1];
		$base_type = 'heterozygous';
		if(length($base1) < length($base2)) { $base_type = 'het_insertion'; }
		if(length($base1) > length($base2)) { $base_type = 'het_deletion'; }
		return ($base1, $base2, $base_type);
	}

	# Ambiguous?
	else {
		die "VCF_struct_determine_bases_and_base_type: Error. Undefined variant: ref $ref_base and cons $cons_base and GT ($GT) = @GT_bases\n";
	}
	#return ($base1, $base2, $base_type);
}

1;
