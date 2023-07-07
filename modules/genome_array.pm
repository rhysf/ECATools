package genomearray;
use strict;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);
use FindBin qw($Bin);
use lib "$Bin/";
use read_VCF_lines;

### rfarrer@broadinstitute.org

sub make_genome_hash_array_from_seq_struct {
	my $seq_struct = $_[0];
	my %genome;
	warn "make_genome_hash_array_from_seq_struct: filling genome array...\n";
	foreach my $id(keys %{$$seq_struct{'seq'}}) {
		my $length = $$seq_struct{'seq_length'}{$id};
		my @array;
		for(my $i=0; $i<=$length; $i++) { $array[$i] = 0; }
		push @{$genome{$id}}, @array;
	}
	return \%genome;
}

sub print_tab_file_from_regions_in_genome_hash {
	my ($fasta, $genome_struct, $find, $outfile) = @_;
	warn "print_tab_file_from_regions_in_genome_hash: finding $find...\n";

	# outfile
	open my $ofh, '>', $outfile or die "Cannot open $outfile : $!";

	my %contigs_seen;

	my $count_found = 0;
	foreach my $id(keys %{$$fasta{'seq'}}) {
		my $length = $$fasta{'seq_length'}{$id};
		my $start;
		my $contiguous = 0;
		#warn "$id length = $length\n";
		for(my $i=0; $i<=$length; $i++) {
			
			# no longer covered - but was on the previous base! - or final base of contig and was/is covered 
			if((($$genome_struct{$id}[$i] ne $find) && ($contiguous eq 1)) || (($i eq $length) && ($contiguous eq 1))) {
				$contiguous = 0;
				$count_found += ($i - $start);

				#warn "genome struct: $id $i ne $find\n";

				# print contig?
				if(!defined $contigs_seen{$id}) {
					$contigs_seen{$id} = 1;
					print $ofh "##$id\n";
				}

				# Print
				my $end = $i;
				print $ofh "$start\t$end\n";

				# new start
				undef $start;
			}

			# Start of coverage
			elsif($$genome_struct{$id}[$i] eq $find) {
				if(!defined $start) { $start = $i; }
				$contiguous = 1;
			}
			# Carry on
			else {}
		}
	}

	close $ofh;

	warn "print_tab_file_from_regions_in_genome_hash: found $count_found\n";
	return 1;
}

sub fill_genome_hash_array_from_vcf {
	my ($genome_hash, $vcf, $settings, $exclude, $include, $min_depth) = @_;

	# sometimes 2 vcf lines for the same position are given (Pilon). Only consider the first one
	my $last_contig_and_position;
	my %variants_found;

	my ($ref_count, $variant_count, $min_depth_filtered, $other_filtered) = (0, 0, 0, 0);
	warn "fill_genome_hash_array_from_vcf: Saving reference positions (1) and variants (2) over genome array from $vcf (setting=$settings, exclude=$exclude, include=$include, min depth=$min_depth)...\n";
	open my $fh, '<', $vcf or die "Cannot open $vcf: $!\n";
	VCF: while(my $line=<$fh>) {
		chomp $line;

		my $VCF_line = vcflines::read_VCF_lines($line);
		my $supercontig = $$VCF_line{'supercontig'};
		my $position = $$VCF_line{'position'};
		my $base_type = $$VCF_line{'base_type0'};

		# Ignore ambigious, or contigs of uninterest
		next VCF if($$VCF_line{'next'} eq 1);
		next VCF if(($exclude !~ m/n/i) && ($supercontig eq $exclude));
		next VCF if(($include !~ m/n/i) && ($supercontig ne $include));

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

		# Min depth filtering
		if((defined $$VCF_line{'DP0'}) && ($$VCF_line{'DP0'} ne '?'))  {
			if($$VCF_line{'DP0'} < $min_depth) {
				$min_depth_filtered++;
				next VCF;
			}
		}

		# heterozygous sites
		if(($settings eq 1) && ($base_type eq 'heterozygous')) {
			$other_filtered++;
			next VCF;
		}
		# indels (snp_multi = E.g. ref: ACTCGTCCTGACT consensus: AACTCGTACTGAC)
		if(($settings =~ m/[12]/) && ($base_type =~ m/insertion|deletion|snp_multi/)) {
			$other_filtered++;
			next VCF;
		}

		# Everything i do want to save
		if($base_type !~ m/heterozygous|snp|insertion|deletion|reference/) {
			$other_filtered++;
			next VCF;
		}

		# save reference bases
		if($base_type =~ m/reference/) {
			$ref_count++;
			$$genome_hash{$supercontig}[$position] = 1;
		} else {
			$variant_count++;
			$variants_found{$base_type}++;
			$$genome_hash{$supercontig}[$position] = 2;
		}
	}

	# summary
	warn "Reference bases:\t$ref_count\n";
	warn "Variant bases:\t$variant_count\n";
	foreach my $type(keys %variants_found) {
		my $count = $variants_found{$type};
		warn "Variant bases (type=$type)\t$count\n";
	}
	warn "Excluded for < min depth:\t$min_depth_filtered\n";
	warn "Excluded for base type:\t$other_filtered\n";
	return $genome_hash;
}

1;