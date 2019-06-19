#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/modules";
use read_FASTA;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -s sequence.fasta > output\n
Output format:   -p Printing options (fasta, nexus, tba, color, fastq, none) [none]
                 -d If printing nexus (dna/protein), if printing TBA (source organism without dashes)  [dna]
		 -z Split fasta (color/base) into x files (n=ignore/all to output, 1=1 sequence per file etc.) [n]\n
Rename:          -a Append word to IDs (n=none) [n]
                 -b Replace IDs for Descriptions (yn) [n]
                 -c If opt_b, which part of Description is wanted (split by space)? [1]
		 -e Drop words from IDs (Separated by comma) E.g. supercontig,_,pilon,1. [n]
		 -q Drop description (y/n) [n]\n
Reorder:         -j Order by None/ignore (n), Biggest to smallest (b), IDs from opt_i (i), Text file from opt_k (k), FASTA file from opt_k (f) [n]
	         -k New order specified in file [n]\n
Modify sequence: -f Truncate x nt from 5' (left) [0]
                 -t Truncate x nt from 3' (right) [0]
	         -m Remove if smaller than this number of characters [1000000000000]
	         -i IDs seperated by comma (n=none) [n]
	         -r Reverse compliment on opt_i (y/n) [n]
		 -n Drop contig id []
		 -l Drop contigs in file []\n
Summarise:       -g Print summary for file (y/n) TURN OFF TO SPEED UP x2 [y]
                 -h Print summary per entry (n=none, s=short (id and length) l=long (id, desc, seq, length etc.)) [n]\n";
our($opt_a, $opt_b, $opt_c, $opt_d, $opt_e, $opt_f, $opt_g, $opt_h, $opt_i, $opt_j, $opt_k, $opt_l, $opt_m, $opt_n, $opt_p, $opt_q, $opt_r, $opt_s, $opt_t, $opt_z);
getopt('abcdefghijklmnpqrstz');
if(!defined $opt_a) { $opt_a = 'n'; }
if(!defined $opt_b) { $opt_b = 'n'; }
if(!defined $opt_c) { $opt_c = 1; }
if(!defined $opt_d) { $opt_d = 'dna'; }
if(!defined $opt_e) { $opt_e = 'n'; }
if(!defined $opt_f) { $opt_f = 0; }
if(!defined $opt_g) { $opt_g = 'y'; }
if(!defined $opt_h) { $opt_h = 'n'; }
if(!defined $opt_i) { $opt_i = 'n'; }
if(!defined $opt_j) { $opt_j = 'n'; }
if(!defined $opt_m) { $opt_m = 1000000000000; }
if(!defined $opt_p) { $opt_p = 'none'; }
if(!defined $opt_q) { $opt_q = 'n'; }
if(!defined $opt_r) { $opt_r = 'n'; }
if(!defined $opt_t) { $opt_t = 0; }
if(!defined $opt_z) { $opt_z = 'n'; }
die $usage unless ($opt_s);
die "Cannot open $opt_s : $!" unless (-e $opt_s);
die "If defined -p, must be fasta, nexus, tba or none: $opt_p\n" if ($opt_p !~ m/fasta|nexus|none|tba|color|fastq/);
die "Setting -g must be y or n\n" if (($opt_g ne 'n') && ($opt_g ne 'y'));
die "Setting -q must be y or n\n" if (($opt_q ne 'n') && ($opt_q ne 'y'));
die "Setting -h must be l, s or n\n" if (($opt_h ne 'n') && ($opt_h ne 'l') && $opt_h ne 's');

# Save sequences
my $fasta = fastafile::fasta_to_struct($opt_s);

### Rename
# Append to ID
if($opt_a ne 'n') { $fasta = fastafile::fasta_struct_append_id($fasta, $opt_a); }
# Replace ID with Description
if($opt_b ne 'n') { $fasta = fastafile::fasta_struct_replace_id_with_desc($fasta, $opt_c); }
# Drop words from ID
if($opt_e ne 'n') { $fasta = fastafile::fasta_struct_remove_words_from_id($fasta, $opt_e); }
# Drop description
if($opt_q eq 'y') { $fasta = fastafile::fasta_struct_remove_description($fasta); }

### Reorder
if($opt_j ne 'n') { $fasta = fastafile::fasta_struct_reorder($fasta, $opt_j, $opt_i, $opt_k); }

### Modify
# Truncate or drop small sequences
if(($opt_f > 0) || ($opt_t > 0) || ($opt_m < 1000000000000)) { $fasta = fastafile::fasta_struct_truncate($fasta, $opt_f, $opt_t, $opt_m); }
# Reverse compliment
if($opt_r ne 'n') { $fasta = fastafile::fasta_struct_reverse_compliment($fasta, $opt_i); } 
# Drop contigs from name or file
if($opt_n) { $fasta = fastafile::fasta_struct_remove_id($fasta, $opt_n); }
if($opt_l) { $fasta = fastafile::fasta_struct_remove_ids_from_text_file($fasta, $opt_l); }

### Summarise and print
my $information;
if(($opt_g ne 'n') || ($opt_h ne 'n')) { $information = fastafile::fasta_summary($opt_s, $opt_g, $opt_h); }
fastafile::fasta_struct_print($fasta, $opt_p, $opt_d, $opt_z);
