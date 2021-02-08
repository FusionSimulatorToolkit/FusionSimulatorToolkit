#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib("$FindBin::Bin/../PerlLib");
use Fasta_reader;

my $num_fusion_genes_select = 200;
my $fasta_file = "fusion_isoforms.db.fasta";


open(my $fh, $fasta_file) or die "Error, cannot open file: $fasta_file";
my $fasta_reader = new Fasta_reader($fasta_file);

while(my $seq_obj = $fasta_reader->next()) {
    

    my $acc = $seq_obj->get_accession();
    my $header = $seq_obj->get_header();
    
    my $sequence = $seq_obj->get_sequence();

    my $left_fusion_seq = $sequence;
    $left_fusion_seq =~ s/[a-z]//g;
    my $right_fusion_seq = $sequence;
    $right_fusion_seq =~ s/[A-Z]//g; 

    
    my $len_left = length($left_fusion_seq);
    my $len_right = length($right_fusion_seq);


    print join("\t", $acc, $len_left, $len_right, $header, $sequence) . "\n"; 

    #print "$acc\n$left_fusion_seq\n$right_fusion_seq\n\n";
}

