#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;
use Nuc_translator;
use Math::CDF;

my $usage = "\n\n\tusage: $0 transcripts.fasta mean_frag_size st_dev\n\n";

my $transcripts_fasta_file = $ARGV[0] or die $usage;
my $mean_frag_size = $ARGV[1] or die $usage;
my $st_dev = $ARGV[2] or die $usage;




my $num_reads = 1000;
my $read_length = 100;

my $sim_reads_outfile_prefix = "sim_short";
open(my $ofh_left_fq, ">${sim_reads_outfile_prefix}_1.fastq") or die $!;
open(my $ofh_right_fq, ">${sim_reads_outfile_prefix}_2.fastq") or die $!;


main: {

    my $fasta_reader = new Fasta_reader($transcripts_fasta_file);

    while (my $seq_obj = $fasta_reader->next()) {
        
        my $acc = $seq_obj->get_accession();
        my $header = $seq_obj->get_header();

        print STDERR "-processing $acc\n";
        
        $acc =~ /^([^\|]+)\|[^\-]+\-\-([^\|]+)/ or die "Error, cannot parse fusion name from acc: $acc";
        my $geneA = $1 or die "Error, didn't parse geneA from $acc";
        my $geneB = $2 or die "Error, didn't parse geneB from $acc";
        
        my $sequence = $seq_obj->get_sequence();
        
        $header =~ /FusedAt:(\d+)/ or die "Error, cannot parse header for FusedAt:";
        my $fusion_pos = $1 or die "Error, couldn't extract fusion position from header: $header";

        for (my $i = 0; $i < $num_reads; $i++) {

            # determine fragment length according to the normal distribution:
            my $frag_length = Math::CDF::qnorm(rand()) * $st_dev + $mean_frag_size;
            $frag_length = int($frag_length);
            
            if ($frag_length < $read_length) {
                $frag_length = $read_length;
            }
            
            # pick a random start point around the fusion breakpoint.
            my $frag_start_pos = $fusion_pos -  int(rand(int($frag_length/2)));
            if ($frag_start_pos < 1) { $frag_start_pos = 1; }

            my $frag_seq = substr($sequence, $frag_start_pos-1, $frag_length);
            # make it strand-specific RF
            
            my $right_read_seq = substr($frag_seq, 0, $read_length); # /2
            my $left_read_seq = substr(&reverse_complement($frag_seq), 0, $read_length); # /1
            
            my $acc = "read^$geneA--$geneB^S$frag_start_pos^F$frag_length^C$i";

            print $ofh_left_fq join("\n", "@" . $acc . "/1",
                                    $left_read_seq,
                                    "+",
                                    'C' x length($left_read_seq)) . "\n";
            
            print $ofh_right_fq join("\n", "@" . $acc . "/2",
                                     $right_read_seq,
                                     "+",
                                     'C' x length($right_read_seq)) . "\n";
            
            
            
            
        }
        
        
    }

    close $ofh_left_fq;
    close $ofh_right_fq;

    print STDERR "\n\n-done\n\n";
    
    
    exit(0);
}


