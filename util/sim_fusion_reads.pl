#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;
use Nuc_translator;
use Math::CDF;


use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);



my $usage = "\n\n\tusage: $0 transcripts.fasta mean_frag_size st_dev\n\n";

my $transcripts_fasta_file;
my $mean_frag_size = 300;
my $st_dev = 50;

my $num_reads = 1000;
my $read_length = 75;
my $sim_reads_outfile_prefix = "sim";

my $usage = <<__EOUSAGE__;

###############################################################
#
#  Required:
#
#  --fusions_fasta <string>      : target transcripts fasta file
#
#  Optional:
#
#  --mean_frag_size <int>       : mean frag size (default: $mean_frag_size)
#
#  --st_dev <int>               : standard deviation (default: $st_dev)
#
#  --num_reads <int>            : number of reads per transcript (default: $num_reads)
#
#  --read_length <int>          : length of reads (default: $read_length)
#
#  --out_prefix <string>        : output prefix for sim reads (default: $sim_reads_outfile_prefix)
#
###############################################################


__EOUSAGE__

    ;

my $help_flag;


&GetOptions( 'h' => \$help_flag,
             
             'fusions_fasta=s' => \$transcripts_fasta_file,
             
             'mean_frag_size=i' => \$mean_frag_size,
             'st_dev=i' => \$st_dev,

             'num_reads=i' => \$num_reads,
             'read_length=i' => \$read_length,

             'out_prefix=s' => \$sim_reads_outfile_prefix,
    );


if ($help_flag) {
    die $usage;
}
unless ($transcripts_fasta_file) {
    die $usage;
}



open(my $ofh_left_fq, ">${sim_reads_outfile_prefix}_1.fastq") or die $!;
open(my $ofh_right_fq, ">${sim_reads_outfile_prefix}_2.fastq") or die $!;


main: {

    my $fasta_reader = new Fasta_reader($transcripts_fasta_file);

    while (my $seq_obj = $fasta_reader->next()) {
        
        my $acc = $seq_obj->get_accession();
        my $header = $seq_obj->get_header();

        print STDERR "-processing $acc\n";

        my ($geneA, $geneB);
        
        if ($acc =~ /^([^\|]+)\|[^\-]+\-\-([^\|]+)/) {
            $geneA = $1;
            $geneB = $2;
        }
        else {
            ($geneA, $geneB) = split(/--/, $acc);
        }
        unless ($geneA && $geneB) {
            die "Error, couldn't parse fusion gene names from $acc";
        }
        
        my $sequence = $seq_obj->get_sequence();
        
        $header =~ /FusedAt:(\d+)/ or die "Error, cannot parse header for FusedAt:";
        my $fusion_pos = $1 or die "Errour, couldn't extract fusion position from header: $header";

        for (my $i = 0; $i < $num_reads; $i++) {

            # determine fragment length asccording to the normal distribution:
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


