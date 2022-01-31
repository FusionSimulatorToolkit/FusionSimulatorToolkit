#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;
use Fastq_reader;

my $usage = "usage: $0 fusions_n_unfused_targets.fasta output_prefix read_length\n\n";

my $targets_fa = $ARGV[0] or die $usage;
my $output_prefix = $ARGV[1] or die $usage;
my $read_len = $ARGV[2] or die $usage;

my $fasta_reader = new Fasta_reader($targets_fa);


my $SPLIT_ANCHOR_REQUIRED = 10;


my $output_left_fq = $output_prefix . ".single_fq";
my $output_stats_file = $output_prefix . ".frag_stats";

open(my $single_fq_ofh, ">$output_left_fq") or die $!;
open(my $stats_ofh, ">$output_stats_file") or die $!;


print $stats_ofh join("\t", 'fusion_name', 'split', 'left_gene', 'counter_left', 'right_gene', 'counter_right') . "\n";

my $counter = 0;
while(1) {

    my $fusion_entry = $fasta_reader->next();
    if (! $fusion_entry) {
        last;
    }
    
    ## examine fusion
    my $fusion_acc = $fusion_entry->get_accession();
    my $split_read_count = &sim_fusion_reads($fusion_entry, $single_fq_ofh, $read_len);
    
    ## examine left unfused
    my $left_orig_entry = $fasta_reader->next();
    my $left_acc = $left_orig_entry->get_accession();
    my $left_unfused_read_count = &sim_unfused_reads($left_orig_entry, $single_fq_ofh, $read_len);
    

    ## examine right unfused
    my $right_orig_entry = $fasta_reader->next();
    my $right_acc = $right_orig_entry->get_accession();
    my $right_unfused_read_count = &sim_unfused_reads($right_orig_entry, $single_fq_ofh, $read_len);
    
    print $stats_ofh join("\t", $fusion_acc, $split_read_count,
        $left_acc, $left_unfused_read_count,
        $right_acc, $right_unfused_read_count) . "\n";

    $counter += 1;
    #if ($counter >= 20) { last; }
}


####
sub sim_fusion_reads {
    my ($fusion_entry, $single_fq_ofh, $read_len) = @_;
    
    my $header = $fusion_entry->get_header();
    $header =~ /FusedAt:(\d+)/;
    my $brkpt_pos = $1 or die "Error, cannot extract breakpoint from header: $header";


    my $count_split = 0;
    
    while ($count_split == 0) {

        my $single_fq = &sim_rnaseq_reads($fusion_entry, $read_len);
        
        my $fq_reader = new Fastq_reader($single_fq);
        
        
        
        while (my $fq_entry = $fq_reader->next()) {
            
            my $core_read_name = $fq_entry->get_core_read_name();
            # @TUBB2B|ENST00000259818.7_396_803_2:0:0_4:0:0_flip0_5/1
            
            my @vals = split(/_/, $core_read_name);
            my $frag_start = $vals[-6];
            my $frag_end = $vals[-5];
            my $flip_info = $vals[-2];
            
            if ($flip_info eq "flip0") {
                $frag_end = $frag_start + $read_len - 1;
            }
            elsif ($flip_info eq "flip1") {
                $frag_start = $frag_end - $read_len + 1;
            }
            else {
                confess "Error, no flip info in $core_read_name";
            }
            
            # check if overlaps breakpoint.
            if ($frag_start <= $brkpt_pos - $SPLIT_ANCHOR_REQUIRED && $frag_end  >=  $brkpt_pos + $SPLIT_ANCHOR_REQUIRED ) {
                
                $count_split += 1;
                
                my $fq_record = $fq_entry->get_fastq_record();
                print $single_fq_ofh $fq_record;
            }
        }
    }
    return($count_split);
}


####
sub sim_unfused_reads {
    my ($seq_entry, $single_fq_ofh, $read_len) = @_;

    # DPP7|ENST00000497375.1 brkpt: 172
    my $header = $seq_entry->get_header();
    $header =~ /brkpt: (\d+)/;
    my $brkpt_pos = $1 or die "Error, cannot extract breakpoint from header: $header";
    
    my ($single_fq) = &sim_rnaseq_reads($seq_entry, $read_len);

    my $fq_reader = new Fastq_reader($single_fq);
    
    my $count_split = 0;
    
    while (my $fq_entry = $fq_reader->next()) {

        my $core_read_name = $fq_entry->get_core_read_name();
        # @TUBB2B|ENST00000259818.7_396_803_2:0:0_4:0:0_flip0_5/1

        my @vals = split(/_/, $core_read_name);
        my $frag_start = $vals[-6];
        my $frag_end = $vals[-5];
        my $flip_info = $vals[-2];

        if ($flip_info eq "flip0") {
            $frag_end = $frag_start + $read_len - 1;
        }
        elsif ($flip_info eq "flip1") {
            $frag_start = $frag_end - $read_len + 1;
        }
        else {
            confess "Error, no flip info in $core_read_name";
        }
        
        # check if overlaps breakpoint.
        if ($frag_start <= $brkpt_pos - $SPLIT_ANCHOR_REQUIRED && $frag_end  >=  $brkpt_pos + $SPLIT_ANCHOR_REQUIRED ) {
          
            $count_split += 1;
                
            my $fq_record = $fq_entry->get_fastq_record();
            print $single_fq_ofh $fq_record;
        }
    }
    

    return($count_split);
}


####
sub sim_rnaseq_reads {
    my ($seq_obj, $read_len) = @_;
    
    my $depth_of_cov = int(rand(91)) + 10; # min 10x cov

    my $header = $seq_obj->get_header();
    my $accession = $seq_obj->get_accession();
    my $sequence = $seq_obj->get_sequence();

    my $num_reads = int ( ( length($sequence) * $depth_of_cov) / $read_len );
    
    my $tmp_fa_file = "tmp.fa";
    open(my $tmp_fa_ofh, ">$tmp_fa_file") or die $!;
    print $tmp_fa_ofh ">$accession\n$sequence\n";
    close $tmp_fa_ofh;

    my $left_fq_file = "tmp.single_fq";
    my $right_fq_file = "tmp.to_delete.fq";

    my $cmd = "wgsim-trans $tmp_fa_file $left_fq_file $right_fq_file -N $num_reads -1 $read_len -2 $read_len -e 0.001 -r 0"; #-e 0 -r 0 -R 0 ";
    &process_cmd($cmd);
    
    unlink($right_fq_file);

    return($left_fq_file);
    
}

####
sub process_cmd {
    my ($cmd) = @_;
    
    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }
        
}
