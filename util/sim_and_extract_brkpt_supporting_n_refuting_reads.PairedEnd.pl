#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;
use Fastq_reader;

my $usage = "usage: $0 fusions_n_unfused_targets.fasta output_prefix read_length\n\n";

my $targets_fa = $ARGV[0] or die $usage;
my $output_prefix = $ARGV[1] or die $usage;
my $read_len = $ARGV[2] or die $usage;

my $fasta_reader = new Fasta_reader($targets_fa);


my $SPLIT_ANCHOR_REQUIRED = 10;
my $FUZZ = 5;

my $output_left_fq = $output_prefix . ".left_fq";
my $output_right_fq = $output_prefix . ".right_fq";

my $output_stats_file = $output_prefix . ".frag_stats";

open(my $left_fq_ofh, ">$output_left_fq") or die $!;
open(my $right_fq_ofh, ">$output_right_fq") or die $!;
open(my $stats_ofh, ">$output_stats_file") or die $!;


print $stats_ofh join("\t", 'fusion_name', 'split', 'spanning', 'left_gene', 'counter_left', 'right_gene', 'counter_right') . "\n";

my $counter = 0;
while(1) {

    my $fusion_entry = $fasta_reader->next();
    if (! $fusion_entry) {
        last;
    }
    
    ## examine fusion
    my $fusion_acc = $fusion_entry->get_accession();
    my ($split_read_count, $span_read_count) = &sim_fusion_reads($fusion_entry, $left_fq_ofh, $right_fq_ofh, $read_len);

    ## examine left unfused
    my $left_orig_entry = $fasta_reader->next();
    my $left_acc = $left_orig_entry->get_accession();
    my $left_unfused_read_count = &sim_unfused_reads($left_orig_entry, $left_fq_ofh, $right_fq_ofh, $read_len, "LEFTcounter");
    

    ## examine right unfused
    my $right_orig_entry = $fasta_reader->next();
    my $right_acc = $right_orig_entry->get_accession();
    my $right_unfused_read_count = &sim_unfused_reads($right_orig_entry, $left_fq_ofh, $right_fq_ofh, $read_len, "RIGHTcounter");
    
    print $stats_ofh join("\t", $fusion_acc, $split_read_count, $span_read_count,
        $left_acc, $left_unfused_read_count,
        $right_acc, $right_unfused_read_count) . "\n";

    $counter += 1;
    #if ($counter >= 20) { last; }
}


####
sub sim_fusion_reads {
    my ($fusion_entry, $left_fq_ofh, $right_fq_ofh, $read_len) = @_;
    

    my $header = $fusion_entry->get_header();
    $header =~ /FusedAt:(\d+)/;
    my $brkpt_pos = $1 or die "Error, cannot extract breakpoint from header: $header";
    

    my $count_span = 0;
    my $count_split = 0;


    while ($count_split == 0) {

        my ($left_fq, $right_fq) = &sim_rnaseq_reads($fusion_entry, $read_len);
        
        my $left_fq_reader = new Fastq_reader($left_fq);
        my $right_fq_reader = new Fastq_reader($right_fq);
        
        
        
        while (my $left_fq_entry = $left_fq_reader->next()) {
            my $right_fq_entry = $right_fq_reader->next();
            
            my $core_read_name = $left_fq_entry->get_core_read_name();
            # @DPP7|ENSG00000176978.9--ELP4|ENSG00000109911.13_18_606_3:0:0_2:0:0_e/1
            
            my @vals = split(/_/, $core_read_name);
            my $frag_start = $vals[-6];
            my $frag_end = $vals[-5];
            
            # check if overlaps breakpoint.
            if ($frag_start > $brkpt_pos || $frag_end < $brkpt_pos) {
                next;
            }
            
            my $lend_read_start = $frag_start;
            my $lend_read_end = $frag_start + $read_len -1;
            
            my $rend_read_start = $frag_end - $read_len + 1;
            my $rend_read_end = $frag_end;
            
            my $read_type = "other";
            
            ### brkpt position:     |
            ##               AGCTAGCTactgactg
            
            
            if ( ($lend_read_start + $SPLIT_ANCHOR_REQUIRED -1 <= $brkpt_pos && $brkpt_pos <= $lend_read_end - $SPLIT_ANCHOR_REQUIRED) 
                 ||
                 ($rend_read_start + $SPLIT_ANCHOR_REQUIRED -1 <= $brkpt_pos && $brkpt_pos <= $rend_read_end - $SPLIT_ANCHOR_REQUIRED) ) {
                
                $count_split += 1;
                $read_type = "split";
                
            }
            
            
            #  LEND span fuzz                         <=======================>--------------------------<=================> 
            #    brkpt                                                       |
            #  REND span fuzz    <=========>--------------------------------<====================>
            #                                                                ^
            #                                                             fuzzy overlap allowed
            
            
            elsif ($lend_read_end <= $brkpt_pos + $FUZZ && $brkpt_pos + 1 - $FUZZ < $rend_read_start) {
                $count_span += 1;
                $read_type = "span";
            }
            
            my $left_fq_record = $left_fq_entry->get_fastq_record();
            $left_fq_record =~ s/\n\+\n/\n\+$read_type\n/;
            $left_fq_record =~ s/_flip/_${read_type}-flip/;
            
            my $right_fq_record = $right_fq_entry->get_fastq_record();
            $right_fq_record =~ s/\n\+\n/\n\+$read_type\n/;
            $right_fq_record =~ s/_flip/_${read_type}-flip/;
            
            print $left_fq_ofh $left_fq_record;
            print $right_fq_ofh $right_fq_record;
        }
        
    }
    
    return($count_split, $count_span);
}


####
sub sim_unfused_reads {
    my ($seq_entry, $left_fq_ofh, $right_fq_ofh, $read_len, $read_type) = @_;

    # DPP7|ENST00000497375.1 brkpt: 172
    my $header = $seq_entry->get_header();
    $header =~ /brkpt: (\d+)/;
    my $brkpt_pos = $1 or die "Error, cannot extract breakpoint from header: $header";
    
    my ($left_fq, $right_fq) = &sim_rnaseq_reads($seq_entry, $read_len);

    my $left_fq_reader = new Fastq_reader($left_fq);
    my $right_fq_reader = new Fastq_reader($right_fq);
    
    my $count_span = 0;
    
    while (my $left_fq_entry = $left_fq_reader->next()) {
        my $right_fq_entry = $right_fq_reader->next();

        my $core_read_name = $left_fq_entry->get_core_read_name();
        # @DPP7|ENSG00000176978.9--ELP4|ENSG00000109911.13_18_606_3:0:0_2:0:0_e/1
        
        my @vals = split(/_/, $core_read_name);
        my $frag_start = $vals[-6];
        my $frag_end = $vals[-5];

        # check if overlaps breakpoint.
        if ($frag_start > $brkpt_pos || $frag_end < $brkpt_pos) {
            next;
        }
    
        $count_span += 1;
        
        my $left_fq_record = $left_fq_entry->get_fastq_record();
        $left_fq_record =~ s/\n\+\n/\n\+$read_type\n/;
        $left_fq_record =~ s/_flip/_${read_type}-flip/;
        
        my $right_fq_record = $right_fq_entry->get_fastq_record();
        $right_fq_record =~ s/\n\+\n/\n\+$read_type\n/;        
        $right_fq_record =~ s/_flip/_${read_type}-flip/;

        print $left_fq_ofh $left_fq_record;
        print $right_fq_ofh $right_fq_record;
    }
    
    return($count_span);
}



####
sub sim_rnaseq_reads {
    my ($seq_obj, $read_len) = @_;
    
    my $depth_of_cov = int(rand(91)) + 10; # min 10x cov

    my $header = $seq_obj->get_header();
    my $accession = $seq_obj->get_accession();
    my $sequence = $seq_obj->get_sequence();

    my $num_reads = int ( ( length($sequence) * $depth_of_cov) / (2 * $read_len ) );
    
    my $tmp_fa_file = "tmp.fa";
    open(my $tmp_fa_ofh, ">$tmp_fa_file") or die $!;
    print $tmp_fa_ofh ">$accession\n$sequence\n";
    close $tmp_fa_ofh;

    my $left_fq_file = "tmp.left_fq";
    my $right_fq_file = "tmp.right.fq";

    my $cmd = "wgsim-trans $tmp_fa_file $left_fq_file $right_fq_file -N $num_reads -1 $read_len -2 $read_len -e 0.001 -r 0"; 
    &process_cmd($cmd);
    
    return($left_fq_file, $right_fq_file);
    
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
