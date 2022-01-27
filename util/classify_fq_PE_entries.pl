#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;
use Fastq_reader;

my $usage = "\n\n\tusage: $0 fusions_n_unfused_targets.fasta left_fq right_fq\n\n";

my $targets_fa = $ARGV[0] or die $usage;
my $left_fq = $ARGV[1] or die $usage;
my $right_fq = $ARGV[2] or die $usage;



my $fasta_reader = new Fasta_reader($targets_fa);


my $SPLIT_ANCHOR_REQUIRED = 10;
my $FUZZ = 5;

my %acc_to_breakpoint;

my $counter = 0;
while(1) {

    my $fusion_entry = $fasta_reader->next();
    if (! $fusion_entry) {
        last;
    }
    
    ## examine fusion
    my $fusion_acc = $fusion_entry->get_accession();
    {
        my $header = $fusion_entry->get_header();
        $header =~ /FusedAt:(\d+)/;
        my $brkpt_pos = $1 or die "Error, cannot extract breakpoint from header: $header";
        $acc_to_breakpoint{$fusion_acc} = $brkpt_pos;
    }

    
    ## examine left unfused
    my $left_orig_entry = $fasta_reader->next();
    my $left_acc = $left_orig_entry->get_accession();
    {
        # DPP7|ENST00000497375.1 brkpt: 172
        my $header = $left_orig_entry->get_header();
        $header =~ /brkpt: (\d+)/;
        my $brkpt_pos = $1 or die "Error, cannot extract breakpoint from header: $header";
        $acc_to_breakpoint{$left_acc} = $brkpt_pos;
    }

    ## examine right unfused
    my $right_orig_entry = $fasta_reader->next();
    my $right_acc = $right_orig_entry->get_accession();
    {
        # DPP7|ENST00000497375.1 brkpt: 172
        my $header = $right_orig_entry->get_header();
        $header =~ /brkpt: (\d+)/;
        my $brkpt_pos = $1 or die "Error, cannot extract breakpoint from header: $header";
        $acc_to_breakpoint{$right_acc} = $brkpt_pos;
    }

    $counter += 1;
    #if ($counter >= 20) { last; }
}


my $left_fq_reader = new Fastq_reader($left_fq);
my $right_fq_reader = new Fastq_reader($right_fq);

while (my $left_fq_entry = $left_fq_reader->next()) {
    my $right_fq_entry = $right_fq_reader->next();
    
    my $core_read_name = $left_fq_entry->get_core_read_name();
    # @DPP7|ENSG00000176978.9--ELP4|ENSG00000109911.13_18_606_3:0:0_2:0:0_e/1
    
    my @vals = split(/_/, $core_read_name);
    
    my $transcript  = $vals[0];
    my $brkpt_pos = $acc_to_breakpoint{$transcript} or die "Error, no brkpt for [$transcript] of $core_read_name";
    
    my $frag_start = $vals[-6];
    my $frag_end = $vals[-5];
    
    my $read_len_left = length($left_fq_entry->get_sequence());
    my $read_len_right = length($right_fq_entry->get_sequence());
    
    
    my $lend_read_start = $frag_start;
    my $lend_read_end = $frag_start + $read_len_left -1;
    
    my $rend_read_start = $frag_end - $read_len_right + 1;
    my $rend_read_end = $frag_end;
    
    
    # check if overlaps breakpoint.
    my $read_type = "other";
    if ($frag_start > $brkpt_pos || $frag_end < $brkpt_pos) {
        $read_type = "No_brkpt_overlap";
    }
    else {
        
        $read_type = "other";
        
        ### brkpt position:     |
        ##               AGCTAGCTactgactg
        
        
        if ( ($lend_read_start + $SPLIT_ANCHOR_REQUIRED -1 <= $brkpt_pos && $brkpt_pos <= $lend_read_end - $SPLIT_ANCHOR_REQUIRED) 
             ||
             ($rend_read_start + $SPLIT_ANCHOR_REQUIRED -1 <= $brkpt_pos && $brkpt_pos <= $rend_read_end - $SPLIT_ANCHOR_REQUIRED) ) {
            
            
            $read_type = "split";
            
        }
        
        
        #  LEND span fuzz                         <=======================>--------------------------<=================> 
        #    brkpt                                                       |
        #  REND span fuzz    <=========>--------------------------------<====================>
        #                                                                ^
        #                                                             fuzzy overlap allowed
        
        elsif ($lend_read_end <= $brkpt_pos + $FUZZ && $brkpt_pos + 1 - $FUZZ < $rend_read_start) {
            
            $read_type = "span";
        }
    }
    
    
    print(join("\t", $transcript, $core_read_name, $brkpt_pos, "$lend_read_start-$lend_read_end", "$rend_read_start-$rend_read_end", $read_type) . "\n");
    
}

exit(0);


