#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\tusage: $0 selected.*.tab \n\n";

my $selected_tab_file = $ARGV[0] or die $usage;


print join("\t", "ensg_pair", "ensg_pair_w_brk", "left_brkpt", "right_brkpt") . "\n";

open(my $fh, $selected_tab_file) or die "Error, cannot open file: $selected_tab_file";
while(<$fh>) {
    chomp;
    my @x = split(/\t/);
    my $ensg_pair = $x[0];
    my $lend_brkpt = $x[1];
    my $rend_brkpt = $x[2];
    
    my $header_info = $x[3];

    my @pts = split(/\s+/, $header_info);

    my $coords_info_lend = $pts[2];
    my $coords_info_rend = $pts[3];

    my $left_brkpt = &parse_brkpt($coords_info_lend, "left");
    my $right_brkpt = &parse_brkpt($coords_info_rend, "right");

    print join("\t", "$ensg_pair\t$ensg_pair^$lend_brkpt^$rend_brkpt", $left_brkpt, $right_brkpt) . "\n";
    

}    


exit(0);

####
sub parse_brkpt {
    my ($coords_info, $fusion_side) = @_;

    if ($coords_info =~ /^(chr\S+):(\d+).*\-(\d+)\[([\+\-])\]$/) {
        my $chrom = $1;
        my $left_side_brkpt = $2;
        my $right_side_brkpt = $3;
        my $orient = $4;
        
        if ($fusion_side eq "left") {
            return("$chrom:$right_side_brkpt:$orient");
        }
        else {
            return("$chrom:$left_side_brkpt:$orient");
        }
    }
    else {
        die "Error, cannot parse coords info: $coords_info";
    }
    
}


