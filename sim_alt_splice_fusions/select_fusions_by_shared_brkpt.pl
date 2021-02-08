#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(shuffle);

my $usage = "usage: $0 fusions.reorg_tab (left|right) num_genes\n";

my $fusions_reorg_tab_file = $ARGV[0] or die $usage;
my $left_or_right = $ARGV[1] or die $usage;
my $num_genes = $ARGV[2] or die $usage;

unless ($left_or_right =~ /^(left|right)$/) {
    die $usage;
}


my %fusion_brkpt_list;

open(my $fh, $fusions_reorg_tab_file) or die $!;
while(my $line = <$fh>) {
    my @vals = split(/\t/, $line);
    my $acc = $vals[0];
    my $brkpt = ($left_or_right eq "left") ? $vals[1] : $vals[2];
    
    push (@{$fusion_brkpt_list{$acc}->{$brkpt}}, \@vals);
}

my @genes = keys %fusion_brkpt_list;

@genes = shuffle(@genes);

my $num_genes_reported = 0;

foreach my $gene (@genes) {
    
    my $gene_brkpt_fusions_href = $fusion_brkpt_list{$gene};
    
    foreach my $gene_brkpt_fusion (shuffle(keys %$gene_brkpt_fusions_href)) {
        
        my @entries = @{$gene_brkpt_fusions_href->{$gene_brkpt_fusion}};
        my $num_entries = scalar(@entries);
    
        if ($num_entries >= 2) {
            my @indices = 0..($num_entries-1);
            @indices = shuffle(@indices);
            print join("\t", @{$entries[$indices[0]]});
            print join("\t", @{$entries[$indices[1]]});
            
            $num_genes_reported += 1;
            last;
        }
        

        
    }
    if ($num_genes_reported >= $num_genes) {
        last;
    }
}


