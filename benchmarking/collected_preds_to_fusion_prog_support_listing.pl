#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\n\tusage: $0 fusion_preds.collected\n\n";

my $preds_file = $ARGV[0] or die $usage;

main: {
    
    my %fusion_to_prog;

    open (my $fh, $preds_file) or die "Error, cannot open file $preds_file";
    while (<$fh>) {
        chomp;
        my ($sample_name, $prog, $fusion_name, $junc_support, $frag_support) = split(/\t/);
        
        $fusion_name = "$sample_name|$fusion_name";
        
        $fusion_to_prog{$fusion_name}->{$prog}++;
    }
    close $fh;

    my @fusion_structs;
    foreach my $fusion_name (keys %fusion_to_prog) {
        my $progs_href = $fusion_to_prog{$fusion_name};
        
        my @prognames = sort keys %$progs_href;
        my $num_progs = scalar(@prognames);
        
        push (@fusion_structs, { fusion_name => $fusion_name,
                                 prognames => \@prognames,
                                 count => $num_progs,
              } );

    }

    @fusion_structs = reverse sort {$a->{count} <=> $b->{count} } @fusion_structs;

    foreach my $fusion_struct (@fusion_structs) {
        print join("\t", $fusion_struct->{fusion_name}, 
                   join(",", @{$fusion_struct->{prognames}}),
                   $fusion_struct->{count},
            ) . "\n";
    }

    exit(0);
}

