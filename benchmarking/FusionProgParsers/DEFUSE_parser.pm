package DEFUSE_parser;

use strict;
use warnings;
use Carp;

sub parse_fusion_result_file {
    my ($defuse_out_file) = @_;

    my @fusions;

    open (my $fh, $defuse_out_file) or die "Error, cannot open file $defuse_out_file";
    my $header = <$fh>;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);

        my $geneA = $x[30];
        my $geneB = $x[31];

        my $junction_count = $x[2]; # splitr_count
        unless ($junction_count =~ /\w/) {
            $junction_count = 0;
        }

        my $spanning_count = $x[56]; # span_count
        unless ($spanning_count =~ /\w/) {
            $spanning_count = 0;
        }

        my $chrA = $x[24];
        my $brkpt_A = $x[37];
        my $chrB = $x[25];
        my $brkpt_B = $x[38];

        my $struct = {

            geneA => $geneA,
            chrA => $chrA,
            coordA => $brkpt_A,

            geneB => $geneB,
            chrB => $chrB,
            coordB => $brkpt_B,

            span_reads => $spanning_count,
            junc_reads => $junction_count,
        };

        push (@fusions, $struct);
    }

    return(@fusions);
}

1; #EOM

