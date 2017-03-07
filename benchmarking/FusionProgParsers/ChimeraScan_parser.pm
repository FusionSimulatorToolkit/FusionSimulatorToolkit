package ChimeraScan_parser;

use strict;
use warnings;
use Carp;

sub parse_fusion_result_file {
    my ($chimeraScan_file) = @_;

    my @fusions;

    open (my $fh, $chimeraScan_file) or die "Error, cannot open file $chimeraScan_file";
    my $header = <$fh>;
    while (<$fh>) {
        if (/^\#/) { next; }
        chomp;
        my @x = split(/\t/);

        my $chrA = $x[0];
        my $chrA_start = $x[1];
        my $chrA_end = $x[2];

        my $chrB = $x[3];
        my $chrB_start = $x[4];
        my $chrB_end = $x[5];

        my $chrA_strand = $x[8];
        my $chrB_strand = $x[9];


        my $geneA = $x[12];
        my $geneB = $x[13];


        my $brkpt_A = ($chrA_strand eq '+') ? $chrA_end : $chrA_start;
        my $brkpt_B = ($chrB_strand eq '+') ? $chrB_end : $chrB_start;


        my $total_frags = $x[16];

        my $junction_count = $x[17];
        my $spanning_count = $total_frags - $junction_count;


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

    close $fh;

    return(@fusions);
}

1; #EOM

