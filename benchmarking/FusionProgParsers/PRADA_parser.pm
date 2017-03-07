package PRADA_parser;

use strict;
use warnings;
use Carp;

sub parse_fusion_result_file {
    my ($prada_file) = @_;

    my @fusions;

    open (my $fh, $prada_file) or die "Error, cannot open file $prada_file";
    my $header = <$fh>;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);

        my $span_count = $x[6];

        my $fusion_info = $x[11];
        # CPNE1:20:34243124_PI3:20:43804502,2

        my @fusion_evidence = split(/\|/, $fusion_info);
        foreach my $f_info (@fusion_evidence) {

            $f_info =~ /^(\S+):([^\:]+):(\d+)_(\S+):([^\:]+):(\d+),(\d+)$/ or die "Error, cannot parse $f_info";

            my $geneA = $1;
            my $chrA = $2;
            my $coordA = $3;

            my $geneB = $4;
            my $chrB = $5;
            my $coordB = $6;

            my $junc_reads = $7;


            my $struct = {
                geneA => $geneA,
                chrA => $chrA,
                coordA => $coordA,

                geneB => $geneB,
                chrB => $chrB,
                coordB => $coordB,

                span_reads => $span_count,
                junc_reads => $junc_reads,
            };

            push (@fusions, $struct);
        }
    }


    close $fh;

    return(@fusions);
}

1; #EOM

