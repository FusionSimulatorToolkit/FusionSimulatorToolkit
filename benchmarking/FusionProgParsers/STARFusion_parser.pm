package STARFusion_parser;

use strict;
use warnings;
use Carp;


sub parse_fusion_result_file {
    my ($starFusion_file) = @_;

    my @fusions;

    open (my $fh, $starFusion_file) or die "Error, cannot open file $starFusion_file";
    while (<$fh>) {
        if (/^\#/) { next; }

        chomp;

        my ($fusion, $junction_reads, $spanning_reads, $splice_type,
            $fusion_gene_A, $chr_coords_A,
            $fusion_gene_B, $chr_coords_B) = split(/\t/);

        my $rest;
        ($fusion_gene_A, $rest) = split(/\^/, $fusion_gene_A);
        ($fusion_gene_B, $rest) = split(/\^/, $fusion_gene_B);

        if ($fusion_gene_A eq $fusion_gene_B) { next; } # no self-fusions

        my ($chrA, $coordA, $orientA) = split(/:/, $chr_coords_A);
        my ($chrB, $coordB, $orientB) = split(/:/, $chr_coords_B);


        my $struct = {
            geneA => $fusion_gene_A,
            chrA => $chrA || ".",
            coordA => $coordA || ".",

            geneB => $fusion_gene_B,
            chrB => $chrB || ".",
            coordB => $coordB || ".",

            span_reads => $spanning_reads,
            junc_reads => $junction_reads,
        };

        push (@fusions, $struct);

    }

    close $fh;

    return(@fusions);
}


1; #EOM

