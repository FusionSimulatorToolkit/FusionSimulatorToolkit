package MapSplice_parser;

use strict;
use warnings;
use Carp;

=mapsplice_format

    Go here: http://www.netlab.uky.edu/p/bioinfo/MapSplice2FusionJunctionFormat

=cut


sub parse_fusion_result_file {
    my ($mapsplice_out_file) = @_;

    my @fusions;

    my $get_unique_gene_list_sref = sub {
        my ($gene_txt) = @_;

        my %genes;
        my @fields = split(/,/, $gene_txt);
        foreach my $gene (@fields) {
            if ($gene) {
                $genes{$gene} = 1;
            }
        }

        my $unique_gene_list = join(",", keys %genes);

        return($unique_gene_list);
    };


    open (my $fh, $mapsplice_out_file) or die "Error, cannot open file $mapsplice_out_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);

        my $geneA =  &$get_unique_gene_list_sref($x[60]);
        my $geneB =  &$get_unique_gene_list_sref($x[61]);

        my $junction_count = $x[4];
        my $spanning_count = $x[27];

        my ($chrA, $chrB) = split(/\~/, $x[0]);
        unless ($chrA =~ /chr/ && $chrB =~ /chr/) {
            confess "Erorr, didn't parse chr vals from $x[0] of $_";
        }

        $chrA =~ s/chr//;
        $chrB =~ s/chr//;

        my $brkpt_A = $x[1];
        my $brkpt_B = $x[2];

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

