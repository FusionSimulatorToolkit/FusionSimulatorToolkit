package FusionCatcher_parser;

use strict;
use warnings;
use Carp;

sub parse_fusion_result_file {
    
    my ($fusionCatcher_file) = @_;
    
    my @fusions;
    
    open (my $fh, $fusionCatcher_file) or die "Error, cannot open file $fusionCatcher_file";
    my $header = <$fh>;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        
        my $geneA = $x[0];
        my $geneB = $x[1];
        
        unless ($geneA =~ /\w/ && $geneB =~ /\w/) { next; } # not scoring fusions not tied to genes.
        
        my $brkpt_A = $x[8];
        my $brkpt_B = $x[9];
        
        my ($chrA, $coordA, $orientA) = split(/:/, $brkpt_A);
        my ($chrB, $coordB, $orientB) = split(/:/, $brkpt_B);
        
        my $spanning_count = $x[4];
        my $junction_count = $x[5];
        
        my $struct = {
            
            geneA => $geneA,
            chrA => $chrA,
            coordA => $coordA,
            
            geneB => $geneB,
            chrB => $chrB,
            coordB => $coordB,
            
            span_reads => $spanning_count,
            junc_reads => $junction_count,
        };
        
        push (@fusions, $struct);
    }
    
    close $fh;
    
    return(@fusions);
}

1; #EOM

