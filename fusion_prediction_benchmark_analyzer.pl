#!/usr/bin/env perl

use strict;
use warnings;
use Set::IntervalTree;


use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);

my $usage = <<__EOUSAGE__;

######################################################################
#
#  --gene_spans <string>      reference annotation gene spans file
#
#  --truth_fusions <string>   file containing a list of the true fusions.
#
#  --fusion_preds <string>    fusion predictions ranked accordingly.
#
#
######################################################################


__EOUSAGE__


    ;



my $help_flag;
my $gene_spans_file;
my $fusion_preds_file;
my $truth_fusions_file;

&GetOptions ( 'h' => \$help_flag,
              'gene_spans=s' => \$gene_spans_file,
              'fusion_preds=s' => \$fusion_preds_file,
              'truth_fusions=s' => \$truth_fusions_file,
    );


unless ($gene_spans_file && $fusion_preds_file && $truth_fusions_file) {
    die $usage;
}


main : {

    my %gene_id_to_coords;
    my %interval_trees = &build_interval_trees($gene_spans_file, \%gene_id_to_coords);
    
    my %TP_fusions = &parse_TP_fusions($truth_fusions_file);

    my %FP_fusions;

    open (my $fh, $fusion_preds_file) or die "Error, cannot open file $fusion_preds_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);

        my $fusion_name = $x[0];
        my ($geneA, $geneB) = split(/--/, $fusion_name);
        
        my @overlapping_genesA = &find_overlapping_genes($geneA, \%interval_trees, \%gene_id_to_coords);

        my @overlapping_genesB = &find_overlapping_genes($geneB, \%interval_trees, \%gene_id_to_coords);
        
        print "geneA: $geneA has overlaps: @overlapping_genesA\n";
        print "geneB: $geneB has overlaps: @overlapping_genesB\n";
    }


    exit(0);
    
    
}

####
sub build_interval_trees {
    my ($gene_spans_file, $gene_id_to_coords_href) = @_;

    my %interval_trees;
    
    open (my $fh, $gene_spans_file) or die "Error, cannot open file $gene_spans_file";
    while (<$fh>) {
        chomp;
        my ($gene_id, $chr, $lend, $rend, $orient, $gene_symbol) = split(/\t/);

        my $itree = $interval_trees{$chr};
        unless (ref $itree) {
            $itree = $interval_trees{$chr} = Set::IntervalTree->new;
        }

        $itree->insert($gene_symbol, $lend, $rend);
    
        $gene_id_to_coords_href->{$gene_symbol} = { chr => $chr,
                                                    lend => $lend,
                                                    rend => $rend };
    }
    close $fh;
    
    return(%interval_trees);
}


####
sub parse_TP_fusions {
    my ($truth_fusions_file) = @_;

    my %TP_fusions;
    
    open (my $fh, $truth_fusions_file) or die $!;
    while (<$fh>) {
        chomp;
        my $fusion = $_;
        if ($fusion =~ /^(\w+)--(\w+)$/) {
            $TP_fusions{$fusion} = 1;
        }
        else {
            die "Error, cannot parse fusion: $fusion as fusionA--fusionB";
        }
    }
    close $fh;

    return(%TP_fusions);
}

####
sub find_overlapping_genes {
    my ($gene, $interval_trees_href, $gene_id_to_coords_href) = @_;

    my $gene_info_struct = $gene_id_to_coords_href->{$gene} or die "Error, not finding gene info for [$gene] ";
    
    my $chr = $gene_info_struct->{chr};
    my $lend = $gene_info_struct->{lend};
    my $rend = $gene_info_struct->{rend};

    my $itree = $interval_trees_href->{$chr} or die "Error, no interval tree stored for chr [$chr]";

    my $overlaps_aref = $itree->fetch($lend, $rend);

    return(@$overlaps_aref);

}

