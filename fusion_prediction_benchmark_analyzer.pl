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
    my %seen_TP;

    open (my $fh, $fusion_preds_file) or die "Error, cannot open file $fusion_preds_file";
    while (<$fh>) {
        chomp;
        if (/^\#/) { next; }
        my @x = split(/\t/);

        my $fusion_name = $x[0];
        
        my ($accuracy_token, $accuracy_explanation);

        if ($seen_TP{$fusion_name}) {
            $accuracy_token = "NA-TP";
            $accuracy_explanation = join("_", qq(already scored $fusion_name as TP));
        }
        elsif ($FP_fusions{$fusion_name}) {
            $accuracy_token = "NA-FP";
            $accuracy_explanation = join("_", qq(already scored $fusion_name as FP));
        }
        elsif ($TP_fusions{$fusion_name}) {
            $accuracy_token = "TP";
            $seen_TP{$fusion_name} = 1;
            $accuracy_explanation = join("_", qq(first encounter of TP $fusion_name));
        }
        else {
            # haven't seen it yet and not a known TP
            # map to overlapping genes on the genome and consider equivalent.
            
            my ($geneA, $geneB) = split(/--/, $fusion_name);
            
            my @overlapping_genesA = &find_overlapping_genes($geneA, \%interval_trees, \%gene_id_to_coords);
            
            my @overlapping_genesB = &find_overlapping_genes($geneB, \%interval_trees, \%gene_id_to_coords);
            
            ## explore combinations between A and B pairs
            my %newly_found_TPs;
            my %existing_TPs;
            my %existing_FPs;
            
            
            foreach my $gA (@overlapping_genesA) {
                foreach my $gB (@overlapping_genesB) {

                    my $candidate_fusion = join("--", $gA, $gB);
                    if ($seen_TP{$candidate_fusion}) {
                        $existing_TPs{$candidate_fusion} = 1;
                    }
                    elsif ($TP_fusions{$candidate_fusion}) {
                        $newly_found_TPs{$candidate_fusion} = 1;
                    }
                    elsif ($FP_fusions{$candidate_fusion}) {
                        $existing_FPs{$candidate_fusion} = 1;
                    }
                }
            }

            if (my @new_TPs = keys %newly_found_TPs) {
                $accuracy_token = "TP";
                foreach my $f (@new_TPs) {
                    $seen_TP{$f} = 1;
                }
                $accuracy_explanation = "chr mapping to first encounter of TP " . join(", ", @new_TPs);
            }
            elsif (my @existing_TPs = keys %existing_TPs) {
                $accuracy_token = "NA-TP";
                $accuracy_explanation = "chr mapping to already scored TP " . join(",", @existing_TPs);
            }
            elsif (my @existing_FPs = keys %existing_FPs) {
                $accuracy_token = "NA-FP";
                $accuracy_explanation = "chr mapping to already scored FP " . join(",", @existing_FPs);
            }
        }

        unless ($accuracy_token) {
            # must be a FP
            $accuracy_token = "FP";
            $accuracy_explanation = "first encounter of FP fusion $fusion_name";
            $FP_fusions{$fusion_name} = 1;
        }

        print join("\t", $accuracy_token, @x, $accuracy_explanation) . "\n";
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
        if ($fusion =~ /^(\S+)--(\S+)$/) {
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
    my ($genes, $interval_trees_href, $gene_id_to_coords_href) = @_;
   

    my @overlapping_genes;
    # allow for comma-separated lists
    foreach my $gene (split(/,/, $genes)) {
        
        my $gene_info_struct = $gene_id_to_coords_href->{$gene} or die "Error, not finding gene info for [$gene] ";
        
        my $chr = $gene_info_struct->{chr};
        my $lend = $gene_info_struct->{lend};
        my $rend = $gene_info_struct->{rend};
        
        my $itree = $interval_trees_href->{$chr} or die "Error, no interval tree stored for chr [$chr]";
        
        my $overlaps_aref = $itree->fetch($lend, $rend);
        
        if ($overlaps_aref && @$overlaps_aref) {
            push (@overlapping_genes, @$overlaps_aref);
        }
    }

    return(@overlapping_genes);

}

