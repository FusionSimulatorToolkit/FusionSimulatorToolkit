#!/usr/bin/env perl

use strict;
use warnings;
use Set::IntervalTree;


use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);

my $usage = <<__EOUSAGE__;

#################################################################################################
#
# Required:
#
#  --progname <string>        name of fusion predictor
#
#  --sample <string>          name of sample being processed
#
#  --gene_spans <string>      reference annotation gene spans file
#
#  --truth_fusions <string>   file containing a list of the true fusions.
#
#  --fusion_preds <string>    fusion predictions ranked accordingly.
#
#
# Optional:
#
#  --allow_reverse_fusion     if true fusion is A--B, allow for it to be reported as B--A
#
#  --allow_paralogs <string>  file containing tab-delimited list of paralog clusters
#                             so if TP_A--TP_B  is a true fusion,
#                               paraA1--paraB2  would be considered an ok proxy and scored as a TP.
#
##################################################################################################


__EOUSAGE__


    ;


my $progname = "";
my $sample = "";

my $help_flag;
my $gene_spans_file;
my $fusion_preds_file;
my $truth_fusions_file;

my $ALLOW_REVERSE_FUSION = 0;

my $ALLOW_PARALOGS = 0;
my $paralogs_file;

&GetOptions ( 'h' => \$help_flag,
              
              'progname=s' => \$progname,
              'sample=s' => \$sample,
              
              'gene_spans=s' => \$gene_spans_file,
              'fusion_preds=s' => \$fusion_preds_file,
              'truth_fusions=s' => \$truth_fusions_file,
              
              'allow_reverse_fusion' => \$ALLOW_REVERSE_FUSION,
              
              'allow_paralogs=s' => \$paralogs_file,
    );


if ($help_flag) { die $usage; }

unless ($progname && $sample
        && 
        $gene_spans_file && $fusion_preds_file && $truth_fusions_file) {

    die $usage;

}

if ($paralogs_file) {
    $ALLOW_PARALOGS = 1;
}


main : {

    my %gene_id_to_coords;
    my %interval_trees = &build_interval_trees($gene_spans_file, \%gene_id_to_coords);
    
    my %TP_fusions = &parse_TP_fusions($truth_fusions_file);

    my %FP_fusions;
    my %seen_TP;

    my %paralog_fusion_to_TP_fusion;
    if ($ALLOW_PARALOGS) {
        %paralog_fusion_to_TP_fusion = &parse_paralogs_integrate_parafusions(\%TP_fusions, $paralogs_file);
    }

    
    # print header
    print join("\t", "#pred_result", "ProgName", "Sample", "FusionName", "J", "S", "explanation") . "\n";

    open (my $fh, $fusion_preds_file) or die "Error, cannot open file $fusion_preds_file";
    while (<$fh>) {
        chomp;
        if (/^\#/) { next; }
        my @x = split(/\t/);
        
        my $fusion_name = $x[0];
        my ($geneA, $geneB) = split(/--/, $fusion_name);
                
        my $reverse_fusion_name = "$geneB--$geneA";
        
        my $using_para_fusion_flag = 0;


        if ($ALLOW_PARALOGS) {
            unless ($TP_fusions{$fusion_name} || ($ALLOW_REVERSE_FUSION && $TP_fusions{$reverse_fusion_name}) ) {
                
                my $adj_fusion;
                
                if (my $TP_fusion = $paralog_fusion_to_TP_fusion{$fusion_name}) {
                    $adj_fusion = $TP_fusion;
                }
                elsif ($ALLOW_REVERSE_FUSION && ($TP_fusion = $paralog_fusion_to_TP_fusion{$reverse_fusion_name}) ) {
                    $adj_fusion = $TP_fusion;
                }
                
                if ($adj_fusion) {
                    # rest important vars here
                    
                    $using_para_fusion_flag = 1;
                    
                    print STDERR "-substituting PARA fusion ($adj_fusion) for TP fusion ($fusion_name)\n";
                    
                    $fusion_name = $adj_fusion;
                    ($geneA, $geneB) = split(/--/, $fusion_name);
                    
                    $reverse_fusion_name = "$geneB--$geneA";
                }
            } 
        }
        

        my ($accuracy_token, $accuracy_explanation);

        ############################
        ## Check for already seen TP
        
        if ($seen_TP{$fusion_name}) {
            $accuracy_token = "NA-TP";
            $accuracy_explanation = "already scored $fusion_name as TP";
        }
        
        elsif ($ALLOW_REVERSE_FUSION && $seen_TP{$reverse_fusion_name}) {
            $accuracy_token = "NA-TP_rev";
            $accuracy_explanation = "already scored $fusion_name (rev) as TP";
        }

        ############################
        ## Check for already seen FP
        
        elsif ($FP_fusions{$fusion_name}) {
            $accuracy_token = "NA-FP";
            $accuracy_explanation = "already scored $fusion_name as FP";
        }

        elsif ($ALLOW_REVERSE_FUSION && $FP_fusions{$reverse_fusion_name}) {
            $accuracy_token = "NA-FP_rev";
            $accuracy_explanation = "already scored $fusion_name (rev) as FP";
        }
        
        ###############################
        ## Check for new TP recognition
        
        elsif ($TP_fusions{$fusion_name}) {
            $accuracy_token = "TP";
            $seen_TP{$fusion_name} = 1;
            $accuracy_explanation = "first encounter of TP $fusion_name";
        }
        
        elsif ($ALLOW_REVERSE_FUSION && $TP_fusions{$reverse_fusion_name}) {
            $accuracy_token = "TP";
            $seen_TP{$fusion_name} = 1;
            $accuracy_explanation = "first encounter of TP $fusion_name (rev)";
        }
        
        else {
            # haven't seen it yet and not a known TP
            # map to overlapping genes on the genome and consider equivalent.
            
            
            
            my @overlapping_genesA = &find_overlapping_genes($geneA, \%interval_trees, \%gene_id_to_coords);
            
            my @overlapping_genesB = &find_overlapping_genes($geneB, \%interval_trees, \%gene_id_to_coords);
            
            if (@overlapping_genesA && @overlapping_genesB) {
                ## explore combinations between A and B pairs
                
                my %newly_found_TPs;
                my %existing_TPs;
                my %existing_FPs;
                
                
                foreach my $gA (@overlapping_genesA) {
                    foreach my $gB (@overlapping_genesB) {
                        
                        my $candidate_fusion = join("--", $gA, $gB);
                        my $reverse_candidate_fusion = join("--", $gB, $gA);
                        
                        ###########################
                        # check for already seen TP
                        if ($seen_TP{$candidate_fusion}) {
                            $existing_TPs{$candidate_fusion} = 1;
                        }

                        elsif ($ALLOW_REVERSE_FUSION && $seen_TP{$reverse_candidate_fusion}) {
                            $existing_TPs{$reverse_candidate_fusion} = 1;
                        }
                        
                        ####################
                        ## Check for new TPs
                        elsif ($TP_fusions{$candidate_fusion}) { 
                            $newly_found_TPs{$candidate_fusion} = 1;
                        }

                        elsif ($ALLOW_REVERSE_FUSION && $TP_fusions{$reverse_candidate_fusion}) {
                            $newly_found_TPs{$reverse_candidate_fusion} = 1;
                        }
                        
                        #############################
                        ## Check for already seen FPs
                        
                        elsif ($FP_fusions{$candidate_fusion}) {
                            $existing_FPs{$candidate_fusion} = 1;
                        }
                        
                        elsif ($ALLOW_REVERSE_FUSION && $FP_fusions{$reverse_candidate_fusion}) {
                            $existing_FPs{$reverse_candidate_fusion} = 1;
                        }
                    }
                }
                
                if (my @new_TPs = keys %newly_found_TPs) {
                    $accuracy_token = "TP";
                    foreach my $f (@new_TPs) {
                        $seen_TP{$f} = 1;
                    }
                    $accuracy_explanation = "chr mapping to first encounter of TP " . join(",", @new_TPs);
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
            
            else {
                ## no overlapping genes found
                my @missing_genes;
                unless (@overlapping_genesA) {
                    push (@missing_genes, $geneA);
                }
                unless (@overlapping_genesB) {
                    push (@missing_genes, $geneB);
                }
                $accuracy_token = "NA-unknown";
                $accuracy_explanation = "cannot find record of " . join(",", @missing_genes);
                
            }
        }
        
        unless ($accuracy_token) {
            # must be a FP
            $accuracy_token = "FP";
            $accuracy_explanation = "first encounter of FP fusion $fusion_name";
            $FP_fusions{$fusion_name} = 1;
        }
        
        
        if ($using_para_fusion_flag) {
            $accuracy_explanation .= " (using_para_fusion)";
        }
        
        $accuracy_explanation =~ s/\s/_/g;
        
        print join("\t", $accuracy_token, $progname, $sample, @x, $accuracy_explanation) . "\n";
    }
    

    ## Report false-negatives (known fusions not predicted)
    
    foreach my $fusion_name (keys %TP_fusions) {
        if (! $seen_TP{$fusion_name}) {
            print join("\t", "FN", $progname, $sample, $fusion_name, 0, 0, "lacking_this_fusion_prediction") . "\n";
        }
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
        
        my $gene_info_struct = $gene_id_to_coords_href->{$gene};
        
        unless ($gene_info_struct) {
            print STDERR "WARNING - not finding a record of gene $gene\n";
            next;
        }
        
        my $chr = $gene_info_struct->{chr};
        my $lend = $gene_info_struct->{lend};
        my $rend = $gene_info_struct->{rend};
        
        my $itree = $interval_trees_href->{$chr} or die "Error, no interval tree stored for chr [$chr]";
        
        my $overlaps_aref = $itree->fetch($lend, $rend);
        
        if ($overlaps_aref && @$overlaps_aref) {
            
            unless ( grep { $_ eq $gene } @$overlaps_aref) {
                die "Error, found overlapping genes for $gene, but doesn't include $gene : { @$overlaps_aref } ";
            }
            push (@overlapping_genes, @$overlaps_aref);
        }
    }

    return(@overlapping_genes);

}

####
sub parse_paralogs_integrate_parafusions {
    my ($TP_fusions_href, $paralogs_file) = @_;
    
    my %gene_to_para_list;
    {
        open (my $fh, $paralogs_file) or die $!;
        while (<$fh>) {
            chomp;
            my @x = split(/\s+/);
            
            foreach my $gene (@x) {
                $gene_to_para_list{$gene} = \@x;
            }
        }
        close $fh;
    }


    
    my %paralog_fusion_to_TP_fusion;

    my @TP_fusions = keys %$TP_fusions_href;

    foreach my $TP_fusion (@TP_fusions) {
        my ($geneA, $geneB) = split(/--/, $TP_fusion);
        
        my @paraA = ($geneA);
        if (my $para_aref = $gene_to_para_list{$geneA}) {
            @paraA = @$para_aref;
        }
        my @paraB = ($geneB);
        if (my $para_aref = $gene_to_para_list{$geneB}) {
            @paraB = @$para_aref;
        }

        foreach my $gA (@paraA) {
            foreach my $gB (@paraB) {

                my $para_fusion = "$gA--$gB";
                
                $paralog_fusion_to_TP_fusion{$para_fusion} = $TP_fusion;
            }
        }
    }
    
    return(%paralog_fusion_to_TP_fusion);

}

