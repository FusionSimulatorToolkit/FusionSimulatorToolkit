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
#  --gene_spans <string>      reference annotation gene spans file
#
#  --truth_fusions <string>   file containing a list of the true fusions.
#
#  --fusion_preds <string>    fusion predictions ranked accordingly.
#
#
# Optional:
#
#  --unsure_fusions <string>   fusions where we're not sure if it's a TP or FP
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



my $help_flag;
my $gene_spans_file;
my $fusion_preds_file;
my $truth_fusions_file;

my $ALLOW_REVERSE_FUSION = 0;

my $ALLOW_PARALOGS = 0;
my $paralogs_file;
my $unsure_fusions_file;

&GetOptions ( 'h' => \$help_flag,
              'gene_spans=s' => \$gene_spans_file,
              'fusion_preds=s' => \$fusion_preds_file,
              'truth_fusions=s' => \$truth_fusions_file,
              'unsure_fusions=s' => \$unsure_fusions_file,
              
              'allow_reverse_fusion' => \$ALLOW_REVERSE_FUSION,
              
              'allow_paralogs=s' => \$paralogs_file,
    );


if ($help_flag) { die $usage; }

unless ($gene_spans_file && $fusion_preds_file && $truth_fusions_file) {
    
    die $usage;

}

if ($paralogs_file) {
    $ALLOW_PARALOGS = 1;
}


main : {

    my %gene_id_to_coords;
    my %interval_trees = &build_interval_trees($gene_spans_file, \%gene_id_to_coords);
    
    my %TP_fusions = &parse_fusion_listing($truth_fusions_file);
    
    my %unsure_fusions;
    if ($unsure_fusions_file) {
        %unsure_fusions = &parse_fusion_listing($unsure_fusions_file);
    }
    
    my %FP_fusions;
    my %seen_TP;

    my %paralog_fusion_to_TP_fusion;
    if ($ALLOW_PARALOGS) {
        %paralog_fusion_to_TP_fusion = &parse_paralogs_integrate_parafusions(\%TP_fusions, $paralogs_file);
    }

    
    my %prog_names;
    # print header
    print join("\t", "#pred_result", "ProgName", "Sample", "FusionName", "J", "S", "explanation") . "\n";

    open (my $fh, $fusion_preds_file) or die "Error, cannot open file $fusion_preds_file";
    while (<$fh>) {
        chomp;
        if (/^\#/) { next; }
        my @x = split(/\t/);
        
        my ($sample, $prog_name, $fusion_name, $J, $S) = @x;
        my ($geneA, $geneB) = split(/--/, $fusion_name);
        
        $prog_names{$prog_name} = 1;
        
        $fusion_name = "$sample|$fusion_name";
        
        my $reverse_fusion_name = "$sample|$geneB--$geneA";
        
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
                    # reset important vars here
                    
                    $using_para_fusion_flag = 1;
                    
                    print STDERR "-substituting PARA fusion ($adj_fusion) for TP fusion ($fusion_name)\n";
                    
                    $fusion_name = $adj_fusion;
                    $fusion_name =~ /\|(\S+)--(\S+)/ or die "Error, cannot decode fusion name: $fusion_name";
                    ($geneA, $geneB) = ($1, $2);
                    
                    $reverse_fusion_name = "$sample|$geneB--$geneA";
                }
            } 
        }
        
        
        my ($accuracy_token, $accuracy_explanation);

        ############################
        ## Check for already seen TP
        
        if ($seen_TP{"$prog_name,$fusion_name"}) {
            $accuracy_token = "NA-TP";
            $accuracy_explanation = "already scored $fusion_name as TP";
        }
        
        elsif ($ALLOW_REVERSE_FUSION && $seen_TP{"$prog_name,$reverse_fusion_name"}) {
            $accuracy_token = "NA-TP_rev";
            $accuracy_explanation = "already scored $fusion_name (rev) as TP";
        }

        ############################
        ## Check for already seen FP
        
        elsif ($FP_fusions{"$prog_name,$fusion_name"}) {
            $accuracy_token = "NA-FP";
            $accuracy_explanation = "already scored $fusion_name as FP";
        }

        elsif ($ALLOW_REVERSE_FUSION && $FP_fusions{"$prog_name,$reverse_fusion_name"}) {
            $accuracy_token = "NA-FP_rev";
            $accuracy_explanation = "already scored $fusion_name (rev) as FP";
        }
        
        ###########
        ## Check to see if we should ignore it
        elsif (%unsure_fusions && $unsure_fusions{$fusion_name}) {
            $accuracy_token = "NA-UNCLASS";
            $accuracy_explanation = "not classifying $fusion_name, in unsure list";
        }
        elsif (%unsure_fusions && $ALLOW_REVERSE_FUSION && $unsure_fusions{$reverse_fusion_name}) {
            $accuracy_token = "NA-UNCLASS_rev";
            $accuracy_explanation = "not classifying $fusion_name (rev), in unsure list";
        }
                
        ###############################
        ## Check for new TP recognition
        
        elsif ($TP_fusions{$fusion_name}) {
            $accuracy_token = "TP";
            $seen_TP{"$prog_name,$fusion_name"} = 1;
            $accuracy_explanation = "first encounter of TP $fusion_name";
        }
        
        elsif ($ALLOW_REVERSE_FUSION && $TP_fusions{$reverse_fusion_name}) {
            $accuracy_token = "TP";
            $seen_TP{"$prog_name,$fusion_name"} = 1;
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
                my %found_unsure;
                
                
                foreach my $gA (@overlapping_genesA) {
                    foreach my $gB (@overlapping_genesB) {
                        
                        my $candidate_fusion = "$sample|" . join("--", $gA, $gB);
                        my $reverse_candidate_fusion = "$sample|" . join("--", $gB, $gA);
                        
                        ###########################
                        # check for already seen TP
                        if ($seen_TP{"$prog_name,$candidate_fusion"}) {
                            $existing_TPs{$candidate_fusion} = 1;
                        }

                        elsif ($ALLOW_REVERSE_FUSION && $seen_TP{"$prog_name,$reverse_candidate_fusion"}) {
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
                        
                        #####################
                        ## Check for unsure status
                        elsif ($unsure_fusions{$candidate_fusion}) {
                            $found_unsure{$candidate_fusion} = 1;
                        }
                        elsif ($ALLOW_REVERSE_FUSION && $unsure_fusions{$reverse_candidate_fusion}) {
                            $found_unsure{$reverse_candidate_fusion} = 1;
                        }
                        
                        #############################
                        ## Check for already seen FPs
                        
                        elsif ($FP_fusions{"$prog_name,$candidate_fusion"}) {
                            $existing_FPs{$candidate_fusion} = 1;
                        }
                        
                        elsif ($ALLOW_REVERSE_FUSION && $FP_fusions{"$prog_name,$reverse_candidate_fusion"}) {
                            $existing_FPs{$reverse_candidate_fusion} = 1;
                        }
                    }
                }
                
                if (my @new_TPs = keys %newly_found_TPs) {
                    $accuracy_token = "TP";
                    foreach my $f (@new_TPs) {
                        $seen_TP{"$prog_name,$f"} = 1;
                    }
                    $accuracy_explanation = "chr mapping to first encounter of TP " . join(",", @new_TPs);
                }
                elsif (my @existing_TPs = keys %existing_TPs) {
                    $accuracy_token = "NA-TP";
                    $accuracy_explanation = "chr mapping to already scored TP " . join(",", @existing_TPs);
                }
                elsif (my @found_unsures = keys %found_unsure) {
                    $accuracy_token = "NA-UNCLASS";
                    $accuracy_explanation = "chr mapping to unsure entry " . join(",", @found_unsures);
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
            $FP_fusions{"$prog_name,$fusion_name"} = 1;
        }
        
        
        if ($using_para_fusion_flag) {
            $accuracy_explanation .= " (using_para_fusion)";
        }
        
        $accuracy_explanation =~ s/\s/_/g;
        
        print join("\t", $accuracy_token, $prog_name, $sample, @x, $accuracy_explanation) . "\n";
    }
    

    ## Report false-negatives (known fusions not predicted)
    
    foreach my $prog_name (keys %prog_names) {
        foreach my $fusion_name (keys %TP_fusions) {
            if (! $seen_TP{"$prog_name,$fusion_name"}) {
                my ($sample_name, $core_fusion_name) = split(/\|/, $fusion_name);
                print join("\t", "FN", $prog_name, $sample_name, $core_fusion_name, 0, 0, "lacking_this_fusion_prediction") . "\n";
            }
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
sub parse_fusion_listing {
    my ($fusions_file) = @_;

    my %fusions;
    
    open (my $fh, $fusions_file) or die $!;
    while (<$fh>) {
        chomp;
        my $fusion = $_;
        if ($fusion =~ /^(\S+)\|(\S+)--(\S+)$/) {
            $fusions{$fusion} = 1;
        }
        else {
            die "Error, cannot parse fusion: $fusion as samplename|fusionA--fusionB";
        }
    }
    close $fh;

    return(%fusions);

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
    my ($orig_fusions_href, $paralogs_file) = @_;
    
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


    
    my %paralog_fusion_to_orig_fusion;

    my @orig_fusions = keys %$orig_fusions_href;

    foreach my $orig_fusion (@orig_fusions) {
        my ($sample, $orig_fusion_name) = split(/\|/, $orig_fusion);
        my ($geneA, $geneB) = split(/--/, $orig_fusion_name);
        
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

                my $para_fusion = "$sample|$gA--$gB";
                
                $paralog_fusion_to_orig_fusion{$para_fusion} = $orig_fusion;
            }
        }
    }
    
    return(%paralog_fusion_to_orig_fusion);

}

