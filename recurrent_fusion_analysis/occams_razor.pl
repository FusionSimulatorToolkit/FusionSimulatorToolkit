#!/usr/bin/env perl

use strict;
use warnings;
use Set::IntervalTree;
use Carp;

use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);

my $usage = <<__EOUSAGE__;

#################################################################################################
#
# Required:
#
#  --fusion_preds <string>    fusion predictions (tab delimited output file where
#                              the first two columns are:   sample_name (tab) fusion_name
#                              
#
#  --gene_spans <string>      reference annotation gene spans file
#
#
# Optional:
#
#
#  --allow_paralogs <string>  file containing tab-delimited list of paralog clusters
#                             so if TP_A--TP_B  is a true fusion,
#                               paraA1--paraB2  would be considered an ok proxy and scored as a TP.
#
##################################################################################################


__EOUSAGE__


    ;


my $help_flag;
my $fusion_preds_file;
my $gene_spans_file;

my $ALLOW_PARALOGS = 0;
my $paralogs_file;

&GetOptions ( 'h' => \$help_flag,
              
              'gene_spans=s' => \$gene_spans_file,
              'fusion_preds=s' => \$fusion_preds_file,
              
              'allow_paralogs=s' => \$paralogs_file,
    );


if ($help_flag) { die $usage; }

unless ($gene_spans_file && $fusion_preds_file) { 

    die $usage;

}

if ($paralogs_file) {
    $ALLOW_PARALOGS = 1;
}


main : {

    my %gene_id_to_coords = &parse_gene_spans($gene_spans_file);
    
    
    my %gene_to_para_list;
    if ($ALLOW_PARALOGS) {
        %gene_to_para_list = &parse_paralog_info($paralogs_file);
    }
    
    
    open (my $fh, $fusion_preds_file) or die "Error, cannot open file $fusion_preds_file";
    my $header = <$fh>;
    chomp $header;
    unless ($header =~ /^\#/) {
        die "Error, table header line must start with '#'";
    }
    
    my %fusions_to_samples = &parse_fusion_preds($fh);
    
    
    ## order fusions according to recurrence and in desc order

    my @prioritized_fusions = &prioritize_recurrent_fusions(\%fusions_to_samples);
    
    # examine fusion gene mappings and update fusion names according to Occam's razor:
    #  if a different pair of gene names match up to the same genomic locus
    #  and it was reported less frequently than a more dominant fusion pair,
    #  assign the dominant fusion pair name to that fusion.

    
    print join("\t", $header, "orig_fusion_name") . "\n";
    
    my %interval_trees;
    my %seen_fusion;
    
    foreach my $fusion (@prioritized_fusions) {

        my ($geneA, $geneB) = split(/--/, $fusion);
        
        my @overlapping_genesA = &find_overlapping_genes($geneA, \%interval_trees, \%gene_id_to_coords);
        
        my @overlapping_genesB = &find_overlapping_genes($geneB, \%interval_trees, \%gene_id_to_coords);
        
        my $use_fusion_name = "";
        
        if (@overlapping_genesA && @overlapping_genesB) {
            ## explore combinations between A and B pairs
                
            foreach my $gA (@overlapping_genesA) {
                foreach my $gB (@overlapping_genesB) {
                    
                    my $candidate_fusion = "$gA--$gB";
                    if ($seen_fusion{$candidate_fusion}) {
                        $use_fusion_name = $candidate_fusion;
                        last;
                    }
                }
                last if $use_fusion_name;
            }
        }
        
        unless ($use_fusion_name) {
            $use_fusion_name = $fusion;
            
            ## add to interval tree for next time
            &add_fusion_genes_to_interval_tree($fusion, \%interval_trees, \%gene_id_to_coords);
            $seen_fusion{$fusion} = 1; # new candidate for others to map to.
        }
        
        ## report output
        &report_fusions($fusions_to_samples{$fusion}, $use_fusion_name);

    }
    
    exit(0);
    
}



####
sub parse_gene_spans {
    my ($gene_spans_file) = @_;
        
    my %gene_id_to_coords;
    
    open (my $fh, $gene_spans_file) or die "Error, cannot open file $gene_spans_file";
    while (<$fh>) {
        chomp;
        my ($gene_id, $chr, $lend, $rend, $orient, $gene_symbol) = split(/\t/);
        
        $gene_id_to_coords{$gene_symbol} = { chr => $chr,
                                             lend => $lend,
                                             rend => $rend };
    }
    close $fh;
    
    return(%gene_id_to_coords);
}

####
sub find_overlapping_genes {
    my ($genes, $interval_trees_href, $gene_id_to_coords_href) = @_;
   

    my @overlapping_genes;
    # allow for comma-separated lists
    foreach my $gene (split(/,/, $genes)) {
        
        my $gene_info_struct = $gene_id_to_coords_href->{$gene};
        
        unless ($gene_info_struct) {
            confess "ERROR - not finding a record of gene $gene\n";
        }
        
        my $chr = $gene_info_struct->{chr};
        my $lend = $gene_info_struct->{lend};
        my $rend = $gene_info_struct->{rend};
        
        my $itree = $interval_trees_href->{$chr};
        unless ($itree) { next; }
        
        my $overlaps_aref = $itree->fetch($lend, $rend);
        
        if ($overlaps_aref && @$overlaps_aref) {
            
            push (@overlapping_genes, @$overlaps_aref);
        }
    }

    return(@overlapping_genes);

}

####
sub parse_paralog_info {
    my ($paralogs_file) = @_;
    
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

    return(%gene_to_para_list);
    

}


####
sub parse_fusion_preds {
    my ($fh) = @_;
        
    
    my %fusions_to_samples;
    
    while (<$fh>) {
        chomp;
        my $line = $_;
        
        if (/^\#/) { next; }
        my @x = split(/\t/);
        
        my $sample_name = $x[0];
        my $fusion_name = $x[1];
        

        push (@{$fusions_to_samples{$fusion_name}->{$sample_name}}, $line);
        
    }
    close $fh;

    return(%fusions_to_samples);
    
}            
 
####
sub prioritize_recurrent_fusions {
    my ($fusions_to_samples_href) = @_;

    my %fusion_sample_counter;

    foreach my $fusion (keys %$fusions_to_samples_href) {
        
        my $num_samples = scalar(keys %{$fusions_to_samples_href->{$fusion}});
        $fusion_sample_counter{$fusion} = $num_samples;
    }

    my @fusions = reverse sort {$fusion_sample_counter{$a} <=> $fusion_sample_counter{$b}} keys %fusion_sample_counter;

    return(@fusions);
}


####
sub add_fusion_genes_to_interval_tree {
    my ($fusion_name, $interval_trees_href, $gene_id_to_coords_href) = @_;
    
    my ($geneA, $geneB) = split(/--/, $fusion_name);
    foreach my $gene ($geneA, $geneB) {
        &add_gene_to_interval_tree($gene, $interval_trees_href, $gene_id_to_coords_href);
    }

    return;
}

####
sub add_gene_to_interval_tree {
    my ($gene, $interval_trees_href, $gene_id_to_coords_href) = @_;

    
    if (exists $interval_trees_href->{__GENE_ADDED}->{$gene}) {
        return;
    }

    my $gene_coords_href = $gene_id_to_coords_href->{$gene};
    my ($chr, $lend, $rend) = ($gene_coords_href->{chr}, $gene_coords_href->{lend}, $gene_coords_href->{rend});
    
    my $itree = $interval_trees_href->{$chr};
    unless (ref $itree) {
        $itree = $interval_trees_href->{$chr} = Set::IntervalTree->new;
    }
    
    $itree->insert($gene, $lend, $rend);
 
    $interval_trees_href->{__GENE_ADDED}->{$gene} = 1;
    
    return;
}

####
sub report_fusions {
    my ($sample_fusions_href, $use_fusion_name) = @_;

    foreach my $sample (keys %$sample_fusions_href) {
        
        my @lines = @{$sample_fusions_href->{$sample}};
        foreach my $line (@lines) {
            my @x = split(/\t/, $line);
            my $fusion_name = $x[1];
            if ($fusion_name ne $use_fusion_name) {
                push (@x, $fusion_name);
                $x[1] = $use_fusion_name;
            }
            else {
                push (@x, "***=***");
            }

            print join("\t", @x) . "\n";
        }
    }

    return;
}
                

