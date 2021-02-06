#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Fasta_reader;
use GFF3_utils;
use Carp;
use Nuc_translator;
use List::Util qw(shuffle);


my $usage = "\n\nusage: $0 gff3_file genome_fasta num_genes min_chimeras_per_gene\n\n";


my $gff3_file = $ARGV[0] or die $usage;
my $fasta_db = $ARGV[1] or die $usage;
my $num_genes = $ARGV[2] or die $usage;
my $min_chimeras_per_gene = $ARGV[3] or die $usage;


my $MIN_CHIMERA_PART_LENGTH = 100;

my $fasta_reader = new Fasta_reader($fasta_db);
my %genome = $fasta_reader->retrieve_all_seqs_hash();

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

my @gene_ids = keys %$gene_obj_indexer_href;
foreach my $gene (values %$gene_obj_indexer_href) {
    $gene->delete_isoforms();
}



my %GENES_USED;

@gene_ids = shuffle(@gene_ids);

my $num_genes = scalar (@gene_ids);


my $num_chimeric_genes_made = 0;
while ($num_chimeras_made < $num_chimeras) {
        
    my $gene_id_left = $gene_ids[ int(rand($num_genes)) ];
    my $gene_id_right = $gene_ids[ int(rand($num_genes)) ];

    if ($gene_id_left =~ /-/) { next; }  # avoid the annotated read-thru transcripts.
    if ($gene_id_right =~ /-/) { next; }
    
    if ($gene_id_left eq $gene_id_right) { next; }

        
    if ($GENES_USED{$gene_id_left} || $GENES_USED{$gene_id_right}) { next; }

    print STDERR "[$num_chimeras_made done]  Testing: $gene_id_left vs. $gene_id_right\n";

    if ($gene_id_left eq $gene_id_right) { next; }

    


    my $gene_obj_left = $gene_obj_indexer_href->{$gene_id_left};
    my $left_chr_seq = $genome{$gene_obj_left->{asmbl_id}};

    unless (&all_consensus_intron_dinucs($gene_obj_left, \$left_chr_seq)) { 
        $GENES_USED{$gene_id_left}++;
        next; 
    }

    my $gene_obj_right = $gene_obj_indexer_href->{$gene_id_right};
    my $right_chr_seq = $genome{$gene_obj_right->{asmbl_id}};
    
    unless (&all_consensus_intron_dinucs($gene_obj_right, \$right_chr_seq)) { 
        $GENES_USED{$gene_id_right}++;
        next; 
    }
    
    
    #########################
    ## Left part of chimera:
    #########################
    
    my @left_exons_sampled = &select_exons($gene_obj_left, "left");
    unless (@left_exons_sampled) { next; }
    
    #print "Gene before: " . $gene_obj_left->toString();


    my @right_exons_sampled = &select_exons($gene_obj_right, "right");
    unless (@right_exons_sampled) { next; }


    # Try all combinations.

    my @chim_entries;
    
    for (my $i = 0; $i <= $#left_exons_sampled; $i++) {
        for (my $j = 0; $j <= $#right_exons_sampled; $j++) {

            my @exons_left = @left_exons_sampled[0..$i];
            my @exons_right = @right_exons_sampled[$j..$#right_exons_sampled];

                
            my $gene_obj_left_copy = $gene_obj_left->clone_gene();
            $gene_obj_left_copy->{mRNA_exon_objs} = \@exons_left;
            $gene_obj_left_copy->refine_gene_object();
            #print "Gene after: " . $gene_obj_left_copy->toString();
            my $left_cdna = $gene_obj_left_copy->create_cDNA_sequence(\$left_chr_seq);
            if (length($left_cdna) < $MIN_CHIMERA_PART_LENGTH) { next; }
    
                
            my $gene_obj_right_copy = $gene_obj_right->clone_gene();
            $gene_obj_right_copy->{mRNA_exon_objs} = \@exons_right;
            $gene_obj_right_copy->refine_gene_object();
            my $right_cdna = $gene_obj_right_copy->create_cDNA_sequence(\$right_chr_seq);
            if (length($right_cdna) < $MIN_CHIMERA_PART_LENGTH) { next; }
            
            
            #############################
            ## Construct the chimera
            #############################
            
            my $part_A_coord_string = &construct_coord_string($gene_obj_left_copy);
            my $part_B_coord_string = &construct_coord_string($gene_obj_right_copy);
            
            my $gene_A = $gene_obj_left_copy->{TU_feat_name};
            my $gene_B = $gene_obj_right_copy->{TU_feat_name};
            
            my $trans_A = $gene_obj_left_copy->{Model_feat_name};
            my $trans_B = $gene_obj_right_copy->{Model_feat_name};
            
            my $chim_entry = ">$gene_A--$gene_B $trans_A--$trans_B $part_A_coord_string $part_B_coord_string FusedAt:" . length($left_cdna) . "\n"
                . uc($left_cdna) . lc($right_cdna) . "\n";
     
            push (@chim_entries, $chim_entry);
        }

        if (scalar(@chim_entries) >= $min_chimeras_per_gene) {
            $num_chimeras_made++;
            
            $GENES_USED{$gene_A}++;
            $GENES_USED{$gene_B}++;

            print join("", @chim_entries);
            
        }
    }
}



exit(0);

####
sub construct_coord_string {
    my ($gene_obj) = @_;

    my $chr = $gene_obj->{asmbl_id};

    my $orient = $gene_obj->get_orientation();

    my @exons = $gene_obj->get_exons();
    
    my @coords;
    foreach my $exon (@exons) {
        my ($end5, $end3) = $exon->get_coords();
        push (@coords, "$end5-$end3");
    }
    
    my $coordstring = "$chr:" . join(",", @coords) . "[$orient]";
 
    return($coordstring);
}



####
sub select_exons {
    my ($gene_obj, $end) = @_;

    my @exons = $gene_obj->get_exons();

    if (scalar @exons == 1) {
        return();
    }

    if ($end eq "left") {
        my $last_index = $#exons -1;
        @exons = @exons[0..$last_index];
        return(@exons);
                
    }
    else {
        my $first_index = 1;
        @exons = @exons[$first_index..$#exons];
        return(@exons);
    }
    
    return(@exons);
}

####
sub all_consensus_intron_dinucs {
    my ($gene_obj, $chr_seq_ref) = @_;

    my @intron_coords = $gene_obj->get_intron_coordinates();

    my $orient = $gene_obj->get_orientation();
    
    my @acceptable_bounds = ("GT-AG", "GC-AG");
    

    foreach my $intron (@intron_coords) {
        my ($lend, $rend) = sort {$a<=>$b} @$intron;

        my $intron_seq = substr($$chr_seq_ref, $lend-1, $rend-$lend + 1);

        if ($orient eq '-') {
            $intron_seq = &reverse_complement($intron_seq);
        }

        my $donor = substr($intron_seq, 0, 2);
        my $acceptor = substr($intron_seq, -2);

        my $splice_pair = uc "$donor-$acceptor";
        
        print STDERR "$intron\t$splice_pair\n";
        
        if (! grep { $_ eq $splice_pair } @acceptable_bounds) {
            return(0);
        }
                
    }
    return(1); # all good
}
