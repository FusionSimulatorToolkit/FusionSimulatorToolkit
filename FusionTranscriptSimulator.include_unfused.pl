#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/PerlLib");
use Gene_obj;
use Fasta_reader;
use GTF_utils;
use Carp;
use Nuc_translator;
use List::Util qw(shuffle);
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use Data::Dumper;

my ($gtf_file, $fasta_db, $num_chimeras, $out_prefix, $help_flag);


if ($ENV{CTAT_GENOME_LIB}) {
    $gtf_file = "$ENV{CTAT_GENOME_LIB}/ref_annot.gtf";
    $fasta_db = "$ENV{CTAT_GENOME_LIB}/ref_genome.fa";
}


my $usage = <<__EOUSAGE;

###################################################################
#
#  ** uses CTAT_GENOME_LIB env var for defaults ** 
#
#    --gtf_file <string>     ref annot gtf file ($gtf_file)
#
#    --ref_genome <string>    ref genome fasta file ($fasta_db)
#
#  reqiored:
#
#   --num_chimeras <int>     number of chimeric transcripts to create
#
#   --out_prefix <string>    output prefix
#
####################################################################

__EOUSAGE

    ;





&GetOptions ( 'h' => \$help_flag,
              'gtf_file=s' => \$gtf_file,
              'ref_genome=s' => \$fasta_db,
              'num_chimeras=i' => \$num_chimeras,
              'out_prefix=s' => \$out_prefix,
    );


if ($help_flag) {
    die $usage;
}

unless ($gtf_file && $fasta_db && $num_chimeras && $out_prefix) {
    die $usage;
}


open(my $out_fasta_ofh, ">$out_prefix.fasta") or die "Error, cannot write to $out_prefix.fasta";
open(my $out_accs_ofh, ">$out_prefix.fusion_list") or die "Error, cannot write to $out_prefix.fusion_list";

my $MIN_CHIMERA_PART_LENGTH = 100;

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GTF_utils::index_GTF_gene_objs($gtf_file, $gene_obj_indexer_href);

my @gene_ids = keys %$gene_obj_indexer_href;
foreach my $gene (values %$gene_obj_indexer_href) {
    $gene->delete_isoforms();
}

@gene_ids = shuffle(@gene_ids);

my $num_genes = scalar (@gene_ids);

## parse genome file
my $fasta_reader = new Fasta_reader($fasta_db);
my %genome = $fasta_reader->retrieve_all_seqs_hash();

################
## make chimeras

my %GENES_USED;

my $num_chimeras_made = 0;
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

    my $gene_obj_right = $gene_obj_indexer_href->{$gene_id_right};
    my $right_chr_seq = $genome{$gene_obj_right->{asmbl_id}};
    
    
    my $min_orig_seqlen = 1000;
    my $left_cdna_orig = $gene_obj_left->create_cDNA_sequence(\$left_chr_seq);
    if (length($left_cdna_orig) < $min_orig_seqlen) {
        next;
    }
    my $right_cdna_orig = $gene_obj_right->create_cDNA_sequence(\$right_chr_seq);
    if (length($right_cdna_orig) < $min_orig_seqlen) {
        next;
    }
    
    unless (&all_consensus_intron_dinucs($gene_obj_left, \$left_chr_seq)) { 
        $GENES_USED{$gene_id_left}++;
        next; 
    }

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
    
    my $gene_obj_left_copy = $gene_obj_left->clone_gene();
    $gene_obj_left_copy->{mRNA_exon_objs} = \@left_exons_sampled;

    $gene_obj_left_copy->refine_gene_object();
    
    #print "Gene after: " . $gene_obj_left_copy->toString();

    
    
    my $left_cdna = $gene_obj_left_copy->create_cDNA_sequence(\$left_chr_seq);
    
    if (length($left_cdna) < $MIN_CHIMERA_PART_LENGTH) { next; }
    
    #print "Left_cdna: $left_cdna\n";

    #############################
    ## Right part of chimera:
    #############################

    my @right_exons_sampled = &select_exons($gene_obj_right, "right");
    unless (@right_exons_sampled) { next; }
    
    my $gene_obj_right_copy = $gene_obj_right->clone_gene();
    $gene_obj_right_copy->{mRNA_exon_objs} = \@right_exons_sampled;
    
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

    my $gene_A_symbol = $gene_obj_left_copy->{com_name};
    $gene_A = $gene_A_symbol . "|" . $gene_A;
    
    my $gene_B_symbol = $gene_obj_right_copy->{com_name};
    $gene_B = $gene_B_symbol . "|" . $gene_B;
    
    
    my $trans_A = $gene_obj_left_copy->{Model_feat_name};
    my $trans_B = $gene_obj_right_copy->{Model_feat_name};

    print $out_fasta_ofh ">$gene_A--$gene_B $trans_A--$trans_B $part_A_coord_string $part_B_coord_string FusedAt:" . length($left_cdna) . "\n"
        . uc($left_cdna) . lc($right_cdna) . "\n";
     
    
    print $out_fasta_ofh ">$gene_A $trans_A brkpt: " . length($left_cdna) . "\n" .
        uc($left_cdna_orig) . "\n";
    
    print $out_fasta_ofh ">$gene_B $trans_B brkpt: " . (length($right_cdna_orig) - length($right_cdna) ) . "\n" .
        uc($right_cdna_orig) . "\n";
    
    #my ($gene_tok_A, $rest) = split(/\|/, $gene_A);
    #my ($gene_tok_B, $rest2) = split(/\|/, $gene_B);
    
    my $gene_tok_A = $gene_A_symbol;
    my $gene_tok_B = $gene_B_symbol;
    
    my $fusion_tok = "$gene_tok_A--$gene_tok_B";
    print $out_accs_ofh "$fusion_tok\n";
    
    $num_chimeras_made++;

    $GENES_USED{$gene_A}++;
    $GENES_USED{$gene_B}++;

    


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
        my $last_index = int(rand($#exons));
        @exons = @exons[0..$last_index];
        
        # require last exon to contain a CDS
        my $last_exon = $exons[$#exons];
        unless ($last_exon->get_CDS_obj()) {
            return();
        }
        
    }
    else {
        my $first_index = 1 + int(rand($#exons));
        @exons = @exons[$first_index..$#exons];
        
        # require first exon to have a CDS:
        my $first_exon = $exons[0];
        unless ($first_exon->get_CDS_obj()) {
            return();
        }
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
