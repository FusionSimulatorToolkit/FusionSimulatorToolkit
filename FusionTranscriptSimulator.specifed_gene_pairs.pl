#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/PerlLib");
use Gene_obj;
use Fasta_reader;
use GFF3_utils;
use Carp;
use Nuc_translator;
use List::Util qw(shuffle);
use Data::Dumper;

my $usage = "\n\nusage: $0 gff3_file genome_fasta gene_pairs_filen\n";


my $gff3_file = $ARGV[0] or die $usage;
my $fasta_db = $ARGV[1] or die $usage;
my $gene_pairs_file = $ARGV[2] or die $usage;


my $MIN_CHIMERA_PART_LENGTH = 25;

my $fasta_reader = new Fasta_reader($fasta_db);
my %genome = $fasta_reader->retrieve_all_seqs_hash();

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

my @gene_ids = keys %$gene_obj_indexer_href;
foreach my $gene (values %$gene_obj_indexer_href) {
    $gene->delete_isoforms();
}

foreach my $gene_id (@gene_ids) {
    
    my ($gene_name, @rest) = split(/\|/, $gene_id);
    $gene_obj_indexer_href->{$gene_name} = $gene_obj_indexer_href->{$gene_id};
}



open(my $fh, $gene_pairs_file) or die $!;
while(<$fh>) {
    chomp;
    my $fusion_pair = $_;
    my ($gene_id_left, $gene_id_right) = split(/--/, $fusion_pair);
    
    my $gene_obj_left = $gene_obj_indexer_href->{$gene_id_left} or die "Error, no gene found based on id [$gene_id_left]";
    my $left_chr_seq = $genome{$gene_obj_left->{asmbl_id}};

    
    my $gene_obj_right = $gene_obj_indexer_href->{$gene_id_right} or die "Error, no gene found based on id [$gene_id_right]";;
    my $right_chr_seq = $genome{$gene_obj_right->{asmbl_id}};
    

    for (1..10) {
        
    
        #########################
        ## Left part of chimera:
        #########################
        
        my @left_exons_sampled = &select_exons($gene_obj_left, "left");
        unless (@left_exons_sampled) { 
            
            print STDERR " - sorry, no exons sampled from left gene: $gene_id_left\n";
            next; 
        }
        
        #print "Gene before: " . $gene_obj_left->toString();
        
        my $gene_obj_left_copy = $gene_obj_left->clone_gene();
        $gene_obj_left_copy->{mRNA_exon_objs} = \@left_exons_sampled;
        
        $gene_obj_left_copy->refine_gene_object();
        
        #print "Gene after: " . $gene_obj_left_copy->toString();
        
        
        
        my $left_cdna = $gene_obj_left_copy->create_cDNA_sequence(\$left_chr_seq);
        
        if (length($left_cdna) < $MIN_CHIMERA_PART_LENGTH) { 
            print STDERR " - sorry, exons sampled from left: $gene_id_left  too short: " . length($left_cdna) . "\n";
            next; 
        }
        
        #print "Left_cdna: $left_cdna\n";
        
        #############################
        ## Right part of chimera:
        #############################
        
        my @right_exons_sampled = &select_exons($gene_obj_right, "right");
        unless (@right_exons_sampled) { 
            print STDERR " - sorry, no exons sampled from right gene: $gene_id_right\n"; 
            next; 
        }
        
        my $gene_obj_right_copy = $gene_obj_right->clone_gene();
        $gene_obj_right_copy->{mRNA_exon_objs} = \@right_exons_sampled;
        
        $gene_obj_right_copy->refine_gene_object();
        
        my $right_cdna = $gene_obj_right_copy->create_cDNA_sequence(\$right_chr_seq);
        if (length($right_cdna) < $MIN_CHIMERA_PART_LENGTH) { 
            print STDERR " - sorry, exons sampled from right: $gene_id_right  too short: " . length($right_cdna) . "\n";
            next; 
        }
        
        
        #############################
        ## Construct the chimera
        #############################
        
        my $part_A_coord_string = &construct_coord_string($gene_obj_left_copy);
        my $part_B_coord_string = &construct_coord_string($gene_obj_right_copy);
        
        my $gene_A = $gene_obj_left_copy->{TU_feat_name};
        my $gene_B = $gene_obj_right_copy->{TU_feat_name};
        
        my $trans_A = $gene_obj_left_copy->{Model_feat_name};
        my $trans_B = $gene_obj_right_copy->{Model_feat_name};
        
        print ">$fusion_pair $gene_A--$gene_B $trans_A--$trans_B $part_A_coord_string $part_B_coord_string FusedAt:" . length($left_cdna) . "\n"
            . uc($left_cdna) . lc($right_cdna) . "\n";
        
        last;
    }

    print STDERR "-done\n";
    
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

    if (scalar @exons > 1) {
        @exons = &exclude_shorty_terminal_exons(@exons);
    }
    
    if (scalar @exons == 1) {
        # split exon
        my $exon = $exons[0];
        my $exon_len = abs($exon->{end5} - $exon->{end3}) + 1;
        my $half_exon_len = int($exon_len/2);
        if ($end eq 'left') {
            if ($exon->{strand} eq '+') {
                $exon->{end3} = $exon->{end3} - $half_exon_len;
            }
            else {
                $exon->{end3} = $exon->{end3} + $half_exon_len;
            }
        }
        else {
            # right side.  Change end5
            if ($exon->{strand} eq '+') {
                $exon->{end5} = $exon->{end5} + $half_exon_len;
            }
            else {
                $exon->{end5} = $exon->{end5} - $half_exon_len;
            }
        }
        return($exon);
    }
    
    if ($end eq "left") {
        my $last_index = int(rand($#exons));
        @exons = @exons[0..$last_index];
                    
    }
    else {
        my $first_index = 1 + int(rand($#exons));
        @exons = @exons[$first_index..$#exons];
        
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


####
sub exclude_shorty_terminal_exons {
    my @exons = @_;
    
    while (@exons) {
        my $first_exon = $exons[0];
        if (abs($first_exon->{end5} - $first_exon->{end3}) < $MIN_CHIMERA_PART_LENGTH) {
            shift @exons;
        }
        else {
            last;
        }
    }
    
    while (@exons) {
        my $last_exon = $exons[$#exons];
        if (abs($last_exon->{end5} - $last_exon->{end3}) < $MIN_CHIMERA_PART_LENGTH) {
            pop @exons;
        }
        else {
            last;
        }

    }
        
    return(@exons);
}

