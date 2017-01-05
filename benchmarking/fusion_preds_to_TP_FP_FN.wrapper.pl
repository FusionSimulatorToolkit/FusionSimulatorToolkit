#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;


use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);

my $usage = <<__EOUSAGE__;



#################################################################################################
#
# Required:
#
#  --truth_fusions <string>   file containing a list of the true fusions.
#
#  --fusion_preds <string>    fusion predictions ranked accordingly.
#
#
# Optional:
#
#  --unsure_fusions <string>  predictions to ignore (as not TP or FP)
#
#  --allow_reverse_fusion     if true fusion is A--B, allow for it to be reported as B--A
#
#  --allow_paralogs         
#
##################################################################################################


__EOUSAGE__


    ;



my $help_flag;
my $fusion_preds_file;
my $truth_fusions_file;
my $unsure_fusions_file;

my $ALLOW_REVERSE_FUSION = 0;

my $ALLOW_PARALOGS = 0;

&GetOptions ( 'h' => \$help_flag,
              
              'fusion_preds=s' => \$fusion_preds_file,
              'truth_fusions=s' => \$truth_fusions_file,
              'unsure_fusions=s' => \$unsure_fusions_file,
              'allow_reverse_fusion' => \$ALLOW_REVERSE_FUSION,
              
              'allow_paralogs' => \$ALLOW_PARALOGS,
    );


if ($help_flag) { die $usage; }

unless ($fusion_preds_file && $truth_fusions_file) {
    die $usage;
}


my $gene_spans_file = join(",", ("/seq/regev_genome_portal/RESOURCES/CTAT_GENOME_LIB/GRCh37_gencode_v19_FL3/ref_annot.gtf.gene_spans",
                                 "/seq/RNASEQ/TOOLS/nFUSE/nfuse-0.2.1/NFUSE_DATA_DIR/Homo_sapiens.GRCh37.62.gtf.gene_spans", # nFuse
                                 "/seq/RNASEQ/TOOLS/ERICSCRIPT/ericscript-0.5.3/lib/data/homo_sapiens/__liftover_to_hg19/ensembl_v83.hg19.genespans.txt", #Ericscript
                                 "/seq/regev_genome_portal/RESOURCES/human/tophat_fusion_resources/__downloads_for_mapping_annots/Homo_sapiens.GRCh37.61.gtf.Chr.gene_spans", #tophatFusion
                                 "/seq/RNASEQ/TOOLS/DEFUSE-data/Homo_sapiens.GRCh37.69.gtf.Chr.gene_spans", #Defuse
                                 "/seq/RNASEQ/TOOLS/FusionHunter/AnnotationFiles_hg19/hg19.ucscKnownGene.gene_spans", # fusionhunter
                           )
    );

#my $paralogs_file = "/seq/regev_genome_portal/RESOURCES/CTAT_GENOME_LIB/GRCh37_gencode_v19_FL3/blastn/nuc_clusters.dat";
my $paralogs_file = "/seq/regev_genome_portal/RESOURCES/CTAT_GENOME_LIB/GRCh37_gencode_v19_FL3/blastn/paralog_clusters.dat";


my $cmd = "$FindBin::Bin/fusion_preds_to_TP_FP_FN.pl --fusion_preds $fusion_preds_file --truth_fusions $truth_fusions_file --gene_spans $gene_spans_file";

if ($unsure_fusions_file) {
    $cmd .= " --unsure_fusions $unsure_fusions_file ";
}
if ($ALLOW_REVERSE_FUSION) {
    $cmd .= " --allow_reverse_fusion ";
}
if ($ALLOW_PARALOGS) {
    $cmd .= " --allow_paralogs $paralogs_file";
}


print STDERR "CMD: $cmd\n";

my $ret = system($cmd);

exit($ret);

