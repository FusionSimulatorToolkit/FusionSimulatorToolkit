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
#  --ignore_fusions <string>  predictions to ignore (as not TP or FP)
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
my $ignore_fusions_file;

my $ALLOW_REVERSE_FUSION = 0;

my $ALLOW_PARALOGS = 0;

&GetOptions ( 'h' => \$help_flag,
              
              'fusion_preds=s' => \$fusion_preds_file,
              'truth_fusions=s' => \$truth_fusions_file,
              'ignore_fusions=s' => \$ignore_fusions_file,
              'allow_reverse_fusion' => \$ALLOW_REVERSE_FUSION,
              
              'allow_paralogs=s' => \$ALLOW_PARALOGS,
    );


if ($help_flag) { die $usage; }

unless ($fusion_preds_file && $truth_fusions_file) {
    die $usage;
}


my $gene_spans_file = "/seq/regev_genome_portal/RESOURCES/CTAT_GENOME_LIB/GRCh37_gencode_v19_FL3/ref_annot.gtf.gene_spans";

my $paralogs_file = "/seq/regev_genome_portal/RESOURCES/CTAT_GENOME_LIB/GRCh37_gencode_v19_FL3/blastn/nuc_clusters.dat";


my $cmd = "$FindBin::Bin/fusion_preds_to_TP_FP_FN.pl --fusion_preds $fusion_preds_file --truth_fusions $truth_fusions_file ";

if ($ignore_fusions_file) {
    $cmd .= " --ignore_fusions $ignore_fusions_file ";
}
if ($ALLOW_REVERSE_FUSION) {
    $cmd .= " --allow_reverse_fusion ";
}
if ($ALLOW_PARALOGS) {
    $cmd .= " --allow_paralogs $paralogs_file";
}

my $ret = system($cmd);

exit($ret);

