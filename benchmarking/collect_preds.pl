#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


my $usage = "usage: $0 fusion_result_file_listing.dat\n\n";

my $fusion_result_file_listing = $ARGV[0] or die $usage;

my $fusion_prog_parser_lib_dir = "$FindBin::Bin/FusionProgParsers";

my %prog_type_to_file_parser = ( 
    'ChimPipe' => 'ChimPipe_parser',
    'ChimeraScan' => 'ChimeraScan_parser',
    'deFuse' => 'DEFUSE_parser',
    'EricScript' => 'EricScript_parser',

    'FusionCatcher' => 'FusionCatcher_parser',
    'FUSION_CATCHER_V0997b' => 'FusionCatcher_parser',
    'FUSION_CATCHER_V0997b_KP' => 'FusionCatcher_KP_parser',
    'FUSION_CATCHER_V0997b_RF1' => 'FusionCatcher_parser',
    'FUSION_CATCHER_V0997b_RF1KP' => 'FusionCatcher_KP_parser',
    'FC_V0997c' => 'FusionCatcher_parser',
    'FC_V0997b' => 'FusionCatcher_parser',
    'FC_V0994e' => 'FusionCatcher_parser',
        
    
    'FusionHunter' => 'FusionHunter_parser',

    'FusionInspector' => 'FusionInspector_parser',
    'starFusionGRCh37_FI' => 'FusionInspector_parser',
    'starFusionGRCh38_FI' => 'FusionInspector_parser',
    
    
    'InFusion' => 'InFusion_parser', 

    'JAFFA-Assembly' => 'JAFFA_parser',
    'JAFFA-Direct' => 'JAFFA_parser',
    'JAFFA-Hybrid' => 'JAFFA_parser',

    'MapSplice' => 'MapSplice_parser',

    'nFuse' => 'NFuse_parser',

    'PRADA' => 'PRADA_parser',

    'SOAP-fuse' => 'SOAPfuse_parser',

    'STAR-Fusion' => 'STARFusion_parser',
    'STAR-Fusion_GR38' => 'STARFusion_parser',
    'STAR_FUSION_GRCh38v24' => 'STARFusion_parser',
    'starFusionGRCh37_gencode_v19_FL4' => 'STARFusion_parser',
    'starFusionGRCh38_gencode_v26' => 'STARFusion_parser',

    
    'TopHat-Fusion' => 'TopHatFusion_parser',
    );

foreach my $module (values %prog_type_to_file_parser) {
    my $module_path = "$fusion_prog_parser_lib_dir/$module.pm";

    require($module_path);

}


main: {


    # print header
    print join("\t", "sample", "prog", "fusion", "J", "S") . "\n";
    
    open(my $fh, $fusion_result_file_listing) or die "Error, cannot open file $fusion_result_file_listing";
    while (<$fh>) {
        chomp;
        my ($sample_name, $prog_name, $result_file) = split(/\t/);


        unless (exists $prog_type_to_file_parser{$prog_name}) {
            die "Error, no parser for prog [$prog_name] ";
        }
        
        my $parser_function = $prog_type_to_file_parser{$prog_name} . "::" . "parse_fusion_result_file";
                
        no strict 'refs';
        my @fusions = &$parser_function($result_file);

        &add_sum_fusions(\@fusions);
        
        @fusions = reverse sort { $a->{sum_frags} <=> $b->{sum_frags} } @fusions;
        
        foreach my $fusion (@fusions) {

            my $fusion_name = join("--", $fusion->{geneA}, $fusion->{geneB});

            my $junc_count = $fusion->{junc_reads};
            my $span_count = $fusion->{span_reads};


            print join("\t", $sample_name, $prog_name, $fusion_name, $junc_count, $span_count) . "\n";
        }
                    
    }
    close $fh;


    exit(0);
}

####
sub add_sum_fusions {
    my ($fusions_aref) = @_;

    foreach my $fusion (@$fusions_aref) {

        $fusion->{sum_frags} = $fusion->{junc_reads} + $fusion->{span_reads};

    }

}
