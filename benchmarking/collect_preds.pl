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
    'FusionHunter' => 'FusionHunter_parser',
    'FusionInspector' => 'FusionInspector_parser',
    'InFusion' => 'InFusion_parser', 
    'JAFFA-Assembly' => 'JAFFA_parser',
    'MapSplice' => 'MapSplice_parser',
    'nFuse' => 'NFuse_parser',
    'PRADA' => 'PRADA_parser',
    'SOAP-fuse' => 'SOAPfuse_parser',
    'STAR-Fusion' => 'STARFusion_parser',
    'TopHat-Fusion' => 'TopHatFusion_parser',
    );

foreach my $module (values %prog_type_to_file_parser) {
    my $module_path = "$fusion_prog_parser_lib_dir/$module.pm";

    require($module_path);

}


main: {

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


