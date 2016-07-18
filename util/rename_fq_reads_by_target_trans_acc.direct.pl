#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../PerlLib";
use Fastq_reader;


my $usage = "\n\n\tusage: $0 combined.transcripts.fasta.RSEM.idx.fa fileA.fastq [fileB.fastq ....]\n\n";

my $idx_fa = $ARGV[0] or die $usage;
shift @ARGV;

my @fq_files = @ARGV;
unless (@fq_files) { die $usage; }

main: {

    my @trans = ('random'); # index zero reserved for random sequence.
    
    open (my $fh, $idx_fa) or die "Error, cannot open file $idx_fa";
    while (<$fh>) {
        if (/^>(\S+)/) {
            my $trans_name = $1;
            push (@trans, $trans_name);
        }
    }
    close $fh;
    
    foreach my $fq_file (@fq_files) {
        my $fastq_reader = new Fastq_reader($fq_file);
        
        my $renamed_fq_file = "$fq_file.renamed.fq";
        open (my $ofh, ">$renamed_fq_file") or die "Error, cannot write to $renamed_fq_file";
        

        while (my $fq_entry = $fastq_reader->next()) {
            my $fq_record = $fq_entry->get_fastq_record();
            my @lines = split(/\n/, $fq_record);
            my $acc_line = $lines[0];
            $acc_line =~ s/^\@//;
            my @pts = split(/_/, $acc_line);
            my $trans_idx = $pts[2];
            my $gene_name = $trans[$trans_idx];
            
            my $new_read_name = "@" . "$gene_name:$acc_line";
            $lines[0] = $new_read_name;

            print $ofh join("\n", @lines) . "\n";
        }

        close $ofh;

    }
            
    exit(0);
    
}

