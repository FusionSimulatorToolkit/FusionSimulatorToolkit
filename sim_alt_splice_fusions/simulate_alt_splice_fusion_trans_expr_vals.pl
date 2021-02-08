#!/usr/bin/env perl

use strict;
use warnings;

use List::Util qw(max);

my $usage = "\n\n\tusage: $0 RSEM.isoforms.results chim_targets.fasta >  target.forSimulation.RSEM.isoforms.results\n\n";

my $template_isoforms_RSEM = $ARGV[0] or die $usage;
my $chim_targets_fasta = $ARGV[1] or die $usage;

my %involved_in_chim;
my %trans_to_gene_id;
{
    open (my $fh, $chim_targets_fasta) or die $!;
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ />(\S+)--(\S+)/) {
            my ($geneA, $geneB) = ($1, $2);
            
            my @x = split(/\s+/, $line);
            my $acc = shift @x;
            $acc =~ s/>//;
            
            $involved_in_chim{$geneA} = 1;
            $involved_in_chim{$geneB} = 1;
            $involved_in_chim{$acc} = 1;
            
        }
        else {
            # regular ref annotation transcript
            my ($trans_id, $gene_id, $gene_symbol) = split(/\s+/, $line);
            $trans_to_gene_id{$trans_id} = $gene_id;
        }
        
    }
    close $fh;
        
}





## Read in existing TPM data:

my %isoform_to_tpm;

my $sum_tpm = 0;
my %gene_to_tpm;

{
    open (my $fh, $template_isoforms_RSEM) or die $!;
    my $header = <$fh>;
    {
        # just to be safe.
        my @x = split(/\t/, $header);
        unless ($x[5] eq "TPM") {
            die "Error, not recognzing rsem data table formatting / header";
        }
    }
    while (<$fh>) {       
        chomp;
        my @x = split(/\t/);
        my $acc = $x[0];
        my $tpm = $x[5];

        $isoform_to_tpm{$acc} = $tpm;
        
        $sum_tpm += $tpm;

        if (my $gene_id = $trans_to_gene_id{$acc}) {
            $gene_to_tpm{$gene_id} += $tpm;
        }
    }
    close $fh;
}


## reset the expression values for fusions and fusion partners

my $MIN_TPM = 1;

{
    open (my $fh, $template_isoforms_RSEM) or die $!;
    my $header = <$fh>;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $acc = $x[0];
        my $tpm = $x[5];
        my $gene_id = $trans_to_gene_id{$acc};
        
        if ($involved_in_chim{$acc} || ($gene_id && exists $involved_in_chim{$gene_id}) ) {
            
            while ($tpm < $MIN_TPM) {
                $tpm = 2**(rand(10));
            }
            
            print STDERR "-chim_related: $acc => tpm: $tpm\n";
        }
        $isoform_to_tpm{$acc} = $tpm;
        
        $sum_tpm += $tpm;
    }
    close $fh;
}

open (my $fh, $template_isoforms_RSEM) or die $!;
my $header = <$fh>;
print $header;
while (<$fh>) {
    chomp;
    my @x = split(/\t/);
    my $trans = $x[0];
    
    my $tpm = &get_TPM($trans, $sum_tpm);
    
    $x[5] = $tpm;

    print join("\t", @x) . "\n";
}

close $fh;

exit(0);


####
sub get_TPM {
    my ($trans, $sum_tpm) = @_;
    
    my $tpm = $isoform_to_tpm{$trans};
    if (defined $tpm) {
        $tpm = $tpm/$sum_tpm * 1e6;
        
        return($tpm);
    }
    else {
        die "ERROR: no tmp found for $trans  ";
    }
    
}

