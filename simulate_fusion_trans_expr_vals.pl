#!/usr/bin/env perl

use strict;
use warnings;

use List::Util qw(max);

my $usage = "\n\n\tusage: $0 RSEM.isoforms.results chim_targets.fasta >  target.forSimulation.RSEM.isoforms.results\n\n";

my $template_isoforms_RSEM = $ARGV[0] or die $usage;
my $chim_targets_fasta = $ARGV[1] or die $usage;

my %involved_in_chim;
{
    open (my $fh, $chim_targets_fasta) or die $!;
    while (<$fh>) {
        chomp;
        if (/>(\S+)--(\S+) (\S+)--(\S+)/) {
            my ($geneA, $geneB, $transA, $transB) = ($1, $2, $3, $4);
            
            my @x = split(/\s+/);
            my $acc = shift @x;
            $acc =~ s/>//;
            
            $involved_in_chim{$transA} = 1;
            $involved_in_chim{$transB} = 1;
            $involved_in_chim{$acc} = 1;
            
        }
    }
    close $fh;
        
}


## Read in existing TPM data:
## reset the expression values for fusions and fusion partners

my %isoform_to_tpm;

my $sum_tpm = 0;
{
    open (my $fh, $template_isoforms_RSEM) or die $!;
    my $header = <$fh>;
    while (<$fh>) {
        
        chomp;
        my @x = split(/\t/);
        my $acc = $x[0];
        my $tpm = $x[5];
        if ($involved_in_chim{$acc}) {
            $tpm = 1 + 2**(rand(10));
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

