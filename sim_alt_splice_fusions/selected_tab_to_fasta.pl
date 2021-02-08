#!/usr/bin/env perl

use strict;
use warnings;

while(<>) {
    chomp;
    my ($acc, $lend_length, $rend_length, $header_info, $sequence) = split(/\t/);
    
    print ">$acc^$lend_length^$rend_length\n$sequence\n";
}


