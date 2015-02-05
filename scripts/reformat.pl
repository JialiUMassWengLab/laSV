#! /usr/bin/perl

use strict;

open (input, "<contigs.fa") or die "Can't open contigs.fa since $!\n";
#open (output, ">contigs1.fa") or die "Can't open contigs1.fa since $!\n";
my $i=1;
while (my $line=<input>) {
    chomp($line);
    $line =~ s/[^\w>]//g;
    if ($line =~ /^>/) {
	my $line1=<input>;
	if (($line1 ne "") && ($line1 !~ /^>/)) {
	    chomp($line1);
	    $line1 =~ s/[^ACGT]//g;
	    if (length($line1) >= 60) {
		print "$line\n$line1\n";
	    }
	}
    }
    elsif ($line ne "") {
	$line =~ s/[^ACGT]//g;
	if (length($line) >= 60) {
	    print ">ContigX$i\n$line\n";
	    $i++;
	}
    }
}
close input;
#close output;
