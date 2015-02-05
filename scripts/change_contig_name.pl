#! /usr/bin/perl

use strict;

open (input, "<$ARGV[0]") or die "Can't open $ARGV[0] since $!\n";
open (output, ">>tmp.fa") or die "Can't open tmp.fa since $!\n";
while (my $line=<input>) {
    if ($line =~ /^>/) {
	my @a=split(/\-/, $line);
	$a[0] =~ s/^>//;
	print output "\>Contig$a[0]\n";
    }
    else {
	print output $line;
    }
}
close input;
close output;

system("mv tmp.fa $ARGV[0]");
