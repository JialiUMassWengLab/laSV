#! /usr/bin/perl

use strict;

my %closeup=();
my %remove=();
open (input, "<$ARGV[0].breakpoints.fa") or die "Can't open $ARGV[0].breakpoints.fa since $!\n";
open (output, ">>tmp.bed") or die "Can't open tmp.bed since $!\n";
while (my $line=<input>) {
    chomp($line);
    if ($line =~ /^>/) {
	$line =~ s/\>//;
	my @a=split(/\:/, $line);
	chop($a[2]); chop($a[4]);
	my $lower1=$a[2]-1;
	my $lower2=$a[4]-1;
	if (($lower1 < 0) || ($lower2 < 0)) {
	    $line =~ s/\:Sample//;
	    $remove{$line}=1;
	    next;
	}

	if (($a[1] eq $a[3]) && (abs($a[2]-$a[4]) < 500000)) {$closeup{$line}=1;}

	print output "$a[1]\t$lower1\t$a[2]\t$line\|BP1\n";
	print output "$a[3]\t$lower2\t$a[4]\t$line\|BP2\n";
    }
}
close input;
close output;

my %family1=(); my %family2=();
my %pos1=(); my %pos2=();
my %flag=();
system("bedtools intersect -a tmp.bed -b $ARGV[1] -wo > intersect.bed");
open (input, "<intersect.bed") or die "Can't open intersect.bed since $!\n";
while (my $line=<input>) {
    chomp($line);
    my @b=split(/\t/, $line);
    my @c=split(/\|/, $b[3]);
    if ($c[1] eq "BP1") {
	$family1{$c[0]}=$b[10];
	$pos1{$c[0]}="$b[4]\:$b[5]\:$b[6]";
    }
    elsif ($c[1] eq "BP2") {
	$family2{$c[0]}=$b[10];
	$pos2{$c[0]}="$b[4]\:$b[5]\:$b[6]";
    }
    if ((abs($b[5]-$b[2]) > 20) && (abs($b[6]-$b[2]) > 19)) {$flag{$c[0]}=1;}
}
close input;

while ((my $key, my $value) = each (%family1)) {
    if (($value eq $family2{$key}) && ($pos1{$key} ne $pos2{$key}) && ($flag{$key} == 1) && (! defined $closeup{$key})) {
	$key =~ s/\:Sample//;
	my @x=split(/\:/, $pos1{$key});
	my @y=split(/\:/, $pos2{$key});
	$remove{$key}=1;
    }
}

system("rm intersect.bed tmp.bed");

open (input, "<$ARGV[0].breakpoints.fa") or die "Can't open $ARGV[0].breakpoints.fa since $!\n";
open (output, ">>$ARGV[0].breakpoints.fa.filt") or die "Can't open $ARGV[0].breakpoints.fa.filt since $!\n";
while (my $line=<input>) {
    chomp($line);
    if ($line =~ /^>/) {
        $line =~ s/\>//;
        $line =~ s/\:Sample//;
        if (! defined $remove{$line}) {
            print output "\>$line\:Sample\n";
            $line=<input>;
            print output "$line";
        }
    }
}
close input;
close output;
system("mv $ARGV[0].breakpoints.fa.filt $ARGV[0].breakpoints.fa");

open (input, "<$ARGV[0].reference1.fa") or die "Can't open $ARGV[0].reference1.fa since $!\n";
open (output, ">>$ARGV[0].reference1.fa.filt") or die "Can't open $ARGV[0].reference1.fa.filt since $!\n";
while (my $line=<input>) {
    chomp($line);
    if ($line =~ /^>/) {
        $line =~ s/\>//;
	$line =~ s/\:Reference1//;
	if (! defined $remove{$line}) {
	    print output "\>$line\:Reference1\n";
            $line=<input>;
            print output "$line";
	}
    }
}
close input;
close output;
system("mv $ARGV[0].reference1.fa.filt $ARGV[0].reference1.fa");

open (input, "<$ARGV[0].reference2.fa") or die "Can't open $ARGV[0].reference2.fa since $!\n";
open (output, ">>$ARGV[0].reference2.fa.filt") or die "Can't open $ARGV[0].reference2.fa.filt since $!\n";
while (my $line=<input>) {
    chomp($line);
    if ($line =~ /^>/) {
	$line =~ s/\>//;
	$line =~ s/\:Reference2//;
	if (! defined $remove{$line}) {
            print output "\>$line\:Reference2\n";
            $line=<input>;
            print output "$line";
	}
    }
}
close input;
close output;
system("mv $ARGV[0].reference2.fa.filt $ARGV[0].reference2.fa");

