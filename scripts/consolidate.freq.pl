#! /usr/bin/perl

use strict;

my %var_num=();
my %ref_num=();
my %var_len=();
my %ref_len=();
my %contigs=();
open (input, "<$ARGV[0]") or die "Can't open $ARGV[0] since $!\n";
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);
    my @b=split(/\:/, $a[0]);
    my $strand1=substr($b[2], length($b[2])-1, 1);
    my $strand2=substr($b[4], length($b[4])-1, 1);
#    chop($b[2]); chop($b[4]);
    my $oppo1=$b[2];
    my $oppo2=$b[4];
    if ($strand1 eq "+") {$oppo1 =~ s/\+/\-/;}
    else {$oppo1 =~ s/\-/\+/;}
    if ($strand2 eq "+") {$oppo2 =~ s/\+/\-/;}
    else {$oppo2 =~ s/\-/\+/;}
#    my $orien="s";
#    if ($strand1 ne $strand2) {$orien="o";}
#    my $id1="$b[1]\:$b[2]\:$b[3]\:$b[4]\:$orien";
#    my $id2="$b[3]\:$b[4]\:$b[1]\:$b[2]\:$orien";    

    my $id1="$b[1]\:$b[2]\:$b[3]\:$b[4]";
    my $id2="$b[3]\:$b[4]\:$b[1]\:$b[2]";    
    my $id3="$b[1]\:$oppo1\:$b[3]\:$oppo2";
    my $id4="$b[3]\:$oppo2\:$b[1]\:$oppo1";
    if (defined $contigs{$id1}) {
	$contigs{$id1} .= "\;$b[0]";
	$var_num{$id1} += $a[1];
	$ref_num{$id1} += $a[3];
    }
    elsif (defined $contigs{$id2}) {
	$contigs{$id2} .= "\;$b[0]";
	$var_num{$id2} += $a[1];
	$ref_num{$id2} += $a[3];
    }
    elsif (defined $contigs{$id3}) {
        $contigs{$id3} .= "\;$b[0]";
        $var_num{$id3} += $a[1];
        $ref_num{$id3} += $a[3];
    }
    elsif (defined $contigs{$id4}) {
        $contigs{$id4} .= "\;$b[0]";
        $var_num{$id4} += $a[1];
        $ref_num{$id4} += $a[3];
    }
    else {
	$contigs{$id1} = $b[0];
        $var_num{$id1} = $a[1];
        $var_len{$id1} = $a[2];
        $ref_num{$id1} = $a[3];
        $ref_len{$id1} = $a[4];
    }
}
close input;

open (output, ">>$ARGV[0].consolidate") or die "Can't open $ARGV[0].consolidate since $!\n";
while ((my $key, my $value) = each (%contigs)) {
    my $sv_sup_norm=sprintf("%.6f", $var_num{$key}/$var_len{$key});
    my $ref_sup_norm=sprintf("%.6f", $ref_num{$key}/$ref_len{$key});
    my $ratio="0.00000";
    if ($sv_sup_norm + $ref_sup_norm > 0) {$ratio=sprintf("%.5f", $sv_sup_norm/($sv_sup_norm + $ref_sup_norm));}

    my @c=split(/\;/, $value);
    print output "$c[0]\:$key\t$var_num{$key}\t$var_len{$key}\t$ref_num{$key}\t$ref_len{$key}\t$ratio\n";
}
close output;

