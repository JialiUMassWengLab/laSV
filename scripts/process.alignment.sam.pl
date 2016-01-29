#! /usr/bin/perl

use strict;
use List::Util qw(sum);

my $last_contig="";
my $last_count=0;
my $alignment="";

open (input, "<$ARGV[0].sam") or die "Can't open $ARGV[0].sam since $!\n";
open (output, ">>$ARGV[0].split.sam") or die "Can't open $ARGV[0].split.sam since $!\n";
while (my $line=<input>) {

    if ($line =~ /^Contig/) {
	chomp($line);
	my @a=split(/\t/, $line);
	
	if ($a[5] ne "*") {
	    if ($a[0] ne $last_contig) {
		if ($last_count > 1) {
		    print output "$alignment";
		}
		$last_contig=$a[0];
		$last_count=1;
		$alignment="$line\n";
	    }
	    else {
		$last_count++;
		$alignment .= "$line\n";
	    }
	}
    }
    else {
	print output "$line";
    }
}	
close input;
close output;

system("samtools view -bS $ARGV[0].split.sam > $ARGV[0].split.bam");
system("bedtools bamtobed -ed -i $ARGV[0].split.bam > tmp1.bed");
system("bedtools bamtobed -cigar -i $ARGV[0].split.bam > tmp2.bed");
my %cigar=();
open (input, "<tmp2.bed") or die "Can't open tmp2.bed since $!\n";
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);
    my $align="$a[0]\:$a[1]\:$a[2]\:$a[3]";
    $cigar{$align}=$a[6];
}
close input;

open (input, "<tmp1.bed") or die "Can't open tmp1.bed since $!\n";
open (output, ">$ARGV[0].discordant.bed") or die "Can't open $ARGV[0].discordant.bed since $!\n";
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);
    my $align="$a[0]\:$a[1]\:$a[2]\:$a[3]";
    my (@cigar_m)=$cigar{$align}=~/(\d+)M/g;
    my (@cigar_i)=$cigar{$align}=~/(\d+)I/g;
    if ($a[4]/sum(@cigar_m,@cigar_i) < 0.04) {
	print output "$line\t$cigar{$align}\n";
    }
}
close input;
close output;

system("rm tmp?.bed");
system("rm $ARGV[0].split.*");

my %filter=();
open (input, "<$ARGV[0].sam") or die "Can't open $ARGV[0].sam since $!\n";
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\s+/, $line);
    my $as=0;
    my $xs=0;
    for my $i (11..$#a) {
        if ($a[$i] =~ /^AS:i:/) {
            $a[$i] =~ s/AS:i://;
            $as=$a[$i];
        }
        elsif ($a[$i] =~ /^XS:i:/) {
            $a[$i] =~ s/XS:i://;
            $xs=$a[$i];
        }
        if (($xs > 0) && ($as-$xs <= 20)) {
            $a[3]--;
            my $id="$a[0]\;$a[2]\;$a[3]";
            $filter{$id}=1;
        }
    }
}
close input;

open (input, "<$ARGV[0].discordant.bed") or die "Can't open $ARGV[0].discordant.bed since $!\n";
open (output, ">$ARGV[0].discordant.filt.bed") or die "Can't open $ARGV[0].discordant.filt.bed since $!\n";
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);
    my $id="$a[3]\;$a[0]\;$a[1]";
    if (! defined $filter{$id}) {
	print output "$line\n";
    }
}
close input;
close output;

#system("bedtools intersect -a $ARGV[0].discordant.filt.bed -b $ARGV[1] -wo -f 0.8 > TEoverlap.bed");
#system("awk -F \"\\t\" '{OFS=\"\\t\"; if ((\$3-\$2-\$15 <= 10) && (\$12 <= 10)) print \$1,\$2,\$3,\$4,\$5,\$6,\$7}' TEoverlap.bed > remove.bed");
#system("bedtools subtract -a $ARGV[0].discordant.filt.bed -b remove.bed -f 1 > $ARGV[0].discordant.bed");
system("mv $ARGV[0].discordant.filt.bed $ARGV[0].discordant.bed");
#system("rm $ARGV[0].discordant.filt.bed TEoverlap.bed remove.bed");
