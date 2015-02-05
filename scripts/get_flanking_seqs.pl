#! /usr/bin/perl

use strict;
use List::Util qw(sum);
use Bio::Seq;

my %ends=();
open (input, "<$ARGV[1]") or die "Can't open $ARGV[1] since $!\n";
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);
    $ends{$a[0]}=$a[1];
}
close input;

my %contigs=();
my %coor=();
open (input, "<$ARGV[0].aligned.discordant.bed") or die "Can't open $ARGV[0].aligned.discordant.bed since $!\n";
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);
    $a[6] =~ s/H/S/g;
    my (@cigar_s)=$a[6]=~/(\d+)S/g;
    my (@cigar_m)=$a[6]=~/(\d+)M/g;
    my (@cigar_i)=$a[6]=~/(\d+)I/g;
    my (@cigar_d)=$a[6]=~/(\d+)D/g;
    $contigs{$a[3]}=1;

    my $chr_start=$a[1];
    my $start=0;
    if ($a[5] eq "+") {
        if (($#cigar_s == 0) && ($a[6] =~ /S$/)) {$start=0;}
        else {$start=$cigar_s[0];}
    }
    else {
        if (($#cigar_s == 0) && ($a[6] !~ /S$/)) {$start=0;}
        else {$start=$cigar_s[$#cigar_s];}
    }

    my $end=$start+sum(@cigar_m,@cigar_i);
    my $chr_end=$chr_start+sum(@cigar_m,@cigar_d);

    #coordinates are 0-based
    my $pos="$start\:$end\:$a[0]\:$chr_start\:$chr_end\:$a[5]";
    if (defined $coor{$a[3]}) {
        $coor{$a[3]}.=";$pos";
    }
    else {$coor{$a[3]}=$pos;}
}
close input;

my %seq=();
open (input, "<$ARGV[0].contigs.fa") or die "Can't open $ARGV[0].contigs.fa since $!\n";
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);
    if ($a[0] =~ />Contig/) {
	$a[0] =~ s/>//;
	if ($contigs{$a[0]} == 1) {
	    $seq{$a[0]}=<input>;
	    chomp($seq{$a[0]});
	}
    }
}
close input;
    
open (output1, ">>$ARGV[0].breakpoints.fa") or die "Can't open $ARGV[0]_seqs_for_freq_estimate.fa since $!\n";
open (output2, ">>$ARGV[0].reference1.fa") or die "Can't open $ARGV[0]_seqs_for_freq_estimate.fa since $!\n";
open (output3, ">>$ARGV[0].reference2.fa") or die "Can't open $ARGV[0]_seqs_for_freq_estimate.fa since $!\n";
while ((my $key, my $value) = each (%coor)) {
    my @a=split(/\;/, $value);
    my %frags=();
    foreach my $frag (@a) {
	my @b=split(/\:/, $frag);
	$frags{$b[0]}=$frag;
    }

    my @sorted_keys;
    foreach my $key1 (sort {$a <=> $b} keys %frags) {
	push (@sorted_keys, $key1);
    }

    my $lefts=""; my $rights=""; my $mid="";
    my $refb1=""; my $refb2=""; my $refc1=""; my $refc2="";
    my $seq_start=0; my $seq_end=0;
    for my $i (0..$#sorted_keys-1) {
	my @b=split(/\:/, $frags{$sorted_keys[$i]});
	my @c=split(/\:/, $frags{$sorted_keys[$i+1]});

	if (($b[2] eq "chrU") && ($c[2] eq "chrU")) {next;}
	if (($b[2] eq $c[2]) && ($b[5] eq "+") && ($c[5] eq "+") && ($b[3] < $c[3]) && (abs($c[3]-$b[4]+$b[1]-$c[0]) <= 150)) {next;}
	elsif (($b[2] eq $c[2]) && ($b[5] eq "-") && ($c[5] eq "-") && ($c[3] < $b[3]) && (abs($b[3]-$c[4]+$b[1]-$c[0]) <= 150)) {next;}

	$lefts=""; $rights=""; $mid=""; $refb1=""; $refb2=""; $refc1=""; $refc2="";
	my $bpb=$b[4]; my $bpc=$c[3];
	if ($b[5] eq "-") {$bpb=$b[3];}
	if ($c[5] eq "-") {$bpc=$c[4];}
	#lb & rp are 1-based
	my $lb=$c[0]+1; my $rb=$b[1];
	if ($c[0] >= $b[1]) {$lb=$b[1]; $rb=$c[0]+1;}
	$mid=substr($seq{$key}, $lb-1, $rb-$lb+1);
	
	$rights=substr($seq{$key}, $rb, $c[1]-$rb);
	if (length($rights) < 500) {
	    if ($c[5] eq "+") {
		$seq_start=$c[4];
		$seq_end=$c[4]-length($rights)+500;
	    }
	    else {
		$seq_start=$c[3]+length($rights)-500;
		$seq_end=$c[3];
	    }
	    if ($seq_start < 0) {$seq_start=0;}
	    if ($seq_end > $ends{$c[2]}) {$seq_end=$ends{$c[2]};}
	    if ($seq_end > $seq_start) {
		system("twoBitToFa -seq=$c[2] -start=$seq_start -end=$seq_end $ARGV[2] tmp_seq.fa");
		open (input, "<tmp_seq.fa") or die "Can't open tmp_seq.fa right rights since $!\n";
		my $header=<input>;
		my $curr_seq="";
		while (my $line=<input>) {
		    chomp($line);
		    $line=uc($line);
		    $curr_seq .= $line;
		}
		close input;
		system("rm tmp_seq.fa");
		if ($c[5] eq "-") {
		    my $seq=Bio::Seq->new(-seq=>$curr_seq);
		    $curr_seq=$seq->revcom->seq;
		}
		$rights .= $curr_seq;
	    }
	}
	else {
	    $rights=substr($rights, 0, 500);
	}
	
	$lefts=substr($seq{$key}, $b[0], $lb-$b[0]-1);
	if (length($lefts) < 500) {
	    if ($b[5] eq "+") {
		$seq_start=$b[3]+length($lefts)-500;
		$seq_end=$b[3];
	    }
	    else {
                $seq_start=$b[4];
                $seq_end=$b[4]-length($lefts)+500;
	    }		
	    if ($seq_start < 0) {$seq_start=0;}
	    if ($seq_end > $ends{$b[2]}) {$seq_end=$ends{$b[2]};}
	    if ($seq_end > $seq_start) {
		system("twoBitToFa -seq=$b[2] -start=$seq_start -end=$seq_end $ARGV[2] tmp_seq.fa");
		open (input, "<tmp_seq.fa") or die "Can't open tmp_seq.fa right lefts since $!\n";
		my $header=<input>;
		my $curr_seq="";
		while (my $line=<input>) {
		    chomp($line);
		    $line=uc($line);
		    $curr_seq .= $line;
		}
		close input;
		system("rm tmp_seq.fa");
		if ($b[5] eq "-") {
		    my $seq=Bio::Seq->new(-seq=>$curr_seq);
		    $curr_seq=$seq->revcom->seq;
		}
		$lefts=$curr_seq.$lefts;
	    }
	}
	else {
	    $lefts=substr($lefts, length($lefts)-500, 500);
	}
	
	print output1 ">$key\:$b[2]\:$bpb$b[5]\:$c[2]\:$bpc$c[5]\:Sample\n";
	print output1 "$lefts$mid$rights\n";

	#Reference sequences
	if ($c[0] < $b[1]) {
	    my $tmp1=$mid.$rights;
	    my $tmp2=$lefts.$mid;
	    $refc1=substr($tmp1, 0, 501);
	    $refb1=substr($tmp2, length($tmp2)-501, 501);
	}
	else {
	    $refc1=substr($seq{$key}, $rb-1, 1).$rights;
	    $refb1=$lefts.substr($seq{$key}, $lb-1, 1);
	}
	#refc2
	if ($c[5] eq "+") {
	    $seq_start=$c[3]-500;
	    $seq_end=$c[3];
	}
	else {
	    $seq_start=$c[4];
	    $seq_end=$c[4]+500;
	}
	if ($seq_start < 0) {$seq_start=0;}
	if ($seq_end > $ends{$c[2]}) {$seq_end=$ends{$c[2]};}
	if ($seq_end > $seq_start) {
	    system("twoBitToFa -seq=$c[2] -start=$seq_start -end=$seq_end $ARGV[2] tmp_seq.fa");
	    open (input, "<tmp_seq.fa") or die "Can't open tmp_seq.fa leftr1 since $!\n";
	    my $header=<input>;
	    while (my $line=<input>) {
		chomp($line);
		$line=uc($line);
		$refc2 .= $line;
	    }
	    close input;
	    system("rm tmp_seq.fa");
	    if ($c[5] eq "-") {
		my $seq=Bio::Seq->new(-seq=>$refc2);
		$refc2=$seq->revcom->seq;
	    }
	}
	#refb2
	if ($b[5] eq "+") {
	    $seq_start=$b[4];
	    $seq_end=$b[4]+500;
	}
	else {
	    $seq_start=$b[3]-500;
	    $seq_end=$b[3];
	}
	if ($seq_start < 0) {$seq_start=0;}
	if ($seq_end > $ends{$b[2]}) {$seq_end=$ends{$b[2]};}
	if ($seq_end > $seq_start) {
	    system("twoBitToFa -seq=$b[2] -start=$seq_start -end=$seq_end $ARGV[2] tmp_seq.fa");
	    open (input, "<tmp_seq.fa") or die "Can't open tmp_seq.fa rightr2 since $!\n";
	    my $header=<input>;
	    while (my $line=<input>) {
		chomp($line);
		$line=uc($line);
		$refb2 .= $line;
	    }
	    close input;
	    system("rm tmp_seq.fa");
	    if ($b[5] eq "-") {
		my $seq=Bio::Seq->new(-seq=>$refb2);
		$refb2=$seq->revcom->seq;
	    }
	}
	#print reference
	if ($c[6] ne "T") {
	    print output2 ">$key\:$b[2]\:$bpb$b[5]\:$c[2]\:$bpc$c[5]\:Reference1\n";
	    print output2 "$refc2$refc1\n";
	}
	if ($b[6] ne "T") {
	    print output3 ">$key\:$b[2]\:$bpb$b[5]\:$c[2]\:$bpc$c[5]\:Reference2\n";
	    print output3 "$refb1$refb2\n";
	}

    }
    
}
close output1;
close output2;
close output3;
