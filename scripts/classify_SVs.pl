#! /usr/bin/perl

use strict;
use List::Util qw(sum);
use Bio::Seq;

my %contigs=();
open (input, "<$ARGV[0].breakpoints.freq") or die "Can't open $ARGV[0].breakpoints.freq.filt since $!\n";
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);
    my @b=split(/\:/, $a[0]);
    my $bp="$b[1]\:$b[2]\:$b[3]\:$b[4]\:$a[5]\:$a[1]\:$a[3]";
    if (defined $contigs{$b[0]}) {
	$contigs{$b[0]} .= ";$bp";
    }
    else {
	$contigs{$b[0]}=$bp;
    }
}
close input;

my %seq=();
open (input, "<$ARGV[0].contigs.fa") or die "Can't open $ARGV[0].contigs.fa since $!\n";
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);
    if ($a[0] =~ />Contig/) {
        $a[0] =~ s/>//;
        if (defined $contigs{$a[0]}) {
            $seq{$a[0]}=<input>;
            chomp($seq{$a[0]});
        }
    }
}
close input;

my @p=split(/\//, $ARGV[1]);
my $reference=$p[$#p];
$reference =~ s/.2bit//;
open (output, ">>$ARGV[0].SVs.vcf") or die "Can't open $ARGV[0].SVs.vcf since $!\n";
print output "##fileformat=VCFv4.1\n##reference=$reference\n##assembly=$ARGV[0].contigs.fa\n";
print output "##INFO=<ID=BKPTID,Number=.,Type=String,Description=\"ID of the assembled alternate allele in the assembly file\">\n";
print output "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">\n";
print output "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n";
print output "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n";
print output "##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">\n";
print output "##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">\n";
print output "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n";
print output "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
print output "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">\n";
print output "##INFO=<ID=MATEORI,Number=1,Type=String,Description=\"Orientation of mate breakends\">\n";
print output "##INFO=<ID=EVENT,Number=1,Type=String,Description=\"ID of event associated to breakend\">\n";
#print output "##INFO=<ID=MEINFO,Number=4,Type=String,Description=\"Mobile element info of the form NAME,START,END,POLARITY\">\n";
print output "##ALT=<ID=DEL,Description=\"Deletion\">\n";
print output "##ALT=<ID=INS:ME,Description=\"Insertion of TE element\">\n";
print output "##FILTER=<ID=f0,Description=\"Frequency equals 0\">\n";
print output "##FILTER=<ID=r10,Description=\"Total supporting reads fewer than 10\">\n";
print output "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print output "##FORMAT=<ID=GF,Number=1,Type=Float,Description=\"Genotype frequency\">\n";
print output "##FORMAT=<ID=BKSUP,Number=1,Type=Integer,Description=\"Number of reads supporting the breakends\">\n";
print output "##FORMAT=<ID=REFSUP,Number=1,Type=Integer,Description=\"Number of reads supporting reference\">\n";
print output "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$ARGV[0]\n";

my $bp_num=0;
my $translocation=0; my $insertion=0; my $deletion=0; my $duplication=0; my $inversion=0; my $flip_dup=0;
while ((my $key, my $value) = each (%contigs)) {
    my @x=split(/\;/, $value);

    if ($#x == 0) {
	system("grep $key $ARGV[0].aligned.discordant.bed > tmp");

	my @y=split(/\:/, $x[0]);
	my $line1="";
	my $line2="";
        my $special_event=0;
	open (input, "<tmp") or die "Can't open tmp since $!\n";
	while (my $line=<input>) {
	    chomp($line);
	    my @a=split(/\t/, $line);
	    my $bp1="$a[1]$a[5]";
	    my $bp2="$a[2]$a[5]";
	    if (($y[0] eq $a[0]) && ($line1 eq "") && (($y[1] eq $bp1) || ($y[1] eq $bp2))) {
		$line1=$line;
                if (($y[3] eq $bp1) || ($y[3] eq $bp2)) {$special_event=1;}
	    }
	    elsif (($y[2] eq $a[0]) &&(($y[3] eq $bp1) || ($y[3] eq $bp2))) {
		$line2=$line;
	    }
            elsif (($y[2] eq $a[0]) && (($y[1] eq $bp1) || ($y[1] eq $bp2)) && ($line2 eq "") && ($special_event==1)) {
                $line2=$line;
            }
	}
	close input;

	chomp($line1); chomp($line2);
	my @a=split(/\t/, $line1);
	my @b=split(/\t/, $line2);
	$a[6] =~ s/H/S/g;
	$b[6] =~ s/H/S/g;
	my (@cigar_s1)=$a[6]=~/(\d+)S/g;
	my (@cigar_s2)=$b[6]=~/(\d+)S/g;
	my (@cigar_m1)=$a[6]=~/(\d+)M/g;
	my (@cigar_m2)=$b[6]=~/(\d+)M/g;
	my (@cigar_i1)=$a[6]=~/(\d+)I/g;
	my (@cigar_i2)=$b[6]=~/(\d+)I/g;
	my $start1=0; my $start2=0; my $end1=0; my $end2=0;
        if ($a[5] eq "+") {
            if (($#cigar_s1 == 0) && ($a[6] =~ /S$/)) {$start1=0;}
            else {$start1=$cigar_s1[0];}
        }
        else {
            if (($#cigar_s1 == 0) && ($a[6] !~ /S$/)) {$start1=0;}
            else {$start1=$cigar_s1[$#cigar_s1];}
        }
	if ($b[5] eq "+") {
            if (($#cigar_s2 == 0) && ($b[6] =~ /S$/)) {$start2=0;}
            else {$start2=$cigar_s2[0];}
        }
        else {
            if (($#cigar_s2 == 0) && ($b[6] !~ /S$/)) {$start2=0;}
            else {$start2=$cigar_s2[$#cigar_s2];}
        }
        $end1=$start1+sum(@cigar_m1,@cigar_i1);
        $end2=$start2+sum(@cigar_m2,@cigar_i2);
#	print "$cigar_s1[0]\t$cigar_s2[0]\n";

        my @Bblock=(); my @Cblock=();
        if ($start1 < $start2) {@Bblock=@a; @Cblock=@b;}
        else {
            @Bblock=@b; @Cblock=@a;
            my $tmp=$start1; $start1=$start2; $start2=$tmp;
            $tmp=$end1; $end1=$end2; $end2=$tmp;
        }
        my $bbs=$Bblock[1]; my $bbe=$Bblock[2];
        my $cbs=$Cblock[1]; my $cbe=$Cblock[2];
        if ($Bblock[5] eq "-") {$bbs=$Bblock[2]; $bbe=$Bblock[1];}
        if ($Cblock[5] eq "-") {$cbs=$Cblock[2]; $cbe=$Cblock[1];}

        my $svtype="";
        if ($Bblock[0] ne $Cblock[0]) {
            $svtype="TRANSLOCATION$translocation";
            $translocation++;
        }
        else {
            if ($Bblock[5] ne $Cblock[5]) {
                $svtype="INVERSION$inversion";
                $inversion++;
            }
            else {
                if ((($cbs < $bbe)&&($Bblock[5] eq "+")) || (($cbs > $bbe)&&($Bblock[5] eq "-"))) {
                    if (($start2 > $end1) && ($start2-$end1 > abs($bbe-$cbs))) {
                        $svtype="INSERTION$insertion";
                        $insertion++;
                    }
                    elsif (($end1 > $start2) && ($end1-$start2 > abs($bbe-$cbs))) {
                        $svtype="DELETION$deletion";
                        $deletion++;
                    }
                    else {
                        $svtype="DUPLICATION$duplication";
                        $duplication++;
                    }
                }
                elsif ($Bblock[5] eq "+") {
                    if ($cbs-$bbe > $start2-$end1) {
                        $svtype="DELETION$deletion";
                        $deletion++;
                    }
                    else {
                        $svtype="INSERTION$insertion";
                        $insertion++;
                    }
                }
                else {
                    if ($bbe-$cbs > $start2-$end1) {
                        $svtype="DELETION$deletion";
                        $deletion++;
                    }
                    else {
                        $svtype="INSERTION$insertion";
                        $insertion++;
                    }
                }
            }
        }

        my $filter="PASS";
        if ($y[5] + $y[6] < 10) {$filter="r10";}
        if ($y[4] == 0) {$filter="f0"; next;}

	my $mate_chr=""; my $mate_pos=0;
	my $orient="+";
	if ((($a[0] eq $y[0]) && ("$a[1]$a[5]" eq $y[1])) || (($a[0] eq $y[2]) && ("$a[1]$a[5]" eq $y[3]))) {
	    if (($a[0] eq $y[0]) && ("$a[1]$a[5]" eq $y[1])) {
		$mate_chr=$y[2]; $mate_pos=$y[3];
	    }
	    elsif (($a[0] eq $y[2]) && ("$a[1]$a[5]" eq $y[3])) {
		$mate_chr=$y[0]; $mate_pos=$y[1];
	    }
	    if ($mate_pos =~ /\-$/) {$orient="-";}
	    chop($mate_pos);
	    if (($start1 < $end2)&&($end1 > $start2)) {
		my $cilen=0;
		if (($end2-$start1) < ($end1-$start2)) {$cilen=$end2-$start1+1;}
		else {$cilen=$end1-$start2+1;}
		my $bp_start=$a[1];
		if ($a[5] eq $b[5]) {$mate_pos -= $cilen;}
		else {$mate_pos += $cilen;}

		my $pos0=$a[1]-1; my $pos1=$a[1];
		system("twoBitToFa -seq=$a[0] -start=$pos0 -end=$pos1 $ARGV[1] tmp_seq.fa");
		open (fa, "<tmp_seq.fa") or die "Can't open tmp_seq.fa since $!\n";
		my $ref=<fa>; $ref=<fa>; chomp($ref); $ref=uc($ref);
		close fa;
		system("rm tmp_seq.fa");
		if ($a[5] eq $b[5]) {
		    print output "$a[0]\t$bp_start\tbnd_$bp_num\t$ref\t\]$mate_chr\:$mate_pos\]$ref\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];CIPOS=0,$cilen;\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		}
		elsif ($a[5] eq "+") {
		    print output "$a[0]\t$bp_start\tbnd_$bp_num\t$ref\t\[$mate_chr\:$mate_pos\[$ref\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];CIPOS=0,$cilen;\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		}
		else {
		    $pos0=$mate_pos-1; $pos1=$mate_pos;
		    system("twoBitToFa -seq=$mate_chr -start=$pos0 -end=$pos1 $ARGV[1] tmp_seq.fa");
		    open (fa, "<tmp_seq.fa") or die "Can't open tmp_seq.fa since $!\n";
		    my $ref=<fa>; $ref=<fa>; chomp($ref); $ref=uc($ref);
		    close fa;
		    system("rm tmp_seq.fa");
		    print output "$mate_chr\t$mate_pos\tbnd_$bp_num\t$ref\t\[$a[0]\:$bp_start\[$ref\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];CIPOS=0,$cilen;\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		}
		$bp_num++;
	    }
	    elsif (($start1 >= $end2) || ($start2 >= $end1)) {
                my $cilen=0;
                if ($start1 > $end2) {$cilen=$start1-$end2+1;}
                else {$cilen=$start2-$end1+1;}

		my $insert="";
		if ($start2 > $end1) {$insert=substr($seq{$a[3]}, $end1, $cilen-1);}
		elsif ($start1 > $end2) {$insert=substr($seq{$a[3]}, $end2, $cilen-1);}
		if (($orient="-") && ($insert ne "")) {my $tmp_seq=Bio::Seq->new(-seq=>$insert);$insert=$tmp_seq->revcom->seq;}

		my $bp_start=$a[1];
		my $pos0=$a[1]-1; my $pos1=$a[1];
		system("twoBitToFa -seq=$a[0] -start=$pos0 -end=$pos1 $ARGV[1] tmp_seq.fa");
		open (fa, "<tmp_seq.fa") or die "Can't open tmp_seq.fa since $!\n";
		my $ref=<fa>; $ref=<fa>; chomp($ref); $ref=uc($ref);
		close fa;
		system("rm tmp_seq.fa");
		if ($a[5] eq $b[5]) {
		    print output "$a[0]\t$bp_start\tbnd_$bp_num\t$ref\t\]$mate_chr\:$mate_pos\]$insert$ref\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		}
		elsif ($a[5] eq "+") {
		    print output "$a[0]\t$bp_start\tbnd_$bp_num\t$ref\t\[$mate_chr\:$mate_pos\[$insert$ref\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		}
		else {
                    $pos0=$mate_pos-1; $pos1=$mate_pos;
                    system("twoBitToFa -seq=$mate_chr -start=$pos0 -end=$pos1 $ARGV[1] tmp_seq.fa");
                    open (fa, "<tmp_seq.fa") or die "Can't open tmp_seq.fa since $!\n";
                    my $ref=<fa>; $ref=<fa>; chomp($ref); $ref=uc($ref);
                    close fa;
                    system("rm tmp_seq.fa");
                    print output "$mate_chr\t$mate_pos\tbnd_$bp_num\t$ref\t\[$a[0]\:$bp_start\[$insert$ref\t.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		}
		$bp_num++;
	    }
	}

	elsif ((($a[0] eq $y[0]) && ("$a[2]$a[5]" eq $y[1])) || (($a[0] eq $y[2]) && ("$a[2]$a[5]" eq $y[3]))) {
	    if (($a[0] eq $y[0]) && ("$a[2]$a[5]" eq $y[1])) {
		$mate_chr=$y[2]; $mate_pos=$y[3];
	    }
	    elsif (($a[0] eq $y[2]) && ("$a[2]$a[5]" eq $y[3])) {
		$mate_chr=$y[0]; $mate_pos=$y[1];
	    }
            if ($mate_pos =~ /\-$/) {$orient="-";}
	    chop($mate_pos);
	    if (($start1 < $end2)&&($end1 > $start2)) {
                my $cilen=0;
                if (($end2-$start1) < ($end1-$start2)) {$cilen=$end2-$start1+1;}
                else {$cilen=$end1-$start2+1;}
		my $bp_start=$a[2];
		if ($a[5] eq $b[5]) {$mate_pos += $cilen;}
		else {$mate_pos -= $cilen;}

		my $pos0=$a[2]-1; my $pos1=$a[2];
		system("twoBitToFa -seq=$a[0] -start=$pos0 -end=$pos1 $ARGV[1] tmp_seq.fa");
		open (fa, "<tmp_seq.fa") or die "Can't open tmp_seq.fa since $!\n";
		my $ref=<fa>; $ref=<fa>; chomp($ref); $ref=uc($ref);
		close fa;
		system("rm tmp_seq.fa");
		if ($a[5] eq $b[5]) {
		    print output "$a[0]\t$bp_start\tbnd_$bp_num\t$ref\t$ref\[$mate_chr\:$mate_pos\[\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];CIPOS=-$cilen,0;\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		}
		elsif ($a[5] eq "+") {
		    print output "$a[0]\t$bp_start\tbnd_$bp_num\t$ref\t$ref\]$mate_chr\:$mate_pos\]\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];CIPOS=-$cilen,0;\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		}
		else {
                    $pos0=$mate_pos-1; $pos1=$mate_pos;
                    system("twoBitToFa -seq=$mate_chr -start=$pos0 -end=$pos1 $ARGV[1] tmp_seq.fa");
                    open (fa, "<tmp_seq.fa") or die "Can't open tmp_seq.fa since $!\n";
                    my $ref=<fa>; $ref=<fa>; chomp($ref); $ref=uc($ref);
                    close fa;
                    system("rm tmp_seq.fa");
                    print output "$mate_chr\t$mate_pos\tbnd_$bp_num\t$ref\t$ref\]$a[0]\:$bp_start\]\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];CIPOS=-$cilen,0;\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		}
		$bp_num++;
	    }
	    elsif (($start1 >= $end2) || ($start2 >= $end1)) {
                my $cilen=0;
                if ($start1 > $end2) {$cilen=$start1-$end2+1;}
                else {$cilen=$start2-$end1+1;}

		my $insert="";
		if ($start2 > $end1) {$insert=substr($seq{$a[3]}, $end1, $cilen-1);}
		elsif ($start1 > $end2) {$insert=substr($seq{$a[3]}, $end2, $cilen-1);}
		if (($orient="-") && ($insert ne "")) {my $tmp_seq=Bio::Seq->new(-seq=>$insert);$insert=$tmp_seq->revcom->seq;}

		my $bp_start=$a[2];
		my $pos0=$a[2]-1; my $pos1=$a[2];
		system("twoBitToFa -seq=$a[0] -start=$pos0 -end=$pos1 $ARGV[1] tmp_seq.fa");
		open (fa, "<tmp_seq.fa") or die "Can't open tmp_seq.fa since $!\n";
		my $ref=<fa>; $ref=<fa>; chomp($ref); $ref=uc($ref);
		close fa;
		system("rm tmp_seq.fa");
		if ($a[5] eq $b[5]) {
		    print output "$a[0]\t$bp_start\tbnd_$bp_num\t$ref\t$ref$insert\[$mate_chr\:$mate_pos\[\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		}
		elsif ($a[5] eq "+") {
		    print output "$a[0]\t$bp_start\tbnd_$bp_num\t$ref\t$ref$insert\]$mate_chr\:$mate_pos\]\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		}
		else {
                    $pos0=$mate_pos-1; $pos1=$mate_pos;
                    system("twoBitToFa -seq=$mate_chr -start=$pos0 -end=$pos1 $ARGV[1] tmp_seq.fa");
                    open (fa, "<tmp_seq.fa") or die "Can't open tmp_seq.fa since $!\n";
                    my $ref=<fa>; $ref=<fa>; chomp($ref); $ref=uc($ref);
                    close fa;
                    system("rm tmp_seq.fa");
                    print output "$mate_chr\t$mate_pos\tbnd_$bp_num\t$ref\t$ref$insert\]$a[0]\:$bp_start\]\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		}
		$bp_num++;
	    }
	}
#	else {print "$key\t$value\t$a[0]\t$a[1]\t$a[2]\t$a[3]\t$b[0]\t$b[1]\t$b[2]\t$b[3]\n";}

	system("rm tmp");
    }
}
close output;
 
open (output, ">>$ARGV[0].SVs.vcf") or die "Can't open $ARGV[0].SVs.vcf since $!\n";
while ((my $key, my $value) = each (%contigs)) {
    my @x=split(/\;/, $value);

    if ($#x >= 1) {
	system("grep $key $ARGV[0].aligned.discordant.bed >> tmp");
	for my $i (0..$#x) {
	    my @y=split(/\:/, $x[$i]);
	    my $line1="";
	    my $line2="";
	    my $special_event=0;
	    open (input, "<tmp") or die "Can't open tmp since $!\n";
	    while (my $line=<input>) {
		chomp($line);
		my @a=split(/\t/, $line);
		my $bp1="$a[1]$a[5]";
		my $bp2="$a[2]$a[5]";
		if (($y[0] eq $a[0]) && ($line1 eq "") && (($y[1] eq $bp1) || ($y[1] eq $bp2))) {
		    $line1=$line;
		    if (($y[3] eq $bp1) || ($y[3] eq $bp2)) {$special_event=1;}
		}
		elsif (($y[2] eq $a[0]) &&(($y[3] eq $bp1) || ($y[3] eq $bp2))) {
		    $line2=$line;
		}
		elsif (($y[2] eq $a[0]) && (($y[1] eq $bp1) || ($y[1] eq $bp2)) && ($line2 eq "") && ($special_event==1)) {
		    $line2=$line;
		}
	    }
	    close input;

	    chomp($line1); chomp($line2);
	    my @a=split(/\t/, $line1);
	    my @b=split(/\t/, $line2);
	    $a[6] =~ s/H/S/g;
	    $b[6] =~ s/H/S/g;
	    my (@cigar_s1)=$a[6]=~/(\d+)S/g;
	    my (@cigar_s2)=$b[6]=~/(\d+)S/g;
	    my (@cigar_m1)=$a[6]=~/(\d+)M/g;
	    my (@cigar_m2)=$b[6]=~/(\d+)M/g;
	    my (@cigar_i1)=$a[6]=~/(\d+)I/g;
	    my (@cigar_i2)=$b[6]=~/(\d+)I/g;
	    my $start1=0; my $start2=0; my $end1=0; my $end2=0;
	    if ($a[5] eq "+") {
		if (($#cigar_s1 == 0) && ($a[6] =~ /S$/)) {$start1=0;}
		else {$start1=$cigar_s1[0];}
	    }
	    else {
		if (($#cigar_s1 == 0) && ($a[6] !~ /S$/)) {$start1=0;}
		else {$start1=$cigar_s1[$#cigar_s1];}
	    }
	    if ($b[5] eq "+") {
		if (($#cigar_s2 == 0) && ($b[6] =~ /S$/)) {$start2=0;}
		else {$start2=$cigar_s2[0];}
	    }
	    else {
		if (($#cigar_s2 == 0) && ($b[6] !~ /S$/)) {$start2=0;}
		else {$start2=$cigar_s2[$#cigar_s2];}
	    }
	    $end1=$start1+sum(@cigar_m1,@cigar_i1);
	    $end2=$start2+sum(@cigar_m2,@cigar_i2);
	    
	    my @Bblock=(); my @Cblock=();
	    if ($start1 < $start2) {@Bblock=@a; @Cblock=@b;}
	    else {
		@Bblock=@b; @Cblock=@a;
		my $tmp=$start1; $start1=$start2; $start2=$tmp;
		$tmp=$end1; $end1=$end2; $end2=$tmp;
	    }
	    my $bbs=$Bblock[1]; my $bbe=$Bblock[2];
	    my $cbs=$Cblock[1]; my $cbe=$Cblock[2];
	    if ($Bblock[5] eq "-") {$bbs=$Bblock[2]; $bbe=$Bblock[1];}
	    if ($Cblock[5] eq "-") {$cbs=$Cblock[2]; $cbe=$Cblock[1];}
	    
	    my $svtype="";
	    if ($Bblock[0] ne $Cblock[0]) {
		$svtype="TRANSLOCATION$translocation";
		$translocation++;
	    }
	    else {
		if ($Bblock[5] ne $Cblock[5]) {
		    $svtype="INVERSION$inversion";
		    $inversion++;
		}
		else {
		    if ((($cbs < $bbe)&&($Bblock[5] eq "+")) || (($cbs > $bbe)&&($Bblock[5] eq "-"))) {
			if (($start2 > $end1) && ($start2-$end1 > abs($bbe-$cbs))) {
			    $svtype="INSERTION$insertion";
			    $insertion++;
			}
			elsif (($end1 > $start2) && ($end1-$start2 > abs($bbe-$cbs))) {
			    $svtype="DELETION$deletion";
			    $deletion++;
			}
			else {
			    $svtype="DUPLICATION$duplication";
			    $duplication++;
			}
		    }
		    elsif ($Bblock[5] eq "+") {
			if ($cbs-$bbe > $start2-$end1) {
			    $svtype="DELETION$deletion";
			    $deletion++;
			}
			else {
			    $svtype="INSERTION$insertion";
			    $insertion++;
			}
		    }
		    else {
			if ($bbe-$cbs > $start2-$end1) {
			    $svtype="DELETION$deletion";
			    $deletion++;
			}
			else {
			    $svtype="INSERTION$insertion";
			    $insertion++;
			}
		    }
		}
	    }
	    
	    my $filter="PASS";
	    if ($y[5] + $y[6] < 10) {$filter="r10";}
	    if ($y[4] == 0) {$filter="f0"; next;}

	    my $mate_chr=""; my $mate_pos=0;
	    my $orient="+";
	    if ((($a[0] eq $y[0]) && ("$a[1]$a[5]" eq $y[1])) || (($a[0] eq $y[2]) && ("$a[1]$a[5]" eq $y[3]))) {
		if (($a[0] eq $y[0]) && ("$a[1]$a[5]" eq $y[1])) {
		    $mate_chr=$y[2]; $mate_pos=$y[3];
		}
		elsif (($a[0] eq $y[2]) && ("$a[1]$a[5]" eq $y[3])) {
		    $mate_chr=$y[0]; $mate_pos=$y[1];
		}
		if ($mate_pos =~ /\-$/) {$orient="-";}
		chop($mate_pos);
		if (($start1 < $end2)&&($end1 > $start2)) {
		    my $cilen=0;
		    if (($end2-$start1) < ($end1-$start2)) {$cilen=$end2-$start1+1;}
		    else {$cilen=$end1-$start2+1;}
		    my $bp_start=$a[1];
		    if ($b[5] eq "+") {$mate_pos -= $cilen;}
		    else {$mate_pos += $cilen;}

		    my $pos0=$a[1]-1; my $pos1=$a[1];
		    system("twoBitToFa -seq=$a[0] -start=$pos0 -end=$pos1 $ARGV[1] tmp_seq.fa");
		    open (fa, "<tmp_seq.fa") or die "Can't open tmp_seq.fa since $!\n";
		    my $ref=<fa>; $ref=<fa>; chomp($ref); $ref=uc($ref);
		    close fa;
		    system("rm tmp_seq.fa");
		    if ($a[5] eq $b[5]) {
			print output "$a[0]\t$bp_start\tbnd_$bp_num\t$ref\t\]$mate_chr\:$mate_pos\]$ref\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];CIPOS=0,$cilen;\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		    }
		    elsif ($a[5] eq "+") {
			print output "$a[0]\t$bp_start\tbnd_$bp_num\t$ref\t\[$mate_chr\:$mate_pos\[$ref\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];CIPOS=0,$cilen;\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		    }
		    else {
			$pos0=$mate_pos-1; $pos1=$mate_pos;
			system("twoBitToFa -seq=$mate_chr -start=$pos0 -end=$pos1 $ARGV[1] tmp_seq.fa");
			open (fa, "<tmp_seq.fa") or die "Can't open tmp_seq.fa since $!\n";
			my $ref=<fa>; $ref=<fa>; chomp($ref); $ref=uc($ref);
			close fa;
			system("rm tmp_seq.fa");
			print output "$mate_chr\t$mate_pos\tbnd_$bp_num\t$ref\t\[$a[0]\:$bp_start\[$ref\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];CIPOS=0,$cilen;\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		    }
		    $bp_num++;
		}
		elsif (($start1 >= $end2) || ($start2 >= $end1)) {
		    my $cilen=0;
		    if ($start1 > $end2) {$cilen=$start1-$end2+1;}
		    else {$cilen=$start2-$end1+1;}

		    my $insert="";
		    if ($start2 > $end1) {$insert=substr($seq{$a[3]}, $end1, $cilen-1);}
		    elsif ($start1 > $end2) {$insert=substr($seq{$a[3]}, $end2, $cilen-1);}
		    if (($orient="-")  && ($insert ne "")) {my $tmp_seq=Bio::Seq->new(-seq=>$insert);$insert=$tmp_seq->revcom->seq;}

      		    my $bp_start=$a[1];
		    my $pos0=$a[1]-1; my $pos1=$a[1];
		    system("twoBitToFa -seq=$a[0] -start=$pos0 -end=$pos1 $ARGV[1] tmp_seq.fa");
		    open (fa, "<tmp_seq.fa") or die "Can't open tmp_seq.fa since $!\n";
		    my $ref=<fa>; $ref=<fa>; chomp($ref); $ref=uc($ref);
		    close fa;
		    system("rm tmp_seq.fa");
		    if ($a[5] eq $b[5]) {
			print output "$a[0]\t$bp_start\tbnd_$bp_num\t$ref\t\]$mate_chr\:$mate_pos\]$insert$ref\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		    }
		    elsif ($a[5] eq "+") {
			print output "$a[0]\t$bp_start\tbnd_$bp_num\t$ref\t\[$mate_chr\:$mate_pos\[$insert$ref\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		    }
		    else {
			$pos0=$mate_pos-1; $pos1=$mate_pos;
			system("twoBitToFa -seq=$mate_chr -start=$pos0 -end=$pos1 $ARGV[1] tmp_seq.fa");
			open (fa, "<tmp_seq.fa") or die "Can't open tmp_seq.fa since $!\n";
			my $ref=<fa>; $ref=<fa>; chomp($ref); $ref=uc($ref);
			close fa;
			system("rm tmp_seq.fa");
			print output "$mate_chr\t$mate_pos\tbnd_$bp_num\t$ref\t\[$a[0]\:$bp_start\[$insert$ref\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		    }
		    $bp_num++;
		}
	    }

	    if ((($a[0] eq $y[0]) && ("$a[2]$a[5]" eq $y[1])) || (($a[0] eq $y[2]) && ("$a[2]$a[5]" eq $y[3]))) {
		if (($a[0] eq $y[0]) && ("$a[2]$a[5]" eq $y[1])) {
		    $mate_chr=$y[2]; $mate_pos=$y[3];
		}
		elsif (($a[0] eq $y[2]) && ("$a[2]$a[5]" eq $y[3])) {
		    $mate_chr=$y[0]; $mate_pos=$y[1];
		}
		if ($mate_pos =~ /\-$/) {$orient="-";}
		chop($mate_pos);
		if (($start1 < $end2)&&($end1 > $start2)) {
		    my $cilen=0;
		    if (($end2-$start1) < ($end1-$start2)) {$cilen=$end2-$start1+1;}
                    else {$cilen=$end1-$start2+1;}
		    my $bp_start=$a[2];
		    if ($b[5] eq "+") {$mate_pos += $cilen;}
		    else {$mate_pos -= $cilen;}

		    my $pos0=$a[2]-1; my $pos1=$a[2];
		    system("twoBitToFa -seq=$a[0] -start=$pos0 -end=$pos1 $ARGV[1] tmp_seq.fa");
		    open (fa, "<tmp_seq.fa") or die "Can't open tmp_seq.fa since $!\n";
		    my $ref=<fa>; $ref=<fa>; chomp($ref); $ref=uc($ref);
		    close fa;
		    system("rm tmp_seq.fa");
		    if ($a[5] eq $b[5]) {
			print output "$a[0]\t$bp_start\tbnd_$bp_num\t$ref\t$ref\[$mate_chr\:$mate_pos\[\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];CIPOS=-$cilen,0;\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		    }
		    elsif ($a[5] eq "+") {
			print output "$a[0]\t$bp_start\tbnd_$bp_num\t$ref\t$ref\]$mate_chr\:$mate_pos\]\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];CIPOS=-$cilen,0;\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		    }
		    else {
			$pos0=$mate_pos-1; $pos1=$mate_pos;
			system("twoBitToFa -seq=$mate_chr -start=$pos0 -end=$pos1 $ARGV[1] tmp_seq.fa");
			open (fa, "<tmp_seq.fa") or die "Can't open tmp_seq.fa since $!\n";
			my $ref=<fa>; $ref=<fa>; chomp($ref); $ref=uc($ref);
			close fa;
			system("rm tmp_seq.fa");
			print output "$mate_chr\t$mate_pos\tbnd_$bp_num\t$ref\t$ref\]$a[0]\:$bp_start\]\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];CIPOS=-$cilen,0;\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		    }
		    $bp_num++;
		}
                elsif (($start1 >= $end2) || ($start2 >= $end1)) {
		    my $cilen=0;
		    if ($start1 > $end2) {$cilen=$start1-$end2+1;}
		    else {$cilen=$start2-$end1+1;}

                    my $insert="";
                    if ($start2 > $end1) {$insert=substr($seq{$a[3]}, $end1, $cilen-1);}
                    elsif ($start1 > $end2) {$insert=substr($seq{$a[3]}, $end2, $cilen-1);}
		    if (($orient="-")  && ($insert ne "")) {my $tmp_seq=Bio::Seq->new(-seq=>$insert);$insert=$tmp_seq->revcom->seq;}

		    my $bp_start=$a[2];
		    my $pos0=$a[2]-1; my $pos1=$a[2];
		    system("twoBitToFa -seq=$a[0] -start=$pos0 -end=$pos1 $ARGV[1] tmp_seq.fa");
		    open (fa, "<tmp_seq.fa") or die "Can't open tmp_seq.fa since $!\n";
		    my $ref=<fa>; $ref=<fa>; chomp($ref); $ref=uc($ref);
		    close fa;
		    system("rm tmp_seq.fa");
		    if ($a[5] eq $b[5]) {
			print output "$a[0]\t$bp_start\tbnd_$bp_num\t$ref\t$ref$insert\[$mate_chr\:$mate_pos\[\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		    }
		    elsif ($a[5] eq "+") {
			print output "$a[0]\t$bp_start\tbnd_$bp_num\t$ref\t$ref$insert\]$mate_chr\:$mate_pos\]\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		    }
		    else {
			$pos0=$mate_pos-1; $pos1=$mate_pos;
			system("twoBitToFa -seq=$mate_chr -start=$pos0 -end=$pos1 $ARGV[1] tmp_seq.fa");
			open (fa, "<tmp_seq.fa") or die "Can't open tmp_seq.fa since $!\n";
			my $ref=<fa>; $ref=<fa>; chomp($ref); $ref=uc($ref);
			close fa;
			system("rm tmp_seq.fa");
			print output "$mate_chr\t$mate_pos\tbnd_$bp_num\t$ref\t$ref$insert\]$a[0]\:$bp_start\]\t\.\t$filter\tSVTYPE=BND;EVENT=$svtype;BKPTID=$a[3];\tGF:BKSUP:REFSUP\t$y[4]\:$y[5]\:$y[6]\n";
		    }
		    $bp_num++;
		}
	    }

	}	
	system("rm tmp");
    }

}
close output;
