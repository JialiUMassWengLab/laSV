#!/share/bin/perl

use strict;
use Bio::Seq;
use List::Util qw(sum);

die "perl $0 <sam> <output prefix>\n" if @ARGV<1;

#my $insert=$ARGV[1];

my %samp=(); my %ref1=(); my %ref2=();
my %type=();

open (input, "<$ARGV[0].breakpoints.fa") or die "Can't open $ARGV[0].breakpoints.fa since $!\n";
while (my $line=<input>) {
    chomp($line);
    if ($line =~ /^>/) {
	$line =~ s/\>//;
        my @a=split(/\:/, $line);
        my $bp="$a[0]\:$a[1]\:$a[2]\:$a[3]\:$a[4]";
	$samp{$bp}=length(<input>)-1;

	system("grep $a[0] $ARGV[0].aligned.discordant.bed > tmp");
	my $line1="";
	my $line2="";
	open (input1, "<tmp") or die "Can't open tmp since $!\n";
	while (my $line0=<input1>) {
	    chomp($line0);
	    my @x=split(/\t/, $line0);
	    my $bp1="$x[1]$x[5]";
	    my $bp2="$x[2]$x[5]";
	    if (($a[1] eq $x[0]) && ($line1 eq "") && (($a[2] eq $bp1) || ($a[2] eq $bp2))) {
		$line1=$line0;
	    }
	    elsif (($a[3] eq $x[0]) && (($a[4] eq $bp1) || ($a[4] eq $bp2))) {
		$line2=$line0;
	    }
	}
	close input1;
        chomp($line1); chomp($line2);
        my @b=split(/\t/, $line1);
        my @c=split(/\t/, $line2);
	$b[6] =~ s/H/S/g;
	$c[6] =~ s/H/S/g;
        my (@cigar_s1)=$b[6]=~/(\d+)S/g;
        my (@cigar_s2)=$c[6]=~/(\d+)S/g;
        my (@cigar_m1)=$b[6]=~/(\d+)M/g;
        my (@cigar_m2)=$c[6]=~/(\d+)M/g;
        my (@cigar_i1)=$b[6]=~/(\d+)I/g;
        my (@cigar_i2)=$c[6]=~/(\d+)I/g;
        my $start1=0; my $start2=0;
        if ($b[5] eq "+") {
            if (($#cigar_s1 == 0) && ($b[6] =~ /S$/)) {$start1=1;}
            else {$start1=$cigar_s1[0]+1;}
        }
        else {
            if (($#cigar_s1 == 0) && ($b[6] !~ /S$/)) {$start1=1;}
            else {$start1=$cigar_s1[$#cigar_s1]+1;}
        }
        if ($c[5] eq "+") {
            if (($#cigar_s2 == 0) && ($c[6] =~ /S$/)) {$start2=1;}
            else {$start2=$cigar_s2[0]+1;}
        }
        else {
            if (($#cigar_s2 == 0) && ($c[6] !~ /S$/)) {$start2=1;}
            else {$start2=$cigar_s2[$#cigar_s2]+1;}
        }
        my $end1=$start1+sum(@cigar_m1,@cigar_i1)-1;
        my $end2=$start2+sum(@cigar_m2,@cigar_i2)-1;
	close input1;
	system("rm tmp");
	
	if (($end1 < $start2 || $end2 < $start1) && ($a[1] eq $a[3]) && (abs($a[2]-$a[4]) < 100)) {$type{$bp}="I";}
	else {$type{$bp}="O";}
    }
}
close input;

open (input, "<$ARGV[0].reference1.fa") or die "Can't open $ARGV[0].breakpoints.fa since $!\n";
while (my $line=<input>) {
    chomp($line);
    if ($line =~ /^>/) {
        $line =~ s/\>//;
        my @a=split(/\:/, $line);
        my $bp="$a[0]\:$a[1]\:$a[2]\:$a[3]\:$a[4]";
        $ref1{$bp}=length(<input>)-1;
    }
}
close input;
open (input, "<$ARGV[0].reference2.fa") or die "Can't open $ARGV[0].breakpoints.fa since $!\n";
while (my $line=<input>) {
    chomp($line);
    if ($line =~ /^>/) {
        $line =~ s/\>//;
        my @a=split(/\:/, $line);
        my $bp="$a[0]\:$a[1]\:$a[2]\:$a[3]\:$a[4]";
        $ref2{$bp}=length(<input>)-1;
    }
}
close input;

system("samtools view -Xf 0x2 $ARGV[0].reference1.sorted.bam | awk -F \"\\t\" '{OFS=\"\\t\"; if (\$9>0) print \$9}' > distance");
system("Rscript $ARGV[1]/dist_mean_sd.R");
open (input, "<dist_mean_sd.txt") or die "Can't open dist_mean_sd.txt since $!\n";
my $dist_mean=<input>; chomp($dist_mean);
my $dist_sd=<input>; chomp($dist_sd);
close input;
system("rm dist_mean_sd.txt distance");

open (output, ">>$ARGV[0].breakpoints.freq") or die "Can't open $ARGV[0].breakpoints.freq since $!\n";
foreach my $bp (keys %samp) {
    my $ref_sup=0;
    my $sv_sup=0;

    my %ps_samp=(); my %me_samp=(); my %pe_samp=(); my %ms_samp=(); my %mm_samp=(); my %dist_samp=();
    my %ps_ref1=(); my %me_ref1=(); my %pe_ref1=(); my %ms_ref1=(); my %mm_ref1=(); my %dist_ref1=();
    my %ps_ref2=(); my %me_ref2=(); my %pe_ref2=(); my %ms_ref2=(); my %mm_ref2=(); my %dist_ref2=();

#    my @z=split(/\:/, $bp);
#    chop($z[2]); chop($z[4]);
#    if (($type{$bp} eq "O") && ($z[1] eq $z[3]) && (abs($z[4]-$z[2])+$samp{$bp}-1000 <= 3*$dist_sd)) {next;}

    system("samtools view -Xf 0x2 $ARGV[0].breakpoints.sorted.bam $bp\:Sample > temp.sam");
    open in,"temp.sam";
    while(<in>)
    {
	chomp;
	my @f=split/\t/,$_,12;
	## read number 1 or 2
	my ($rnum)=$f[1]=~/(\d)$/;
	
	## NM:Z:* 
	my @z=split(/\s+/, $f[11]);
	my @x=split(/\:/, $z[1]);
	my $nm=$x[2];
	if (! defined $mm_samp{$f[0]}) {$mm_samp{$f[0]}=$nm;}
	else {$mm_samp{$f[0]} += $nm;}
	
	## Coordinate
	my (@cigar_m)=$f[5]=~/(\d+)M/g;
	my (@cigar_d)=$f[5]=~/(\d+)D/g;
	my (@cigar_s)=$f[5]=~/(\d+)S/g;
	my (@cigar_i)=$f[5]=~/(\d+)I/g;
	my $aln_ln=sum(@cigar_m,@cigar_d);

	$f[0] =~ s/\#//;
	if (($nm <= 5) && ($aln_ln >= 40)) {
	    if ($f[1]=~/r/)
	    {
		$ms_samp{$f[0]}=$f[3];
		$me_samp{$f[0]}=$f[3]+$aln_ln-1;
	    }
	    else
	    {
		$ps_samp{$f[0]}=$f[3];
		$pe_samp{$f[0]}=$f[3]+$aln_ln-1;
		$dist_samp{$f[0]}=abs($f[8]);
	    }    	
	    $mm_samp{$f[0]} += sum(@cigar_s);
	}
    }
    close in;
    system("rm temp.sam");

    my $ref_len=0;
    if (defined $ref1{$bp}) {
	system("samtools view -Xf 0x2 $ARGV[0].reference1.sorted.bam $bp\:Reference1 > temp.sam");
	open in,"temp.sam";
	while(<in>)
	{
	    chomp;
	    my @f=split/\t/,$_,12;
	    ## read number 1 or 2
	    my ($rnum)=$f[1]=~/(\d)$/;
	    
	    ## NM:Z:*
	    my @z=split(/\s+/, $f[11]);
	    my @x=split(/\:/, $z[1]);
	    my $nm=$x[2];
	    if (! defined $mm_ref1{$f[0]}) {$mm_ref1{$f[0]}=$nm;}
	    else {$mm_ref1{$f[0]} += $nm;}

	    ## Coordinate
	    my (@cigar_m)=$f[5]=~/(\d+)M/g;
	    my (@cigar_d)=$f[5]=~/(\d+)D/g;
	    my (@cigar_s)=$f[5]=~/(\d+)S/g;
	    my (@cigar_i)=$f[5]=~/(\d+)I/g;
	    my $aln_ln=sum(@cigar_m,@cigar_d);

            $f[0] =~ s/\#//;
	    if (($nm <= 5) && ($aln_ln >= 40)) {
		if ($f[1]=~/r/)
		{
		    $ms_ref1{$f[0]}=$f[3];
		    $me_ref1{$f[0]}=$f[3]+$aln_ln-1;
		}
		else
		{
		    $ps_ref1{$f[0]}=$f[3];
		    $pe_ref1{$f[0]}=$f[3]+$aln_ln-1;
                    $dist_ref1{$f[0]}=abs($f[8]);
		}
		$mm_ref1{$f[0]} += sum(@cigar_s);
	    }
	}
	close in;
	
	$ref_len += $ref1{$bp};
	system("rm temp.sam");
    }

    if (defined $ref2{$bp}) {
	system("samtools view -Xf 0x2 $ARGV[0].reference2.sorted.bam $bp\:Reference2 > temp.sam");
	open in,"temp.sam";
	while(<in>)
	{
	    chomp;
	    my @f=split/\t/,$_,12;
	    ## read number 1 or 2
	    my ($rnum)=$f[1]=~/(\d)$/;
	    
	    ## NM:Z:*
            my @z=split(/\s+/, $f[11]);
            my @x=split(/\:/, $z[1]);
            my $nm=$x[2];
	    if (! defined $mm_ref2{$f[0]}) {$mm_ref2{$f[0]}=$nm;}
	    else {$mm_ref2{$f[0]} += $nm;}
	    
	    ## Coordinate
	    my (@cigar_m)=$f[5]=~/(\d+)M/g;
	    my (@cigar_d)=$f[5]=~/(\d+)D/g;
	    my (@cigar_s)=$f[5]=~/(\d+)S/g;
	    my (@cigar_i)=$f[5]=~/(\d+)I/g;
	    my $aln_ln=sum(@cigar_m,@cigar_d);

            $f[0] =~ s/\#//;
	    if (($nm <= 5) && ($aln_ln >= 40)) {
		if ($f[1]=~/r/)
		{
		    $ms_ref2{$f[0]}=$f[3];
		    $me_ref2{$f[0]}=$f[3]+$aln_ln-1;
		}
		else
		{
		    $ps_ref2{$f[0]} = $f[3];
		    $pe_ref2{$f[0]} = $f[3]+$aln_ln-1;
                    $dist_ref2{$f[0]}=abs($f[8]);
		}
		$mm_ref2{$f[0]} += sum(@cigar_s);
	    }
	}
	close in;

	$ref_len += $ref2{$bp};
	system("rm temp.sam");
    }

    
    my %samp_reads=();
    foreach my $id (keys %ps_samp)
    {
	if ((defined $me_samp{$id}) && ($mm_samp{$id} <= 3) && (abs($dist_samp{$id}-$dist_mean)/$dist_sd <= 2.5)) {
	    my $pt2=$samp{$bp}-500;

	    if ($type{$bp} eq "O") {
		my $ifstradle1 = $ps_samp{$id}-491;
		my $ifstradle2 = $me_samp{$id}-$pt2-10;
		if (($ifstradle1 < 0) && ($ifstradle2 > 0)) {
		    $samp_reads{$id}=1;
		    $sv_sup++;
		    if (($pe_samp{$id}-511 > 0) || ($ms_samp{$id}-$pt2+10 < 0)) {$samp_reads{$id}++;}
		}
	    }
	    elsif ($type{$bp} eq "I") {
		my $ifstradle1 = ($ps_samp{$id}-491)*($me_samp{$id}-511);
		my $ifstradle2 = ($ps_samp{$id}-$pt2+10)*($me_samp{$id}-$pt2-10);
		if (($ifstradle1 < 0) || ($ifstradle2 < 0)) {
		    $samp_reads{$id}=1;
		    $sv_sup++;
		    my $ifspan1 = ($pe_samp{$id}-511)*($pe_samp{$id}-$pt2+10);
		    my $ifspan2 = ($ms_samp{$id}-511)*($ms_samp{$id}-$pt2+10);
		    if (($ifspan1 < 0) || ($ifspan2 < 0)) {$samp_reads{$id}++;}
		}
	    }
	}
    }
    open (output1, ">>undecided_reads") or die "Can't open undecided_reads since $!\n";
    foreach my $id (keys %ps_ref1)
    {
	if ((defined $me_ref1{$id}) && ($mm_ref1{$id} <= 5) && (abs($dist_ref1{$id}-$dist_mean)/$dist_sd <= 2.5)) {
            my $ifstradle = ($ps_ref1{$id}-491)*($me_ref1{$id}-511);
            if ($ifstradle < 0) {
		my $ifspan1 = ($ps_ref1{$id}-491)*($pe_ref1{$id}-511);
		my $ifspan2 = ($ms_ref1{$id}-491)*($me_ref1{$id}-511);
		if ((! defined $samp_reads{$id}) || ($samp_reads{$id} == 0)) {
		    $ref_sup++;
		}
		elsif ((($ifspan1 < 0)||($ifspan2 < 0)) && ($samp_reads{$id} == 1)) {
		    $ref_sup++; 
		    $sv_sup--;
		    $samp_reads{$id}=0;
		}
		elsif (($ifspan1 >= 0) && ($ifspan2 >= 0) && ($samp_reads{$id} == 2)) {;}
		else {
		    print output1 "$id\t$dist_samp{$id}\t$mm_samp{$id}\t$dist_ref1{$id}\t$mm_ref1{$id}\n";
		}
            }
        }
    }
    foreach my $id (keys %ps_ref2)
    {
        if ((defined $me_ref2{$id}) && ($mm_ref2{$id} <= 5) && (abs($dist_ref2{$id}-$dist_mean)/$dist_sd <= 2.5)) {
            my $ifstradle = ($ps_ref2{$id}-491)*($me_ref2{$id}-511);
            if ($ifstradle < 0) {
                my $ifspan1 = ($ps_ref2{$id}-491)*($pe_ref2{$id}-511);
                my $ifspan2 = ($ms_ref2{$id}-491)*($me_ref2{$id}-511);
		if ((! defined $samp_reads{$id}) || ($samp_reads{$id} == 0)) {
		    $ref_sup++;
		}
                elsif ((($ifspan1 < 0)||($ifspan2 < 0)) && ($samp_reads{$id} == 1)) {
                    $ref_sup++;
                    $sv_sup--;
		    $samp_reads{$id}=0;
                }
                elsif (($ifspan1 >= 0) && ($ifspan2 >= 0) && ($samp_reads{$id} == 2)) {;}
		else {
		    print output1 "$id\t$dist_samp{$id}\t$mm_samp{$id}\t$dist_ref2{$id}\t$mm_ref2{$id}\n";
		}
            }
        }
    }
    close output1;

    if (-s "undecided_reads") {
	system("Rscript $ARGV[1]/decide_reads.R $dist_mean $dist_sd");
	if (-s "reads_change.txt") {
	    open (input, "<reads_change.txt") or die "Can't open reads_change.txt since $!\n";
	    while (my $line=<input>) {
		chomp($line);
		if ($samp_reads{$line} > 0) {
		    $sv_sup--;
		    $samp_reads{$line}=0;
		}
		$ref_sup++;
	    }
	    close input;
	    system("rm reads_change.txt");
	}

	if (-s "still_undecided.txt") {
	    open (input, "<still_undecided.txt") or die "Can't open still_undecided.txt since $!\n";
	    while (my $line=<input>) {
		chomp($line);
		if ($samp_reads{$line} > 0) {
		    $sv_sup--;
		    $samp_reads{$line}=0;
		}
	    }
	    close input;
	    system("rm still_undecided.txt");
	}
    }
    system("rm undecided_reads");

    if ($ref_len > 0) {
	if ($type{$bp} eq "O") {
	    if ($samp{$bp}-1000 > $dist_mean) {next;}
	    $samp{$bp}=2*$dist_mean+1000-$samp{$bp};
	}
	elsif ($samp{$bp}-1000 > 2*$dist_mean) {$samp{$bp}=4*$dist_mean;}
        else {$samp{$bp}=2*$dist_mean+$samp{$bp}-1000;}

	if ((defined $ref1{$bp}) && (defined $ref2{$bp})) {$ref_len=4*$dist_mean;}
	else {$ref_len=2*$dist_mean;}

	my $sv_sup_norm=sprintf("%.6f", $sv_sup/$samp{$bp});
	my $ref_sup_norm=sprintf("%.6f", $ref_sup/$ref_len);
	my $ratio="0.00000";
	if ($sv_sup_norm + $ref_sup_norm > 0) {$ratio=sprintf("%.5f", $sv_sup_norm/($sv_sup_norm + $ref_sup_norm));}
	
	print output "$bp\t$sv_sup\t$samp{$bp}\t$ref_sup\t$ref_len\t$ratio\n";
    }

}
close output;
