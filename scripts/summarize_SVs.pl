#! /usr/bin/perl

use strict;

my @files=<*.SVs.vcf>;
foreach my $file (@files) {
    my $title=$file;
    $title =~ s/vcf/summary/;
    open (input, "<$file") or die "Can't open $file since $!\n";
    open (output, ">>$title") or die "Can't open $title since $!\n";
    while (my $line=<input>) {
	if ($line !~ /^\#/) {
	    chomp($line);
	    my @a=split(/\t/, $line);
	    my @z=split(/\:/, $a[9]);
	    
	    if (($z[1] >= 4) && ($a[7] =~ /EVENT/) && ($a[6] eq "PASS")) {
		my @b=split(/\;/, $a[7]);
		$b[1] =~ s/EVENT=//;
		if ($b[1] eq "") {$b[1]="Unknown";}
		
		my $mate="";
		my $insert="";
		if ($a[4] =~ /\[/) {
		    my @c=split(/\[/, $a[4]);
		    for my $i (0..$#c) {
			if ($c[$i] =~ /\:/) {$mate=$c[$i];}
			elsif ($c[$i] =~ /[ACGT]/) {$insert=$c[$i];}
		    }
		}
		elsif ($a[4] =~ /\]/) {
		    my @c=split(/\]/, $a[4]);
		    for my $i (0..$#c) {
			if ($c[$i] =~ /\:/) {$mate=$c[$i];}
			elsif ($c[$i] =~ /[ACGT]/) {$insert=$c[$i];}
		    }
		}
		
		my @d=split(/\:/, $mate);
		my $dist="NA"; my $lower=0; my $upper=0;
		if ($d[0] eq $a[0]) {
		    $dist=abs($d[1]-$a[1]);
		    if ($b[1] =~ /INSERTION/) {$dist=length($insert);}
		    if ($a[1] <= $d[1]) {$lower=$a[1]-1; $upper=$d[1];}
		    else {$lower=$d[1]-1; $upper=$a[1];}
		}
		else {
		    $lower=$a[1]-1000; $upper=$a[1]+1000;
		    if ($lower < 0) {$lower=0;}
		}
		
		print output "$a[0]\t$lower\t$upper\t$a[2]\t$b[1]\t$dist";
		
		for my $i (9..$#a) {
		    my @e=split(/\:/, $a[$i]);
		    if ($b[1] =~ /DUPLICATION/) {
			my $title1=$file;
			$title1 =~ s/SVs.vcf/breakpoints.freq/;
			$b[2] =~ s/BKPTID=//;
			system("grep $b[2] $title1 > tmp");
			my $new_freq=0;
			open (input1, "<tmp") or die "Can't open tmp since $!\n";
			while (my $line1=<input1>) {
			    chomp($line1);
			    my @q=split(/\t/, $line1);
			    if (($q[1] == $e[1]) && ($q[3] == $e[2]) && ($q[5] == $e[0])) {
				if ($q[3] == 0) {
				    $new_freq="1.0000"; 
				    next;
				}
				$new_freq=sprintf("%.4f", ($q[1]/$q[2])/($q[3]/$q[4]));
				if ($new_freq > 1) {$new_freq="1.0000";}
				last;
			    }
			}
			close input1;
			system("rm tmp");
			print output "\t$new_freq\t$e[1]";
		    }
		    else {
			print output "\t$e[0]\t$e[1]";
		    }
		}
		
		print output "\n";
		
		if ($d[0] ne $a[0]) {
		    $lower=$d[1]-1000; $upper=$d[1]+1000;
		    if ($lower < 0) {$lower=0;}
		    print output "$d[0]\t$lower\t$upper\t$a[2]\t$b[1]\t$dist";
		    
		    for my $i (9..$#a) {
			my @e=split(/\:/, $a[$i]);
			print output "\t$e[0]\t$e[1]";
		    }
		    
		    print output "\n";
		}
		
	    }
	}
    }
    close input;
    close output;
}
	   
