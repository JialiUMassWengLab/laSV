#! /usr/bin/perl

use strict;

my %branchlen=();

my $last_read="";
my $kmers="";
my $past="";
my $outfile=$ARGV[0].".read.span";
#$outfile =~ s/map.to.branch.sam/read.span/;

open (output, ">>$outfile") or die "Can't open $outfile since $!\n";
while (my $line=<STDIN>) {
    chomp($line);
    my @a=split(/\t/, $line);

    if ($a[0] eq "\@SQ") {
	$a[1] =~ s/SN\:b\_//;
	$a[2] =~ s/LN\://;
	$branchlen{$a[1]}=$a[2];
    }

    elsif ($a[5] ne "*") {
	$a[2] =~ s/b_//;
	my $match=0;
	my $edit_distance=100;
	my $clip1=0;
	my $clip2=0;
	
	my @b=split(/M/,$a[5]);
	if ($b[0] =~ /S/) {
	    my @c=split(/S/,$b[0]);
	    $clip1=1;
	    $match=$c[1];
	}
	elsif ($b[0] =~ /H/) {
            my @c=split(/H/,$b[0]);
            $clip1=1;
            $match=$c[1];
	}
	else {$match=$b[0];}
	if (($b[1] =~ /S/) || ($b[1] =~ /H/)) {$clip2=1;}
	
	my @d=split(/\:/, $a[11]);
	$edit_distance=$d[2];

#	if (($match < $ARGV[1]) or (($clip1==1) and ($a[3] != 1)) or (($clip2==1) and ($a[3]+$match-1 != $branchlen{$a[2]})) or ($edit_distance > 0)) {}
	if (($match < $ARGV[1]) or ($edit_distance > 0)) {}
	
	else {
	    my $curr_read=$a[0];	
            $curr_read =~ s/\/1$//;
            $curr_read =~ s/\/2$//;
	    if ($curr_read eq $last_read) {
		if ($past !~ /$a[2]/) {
		    $kmers=$kmers."$a[2]";
		    $past=$past.",$a[2]";
		}
	    }
	    else {
		if (($last_read ne "") && (length($kmers) > $ARGV[1])) {
#		    print output "$last_read $kmers\n";
		    print output "$kmers\n";
		}
		$last_read=$curr_read;
		$kmers=$a[2];
		$past=$a[2];
	    }
	}

    }

}
close input;
close output;

