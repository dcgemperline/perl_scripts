#!/usr/bin/perl -w
use strict;

my $myfile;
my @lines;
$myfile ="FBXCDS_testing.txt";

open FILE, "<$myfile" or die $!;
while( <FILE> ) {
        next if /^(\s)*$/;  # skip blank lines
        chomp;              # remove trailing newline characters
        push @lines, $_;    # push the data line onto the array
    }
close FILE;

my $count=0;
foreach(@lines){
	if ($_ =~ m/^\>/){
	$count++;
	}
}

print $count;