#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;

my $db1; #TAIR10
$db1 = Bio::DB::Fasta->new('TAIR10_chr_all.fasta');
my $db2; #bac sequences in fasta format
$db2 = Bio::DB::Fasta->new('all_fbox_bac_at_check_positions_.fasta');
my @ids=$db2->get_all_ids;

my $bacfile = "all_fbox_bac_at_check_positions_.fasta";
my @lines;

open FILE, "<$bacfile" or die $!;
while( <FILE> ) {
        next if /^(\s)*$/;  # skip blank lines
        chomp;              # remove trailing newline characters
        push @lines, $_;    # push the data line onto the array
    }
close FILE;

my %lookuphash;
my @header;
my @header2;
foreach(@lines){
	if ($_ =~ m/^\>/){
		@header = (split(/\| /,$'));
		my $key = $header[0];
		$key =~ s/^\s+|\s+$//g; # remove whitespace
		@header2=split(/-/,$header[1]);
		my @data = ($header2[0], $header2[1], $header2[2], $header2[3]);
		$lookuphash {$key} = [@data];
	}
}
my $count=0;
foreach (sort(@ids)){
	 my $name=$_;
	 my $chromasome =$lookuphash{$_}[0];
	 my $strand =$lookuphash{$_}[1];
	 my $start;
	 my $end;
	 if ($strand eq "plus"){
	 	$start =$lookuphash{$_}[2];
	 	$end =$lookuphash{$_}[3];
	 }
	 if ($strand eq "minus"){
	 	$start =$lookuphash{$_}[3];
	 	$end =$lookuphash{$_}[2];
	 }
	my $object =$db2->get_Seq_by_id($_);
	my $testseq = $object->seq;
	my $testseq2 = getSeq($chromasome, $start, $end);
	if ($testseq eq $testseq2){
		$count++
	}
}

print $count;



sub getSeq{
	##Get Sequence from Chromasome
	my $chromosome = $_[0];
	my $start = $_[1];
	my $end = $_[2];
	my $strand=$_[3];
	my $myseq =$db1->seq($chromosome, int($start) => int($end)); ## net to fix int for $chromasome (mitochondira and chloroplast are not ints and are viable input)
	return $myseq;
}
