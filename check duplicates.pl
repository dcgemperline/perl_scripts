#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;

my $geno; #TAIR10 For genomic Sequences
$geno = Bio::DB::Fasta->new('all_fbox_bac_at_check_positions_.fasta');
my $cds; #CDS from Zhihua
$cds = Bio::DB::Fasta->new('at_cds_list.fasta');
my @ids=$geno->get_all_ids;
my @ids2=$cds->get_all_ids;

my $myfile = "FBXoutput_testing.txt";
my @lines;

open FILE, "<$myfile" or die $!;
while( <FILE> ) {
        #next if /^(\s)*$/;  # skip blank lines
        chomp;              # remove trailing newline characters
        push @lines, $_;    # push the data line onto the array
    }
close FILE;

my %lookuphash;
my @header;
my $indx;
my $count;
my $genomiccount;
populatehash();
checkCDS();
checkGenomic();

sub populatehash{
	for ($indx=0; $indx<scalar(@lines); $indx++){
	#for ($indx=0; $indx<2; $indx++){
		if ($lines[$indx] =~ m/^\>/){
			@header = (split(/\>|\||\#/,$lines[$indx]));
			my $namekey = $header[1];
			my $exonnumber = $header[5];
			my $preCDSline=$lines[$indx+2*$exonnumber+1];
			my $preGenomic=$lines[$indx+2*$exonnumber];
			my @CDSlinearray =split(/CDS\t/, $preCDSline);
			my @Genomiclinearray=split(/\t/,$preGenomic);
			my $CDSline =$CDSlinearray[1];
			my $Genomicline = $Genomiclinearray[2];
			my @data = ($CDSline,$Genomicline);
			$lookuphash {$namekey} = [@data];
		}
	}
}

sub checkCDS {
	my $count=0;
	foreach(@ids2){
	my $obj=$cds->get_Seq_by_id($_);
	my $seq =$obj->seq;
	my $testseq = $lookuphash{$_}[0];
		if($seq eq $testseq) {
		$count++
		}	
	}
	print $count."\n";
}

sub checkGenomic {
	my %genomiclookuphash;
	my @lines2;
	my $indx;
	my $bac_seq_file ='all_fbox_bac_at_check_positions_.fasta';
	open FILE, "<$bac_seq_file" or die $!;
	while( <FILE> ) {
        #next if /^(\s)*$/;  # skip blank lines
        chomp;              # remove trailing newline characters
        push @lines2, $_;    # push the data line onto the array
    }
	close FILE;
	for ($indx=0; $indx<scalar(@lines2); $indx++){
		if ($lines2[$indx] =~ m/^\>/){
			@header = (split(/\>|\||\#/,$lines2[$indx]));
			my $key = $header[1];
			$key =~ s/^\s+|\s+$//g; # remove whitespace
			my $genomicsequence = $lines2[$indx+1];
			$genomiclookuphash {$key} = $genomicsequence;
		}
	}
	
	foreach (sort(@ids)){
		my $refseq;
		my $queryseq;
		$refseq = $genomiclookuphash{$_};
		$queryseq = $lookuphash{$_}[1];
		if ($queryseq eq $refseq){
			$genomiccount++;
		}

	}
print $genomiccount;
}