#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

my $directory = 'Phytozome_v8';
opendir DIR, $directory or die "$directory does not exist";
my @files=readdir(DIR);
@files = sort(@files);
my @fastafiles;
foreach(@files){
	if ($_ =~ /\.fa$/){
	push(@fastafiles, $_)	
	}
}
chdir($directory);

my %HbYX_sequences;

foreach(@fastafiles){
	my $inputfile=$_;
	my $cleanedfile = cleanfasta($inputfile);
	HbyX_Finder($cleanedfile);
	Output_printer($cleanedfile);
	%HbYX_sequences=();
}



sub filetoarray{
	open FILE, "<$_[0]" or die $!;
	my @lines =<FILE>;
	close FILE;
	return @lines;
}

sub cleanfasta{
#This removes Newline Chars from Fasta files in the sequences (ex 10nt per line to total sequence), makes them easier to handle
#input file
my $file1 ="$_[0]";
my $workingdirectory =cwd();
open FILE, "<$file1" or die 'Cant Open file from input file list at'.$workingdirectory;
#output file
chdir('Cleaned_Fasta');
my $cleanfile ="formatted_".$file1;
open FILE2, ">$cleanfile";
while (<FILE>) {
	next if /^#/;
	die unless /^>/;
	chomp;
	my $header= $_;
	my $seq = '';
	while (<FILE>) {
		next if /^#/;
		last if /^>/;
		chomp;
		$seq .= $_;
	}
	print FILE2 "$header\n$seq\n";
	redo if /^>/;
	last if eof;
}
close FILE;
chdir ($workingdirectory);
return $cleanfile;
}

sub HbyX_Finder{
	my $workingdir=cwd();
	chdir('Cleaned_Fasta');
	my @protein_fasta_file_array = filetoarray($_[0]);
	chdir($workingdir);
	for (my $index=0; $index<(my $size=@protein_fasta_file_array); $index++){
		if ($protein_fasta_file_array[$index] =~ m/^\>/){
			my @header = (split(/\|/,$'));
			my $key = $header[0];
			#my $description = $header[2];
			my $sequence = $protein_fasta_file_array[$index+1];
			$key =~ s/^\s+|\s+$//g; # remove whitespace
			chomp($sequence);
			if ($sequence =~ /\*$/){
				chop($sequence);
			}
			if ($sequence =~ /[WYF]..[LIV]/){
				#my @sequence_description_array = ($sequence,$description);
				$HbYX_sequences{$key}=$sequence;
			}
				
		}
		
	}
	#loop through and store ID,  description, and sequence as strings
	#if regex of the sequence contains the HbYX Motif add it to the array containing the sequence and description
	#based on TAIR10 protein sequence descriptions in this format
	# splitstring based on | with the first element being an ID [0], the third being the description [2]
}


sub Output_printer{
	my $directory = cwd();
	chdir('Output') or die "Failed to change to Output Directory";
	open (OUTPUTFILE, '>AIM_'.$_[0]);
	foreach my $key (keys %HbYX_sequences){
		my $id =$key;
		my $sequence = $HbYX_sequences{$key};
		print OUTPUTFILE "\>".$id."\|".$_[0]."\n";
		print OUTPUTFILE $sequence."\n";
	}
	close OUTPUTFILE;
	chdir($directory);
	
}



