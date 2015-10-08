#!/usr/bin/perl -w
use strict;
##blat command that I ran
# blat TAIR10_chr_all.fasta  at_cds_list.fasta -out=psl new_output1.psl -tileSize=12 -maxGap=0 -minIdentity=100 -noHead 
# blat all_fbox_bac_at_check_positions_.fasta at_cds_list.fasta -out=psl new_output1.psl -tileSize=12 -maxGap=0 -minIdentity=100 -noHead
#Dependencies
use Bio::DB::Fasta;
use Bio::SeqIO;
use Try::Tiny;

#class BlatLine
package BlatLine;

# Declare Database before use
my $db;
my %lookuphash;

#constructor
sub new {
    my ($class) = @_;
    my $self = {
        _name 			=> undef,
        _strand  		=> undef,
        _chNum      	=> undef,
        _chStart   		=> undef,
        _chEnd 			=> undef,
        _blkCt  		=> undef,
        _blkSz      	=> undef,
        _qstart   		=> undef,
        _tstart   		=> undef,
        _exon	    	=> undef,
        _numexon		=> undef,
        _intron	   		=> undef,
        _matchsize		=> undef,
        _genomicstart	=> undef,
        _genomicend		=> undef
        
    };
    bless $self, $class;
    return $self;
}

#accessor method for Blatline name
sub name {
    my ( $self, $name ) = @_;
    $self->{_name} = $name if defined($name);
    return $self->{_name};
}
#accessor method for Blatline sstartmathtrand
sub strand {
    my ( $self, $strand ) = @_;
    $self->{_strand} = $strand if defined($strand);
    return $self->{_strand};
}
#accessor method for Blatline chNum
sub chNum {
    my ( $self, $chNum ) = @_;
    $self->{_chNum} = $chNum if defined($chNum);
    return $self->{_chNum};
}
#accessor method for Blatline chStart
sub chStart {
    my ( $self, $chStart ) = @_;
    $self->{_chStart} = $chStart if defined($chStart);
    return $self->{_chStart};
}
#accessor method for Blatline chEnd
sub chEnd {
    my ( $self, $chEnd ) = @_;
    $self->{_chEnd} = $chEnd if defined($chEnd);
    return $self->{_chEnd};
}
#accessor method for Blatline blkCt
sub blkCt {
    my ( $self, $blkCt ) = @_;
    $self->{_blkCt} = $blkCt if defined($blkCt);
    return $self->{_blkCt};
}
#accessor method for Blatline blkSz
sub blkSz {
	my ( $self, $blkSz ) = @_;
	my $strand = $_[2];
    my @splits = split(/,/,$blkSz);
    my $splitslength = scalar(@splits);
    
    if ($strand eq "-"){
    	@splits=reverse(@splits)
    }
    
    for (my $indx = 0; $indx < $splitslength; $indx++)
    {
    		$self->{_blkSz}[$indx] = $splits[$indx];
    }
    
    
    return $self->{_blkSz};
}
sub blkSzAt{
	my ( $self,$indx) = @_;
	return $self->{_blkSz}[$indx];
}
#accessor method for Blatline qstart
sub qstart {
    my ( $self, $qstart ) = @_;
    $self->{_qstart} = $qstart if defined($qstart);
    return $self->{_qstart};
}
#accessor method for Blatline tstart
sub tstart {
    my ( $self, $tstart ) = @_;
    my $strand =$_[2];
    my @splits = split(/,/,$tstart);
    pop(@splits); #newline causes issues
 

    my $splitslength = scalar(@splits);
    
    for (my $indx = 0; $indx < $splitslength; $indx++)
    {
    	$splits[$indx] = $splits[$indx] + 1+;
    }
    
    
    if ($strand eq "-"){
    	@splits=reverse(@splits);
    }
    
        for (my $indx = 0; $indx < $splitslength; $indx++)
    {
    	$self->{_tstart}[$indx] = $splits[$indx];
    }
    
    #$self->{_tstart} = $tstart if defined($tstart);
    

    
    return $self->{_tstart};
}
sub tstartAt{
	my ( $self,$indx) = @_;
	return $self->{_tstart}[$indx];
}
#accessor method for Blatline numexon
sub numexon {
    my ( $self, $numexon ) = @_;
    $self->{_numexon} = $numexon if defined($numexon);
    return $self->{_numexon};
}

#accessor method for Blatline genomicstart
sub genomicstart {
    my ( $self, $genomicstart ) = @_;
    $self->{_genomicstart} = $genomicstart if defined($genomicstart);
    return $self->{_genomicstart};
}

#accessor method for Blatline genomicend
sub genomicend {
    my ( $self, $genomicstart ) = @_;
    $self->{_genomicend} = $genomicstart if defined($genomicstart);
    return $self->{_genomicend};
}

#accessor method for Blatline exonAdd
sub exonAdd {
    my ( $self, $exon ) = @_;
    my $getridofbullshitindirection = $self->{_exon};
    my $nextIndx;
    if(!defined($getridofbullshitindirection))
    {
    	$nextIndx=0;
    }
    else
    {
    	$nextIndx = $self->{_numexon};
    }
       
    $self->{_exon}[$nextIndx] = $exon if defined($exon);
    
    $self->{_numexon} = $nextIndx + 1;
     
    return $self->{_exon};
}
#accessor method for intronAdd
sub intronAdd {
    my ( $self, $intron ) = @_;
    my $getridofbullshitindirection = $self->{_exon};
    my $nextIndx;
    if(!defined($getridofbullshitindirection))
    {
    	$nextIndx=0;
    }
    else
    {
    	$nextIndx = $self->{_numexon}-1; #we will have added exactly 1 additional exon before we save the intron
    }
    
    $self->{_intron}[$nextIndx] = $intron if defined($intron);
    
    return $self->{_intron};
}
#accessor method for Blatline intron
sub intron {
    my $self;
    return $self->{_intron};
}
#accessor method for Blatline exon
sub exon {
    my $self;
    return $self->{_exon};
}

sub matchsize {
    my ( $self, $matchsize ) = @_;
    $self->{_matchsize} = $matchsize if defined($matchsize);
    return $self->{_matchsize};
}

#accessor method for print
sub print {
	my $outputfilename;
	my $garbagefilename;
	$outputfilename = "FBXoutput_testing.txt";
	$garbagefilename ="garbage.txt";
	open (MYFILE, '>>'.$outputfilename);
	open (MYFILE2,'>>'.$garbagefilename);
    my ($self) = @_;
	my $lookupstart;
	my $lookupend;
	my $lookupname;
	$lookupstart=get_reference_start($self->name);
	$lookupend=get_reference_stop($self->name);	
	if ($self->genomicstart == $lookupstart && $self->genomicend == $lookupend){
		printf MYFILE (">%s|(%s)|Chr:%s|Exon#%s\n",$self->name, $self->strand, $self->chNum, $self->numexon);
    	$self->printExons();
    	$self->printIntrons();
    	$self->printAll();
    	$self->printmycodingsequence();
		}
	else {
		printf MYFILE2 (">%s|(%s)|Chr:%s\n",$self->name, $self->strand, $self->chNum);
		printf MYFILE2 ("Start:%s - End:%s\n LookupStart:%s -LookupEnd:%s\n", $self->genomicstart,$self->genomicend,$lookupstart,$lookupend);
	}
	

}

sub get_reference_start{
	my $name =$_[0];
	my $chromasome =$lookuphash{$name}[0];
	my $strand =$lookuphash{$name}[1];
	my $start;
	#if ($strand eq "plus"){
	#	$start =$lookuphash{$name}[2];
	#}
	#if ($strand eq "minus"){
	#	$start =$lookuphash{$name}[3];
	#}
	$start =$lookuphash{$name}[2];
	return $start;	
}

sub get_reference_stop{

	my $name =$_[0];
	my $chromasome =$lookuphash{$name}[0];
	my $strand =$lookuphash{$name}[1];
	my $stop;
	#if ($strand eq "plus"){
	#	$stop =$lookuphash{$name}[3];
	#}
	#if ($strand eq "minus"){
	#	$stop =$lookuphash{$name}[2];
	#}
	$stop=$lookuphash{$name}[3];
	return $stop;
}

#method for printing exons
sub printExons
{
    my ($self) = @_;
	my $idx;
    #printf MYFILE ("  Exons:\n");
    my $size = $self->{_numexon};
    my $exonstart;
    my $exonend;
    my $count=1;
    for($idx=0; $idx< $size; $idx++)
	{
		$exonstart=$self->{_tstart}[$idx];
		$exonend=$self->{_tstart}[$idx] + $self->{_blkSz}[$idx];
		printf MYFILE ("Exon%s\t(%s,%s)\t",$count, $exonstart,$exonend-1);
		printf MYFILE ("%s",$self->{_exon}[$idx]);
		printf MYFILE ("\n");
		$count++;
	}
}

sub printmycodingsequence{
	my ($self) = @_;
	my $idx;
    my $size = $self->{_numexon};
    my $exonstart;
    my $exonend;
    printf MYFILE ("CDS\t");
    for($idx=0; $idx< $size; $idx++)
	{
		$exonstart=$self->{_tstart}[$idx];
		$exonend=$self->{_tstart}[$idx] + $self->{_blkSz}[$idx];
		printf MYFILE ("%s",$self->{_exon}[$idx]);
	}
	printf MYFILE ("\n");
	
}

#method for printing introns
sub printIntrons
{
    my ($self) = @_;
	my $idx;
	my $intronstart;
	my $intronend;
    #printf MYFILE ("  Introns:\n");
    my $size = $self->{_numexon}-1;
    my $strand = $self->strand;
    my $count=1;
    for($idx=0; $idx< $size; $idx++)
	{
		if ($strand eq '+'){
			$intronstart=$self->{_tstart}[$idx]+$self->{_blkSz}[$idx]; # +1 OFSET so you dont grab last nucleotide
			$intronend=$self->{_tstart}[$idx+1] - 1;	
		}
		if ($strand eq '-'){
			$intronstart=$self->{_tstart}[$idx+1]+$self->{_blkSz}[$idx+1]; # +1 OFSET so you dont grab last nucleotide
			$intronend=$self->{_tstart}[$idx] - 1;
		}

		printf MYFILE ("Intron%s\t(%s,%s)\t",$count, $intronstart,$intronend);
		printf MYFILE ("%s",$self->{_intron}[$idx]);
		printf MYFILE ("\n");
	}
	if($size==0)
	{
		#printf MYFILE ("   none\n");
	}
}
#method for printing genomic sequence
sub printAll
{
    my ($self) = @_;
	my $idx;
	my $size;
    my $genomicstart;
    my $genomicend;
    my $strand;
    $strand =$self->strand;
    $size = $self->{_numexon}-1;
    $genomicstart=$self->genomicstart();
    $genomicend=$self->genomicend();
  	#my $chromnum=$self->chNum;
    #$genomicstart=$self->{_tstart}[0];
    #$genomicend=$self->{_tstart}[$size] + $self->{_blkSz}[$size]-1;
    printf MYFILE ("Genomic\t (%s,%s) \t", $genomicstart,$genomicend);
    #my $genomicsequence =getSeq($chromnum, $genomicstart,$genomicend, $strand);
    #printf MYFILE ("%s\n",$genomicsequence);
    # need to get the starts% and stops of the genomic sequence here, shouldnt be too bad
    #
    printf MYFILE ("%s",$self->{_exon}[0]);
    for($idx=0; $idx<$size; $idx++)
	{
		printf MYFILE ("%s%s",$self->{_intron}[$idx],$self->{_exon}[$idx+1]);
	}
	printf MYFILE ("\n");
}

1;

##Opens psl file
###########################
###### MAIN FUNCTION ######
###########################

my @lines = openpslfile("new_output1.psl");
my @BlatLines = parseLines(@lines);
LoadReference();
LoadDatabase("TAIR10_chr_all.fasta");
BlatlineMath(@BlatLines);

#need a subroutine to double check positions of each gene

PrintBLatLines(@BlatLines);

##########################
######## OTHERS ##########
##########################
sub openpslfile{
	open FILE, "<$_[0]" or die $!;
	my @lines =<FILE>;
	close FILE;
	return @lines;
}
sub LoadReference{
	my $bacfile = "all_fbox_bac_at_check_positions_.fasta"; #all_fbox_bac_at_check_positions_.fasta
	my @lines;

	open FILE, "<$bacfile" or die $!;
	while( <FILE> ) {
        next if /^(\s)*$/;  # skip blank lines
        chomp;              # remove trailing newline characters
        push @lines, $_;    # push the data line onto the array
    }
	close FILE;
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
	
	
}

sub parseLines{
	my @input = @_;
	my @BlatLines;
	my @formatted;
	my $count;
	my $myvar;
	#printf("input:%s \n\n", scalar(@input) );
	for($count=0; $count<scalar(@input); $count++){	
		@formatted = split (/\t/,$input[$count]);
		my $line = new BlatLine();
		$line->name($formatted[9]);
		my $offset =get_reference_start($formatted[9]);
		$line->strand($formatted[8]);
		$line->chNum($formatted[13]);
		$line->chStart($formatted[14]);
		$line->chEnd($formatted[15]);
		$line->blkCt($formatted[17]);
		$line->blkSz($formatted[18],$formatted[8]);
		$line->qstart($formatted[19]);
		$line->tstart($formatted[20],$formatted[8]);
		$line->matchsize($formatted[12]-$formatted[11]);
		$line->genomicstart($formatted[15]+1+$offset);
		$line->genomicend($formatted[16]+$offset);
		
		$BlatLines[$count]=$line;
		

		# 9 is name
		# 8 is plus or minus
		# 13 is Chr #
		# 14 is Chr Start
		# 15 is Chr End
		# 17 is Block Counts (number of blocks matched)
		# 18 is Block Sizes (or block lengths)
		# 19 is qStarts
		# 20 is tStarts
	}
	# Debuggning BLATLINE OBJECTS
	#printf("BlatLines:%s \n\n", scalar(@BlatLines) );
	#for($count=1; $count<=scalar(@BlatLines)&&$count<10; $count++)
	#{
	#	$BlatLines[$count]->print()
	#}
	return @BlatLines;	
}


sub PrintBLatLines{
	#printf("BlatLines:%s \n\n", scalar(@_) );
	my $count;
	#for($count=0; $count<=scalar(@_)&&$count<10; $count++) #debugging 
	for($count=0; $count<=scalar(@_)-1; $count++) #printall
	{	
		$_[$count]->print();
	}
}


## Goes through Each Blatline object and determines number of exons and introns and fixes the positioning
sub BlatlineMath{
	
	my @BlatLines = @_;
	my $lineidx;	
	for($lineidx = 0; $lineidx < scalar(@BlatLines); $lineidx++)
	#for($lineidx = 0; $lineidx < 10; $lineidx++) #used for debugging, only run through the first ten blatlines
	{
		my $ctidx;
		my $blatLine = $BlatLines[$lineidx];
		my $lines = $blatLine->blkCt;
		my $strand = $blatLine->strand;
		my $objname =$blatLine->name;
		################################################### This code looks at number of exons and only adds introns if there is more than one exon
		my $prevexon = getSeq($blatLine->chNum, $blatLine->tstartAt(0), $blatLine->tstartAt(0)+$blatLine->blkSzAt(0)-1, $strand);
		$blatLine->exonAdd($prevexon);
		for($ctidx = 1; $ctidx < $lines; $ctidx++)
		{
			#exon position
			#intron position
			my $exonpos_start=$blatLine->tstartAt($ctidx); #+1 OFSET so you dont grab last nucleotide
			my $exonpos_end=$blatLine->tstartAt($ctidx)+$blatLine->blkSzAt($ctidx)-1;
			my $intronpos_start;
			my $intronpos_end;
			if ($strand eq '+'){
				$intronpos_start=$blatLine->tstartAt($ctidx-1)+$blatLine->blkSzAt($ctidx-1);#+1 OFSET so you dont grab last nucleotide
				$intronpos_end=$blatLine->tstartAt($ctidx)-1;
			}
			if ($strand eq '-'){
				$intronpos_start=$blatLine->tstartAt($ctidx)+$blatLine->blkSzAt($ctidx);#+1 OFSET so you dont grab last nucleotide
				$intronpos_end=$blatLine->tstartAt($ctidx-1)-1;
			}			
			my $curexon = getSeq($blatLine->chNum, $exonpos_start, $exonpos_end, $strand); 
			my $intron = getSeq($blatLine->chNum,$intronpos_start,$intronpos_end, $strand);
			$blatLine->intronAdd($intron);
			$blatLine->exonAdd($curexon);
			$prevexon = $curexon;
		}
		###################################################
		### Return the exon and intron to the proper blatline object (using the index of what element we are processing)
		$BlatLines[$lineidx] = $blatLine;
		
	}
	return @BlatLines;
}

sub LoadDatabase{
	$db = Bio::DB::Fasta->new($_[0]);
}

sub getSeq{
	##Get Sequence from Chromasome
	my $chromosome = $_[0];
	my $start = $_[1];
	my $end = $_[2];
	my $strand=$_[3];
	my $myseq;
	if ($strand eq "+"){
		$myseq =$db->seq($chromosome, $start => $end); ## net to fix int for $chromasome (mitochondira and chloroplast are not ints and are viable input)
	}
	if ($strand eq "-"){
		$myseq =$db->seq($chromosome, $end => $start); ## net to fix int for $chromasome (mitochondira and chloroplast are not ints and are viable input)
	}
	return $myseq;
}

