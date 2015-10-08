#!/usr/bin/perl -w
use strict;
##blat command that I ran

# blat TAIR9_chr_all.fasta  at_cds_list.fasta -out=psl new_output1.psl -maxGap=0 -minIdentity=100 -noHead

#Dependencies
use Bio::DB::Fasta;
use Bio::SeqIO;
use Try::Tiny;

#class BlatLine
package BlatLine;

# Declare Database before use
my $db;

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
        _intron	   		=> undef
        
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
    my @splits = split(/,/,$blkSz);
    my $splitslength = scalar(@splits);
    
    for (my $indx = 0; $indx < $splitslength; $indx++) # There is aproblem here we are not subtracting from the last blocksize, we are still picking up the last blocksize, check on the desktop to see if this is correct, I beleive it is on the desktop
    {
    	if($indx == $splitslength-1)
    	{
    		$self->{_blkSz}[$indx] = $splits[$indx] - 1;
    	}
    	else
    	{
    		$self->{_blkSz}[$indx] = $splits[$indx];
    	}
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
    my @splits = split(/,/,$tstart);
    pop(@splits); #newline causes issues
    my $splitslength = scalar(@splits);
    
    for (my $indx = 0; $indx < $splitslength; $indx++)
    {
    	$self->{_tstart}[$indx] = $splits[$indx] + 1;
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

#accessor method for print
sub print {
	my $outputfilename;
	$outputfilename = "output_test.txt";
	open (MYFILE, '>>'.$outputfilename);
    my ($self) = @_;

    #print Blatline info
    printf MYFILE (">%s|(%s)|Chr:%s\n",$self->name, $self->strand, $self->chNum);
    #printf MYFILE ("  Chromosome Start-End (Number): %s - %s ( %s )\n",$self->chStart, $self->chEnd,$self->chNum );
    #printf MYFILE ("  Block Count: %s\n",$self->blkCt);
    #printf("  Block Size: %s\n",$self->blkSz);
    #printf MYFILE ("  qStart: %s\n",$self->qstart);
    #printf("  tStart: %s",$self->tstart);
    $self->printExons();
    $self->printIntrons();
    $self->printAll();
    printf("\n\n");
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
    for($idx=0; $idx< $size; $idx++)
	{
		$exonstart=$self->{_tstart}[$idx];
		$exonend=$self->{_tstart}[$idx] + $self->{_blkSz}[$idx];
		
		printf MYFILE ("Exon\t(%s,%s)\t",$exonstart,$exonend);
		printf MYFILE ("%s",$self->{_exon}[$idx]);
		printf MYFILE ("\n");
	}
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
    for($idx=0; $idx< $size; $idx++)
	{
		$intronstart=$self->{_tstart}[$idx]+$self->{_blkSz}[$idx] + 1; # +1 OFSET so you dont grab last nucleotide
		$intronend=$self->{_tstart}[$idx+1] - 1;
		printf MYFILE ("Intron\t(%s,%s)\t",$intronstart,$intronend);
		printf MYFILE ("%s",$self->{_intron}[$idx]);
		printf MYFILE ("\n");
	}
	if($size==0)
	{
		printf MYFILE ("   none\n");
	}
}
#method for printing genomic sequence
sub printAll
{
    my ($self) = @_;
	my $idx;
    
    printf MYFILE ("Genomic_Sequence\t");
    # need to get the starts and stops of the genomic sequence here, shouldnt be too bad
    #
    printf MYFILE ("%s",$self->{_exon}[0]);
    
    my $size = $self->{_numexon}-1;
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

my @lines = openpslfile("new_output.psl");
my @BlatLines = parseLines(@lines);
LoadDatabase("TAIR9_chr_all.fasta");
BlatlineMath(@BlatLines);
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


sub parseLines{
	my @input = @_;
	my @BlatLines;
	my @formatted;
	my $count;
	my $myvar;
	#printf("input:%s \n\n", scalar(@input) );
	for($count=0; $count<scalar(@input); $count++){	
		@formatted = split (/\t/,$input[$count]);
		my $line = new BlatLine() ;
		$line->name($formatted[9]);
		$line->strand($formatted[8]);
		$line->chNum($formatted[13]);
		$line->chStart($formatted[14]);
		$line->chEnd($formatted[15]);
		$line->blkCt($formatted[17]);
		$line->blkSz($formatted[18]);
		$line->qstart($formatted[19]);
		$line->tstart($formatted[20]);
		
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
	for($count=0; $count<=scalar(@_)&&$count<10; $count++) #debugging 
	#for($count=0; $count<=scalar(@_)-1; $count++) #printall
	{
		$_[$count]->print();
	}
}


## Goes through Each Blatline object and determines number of exons and introns and fixes the positioning
sub BlatlineMath{
	
	my @BlatLines = @_;
	my $lineidx;
	#for($lineidx = 0; $lineidx < scalar(@BlatLines); $lineidx++)
	for($lineidx = 0; $lineidx < 10; $lineidx++) #used for debugging, only run through the first ten blatlines
	{
		my $ctidx;
		my $blatLine = $BlatLines[$lineidx];
		my $lines = $blatLine->blkCt;
		#my $startstr = $blatLine->tstart;
		#my $sizestr = $blatLine->blkSz;
		my $strand = $blatLine->strand;
		#my @starts = fixstarts($startstr);
		#my @sizes= fixsizes($sizestr);
		# might want to write these @starts and @sizes out to the blatline object
		
		################################################### This code looks at number of exons and only adds introns if there is more than one exon 
		#need to get intron and exon positions at each step and store them into an object for each blatline
		my $prevexon = getSeq($blatLine->chNum, $blatLine->tstartAt(0), $blatLine->tstartAt(0)+$blatLine->blkSzAt(0));
		$blatLine->exonAdd($prevexon);
		for($ctidx = 1; $ctidx < $lines; $ctidx++)
		{
			#exon position
			#intron position
			my $exonpos_start=$blatLine->tstartAt($ctidx); #+1 OFSET so you dont grab last nucleotide
			my $exonpos_end=$blatLine->tstartAt($ctidx)+$blatLine->blkSzAt($ctidx);
			
			my $intronpos_start=$blatLine->tstartAt($ctidx-1)+$blatLine->blkSzAt($ctidx-1)+1;#+1 OFSET so you dont grab last nucleotide
			my $intronpos_end=$blatLine->tstartAt($ctidx)-1;
			
			my $curexon = getSeq($blatLine->chNum, $exonpos_start, $exonpos_end); 
			my $intron = getSeq($blatLine->chNum,$intronpos_start,$intronpos_end);
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
	my $myseq =$db->seq($chromosome, int($start) => int($end)); ## net to fix int for $chromasome (mitochondira and chloroplast are not ints and are viable input)
	return $myseq;
}

sub Output {
	## Generates Output in a Friendly File Format
	
	
	
	
}
