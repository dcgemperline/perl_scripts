#!/usr/bin/perl

##blat command that I ran

# blat TAIR9_chr_all.fasta  at_cds_list.fasta -out=psl new_output1.psl -maxGap=3 -minIdentity=95 -noHead

#Dependencies
use Bio::DB::Fasta;
use Bio::SeqIO;

#class BlatLine
package BlatLine;


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
        _fixedstartpos	=> undef,
        _fixedsizepos	=> undef
        
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
    $self->{_blkSz} = $blkSz if defined($blkSz);
    return $self->{_blkSz};
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
    $self->{_tstart} = $tstart if defined($tstart);
    return $self->{_tstart};
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
    return $self->{_intron};
}
#accessor method for Blatline exon
sub exon {
    return $self->{_exon};
}

#accessor method for fixedpos
sub fixedstartpos {
    my ( $self, $tstart ) = @_;
    $self->{_fixedstartpos} = $fixedstartpos if defined($fixedstartpos);
    return $self->{_fixedstartpos};
}

#accessor method for fixedsize
sub fixedsizepos {
    my ( $self, $tstart ) = @_;
    $self->{_fixedsizepos} = $fixedsizepos if defined($fixedsizepos);
    return $self->{_fixedsizepos};
}

#accessor method for print
sub print {
    my ($self) = @_;

    #print Blatline info
    printf("Name (strand):%s (%s)\n",$self->name, $self->strand);
    printf("  Chromosome Start-End (Number): %s - %s ( %s )\n",$self->chStart, $self->chEnd,$self->chNum );
    printf("  Block Count: %s\n",$self->blkCt);
    printf("  Block Size: %s\n",$self->blkSz);
    printf("  qStart: %s\n",$self->qstart);
    printf("  tStart: %s",$self->tstart);
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
    printf("  Exons:\n");
    my $size = $self->{_numexon};
    for($idx=0; $idx< $size; $idx++)
	{
		printf("   %s",$self->{_exon}[$idx]);
		printf("\n");
	}
}

#method for printing introns
sub printIntrons
{
    my ($self) = @_;
	my $idx;
    printf("  Intron:\n");
    my $size = $self->{_numexon}-1;
    for($idx=0; $idx< $size; $idx++)
	{
		printf("   %s",$self->{_intron}[$idx]);
		printf("\n");
	}
	if($size==0)
	{
		printf("   none\n");
	}
}
#method for printing genomic sequence
sub printAll
{
    my ($self) = @_;
	my $idx;
    
    printf("  Genomic Sequence:\n");
    printf("   %s",$self->{_exon}[0]);
    
    my $size = $self->{_numexon}-1;
    for($idx=0; $idx<$size; $idx++)
	{
		printf("%s%s",$self->{_intron}[$idx],$self->{_exon}[$idx+1]);
	}
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
	printf("BlatLines:%s \n\n", scalar(@_) );
	#for($count=0; $count<=scalar(@_)&&$count<10; $count++) #debugging 
	for($count=0; $count<=scalar(@_)-1; $count++) #printall
	{
		$_[$count]->print();
	}
}

## Correcting the Index of Sequence Position Because Blat uses 0 based index;
sub fixstarts{
	my @starts = split (/,/,$_[0]);
	pop (@starts);
	my $length =scalar(@starts);
	my @newstarts;
	my $incstart;
	my $getstart;
	for my $mathidx (0.. $length-1)
	{
		$getstart = $starts[$mathidx];
		$incstart=$getstart +1;
		$newstarts[$mathidx]=$incstart;
	}
	return @newstarts;
	#additionally need to fix the last exon length to -1 so that the offset will be correct
	#might want to rewrite this as a subfunction to deal with chromosome coordintes instead of blat coordinates for plus and minus strands before calling getseq
	
}
## fixes last element in sizes to be in parity with the exon start index
sub fixsizes{
	my @sizes = split (/,/,$_[0]);
	my $sizelength = scalar(@sizes);
	my @newsizes;
	#To get the proper lengths for the last exon you need to change size -1
	for my $sizeidx (0.. $sizelength-1)
	{
		if ($sizeidx!=$sizelength-1){ ## as long as this isnt the last size element just populate newsizes
			$getsize = $sizes[$sizeidx];
			$incstart=$getsize;
			$newsizes[$sizeidx]=$incstart; 
			}
		# if it is the last size then decrement by 1 to get proper exon length
		else 
		{
		$getsize = $sizes[$sizeidx];
		$incstart=$getsize-1;
		$newsizes[$sizeidx]=$incstart;
		}
		 
	}
	return @newsizes;
}

## Goes through Each Blatline object and determines number of exons and introns and fixes the positioning
sub BlatlineMath{
	
	my @BlatLines = @_;
	my $lineidx;
	for($lineidx = 0; $lineidx < scalar(@BlatLines); $lineidx++)
	#for($lineidx = 0; $lineidx < 10; $lineidx++) used for debugging, only run through the first ten blatlines
	{
		my $ctidx;
		my $blatLine = $BlatLines[$lineidx];
		my $lines = $blatLine->blkCt;
		my $startstr = $blatLine->tstart;
		my $sizestr = $blatLine->blkSz;
		my $strand = $blatLine->strand;
		my @starts = fixstarts($startstr);
		my @sizes= fixsizes($sizestr);
		# might want to write these @starts and @sizes out to the blatline object
		
		################################################### This code looks at number of exons and only adds introns if there is more than one exon 
		#need to get intron and exon positions at each step and store them into an object for each blatline
		my $prevexon = getSeq($blatLine->chNum, $starts[0], $starts[0]+$sizes[0]);
		$blatLine->exonAdd($prevexon);
		for($ctidx = 1; $ctidx < $lines; $ctidx++)
		{
			#exon position
			#intron position
			my $exonpos_start=$starts[$ctidx]+1; #+1 OFSET so you dont grab last nucleotide
			my $exonpos_end=$starts[$ctidx]+$sizes[$ctidx];
			my $intronpos_start=$starts[$ctidx-1]+$sizes[$ctidx-1]+1;#+1 OFSET so you dont grab last nucleotide
			my $intronpos_end=$starts[$ctidx];
			my $curexon = getSeq($blatLine->chNum, $exonpos_start, $exonpos_end); 
			my $intron = getSeq($blatLine->chNum,$intronpos_start,$intronpos_end);
			$blatLine->intronAdd($intron);
			$blatLine->exonAdd($curexon);
			
			$prevexon = $curexon;
		}
		###################################################
		### Return the exon and intron to the proper blatline object (using the index of what element we are processing)
		$blatLine->fixedstartpos(@starts);
		$blatLine->fixedsizepos(@sizes);
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