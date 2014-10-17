#**********************USAGE****************************************#
# prompt> perl RIC.pl <sequence file> PERM_MODE_YES_OR_NO NUM_PERMS #
#                                                                   #
#*******************************************************************#

use strict;
use File::Basename;

my @model;
my %scores;
my @seqs;
my $numseqs;
my @colcounts;
my @ricseq;
my $rstype;
my $outfile;
my $searchfile;
my $tempseqfile;
##my $revcomp_stat;
#my $revcomp_arg;
my $perm_arg;
my $perm_stat;
my @exts = qw(.txt .fa .fasta);
my $num_perms=1;

#Read in the sequence(s) to be searched for RIC scoring.
$searchfile = $ARGV[0];

if ($searchfile eq ""){
    print "Enter the name of the file containing the sequences to search: ";
    $searchfile = <STDIN>;
}


$perm_arg=$ARGV[1];
if(!( ($perm_arg=~m/^yes$/i) || ($perm_arg=~m/^no$/i) ) )
	{
	die "Error, expect permutation argument either 'yes' or 'no' !\n";
	}
else
	{
	if($perm_arg=~m/^yes/i)
		{
		$perm_stat=1;
		}
	else
		{
		$perm_stat=0;
		}
	}


$num_perms=$ARGV[2];
if(!($num_perms=~m/^[0-9]+$/))
	{
	die "Invalid 'num_perms' value $num_perms !\n";
	}
if($num_perms<1)
	{
	die "Invalid 'num_perms' value  $num_perms . Must be a positive integer!\n";
	}




#my $h="hello";
#my $ph=&PERMUTE_STR($h);
#my %hchello=&getCharCountHash($h);
#for my $k  (keys(%hchello))
#	{
#	print "$k is a key in the map its value is  ".$hchello{$k}."\n";
#	}
#my $pcRes=verifyPermutation($h,$ph);
#print "The pcres is $pcRes\n";
#$ph="hell";
#$pcRes=verifyPermutation($h,$ph);
#print "The pcres is $pcRes\n";
#die "perm_mode=".$perm_stat." and num_perms=".$num_perms." AND PH=$ph\n";



$searchfile =~ s/\n//g;
$outfile = $searchfile;

#Calculate Scores for 12RSS
$outfile = basename($searchfile,@exts).".12RSS.scores";
&GETMODEL("MM12.model");
$numseqs = 0;
&GETSEQS("MM12RSS.fasta");
&INITMODEL;
&RICSCORE(-45);

#Calculate Scores for 23RSS
@seqs=();
$outfile = basename($searchfile,@exts).".23RSS.scores";
&GETMODEL("MM23.model");
$numseqs = 0;
&GETSEQS("MM23RSS.fasta");
&INITMODEL;
&RICSCORE(-65);


#________end of program________

#http://bioinfo2.ugr.es/documentation/Perl_Cookbook/ch04_18.htm
# fisher_yates_shuffle( \@array ) : generate a random permutation
# of @array in place
sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}



sub PERMUTE_STR {
	my $inStr=$_[0];
	#print "instr is $inStr\n";
	#FROM http://perlmeme.org/faqs/manipulating_text/string_characters.html
	my @chars = map substr( $inStr, $_, 1), 0 .. length($inStr) -1;
	#for(my $c=0;$c<scalar(@chars);$c++)
	#	{
	#	print "chars[c] = ".$chars[$c]."\n";
	#	}
	fisher_yates_shuffle( \@chars );    # permutes @array in place
	#print "....permute....\n";
	#for(my $c=0;$c<scalar(@chars);$c++)
	#	{
	#	print "chars[c] = ".$chars[$c]."\n";
	#	}
	my $retStr=join('',@chars);
	#print "The retStr is $retStr\n";
	return $retStr;
	}



sub getCharCountHash {
	my $inStr=$_[0];
	#print "Input to getCharCountHash is $inStr\n";
	my @chars = map substr( $inStr, $_, 1), 0 .. length($inStr) -1;
	my %hash_ref ;
	for(my $c=0;$c<scalar(@chars);$c++)
		{
		#print "In init looking at ".$chars[$c]."...\n";
		if( exists  $hash_ref{ $chars[$c] })
			{
			$hash_ref{$chars[$c]}=$hash_ref{$chars[$c]}+1;
			}
		else
			{
			$hash_ref{$chars[$c]}=1;
			}
		#print "The value is now ".$hash_ref{$chars[$c]}."\n";
		}
	return %hash_ref
	}


sub verifyPermutation {
	#VERIFY PERMUTATION BOTH-WAYS
	my $p1=$_[0];
	my $p2=$_[1];
	if($p1 eq $p2)
		{
		#the two strings shouldn't be equal, but one should be a permutation of the other!
		return 0;
		}
	#print "Verified that $p1  not eq $p2\n";
	my $firstWay=verifyPermutationOW($p1,$p2);
	my $secondWay=verifyPermutationOW($p2,$p1);
	if($firstWay==1 && $secondWay==1)
		{
		return 1;
		}
	else
		{
		return 0;
		}
	}


sub verifyPermutationOW {
	#VERIFY PERMUTATION ONE-WAY
	my $p1=$_[0];
	my $p2=$_[1];
	#print "In VP p1 is $p1 and p2 is $p2\n";
	my %m1=&getCharCountHash($p1);
	my %m2=&getCharCountHash($p2);
	my $nk1=scalar(keys(%m1));
	my $nk2=scalar(keys(%m2));

	#print "\n\n\n";
	#for my $k  (keys(%m1))
	#	{
	#	print "$k is a key in the m1 map its value is  ".$m1{$k}."\n";
	#	}
	#print "\n\n\n";

	if($nk1!=$nk2)
		{
		#mismatch in size of key sets!
		#print "mismatch in key sizes!\n";
		return 0;
		}

	for my $k  (keys(%m1))
		{
		my $m1val=$m1{$k};
		if(exists  $m2{$k})
			{
			my $m2val=$m2{$k};
			if($m2val!=$m1val)
				{
				#value mismatch!
				return 0;
				}
			}
		else
			{
			#doesn't exist!
			return 0;
			}
		}
	return 1;
	}




sub GETMODEL{
    my ($modelfile);
    my ($tempmodel);
    my ($tempsum);
    my (@tempcols);
    my ($collength);
    my($i);
    if($_[0] eq ""){
	print "Models must be provided in a file with the filename given at runtime, or you can enter the model here.\n";
	print "Models must be in the form (x,y,z)(a,b)...etc.\n";
	print "Enter a model here: ";
	$tempmodel = <STDIN>;
    } else {
	$modelfile = $_[0];
	open(FIN,"<$modelfile");
	print "Opening $modelfile\n";
	$tempmodel = <FIN>;
	close(FIN);
    }
    
    $tempmodel =~ s/\n//g;
    $tempmodel =~ s/\r//g;
    $tempsum = 0;
    @model = split(/\)/,$tempmodel);
    foreach (@model){
	$_ =~ s/\(//g;
	if(/,/){
	    @tempcols = split(/,/,$_);
	    foreach (@tempcols){
		$tempsum += $_;
	    }
	} else {
	    $tempsum += $_;
	}
    }

    if(basename($modelfile,".model") eq "MM12"){
	$rstype = 12;
	$collength = 28;
    } else {
	$rstype = 23;
	$collength = 39;
    }
    
    for($i=0;$i<$collength;$i++){
	$colcounts[$i] = 0;
    }
}

sub GETSEQS{
    my($tempfile);
    $tempfile = $_[0];
    open(FIN,"<$tempfile")
	or die "Cannot open $tempfile. Program will terminate.\n";
    while(<FIN>){
	if(/^>/){
	    next;
	}
	$_ =~ s/\s//g;
	$_ =~ s/\n//g;
	$_ =~ s/\r//g;
	$_ =~ tr/[ACGTN]/[acgtn]/;
	push(@seqs,$_);
	$numseqs += 1;
    }
    close(FIN);
    print "$numseqs sequences read in.\n";
}

sub INITMODEL{
    (my @positions);
    (my $tempval);
    (my $m);
    (my $s);
    (my $x);
    (my $y);
    (my $pos);
    (my $position);
    (my @sequence);
    (my $comb);
    (my $key);
    (my $value);
    print "Initializing model for ".$rstype."rs's . . .\n";
    foreach $m (@model){
	@positions = split(/,/,$m);
	@positions = sort {$x<=>$y}@positions;
	foreach $s (@seqs){
	    $comb = "";
	    @sequence = split(//,$s);
	    foreach $pos (@positions){
		$comb .= $sequence[$pos-1];
		$position = $pos;
	    }
	    $comb = $position.$comb;
	    if($comb =~ /n|\./){
		$colcounts[$position-1] -= 1;
	    }
	    if(exists $scores{$comb}){
		$scores{$comb} = $scores{$comb} + 1;
	    } else {
		$scores{$comb} = 1;
	    }
	}
    }
    print "Initialization complete.\n";
}





sub RICSCORE{
    (my $a);
    (my $x);
    (my $y);
    (my $m);
    (my $q);
    (my $pos);
    (my $lastpos);
    (my $curr);
    (my $class);
    (my $counter);
    (my $lookup);    
    (my $prob);
    (my $totalprod);
    (my $phony);
    (my $end);
    (my @positions);
    (my $searchspace);
    (my $conserved);
    (my $fileempty);
    (my $lcount);
    (my $tempin);
    (my $fastaheader);
    $a = 2;

    if($rstype == 12){
		$searchspace = 28;
    	} else {
		$searchspace = 39;
    }
    
    for($counter=0;$counter<$searchspace;$counter++){
		$colcounts[$counter] += $numseqs;
    }
    
    print "Computing RIC scores for $searchfile.\nSaving to $outfile. . .\n";
    my $permID=0;
    if($perm_stat==0) 
	{
		$num_perms=1;
	    open(FOUT,">$outfile");
	}
    else
	{

	    $outfile.=".gz";
	    open(FOUT,"|gzip -c >> $outfile");
	}

    
    for(my $permID=1;$permID<=$num_perms;$permID++)
	{
	    my $pct_cmplt=(($permID*1.0)/($num_perms*1.0))*100.0;
	    print STDERR "Running with permutation ID=".$permID." of ".$num_perms." ; pct = ".$pct_cmplt."\n";

	    open(ISEQ, "<$searchfile");
	    my $tempPath="/dev/shm/$outfile.temp.txt";
	    #my $numTempPathIO=0;

	    while($tempin = <ISEQ>){
		
			if (-e $tempPath){
				unlink($tempPath);
			}
			open O, ">$tempPath";
			$fastaheader = $tempin;
			$fastaheader =~ s/\n//g;
			$fastaheader =~ s/\r//g;
		unless($tempin =~ /^>/){
				print "Invalid input file. Input must be in FASTA format.\n";
				exit;
			}
			$tempin = <ISEQ>;
			print O "$tempin";
			close (O);
		
			open(FIN,"<$tempPath")
			or die "Cannot open searchfile $tempPath.\n";
		#$numTempPathIO++;
	    	$fileempty = 0;
			$phony = 0 ;
			$tempin = <FIN>;
			$tempin =~ s/\n//g;
			$tempin =~ s/\r//g;
			$tempin =~ tr/[A-Z]/[a-z]/;
			if($perm_stat==1)
				{
				#permute here
				my $permuted=&PERMUTE_STR($tempin);
				#my $pcRes=verifyPermutation($permuted,$tempin);
				#if($pcRes==0)
				#	{
				#	die "Apperently string 1 '$tempin' didn't get correctly permuted to '$permuted' !\n";
				#	}
				$tempin=$permuted;
				}
			@ricseq = split(//,$tempin);

			while($fileempty==0){
				$totalprod = 0;
				if(scalar(@ricseq) < $searchspace){last;}
				$phony += 1;
				$conserved = $ricseq[0].$ricseq[1];
				if($conserved eq "ca"){
				   
				    foreach $m (@model){
						$lookup = "";
						@positions = split(/,/,$m);
						@positions = sort {$x<=>$y}@positions;
						$class = scalar(@positions);
						foreach $pos (@positions){
						    $curr = $ricseq[$pos-1];
						    $lookup .= $curr;
						    $lastpos = $pos;
						}
						$lookup = $lastpos.$lookup;
						if(exists $scores{$lookup}){
						    $q = $scores{$lookup};
						} else {
							$q = 0;
						}
						$prob = ($q + ($a/(4**$class)))/($colcounts[$lastpos-1] + $a);
						$prob = log($prob);
						$totalprod = $totalprod + $prob;
				    } #for each model
				    
				    $end = $phony+$searchspace-1;

				    if($totalprod >= $_[0]){
					    print FOUT "$fastaheader\t";
					    print FOUT "$phony\t$end\t";
					    for($counter=0;$counter<$searchspace;$counter++){
							print FOUT "".$ricseq[$counter];
					    }
					    print FOUT "\t$totalprod";
					    if($perm_stat==1)
							{
							print FOUT "\t$permID";
							}
					    print FOUT "\n";
				    } #if total prod >= 

				}
				shift(@ricseq);
			} #end of while($fileempty==0)
		    
		if (-e $tempPath){
		   unlink($tempPath);
		}
	
		}#end of while(loopthroughtempseq)
	    close(ISEQ);








	} # END perm ID loop

    close(FOUT);

    #print "numTempPathIO is $numTempPathIO\n";
}








