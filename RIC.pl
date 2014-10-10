#**********************USAGE**********************#
# prompt> perl RIC.pl <sequence file>             #
#                                                 #
#*************************************************#

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
my @exts = qw(.txt .fa .fasta);

#Read in the sequence(s) to be searched for RIC scoring.
$searchfile = $ARGV[0];
if ($searchfile eq ""){
    print "Enter the name of the file containing the sequences to search: ";
    $searchfile = <STDIN>;
}
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
$outfile = basename($searchfile,@exts).".23RSS.scores";
&GETMODEL("MM23.model");
$numseqs = 0;
&GETSEQS("MM23RSS.fasta");
&INITMODEL;
&RICSCORE(-65);


#________end of program________


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
    
    open(FOUT,">$outfile");
    open(ISEQ, "<$searchfile");
    my $tempPath="/dev/shm/$outfile.temp.txt"
    while($tempin = <ISEQ>){
		
		if (-e $tempPath){
			unlink($tempPath);
		}
		open O, ">$outfile.temp.txt";
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
    	$fileempty = 0;
		$phony = 0 ;
		$tempin = <FIN>;
		$tempin =~ s/\n//g;
		$tempin =~ s/\r//g;
		$tempin =~ tr/[A-Z]/[a-z]/;
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
			    }
			    
			    $end = $phony+$searchspace-1;

			    if($totalprod >= $_[0]){
				    print FOUT "$fastaheader\t";
				    print FOUT "$phony\t$end\t";
				    for($counter=0;$counter<$searchspace;$counter++){
						print FOUT $ricseq[$counter];
				    }
				    print FOUT "\t$totalprod";
				    print FOUT "\n";
			    }

			}
			shift(@ricseq);
		}#end of while($fileempty==0)
	    
        if (-e $tempPath){
           unlink($tempPath);
        }
	
	}#end of while(loopthroughtempseq)
    close(I);
    close(FOUT);
}
