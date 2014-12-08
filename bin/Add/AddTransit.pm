
package AddTransit;

# This package reads in fasta file and runs TargetP to find mitochondrial transit peptides and signal peptides for proteins
#       secreted by the classical pathway (using SignalP) - the cleave site is based on HMM
#
# It also find the length of the sequence - which is added to the description
# For proteins with a transit/signal peptide this routine
# returns the IPI entry with an additional N-terminal peptide that is added to the C-terminal after a J
# The J is introduced as a theoretical AA - and the idea is to set trypsin up to cleave also before and after J
# For proteins without a N-terminal signal peptide, only the length is added to the header
#
# Example of a header after introduction or a new N-terminal peptide - at the C-terminus
# >IPI00000024.1|SWISS-PROT:Q08174-1  lng=1026 Tax_Id=9606 Isoform 1 of Protocadherin-1 precursor M[1028,T,36,T]
# In this case the protein is 1026 residues long and have a peptide added starting at position 1028 (residue T) while
# this corresponds to an original position starting at 36
# Again this means that the transit peptide is found in position 1-35
#
# Alternatively for predicted signal peptide
# Example of a header after introduction or a new N-terminal peptide - at the C-terminus
# >IPI00000111.6|SWISS-PROT:Q7L1S5-1|ENSEMBL:ENSP00000284224|REFSEQ:NP_113610|H-INV:HIT000252532|VEGA:OTTHU
# MP00000073066 Tax_Id=9606 GalNAc-4-sulfotransferase 2  lng=443 # SP[445,L,25,L]#
# In this case the original protein is 443 residues long and have a peptide added to start at position 445
# (residue L), which corresponds to the original position 25
# Again this means that the signal peptide is found in position 1-24
# The # indicate that the sequence has been modified
# A second # shows that the description has not been truncated - but is not added until later
#
# Example of sequence:
#MQPSEMVMNPKQVFLSVLIFGVAGLLLFMYLQVWIEEQHTGRVEKRREQKVTSGWGPVKYLRPVPRIMSTEKIQEHITNQNPKFHMPEDVREKKENLLLNSERST
#RLLTKTSHSQGGDQALSKSTGSPTEKLIEKRQGAKTVFNKFSNMNWPVDIHPLNKSLVKDNKWKKTEETQEKRRSFLQEFCKKYGGVSHHQSHLFHTVSRIYVED
#KHKILYCEVPKAGCSNWKRILMVLNGLASSAYNISHNAVHYGKHLKKLDSFDLKGIYTRLNTYTKAVFVRDPMERLVSAFRDKFEHPNSYYHPVFGKAIIKKYRP
#NACEEALINGSGVKFKEFIHYLLDSHRPVGMDIHWEKVSKLCYPCLINYDFVGKFETLEEDANYFLQMIGAPKELKFPNFKDRHSSDERTNAQVVRQYLKDLTRT
#ERQLIYDFYYLDYLMFNYTTPFLJLLLFMYLQVWIEEQHTGRVEK
# Note the J and the following sequence VPKFDQNLDTK
#

##################################################
### - Lulu Zheng - March 2013 #####
##################################################

{
    # Global variable to keep count of existing objects
    my $_count = 0;
	# Manage the count of existing objects
    sub get_count {
        $_count;
    }
    sub _incr_count {
        ++$_count;
    }
    sub _decr_count {
        --$_count;
    }
}

sub new{
	my ($class, %args) = @_;
	my $self = bless{
		_file => $args{-file},
		_targetp => $args{-targetp} || 'N',
	},$class;
	
	$class->_incr_count();
	return $self;
}

sub AddTransit{
	my ($self) = @_;
	
	# Clean up tmp files from last run
	unlink "../tmp/trfasta.tmp" if -e "../tmp/trfasta.tmp";

	# It is important that the tmp dir is setup for targetp to be the same as for SAPshodt
	# (see installation notes)
	unlink glob "../tmp/target*/*" || 1;
	rmdir glob <../tmp/target*> || 1;
	# (Important since TargetP can skip sequences if these temporary directories
	#       e.g. after an interrupted run has not been removed)

	open FSA, $self->{_file} or die; # we reformat fasta to make it more easy to handle
	$preseq = "";
	while (<FSA>) {
		X: if (/^>([\w\-]+)(\s.*)$/) {
				$pretag=$1;
				$description=$2;
				while ($_ = <FSA> and !/^>/) {
					chomp;
					$preseq.=$_;
				}
				open(MYNEWFILE, ">>../tmp/trfasta.tmp"); #open for write
				print MYNEWFILE ">$pretag$description\t$preseq\n"; #printing the sequence in one line
					$pretag="";
					$preseq="";
					$description="";
				if (/^>/) {
					goto X;
				}
			}
        close(MYNEWFILE); #*** Close the file ***
	}
	close FSA;

	#Now read the formatted fasta file
	open(MYFASFILE, "../tmp/trfasta.tmp"); #open to read the nonfasta lines
	my @filenames = ();
	my $outfile;
	if($self->{_file} =~ /\//){ @filenames = split /\//,$self->{_file}; $outfile = "tr_sp_".$filenames[@filenames-1];}
	else{ $outfile = "tr_sp_".$self->{_file}; }
	open $trfile,">../tmp/$outfile";
	while (<MYFASFILE>) {
		if (/^>([\w\-]+)(\s.*)\s(OS=Homo sapiens.*)\t(.*)$/) {
			$tag=$1;
			# We split the description in two to be able to introduce the length at a early stage
			# The original length is namely the easy way to spot whether any of the novel peptides
			# has been identified in results from search engine
			$desc1=$2;
			$desc2=$3;
			$seq=$4;
			$lng=(length $seq);                                          
		}
		# we do not analyse all those shorter than 60
		# But we take N-terminal sequences up to 300 to fetch cases where a tryptic peptide with one miscleave
		#  is longer
		if ($seq =~ /^(.{60,300}).*$/) {                        
			$ntermseq=$1;                        
			&predict_target_and_output;                        
		}
		else {
		# else it is shorter than 60 AA - and hence not analysed
		# Now we print the complete original entry with only the addition of the sequence length to the header
			print $trfile ">$tag$desc1 $desc2 lng=$lng\n$seq\n";
		}
		$tag="";
		$desc1="";
		$desc2="";
		$seq="";
		$lng="";
		$to="";
		$newstart="";
		$oldfirstAA="";
		$newfirstAA="";
	}
	close(MYFASFILE); #*** Close the file ***
	close($trfile);

	unlink glob "../tmp/*tmp";  #cleaning up the last tmp files

	
	sub   predict_target_and_output {
	# This routine:
	# 1)makes the input file to TargetP
	# 2)run TargetP
	# 3)parses the TargetP output
	# 4)create novel peptides for proteins having transit or signal peptides
	# 5)print out the sequences (and if existing: novel peptides and extra description)
        # First the sequence is put into a file
         open(MYINFILE, ">../tmp/targetpinput.tmp"); #open for write
         print MYINFILE ">$tag\n$ntermseq\n";
         close(MYINFILE); #*** Close the file ***

         # Run targetp and catch output  - this way only running one sequence at a time (room for optimization)
         # Cut of for mitochondria set to 0.65 to give a specificity of minimum 0.90 (p 1009 in Emanuelsesson et al. 2000 J mol Biol 300, 1005-1016)
         # Note, that only reliability class 1 and 2 are included anyway (hence no cut of for signal peptides)
		 if($self->{_targetp} eq 'Y'){
			system "../bin/targetp-1.1/targetp -N -c -t 0.65 ../tmp/targetpinput.tmp > ../tmp/targetpres.tmp";
		 }
		 else{
			print $trfile ">$tag$desc1 $desc2  lng=$lng\n$seq\n";
			return;
		 }
		 
         # Open the output file, divide into bad prediction (print entry without changes) and  good preditions
         #      for the good prediction the novel N-terminal tryptic peptide is created
         open(MYRESFILE, "../tmp/targetpres.tmp") or print $trfile ">$tag$desc1 $desc2  lng=$lng\n$seq\n";  #We print the record if there is no outputfile
                                                                                           # from TargetP, which happens when TargetP forget
                                                                                           # to clean up its temporary files
               while (<MYRESFILE>) {
                    #matching  Name              Len            mTP     SP  other  Loc  RC  TPlen
                    #matching  IPI00000012.4     155          0.058  0.979  0.012   S    1      -
                    if (/^([\w\-]+)\s+\d+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\d+)\s+(\S+)/) {   #pick the relevant results
                     $tag=$1;
                     $prediction=$2;
                     $reliability=$3;
                     $cleavesite=$4;

                     #Print those with low reliability with no cleavesite or without mitochondria or signal peptide prediction
                        if(($reliability >= 3) || ($cleavesite eq "-") || ($prediction eq "-") || ($prediction eq "*"))  {
                        # Here we print the output for sequences that went through a transit peptide prediction but were negative
                        # we add nothing except the length
							print $trfile ">$tag$desc1 $desc2 lng=$lng\n$seq\n"; #output original entry for sequences without transit peptide
						}
						else {
                        # now we need to pull out the tryptic peptide at the cleavesite position from $seq
                        # We go for the first tryptic peptide with one miscleavage after the transit/signal peptide
							if ($seq =~ /^.{$cleavesite}(.*?[K|R])([^P].*?[K|R])[^P].*$/) {
								$newntermpep=$1;
								$newntermpep2=$2;

                                # Here we catch cases where the second peptide contain more than one K/R
                                #  which can happen in the regular expressio above if first character in second peptide is K/R
								if ($newntermpep2 =~ /^([K|R])[^P]/) {
                                $corrected=$1;
                                $newntermpep2=$corrected;
                                }

                        #Now we join the two tryptic peptides to create the full miscleaved peptide
                        $newntermseq = $newntermpep .  $newntermpep2;

                        # We open a file for making a register of novel peptides
                        open(MYPEPREGFILE, ">>../tmp/PeptideRegister.txt");

                        #We print the right tag for [TR] mitochondria transit or [SP] signal peptide
                               # Here we print the output for sequences that went through a transit peptide prediction

                                # The FROM position is the $newntermseq (position after cleavage site)
                                $from=$cleavesite + 1;
                                # The TO positon is the length of the peptide plus start position
                                $to=(length $newntermseq)+ $cleavesite;   #May become usefull in later version
                                #The original start position of the N-terminal is 1
                                # The new start position of the N-terminal is
                                $newstart = (length($seq) +2); # second position after $seq because we add a J as separater
                                # The new start AA is
                                $newfirstAA  = substr($newntermseq,0,1);
                                
                        if ($prediction eq "M") {
                                #Note that the first AA does not change hence both are newfirstAA
                                print $trfile ">$tag$desc1 $desc2 lng=$lng # M[$newstart,$newfirstAA,$from,$newfirstAA]\n";
                                # Here we print to the file of all the introduced peptides
                                print MYPEPREGFILE "$tag $newntermpep $newntermseq $from $to transit $newstart $newfirstAA $from\n";
                                }
                        if ($prediction eq "S") {
                                print $trfile ">$tag$desc1 $desc2 lng=$lng # SP[$newstart,$newfirstAA,$from,$newfirstAA]\n";
                                print MYPEPREGFILE "$tag $newntermpep $newntermseq $from $to Signal $newstart $newfirstAA $from\n";
                                }
                                # For those that had a transit peptide the new terminal peptide is included at the C-terminus
                                # We use the theoretic Amino acid J and require trypsin to cleave before and after J
                                # We add the original position of the new N-term peptide that is added
                                print $trfile "$seq";
                                print $trfile "J$newntermseq\n";

                        close(MYPEPREGFILE); # close peptide register file

                        }
                        else{
							# if we can not find a tryptic peptide after the new N-terminus we print record without change
							print $trfile ">$tag$desc1 $desc2 lng=$lng\n$seq\n";
                        }

                     }
                     }
               }
         close(MYRESFILE); #*** Close the file ***

            #Make sure to reset
              $lng="";
              $ntermseq="";
              $cleavesite="";
              $newntermseq="";
              $prediction="";
              $to="";
              $from="";
              $newntermseq="";
              $reliability="";
              $pretag="";
              $newstart="";
              $oldfirstAA="";
              $newfirstAA="";
	}

}

sub DESTROY {
    my($self) = @_;
    $self->_decr_count();
}

1;


__END__
