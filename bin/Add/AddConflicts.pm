package AddConflicts;

# This package reads in fasta file and uses a preformatted Swiss-Prot conflict data file - parsed from Swiss-Prot
# This package is part of the framework creating a MS friendly version of uniprot fasta sequence database
#
# It creates tryptic peptides around the given CONFLICT sites and append the peptide at the C-terminal of the uniprot entry
#
# This script also introduces decoy fasta entries where each sequence entry is reversed and the identifier changed from uniprot to REV
#
# Peptides are appended with J as a divider
# The J is introduced as a theoretical AA - and the idea is to set trypsin up to cleave also before and after J
#
# Example of a header after introduction of a CONFLICT annotation, appended at the C-terminus (newline introduced here)
# In this example there is already peptides added one snp peptide in previous step
#>O60925 Tax_Id=9606 Prefoldin subunit 1 lng=122 #SNP[129,G,104,S]CON[142,V,86,I]#
#
# The # indicate that the sequence has been modified
# The second # shows that the description has not been truncated - and is added here
#
# A conflict peptide has been introduces and appended, so that the original position 86 (residue I)
# now is found in position 142 and changed to residue V
#
# MAAPVDLELKKAFTELQAKVIDTQQKVKLADIQIEQLNRTKKHAHLTDTEIMTLVDETNMYEGVGRMFILQSKEAIHSQ
# LLEKQKIAEEKIKELEQKKSYLERSVKEAEDNIREMLMARRAQJSYLERGVKEAEDNIRJQKVAEEKIK
# Note the J and the following sequence QKVAEEKIK which
# makes up the new tryptic peptide with two flanking peptides
#
# - Lulu Zheng  - March 2013

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
		_infile => $args{-infile},
		_outfile => $args{-outfile},
},$class;

$class->_incr_count();
	return $self;
}

sub AddConflicts{
	my ($self) = @_;

	my @filenames = ();
	my $decoyfile;
	if($self->{_outfile} =~ /\//){ 
		@filenames = split /\//,$self->{_outfile};  
		$filenames[@filenames-1] = "decoy_".$filenames[@filenames-1];
		$decoyfile = join("/",@filenames);
	}
	else{ $decoyfile = "decoy_".$self->{_outfile}; }

	# Cleaning up from last run
	unlink "../tmp/fastaconflict.tmp" if -e "../tmp/fastaconflict.tmp";  #cleaning up the  tmp files
	unlink "$self->{_outfile}" if -e $self->{_outfile}; #cleaning up the database from last run
	unlink "$decoyfile" if -e $decoyfile; #cleaning up the database from last run
	unlink "../tmp/unused_conflicts.txt" if -e "../tmp/unused_conflicts.txt";  #cleaning up the wrong file from last run

	open FSA, $self->{_infile} or die; # we reformat fasta to make it more easy to handle
	while (<FSA>) {
		X: if (/^>([\w\-]+)(\s.*)$/) {
				$pretag=$1;
				$description=$2;
				while ($_ = <FSA> and !/^>/) {
					chomp;
					$preseq.=$_;
				}
				open(MYNEWFILE, ">>../tmp/fastaconflict.tmp"); #open for write
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


	#Now read the formatted fasta file and do the main routine
	open(MYFASFILE, "../tmp/fastaconflict.tmp"); #open to read the nonfasta lines
	while (<MYFASFILE>) {
		if (/^>([\w\-]+)(\s.*)\t(.*)$/) {
			$tag=$1;
			$desc=$2;
			$seq=$3;
		}
		@snppeps = ();
		@snpdescr = ();
		#Jump to Subroutine to identify and create the conflict peptides
		&create_conflictpeptide_and_output;

		# We only print when there no longer are more hits to uniprot in Conflict annotation file

		#First join the results to avoid the spacecharacters
		$outpep = join "", @snppeps;
		$outdesc = join "", @snpdescr;

		#We add the tag "SequenceChanged" if conflict peptides are added - and the tag is not already added
		# We simply check whether there is a J in the appended sequence but not in the original sequence
		$seqchange="";
		if ( ($seq !~ /J/ ) && ($outpep =~ /J/ ) ) {
			$seqchange = " # ";
		}
		# A tag to tag end of description in all modified sequences
		$seqchange2="";
		if ( ($seq =~ /J/ ) || ($outpep =~ /J/ ) ) {
			$seqchange2 = " #";
		}


		#Then concatenate the outpep with the seq
		$completeseq = $seq . $outpep;
		#Then split it into 60 residue chunks like the uniprot format
		@outseq = $completeseq =~ /.{1,60}/g;

		#now print the new uniprot fasta file entry
		# Only contains new description and peptides if conflict peptides exists - else original entry preserved

		open(MYDBFILE, ">>$self->{_outfile}"); #open for write
		print MYDBFILE ">$tag$desc$seqchange$outdesc$seqchange2\n"; #output entries
		#also print sequence in chunks
		print MYDBFILE "$_\n" foreach (@outseq); 
		close(MYDBFILE); #*** Close the file ***

		# Now we print an uniprot with decoy also the decoy database entries where all sequences are reversed
		#  and REV is addded to the uniprot identifier
		open(MYDECOYFILE, ">>$decoyfile"); #open for write
		print MYDECOYFILE ">$tag$desc$seqchange$outdesc$seqchange2\n"; #printing the ordinary entry
		#also print sequence in chunks
		print  MYDECOYFILE "$_\n" foreach (@outseq);

		#Create reverse sequence
		$reverseseq = reverse($completeseq);
		#Then split it into 60 residue chunks like the uniprot format
		@outrevseq = $reverseseq =~ /.{1,60}/g;

		#now print the new decoy uniprot fasta file entry
		print MYDECOYFILE ">Rev:$tag$desc$seqchange$outdesc$seqchange2\n"; #printing the decoy entry
		#also print sequence in chunks
		foreach (@outrevseq) {
			print MYDECOYFILE "$_\n";
		}
		close(MYDECOYFILE); #*** Close the file ***



		#Now make sure to reset before next entry
		$seqchange="";
		$tag="";
		$desc="";
		$seq="";
		$to="";
		$uniprottag="";
		$position1="";
		$position2="";
		$residueold="";
		$residuenew="";
		$frontleader="";
		$hit="";
		$frontseq="";
		$postseq="";
		$pep4final="";
		$pep2final="";
		$postseqfull="";
		$frontseqfull="";
		$to="";
		$from="";
		$pep1="";
		$pep2="";
		$pep3="";
		$pep4="";
		@snppeps=();
		@snpdescr=();
		@outrevseq=();
		@outseq=();
		$residue1="";
		$residue2="";
		$completeseq="";
		$seqchange="";
		$seqchange2="";

	}
	close(MYFASFILE); #*** Close the file ***

	# This subroutine read the conflict annotation file and find relevant conflicts for a given uniprot entry
	# it then checks whether the report conflict is in accordance with the existing residue
	# if so it will create a new tryptic peptide around the conflict and return to output it
	# if these is no relevant data nothing is return and the original uniprot entry is output
	sub   create_conflictpeptide_and_output {

		#Now open conflict annotation file and read uniprot_ID position of conflict and what change: from AA to AA
		open(MYCONFLICTRESFILE, "../tmp/swiss_conflicts.txt") or return;
		while (<MYCONFLICTRESFILE>) { 
		$flag=0;
			if(/^([\w\-]+)\s+\w+\s+(\d+)\s+(\d+)\s+Missing\s+.*$/){
				$flag = 1;
				$uniprottag=$1;
				$position1=$2;
				$position2=$3;
				$residue1=substr($seq,$position1-1,$position2-$position1+1);
				$residue2="";
			}
			if (/^([\w\-]+)\s+\w+\s+(\d+)\s+(\d+)\s+(\w+)\s+\S+\s+(\w+)\s+.*$/) {   #pick the relevant results
				$flag = 1;
				$uniprottag=$1;
				$position1=$2;
				$position2=$3;
				$residue1=$4;
				$residue2=$5;
			}
			if($flag == 1){
				#Only use conflict data relevant for the given uniprot entry
				if ($uniprottag eq $tag){
				# now we need to pull out the tryptic peptide at this position from $seq
				# We go for the first tryptic peptide with one miscleavage before and after the conflict position
				# Here we check whether the residue at the known position is correct
				$frontleader=$position1-1;
				$change_len = $position2 - $position1 + 1;
				if ($seq =~ /^(.{$frontleader})(.{$change_len})([^J]*)/){
					$frontseq=$1;
					$hit=$2;
					$postseq=$3;

					# We test both whether residue1 or 2 is the correct match
					if ($hit eq $residue1){
						$residueold=$residue1;
						$residuenew=$residue2;
					}
					if($hit eq $residue2){
						$residueold=$residue2;
						$residuenew=$residue1;
					}

					#Now we do the actualcheck
					if ($hit eq $residueold){
					#Then we can substitute for new residue and create a new peptide
					#We do it in two steps producing the header and trailer
					#Here we take the two last tryptic peptides in the front seq
					# We add the residuenew to be able to check whether it is a Proline
						$frontseqfull= $frontseq . $residuenew;
						if ($frontseqfull =~ /.*[K|R]([^P].*?[K|R])([^P].*?)$/){
							$front_ok="OK";
							$pep1=$1;
							$pep2pre=$2;
							#Now we chop of the last character of pep2 which was the residuenew (also start of pep3)
							$pep2= substr($pep2pre, 0, -$change_len);    #chop of the last characters
						}
						else{
							$front_ok="OK";
							$pep1="";
							$pep2pre=$frontseq;
						}

						# We need to deal with cases where pep1 is only two K/R residues
						#Not really nessesary but will look better
						if ($pep1 =~/^[K|R]([K|R])$/){
							#Then pep1 is the single residue K/R peptide
							$pep1=$1;
						}

						#Here we take the two first tryptic peptides in the leader seq
						# We add the snp position to be able to deal with cases where the new position is R or K (not followed by P)
						# This way we can merge peptide 2 and 3 and be sure we got a complete trytic peptide
						$postseqfull= $residuenew . $postseq;
						if ($postseqfull =~ /^(.*?[K|R])([^P].*?[K|R])[^P].*/){
							$post_ok="OK";
							$pep3=$1;
							$pep4=$2;
						}
						else{
							$post_ok="OK";
							$pep3="";
							$pep4=$postseqfull;
						}
						# We need to deal with those peptides that had a K or R in the position above where it was NOT P
						# I.e. Pep2 could start with a K/R and should not be longer if not followed by P
						# which is because we are forced to demand one AA before the K/R (Namely the NOT P demand)
						# here then - pep2 should be split into two peptides in this case - there become pep1 and pep2

						if ($pep2 =~/^([K|R])[^P]/){
							#Then pep1 is the single residue K/R peptide
							$pep1=$1;
							#Then pep2 is the rest of the peptide pep2
							$pep2final=substr($pep2, 1);    #chop all but first character
							#Now change pep2 to its final string
							$pep2=$pep2final;
						}

						# Also we need to deal with cases where pep2 after removal og new residue is
						# only one residue long namely K or R
						if ($pep2 =~/^([K|R])$/){
							#Then pep1 is the single residue K/R peptide
							$pep1=$1;
							#Then pep2 is nothing - since there is cleaved after this residue
							$pep2="";
						}

						# I.e. Pep4 may start with a K/R and should not be longer if not followed by P
						# In this case pep4 should only be the first K/R residue
						if ($pep4 =~/^([K|R])[^P]/){
							#Then pep4 is the single residue K/R peptide
							$pep4=$1;
						}

						# We add the original position of the new N-term peptide that is added
						# The FROM position is calculated
						$from = $position1 - (length($pep1)+length ($pep2));

						# The TO position is the length of the peptide plus start position
						$to = $position1 + (length($pep3)+length ($pep4)-1);

						#Now we collect the results into arrays, so that they can be printed with uniprot entry once we have all snp results
						# We only accept the results if we got all 4 peptides - hence we wont include cases where we are missing parts of peptide
						# which will occur in cases where the SNP residue is very close to the N or C terminus.

						if ($front_ok || $post_ok){
							# We use the theoretic Amino acid J and require trypsin to cleave before and after J
							push(@snppeps, "J");     # We add the J here to get the first snp position correct

							# A little sidetrack finding the new position of the added snp
							# The new position of the substituted residue in the snp
							$lengseq = length($seq);
							$pepsofar = join "", @snppeps;
							$peplengthsofar = length($pepsofar);
							$PosInNewPeptide = (length($pep1) +length($pep2));  #Length of new peptide in front of con

							#First position after seq and so far added peptides
							# + length of to be added peptide (but only part in front of conflict
							$newposition =  $lengseq + $peplengthsofar + $PosInNewPeptide +1;

							#We push the snp description into an array
							push(@snpdescr, "CON[");
							push(@snpdescr, $newposition);
							push(@snpdescr, ",");
							push(@snpdescr, $residuenew);
							push(@snpdescr, ",");
							push(@snpdescr, $position1);  # This is the position in the original sequence
							push(@snpdescr, ",");
							push(@snpdescr, $residueold);
							push(@snpdescr, "]");

							#We push the snp peptides into an array
							push(@snppeps,$pep1);
							push(@snppeps,$pep2);
							push(@snppeps,$pep3);
							push(@snppeps,$pep4);


							# while working with one conflict at a time we open the file for contining registering novel peptides
							open(MYPEPREGFILE, ">>../tmp/PeptideRegister.txt");
							print MYPEPREGFILE "$tag $pep2$pep3 $pep1$pep2$pep3$pep4 $from $to Conflict $newposition $residuenew $position1 $residueold\n";
							close(MYPEPREGFILE); # close peptide register file

						}

					#Now empty variables
					$pep4final="";
					$pep2final="";
					$postseqfull="";
					$pep2final="";
					$frontseqfull="";
					$pep2pre="";
					$pep4final="";
					$pepnew1="";
					$to="";
					$from="";
					$pep1="";
					$pep2="";
					$pep3="";
					$pep4="";
					$post_ok="";
					$front_ok="";
					$lengseq = "";
					$peplengthsofar = "";
					$newposition = "";
					}
					else{   #If the residue position is wrong we send it to a file
						open(MYBADFILE, ">>../tmp/unused_conflicts.txt"); #open for write, overwrite
						print MYBADFILE ">$uniprottag Position: $position1 The found residue $hit\t The expected $residueold\n"; #printing the sequence in one line
						close(MYBADFILE); # *** Close the file ***
					}

				}

				}
			}
		}
		close(MYCONFLICTRESFILE); #*** Close the file ***
	}
	unlink "../tmp/fastaconflict.tmp" if -e "../tmp/fastaconflict.tmp";
}

sub DESTROY {
	my($self) = @_;
	$self->_decr_count();
}

1;
  
  
__END__
  