
package AddSAPs;

# This package reads in fasta file and uses a SNPs data file from Emsembl BioMart
# This package is part of the framework creating a MS friendly version of uniprot fasta sequence database
#
# This package creates tryptic peptides around the given SAP sites and append the peptide at the C-terminal of the uniprot entry
# Peptides are appended with J as a divider
# The J is introduced as a theoretical AA - and the idea is to set trypsin up to cleave also before and after J
# For proteins without a N-terminal signal peptide, only the length is added to the header
#
# Example of a header after introduction of two SNP peptides, appended at the C-terminus (Note: newline introduced here)
#>Q14773 Intercellular adhesion molecule 4 OS=Homo sapiens GN=ICAM4 PE=2 SV=1 lng=271 # SNP[298,L,208,V]SNP[331,L,48,F]
# The # indicate that the sequence has been modified
# A second # shows that the description has not been truncated - but is not added until later

# In this case the original protein is 271 residues long and have two SNP peptide added that
# are appended so that the alternative residue is found in position 298 and 331, respectively
# They have the original positions 208 (residue V) and 48 (residue F)
# To ensure that peptides with one miscleavage are supported - an additional peptide is included both up- and downstream
# of the peptide containing the snp.
#
#Example of uniprot entry containing:
#  -two SNP peptides appended
# Note the J's used to separate extra peptides
#MGSLFPLSLLFFLAAAYPGVGSALGRRTKRAQSPKGSPLAPSGTSVPFWVRMSPEFVAVQPGKSVQLNCSNSCPQPQNS
#SLRTPLRQGKTLRGPGWVSYQLLDVRAWSSLAHCLVTCAGKTRWATSRITAYKPPHSVILEPPVLKGRKYTLRCHVTQV
#FPVGYLVVTLRHGSRVIYSESLERFTGLDLANVTLTYEFAAGPRDFWQPVICHARLNLDGLVVRNSSAPITLMLAWSPA
#PTALASGSIAALVGILLTVGAAYLCKCLAMKSQAJFTGLDLANVTLTYEFAAGPRDFWQPLICHARLNLDGLVVRJAQS
#PKGSPLAPSGTSVPLWVRMSPEFVAVQPGK
#
##################################################
### - Lulu Zheng - March 2013 #####
##################################################

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;

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
		_sap => $args{-sap},
		_stop_gain => $args{-stop_gain},
		#_problematic => $args{-problematic} || 0,
	},$class;
	
	$class->_incr_count();
	return $self;
}

sub AddSAPs{
	my ($self) = @_;
	
	unlink "../tmp/bigseqs.fasta" if -e "../tmp/bigseqs.fasta";
	
	my %pep_hash = ();
	my $str = Bio::SeqIO->new(-file=> $self->{_file},-format => "fasta");
	while(my $seq = $str->next_seq()){
		my $id = $seq->display_id();
		if($seq->display_id() =~ /\|/){
			my @ids = split /\|/,$seq->display_id();
			$id = $ids[1];
		}
		$pep_hash{$id} = $seq;
	}

	my $seq_obj = Bio::Seq->new(-id => "00001",-seq => "aaaaa");
	my @filenames = ();
	my $outfile;
	if($self->{_file} =~ /\//){ @filenames = split /\//,$self->{_file};$outfile = "../tmp/sap_".$filenames[@filenames-1];}
	else{ $outfile = "../tmp/sap_".$self->{_file}; }
	my $out = Bio::SeqIO->new(-file => ">$outfile",-format => "fasta");

	my $snp = $self->{_sap};
	my $stop_gain = $self->{_stop_gain};
	

	foreach my $protein(keys %{$snp}){
		next unless exists $pep_hash{$protein};
		
		my @locs = sort {$a <=> $b} keys %{$snp->{$protein}};
		my $string = $pep_hash{$protein}->seq();
		my $pri_string = $string;
		$pri_string =~ s/J.*//;
		my $desc = $pep_hash{$protein}->desc();
		
		my @substrs = ();
		my @positions = ();
		
		while($pri_string =~ /([K|R])/g){
			next if substr($string,$-[1]+1,1) eq "P";
			push @positions, $-[1] + 1;
			push @substrs,$`.$1;
		}
		push @positions,length($pri_string);
		push @substrs,$pri_string;
		
		for(my $i = 1; $i < @positions; $i++){
			$substrs[$i] = substr($substrs[$i],$positions[$i-1], $positions[$i] - $positions[$i-1]);
		}
		
		my $trypsin_pep = "";
		my $trypsin_desc = "";
		
		for(my $i = 0; $i < @positions; $i++){
			my @hits = ();
			for(my $j = 0; $j < @locs; $j++){
				my @mutations = split /\//,$snp->{$protein}->{$locs[$j]};
				if($i == 0){	#### for the first trypsin peptide in the primary sequence
					if($locs[$j] >= 1 && $locs[$j] + length($mutations[0]) - 1 <= $positions[$i]){
						push @hits,$locs[$j];
					}
				}
				else{
					if( $locs[$j] >= $positions[$i-1] + 1 && $locs[$j] + length($mutations[0]) - 1 <= $positions[$i]){
						push @hits,$locs[$j] - $positions[$i-1];
					}
				}
			}

			next unless @hits;

			my @peps = ();
			my @descs = ();
			my $first = substr($substrs[$i],0,$hits[0] - 1);
			push @peps,$first;
			push @descs,"";
			my $m;
			my @mutations;
			my $num = 0;
			for($m = 0;$m < @hits; $m++){
				@mutations = ();
				if($i == 0){ @mutations = split /\//,$snp->{$protein}->{$hits[$m]}; }
				else{ 
					my $pos = $hits[$m] + $positions[$i-1];
					@mutations = split /\//,$snp->{$protein}->{$pos}; 
				}
				$num += @mutations - 1;		#### to count the number of SAPs in each tryptic peptide
			}
			
			#-------------------------------------------------------------------------------
			# for the problematic peptides, just use SAPs to subsititute the original AA to generate SAP heterogeneity sequences
			if($num >= 21){
				@peps = ();
				@descs = ();
				for($m = 0;$m < @hits; $m++){
					@mutations = ();
					if($i == 0){ @mutations = split /\//,$snp->{$protein}->{$hits[$m]}; }
					else{ 
						my $pos = $hits[$m] + $positions[$i-1];
						@mutations = split /\//,$snp->{$protein}->{$pos}; 
					}
					for(my $k = 1; $k < @mutations; $k++){
						push @peps,substr($substrs[$i],0,$hits[$m] - 1).$mutations[$k].substr($substrs[$i],$hits[$m]);
						push @descs,$mutations[0].":".$hits[$m].":".$mutations[$k] 
					}
				}
			}
			#-------------------------------------------------------------------------------
			# for the normal peptides, permutate all SAPs to generate SAP heterogeneity sequences 
			else{
				for($m = 0;$m < @hits; $m++){
					@mutations = ();
					if($i == 0){ @mutations = split /\//,$snp->{$protein}->{$hits[$m]}; }
					else{ 
						my $pos = $hits[$m] + $positions[$i-1];
						@mutations = split /\//,$snp->{$protein}->{$pos}; 
					}
					if($m != 0){ 
						for(my $n = 0; $n < @peps; $n++){
							my $substr_len = length($mutations[0]);
							#print $mutations[0]."\t".$substrs[$i]."\tlen:".length($substrs[$i])."\tpos1:".$hits[$m-1]."\tpos2:".$hits[$m]."\n";
							$peps[$n] .= substr($substrs[$i],$hits[$m-1]+$substr_len-1,$hits[$m] - $hits[$m-1] - $substr_len);
						}
					}
					my @tmp = @peps;
					my @tmp2 = @descs;
					@peps = ();
					@descs = ();
				
					for(my $n = 0; $n < @tmp; $n++){
						for(my $k = 0; $k < @mutations; $k++){
							push @peps,$tmp[$n].$mutations[$k];
							if($k == 0){ push @descs,$tmp2[$n]; }
							else{
								push @descs,$tmp2[$n].$mutations[0].":".$hits[$m].":".$mutations[$k] if $tmp2[$n] eq "";
								push @descs,$tmp2[$n]."|".$mutations[0].":".$hits[$m].":".$mutations[$k] if $tmp2[$n] ne "";
							}
						}
					}
					@tmp = ();
					@tmp2 = ();
				}
				##### to join the last part
				for(my $j = 0; $j < @peps; $j++){
					my $substr_len = length($mutations[0]);##(length($snp{$hits[$m-1]})-1)/2;
					next if $hits[$m-1]+$substr_len-1 >= length($substrs[$i]);
					$peps[$j] .= substr($substrs[$i],$hits[$m-1]+$substr_len-1);
				}
			}
			
			my $desc_str = "";
			for(my $j = 0 ; $j < @peps; $j++){
				next if $peps[$j] eq $substrs[$i];	# to remove the identical peptide as the primary.
				#print $peps[$j],"\n";
				my @tmp_desc = ();
				if($descs[$j] =~ /\|/){
					my @tmp = split /\|/,$descs[$j];
					push @tmp_desc,@tmp;	
				}
				else{ push @tmp_desc,$descs[$j]; }
				
				if(@positions == 1){ 
					#push @trypsin_peps, $peps[$j];
					$trypsin_pep .= "J".$peps[$j];
					foreach my $desc(@tmp_desc){
						my @tmp = split /\:/,$desc;
						my $pri_pos = $tmp[1];
						my $new_pos = length($string) + length($trypsin_pep) + $tmp[1] + 1;		##### the last "1" means the lenght of "J"
						$desc_str .= "SNP[".$new_pos.",".$tmp[2].",".$pri_pos.",".$tmp[0]."]";
					}
					#push @trypsin_descs,$desc_str;
				}
				elsif(@positions == 2){
					if($i == 0){ 
						foreach my $desc(@tmp_desc){
							my @tmp = split /\:/,$desc;
							my $pri_pos = $tmp[1];
							my $new_pos = length($string) + length($trypsin_pep) + $tmp[1] + 1;
							$desc_str .= "SNP[".$new_pos.",".$tmp[2].",".$pri_pos.",".$tmp[0]."]";
						}
						#push @trypsin_descs,$desc_str;
						#$trypsin_desc .= $desc_str if $j == @peps - 1;
						
						if($peps[$j] =~ /.*?[K|R][^P].*?[R|K]$/){ $trypsin_pep .= "J".$peps[$j]; }
						else{ $trypsin_pep .= "J".$peps[$j].$substrs[$i+1];}
					}
					else{ 
						if($peps[$j] =~ /^[K|R][^P]/){ 
							#push @trypsin_peps,$peps[$j];
							foreach my $desc(@tmp_desc){
								my @tmp = split /\:/,$desc;
								my $pri_pos = $tmp[1] + $positions[$i-1];
								my $new_pos = length($string) + length($trypsin_pep) + $tmp[1] + 1;
								$desc_str .= "SNP[".$new_pos.",".$tmp[2].",".$pri_pos.",".$tmp[0]."]";
							}
							#push @trypsin_descs,$desc_str;
							#$trypsin_desc .= $desc_str if $j == @peps - 1;
							$trypsin_pep .= "J".$peps[$j];
						}
						else{ 
							#push @trypsin_peps,$substrs[$i-1].$peps[$j];
							foreach my $desc(@tmp_desc){
								my @tmp = split /\:/,$desc;
								my $pri_pos = $tmp[1] + $positions[$i-1];
								my $new_pos = length($string) + length($trypsin_pep) + length($substrs[$i-1]) + $tmp[1] + 1;
								$desc_str .= "SNP[".$new_pos.",".$tmp[2].",".$pri_pos.",".$tmp[0]."]";
							}
							#push @trypsin_descs,$desc_str;
							#$trypsin_desc .= $desc_str if $j == @peps - 1;
							$trypsin_pep .= "J".$substrs[$i-1].$peps[$j];
						}
					}
				}
				else{
					if($i == 0){
						foreach my $desc(@tmp_desc){
								my @tmp = split /\:/,$desc;
								my $pri_pos = $tmp[1];
								my $new_pos = length($string) + length($trypsin_pep) + $tmp[1] + 1;
								$desc_str .= "SNP[".$new_pos.",".$tmp[2].",".$pri_pos.",".$tmp[0]."]";
						}
						#$trypsin_desc .= $desc_str if $j == @peps - 1;
						if($peps[$j] =~ /.*?[K|R][^P].*?[R|K]$/) { $trypsin_pep .= "J".$peps[$j]; }
						elsif($peps[$j] =~ /[K|R]$/){ $trypsin_pep .= "J".$peps[$j].$substrs[$i+1];}
						else{ $trypsin_pep .= "J".$peps[$j].$substrs[$i+1].$substrs[$i+2];}
					}
					
					elsif($i == @positions - 1){
						if($peps[$j] =~ /^P/){ 
							foreach my $desc(@tmp_desc){
								my @tmp = split /\:/,$desc;
								my $pri_pos = $tmp[1] + $positions[$i-1];
								my $new_pos = length($string) + length($trypsin_pep) + length($substrs[$i-2].$substrs[$i-1]) + $tmp[1] + 1;
								$desc_str .= "SNP[".$new_pos.",".$tmp[2].",".$pri_pos.",".$tmp[0]."]";
							}
							#$trypsin_desc .= $desc_str if $j == @peps - 1;
							$trypsin_pep .= "J".$substrs[$i-2].$substrs[$i-1].$peps[$j];
							#push @trypsin_peps,$substrs[$i-2].$substrs[$i-1].$peps[$j]; 
						}
						elsif($peps[$j] =~ /^[K|R][^P]/){ 
							#push @trypsin_peps,$peps[$j]; 
							foreach my $desc(@tmp_desc){
								my @tmp = split /\:/,$desc;
								my $pri_pos = $tmp[1] + $positions[$i-1];
								my $new_pos = length($string) + length($trypsin_pep) + $tmp[1] + 1;
								$desc_str .= "SNP[".$new_pos.",".$tmp[2].",".$pri_pos.",".$tmp[0]."]";
							}
							#$trypsin_desc .= $desc_str if $j == @peps - 1;
							$trypsin_pep .= "J".$peps[$j];
						}
						else{ 
							#push @trypsin_peps, $substrs[$i-1].$peps[$j];
							foreach my $desc(@tmp_desc){
								my @tmp = split /\:/,$desc;
								my $pri_pos = $tmp[1] + $positions[$i-1];
								my $new_pos = length($string) + length($trypsin_pep) + length($substrs[$i-1]) + $tmp[1] + 1;
								$desc_str .= "SNP[".$new_pos.",".$tmp[2].",".$pri_pos.",".$tmp[0]."]";
							}
							#$trypsin_desc .= $desc_str if $j == @peps - 1;
							$trypsin_pep .= "J".$substrs[$i-1].$peps[$j];
						}
					}
					else{
						my $new_pep;
						if($peps[$j] =~ /^[K|R][^P]/){ 
							foreach my $desc(@tmp_desc){
								my @tmp = split /\:/,$desc;
								my $pri_pos = $tmp[1] + $positions[$i-1];
								my $new_pos = length($string) + length($trypsin_pep) + $tmp[1] + 1;
								$desc_str .= "SNP[".$new_pos.",".$tmp[2].",".$pri_pos.",".$tmp[0]."]";
							}
							#$trypsin_desc .= $desc_str if $j == @peps - 1;
							
							$new_pep = $peps[$j].$substrs[$i+1];
							if($peps[$j] =~ /^[K|R][^P].*?[^K,R]$/){
								if($i+2 <= @positions-1){ $trypsin_pep .= "J".$new_pep.$substrs[$i+2];}
								else{ $trypsin_pep .= "J".$new_pep;}
							}
							else{ $trypsin_pep .= "J".$new_pep;}
						}
						elsif($peps[$j] =~ /.*?[K|R][^P].*?[R|K]$/){
							$new_pep = $substrs[$i-1].$peps[$j];
							if($peps[$j] =~ /^P.*?[K|R][^P].*?[K|R]$/){
								if($i-2 >= 0){
									#push @trypsin_peps,$substrs[$i-2].$new_pep;
									foreach my $desc(@tmp_desc){
										my @tmp = split /\:/,$desc;
										my $pri_pos = $tmp[1] + $positions[$i-1];
										my $new_pos = length($string) + length($trypsin_pep) + length($substrs[$i-2].$substrs[$i-1]) + $tmp[1] + 1;
										$desc_str .= "SNP[".$new_pos.",".$tmp[2].",".$pri_pos.",".$tmp[0]."]";
									}
									#$trypsin_desc .= $desc_str if $j == @peps - 1;
									$trypsin_pep .= "J".$substrs[$i-2].$new_pep;
								}
								else{ 
									#push @trypsin_peps,$new_pep;
									foreach my $desc(@tmp_desc){
										my @tmp = split /\:/,$desc;
										my $pri_pos = $tmp[1] + $positions[$i-1];
										my $new_pos = length($string) + length($trypsin_pep) + length($substrs[$i-1]) + $tmp[1] + 1;
										$desc_str .= "SNP[".$new_pos.",".$tmp[2].",".$pri_pos.",".$tmp[0]."]";
									}
									#$trypsin_desc .= $desc_str if $j == @peps - 1;
									$trypsin_pep .= "J".$new_pep;
								}
							}
						}
						else{ 
							$new_pep = $substrs[$i-1].$peps[$j].$substrs[$i+1]; 
							if($new_pep =~ /^[^P].*?[R|K][^P].*?[R|K]$/){ 
								#push @trypsin_peps,$new_pep;
								foreach my $desc(@tmp_desc){
										my @tmp = split /\:/,$desc;
										my $pri_pos = $tmp[1] + $positions[$i-1];
										my $new_pos = length($string) + length($trypsin_pep) + length($substrs[$i-1]) + $tmp[1] + 1;
										$desc_str .= "SNP[".$new_pos.",".$tmp[2].",".$pri_pos.",".$tmp[0]."]";
								}
								#$trypsin_desc .= $desc_str if $j == @peps - 1;
								$trypsin_pep .= "J".$new_pep;
							}
							elsif($new_pep =~ /^P.*?[R|K][^P].*?[R|K]$/){
								if($i-2 >= 0){ 
									#push @trypsin_peps,$substrs[$i-2].$new_pep;
									foreach my $desc(@tmp_desc){
										my @tmp = split /\:/,$desc;
										my $pri_pos = $tmp[1] + $positions[$i-1];
										my $new_pos = length($string) + length($trypsin_pep) + length($substrs[$i-2].$substrs[$i-1]) + $tmp[1] + 1;
										$desc_str .= "SNP[".$new_pos.",".$tmp[2].",".$pri_pos.",".$tmp[0]."]";
									}
									#$trypsin_desc .= $desc_str if $j == @peps - 1;
									$trypsin_pep .= "J".$substrs[$i-2].$new_pep;
								}
								else{ #push @trypsin_peps,$new_pep;
									foreach my $desc(@tmp_desc){
										my @tmp = split /\:/,$desc;
										my $pri_pos = $tmp[1] + $positions[$i-1];
										my $new_pos = length($string) + length($trypsin_pep) + length($substrs[$i-1]) + $tmp[1] + 1;
										$desc_str .= "SNP[".$new_pos.",".$tmp[2].",".$pri_pos.",".$tmp[0]."]";
									}
									#$trypsin_desc .= $desc_str if $j == @peps - 1;
									$trypsin_pep .= "J".$new_pep;
								}
							}
							elsif($new_pep =~ /^[^P].*?[^K,R][^P].*?[R|K]$/){
								foreach my $desc(@tmp_desc){
									my @tmp = split /\:/,$desc;
									my $pri_pos = $tmp[1] + $positions[$i-1];
									my $new_pos = length($string) + length($trypsin_pep) + length($substrs[$i-1]) + $tmp[1] + 1;
									$desc_str .= "SNP[".$new_pos.",".$tmp[2].",".$pri_pos.",".$tmp[0]."]";
								}
								#$trypsin_desc .= $desc_str if $j == @peps - 1;
								if($i+2 <= @positions-1){ $trypsin_pep .= "J".$new_pep.$substrs[$i+2];}
								else{ $trypsin_pep .= "J".$new_pep;}
							}
						}
					}
				}
			}
			@peps = ();
			@descs = ();
			$trypsin_desc .= $desc_str;
		}
		
		my $id = $protein;
		my $change = "";
		if($string !~ /J/ && $trypsin_pep =~ /J/){ $change = " # "; }
		if(length($trypsin_pep) > 10000000){ 
			open OUT,">>../tmp/bigseqs.fasta"; 
			print OUT ">$id $desc$change$trypsin_desc\n$string$trypsin_pep\n"; 
			close OUT;
		}
		else{
			$seq_obj->display_id($id);
			$seq_obj->desc($desc.$change.$trypsin_desc);
			$seq_obj->seq($string.$trypsin_pep);
			$out->write_seq($seq_obj);
		}

		
		my $count = 1;
		foreach my $loc(sort {$a <=> $b} keys %{$stop_gain->{$protein}}){
			my $string = $pep_hash{$protein}->subseq(1,$loc-1);
			my $id = $protein."_".$count;
			my $desc = "truncated isform ".$count.": ".$pep_hash{$protein}->desc()." # Variation: $loc".$stop_gain->{$protein}->{$loc}." #";
			$seq_obj->display_id($id);
			$seq_obj->desc($desc);
			$seq_obj->seq($string);
			$out->write_seq($seq_obj);
			$count++;
		}
	}

	foreach my $protein(keys %pep_hash){
		unless(exists $snp->{$protein}){ $out->write_seq($pep_hash{$protein}); }
	}	
	
	if (-e "../tmp/bigseqs.fasta"){
		open PRI,">>$outfile";
		open BIG,"../tmp/bigseqs.fasta";
		while(my $line = <BIG>){ print PRI $line;}
		close PRI;
		close BIG;
	}
	unlink "../tmp/bigseqs.fasta" if -e "../tmp/bigseqs.fasta";
}

sub DESTROY {
    my($self) = @_;
    $self->_decr_count();
}

1;
