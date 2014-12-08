#!/usr/bin/perl

# SAPshodt.pl March 2013 $

use strict;
use warnings;

use lib "Bio";
use lib "Add";
use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;

use Add::Check;
use Add::CreateSwissConflicts;
use Add::AddTransit;
use Add::AddSAPs;
use Add::AddConflicts;

#-------------------------------------------------------------------------------

# get the user options
my ( $help, $fasta, $sapfile, $spfile, $inc, $targetp, $outfile );
GetOptions( 'help'         => \$help,
            'fasta=s'      => \$fasta,
			'sapfile=s'		=> \$sapfile,
			'targetp=s'		=> \$targetp,
			'spfile=s'      => \$spfile,
            'inc=s'        => \$inc,
            'outfile=s'    => \$outfile,
            'h'            => \$help,
);

help() if $help;
help() unless ( $sapfile and $fasta ); # required options


#-------------------------------------------------------------------------------

# check the input parameters

unless ( defined $sapfile and defined $fasta ){
	print "FATAL: must specify both -snp and -fasta\n";
	exit;
}

unless (-s $sapfile){
	print "FATAL: can't find file $sapfile\n";
	exit;
}

unless (-s $fasta){
	print "FATAL: can't find file $fasta\n";
	exit;
}

if ( $outfile and -s $outfile ){
	print "FATAL: output file $outfile already exists\n";
	exit;
}

#-------------------------------------------------------------------------------
# read the fasta file

my %pep_hash = ();
my $str = Bio::SeqIO->new(-file => $fasta, -format => 'fasta');
while(my $seq = $str->next_seq()){
	my $id = $seq->display_id();
	if($seq->display_id() =~ /\|/){
		my @ids = split /\|/,$seq->display_id();
		$id = $ids[1];
	}
	$pep_hash{$id} = $seq;
}

my %inclusion = ();
if(!defined($inc) or $inc eq "all"){ %inclusion = %pep_hash;}
else{
	open FILE,$inc;
	while(<FILE>){
		chomp;
		my $id = $_;
		$id =~ s/\r$// if $id =~ /\r$/;
		if(/\|/){
			my @ids = split /\|/;
			$id = $ids[1];
		}
		unless(exists $pep_hash{$id}){ die qq(FATAL: "$_" doesn't exist in the "$fasta". Please make sure the versions of protein accession numbers are the same ); }
		$inclusion{$id} = $pep_hash{$id};
	}
	close FILE;
}

undef %pep_hash;
#-------------------------------------------------------------------------------
# extract the SAPs information from $sapfile. The format should be like this (tab delimited):
# H0YG80	4	L/Q	rs143466490	T/A	ENSP00000442112
# H0YG80	7	M/T	TMP_ESP_15_63889611	T/C	ENSP00000442112
# H0YG80	32	L/F	TMP_ESP_15_63889685	C/T	ENSP00000442112
# H0YG80	42	V/M	TMP_ESP_15_63889715	G/A	ENSP00000442112
# H0YG80	44	E/D	rs142301663	G/C	ENSP00000442112
# H0YG80	55	R/C	rs148140490	C/T	ENSP00000442112
# 
# The first three columns are required, which means protein accession, SAP location and SAP. 
# And the last three columns here are only to indicate the source of SAPs
# Note: no header columns are allowed  

my %sap = ();
my %stop_gain = ();

open FILE,$sapfile;
while(<FILE>){
	chomp;
	my @columns = split /\t/;
	
	if($columns[2]){ 
			if($columns[2] =~ /\*/){ 
				next if $columns[1] == 1;
				$stop_gain{$columns[0]}{$columns[1]} = $columns[2]; 
			}
			else{ 
				if(exists $sap{$columns[0]}{$columns[1]}){
					my @mutations1 = split /\//,$sap{$columns[0]}{$columns[1]};
					my @mutations2 = split /\//,$columns[2];
					#unless($mutations1[0] eq $mutations2[0]){ die $columns[0]."\t".$columns[1]."\t".$mutations1[0]."\t".$mutations2[0]."\n"; }
					next unless($mutations1[0] eq $mutations2[0]);
					my $num = $#mutations1;
					my $flag = 0;
					foreach my $mu(@mutations1[1..$num]){
						if($mutations2[1] eq $mu){ $flag = 1; last;}
					}
					next if $flag == 1;
					$sap{$columns[0]}{$columns[1]} .= "/".$mutations2[1];
				} 
				else{
					$sap{$columns[0]}{$columns[1]} = $columns[2];
				}
			}
	}	
}
close FILE;

#-------------------------------------------------------------------------------
# check whether there are some problematic sequences that containing too many SAPs in the trypsic region

my $check_obj = Check->new(-snp => \%sap,-pep => \%inclusion);
my $problematic = $check_obj->Check();
my @keys = keys %{$problematic};
if(@keys){
	open PRO,">../tmp/problematic_proteins.txt";
	print PRO "Accession\tDescription\n";
	foreach my $id(@keys){
		print PRO $id."\t".$problematic->{$id}."\n";
	}
	close PRO;
	if(@keys > 20){
		print "\nThere are more than 20 proteins containing too many SAPs (>20) in their certain tryptic peptide\n";
		print "Please see the proteins list in the file ../tmp/problematic_proteins.txt\n";
	}else{
		print "\nThe following proteins contain too many SAPs (>20) in their certain tryptic peptide: \n\n";
		foreach my $id(@keys){
			print $id."\t".$problematic->{$id}."\n";
		}
		print "\nYou can also see the proteins list in the file ../tmp/problematic_proteins.txt\n";
	}
	print "\nThe SAPs in the corresponding peptide/s won't be permutated.\n";
	print "They're only used to substitute their original amino acids to generate SAP heterogeneity sequences.\n\n";
}
else{
	print "\nAll the input sequences are good. We will permutate all SAPs in each tryptic peptide\n\n";
}

#print "Then, please select whether to generate the SAP heterogeneity sequences for those problematic proteins.[Y/N]: ";
#my $choice = <STDIN>;
#chomp($choice);

#if ($choice eq "Y" or $choice eq "y" or $choice =~ /yes/i){}

#-------------------------------------------------------------------------------
# creat conflict sequences from swissprot annotations if possible

if(defined $spfile){
	my $sp_anno = CreateSwissConflicts->new(-file => $spfile);
	$sp_anno->CreateSwissConflicts();		#### This will produce a file named "swiss_conflicts.txt" in the dir "../tmp"
}else{
	unlink "../tmp/swiss_conflicts.txt" if -e "../tmp/swiss_conflicts.txt";
}

#-------------------------------------------------------------------------------
####  CREATE & APPEND PEPTIDES  ####

# Run TargetP prediction for all input sequences
# This will append new peptides from mitochrondrial transit peptides &  from signal peptides
open TMP,">../tmp/tmpfasta.fasta";
foreach(keys %inclusion){
	print TMP ">$_"." ".$inclusion{$_}->desc()."\n";
	print TMP $inclusion{$_}->seq()."\n";
}
close TMP;

unless(defined $targetp){ $targetp = "N"; }
else{
	if($targetp =~ /^Y$/i or $targetp =~ /^yes$/i ){ $targetp = 'Y'; }
	else{ $targetp = "N"; }
}

my $addtransit = AddTransit->new(-file => "../tmp/tmpfasta.fasta", -targetp => $targetp);
$addtransit->AddTransit();					#### This will produce a file named "tr_sp_tmpfasta.fasta" in the dir "../tmp"
unlink "../tmp/tmpfasta.fasta";

# Use the package AddSAPs.pm to convert  ../tmp/tr_sp_tmpfasta.fasta to include SAP heterogeneity peptides

my $sapfasta = AddSAPs->new(
					-file => "../tmp/tr_sp_tmpfasta.fasta",
					-sap => \%sap,
					-stop_gain => \%stop_gain,
					#-problematic => $choice,
				);
$sapfasta->AddSAPs();						#### This will produce a file named "sap_tr_sp_tmpfasta.fasta" in the dir "../tmp"
unlink "../tmp/tr_sp_tmpfasta.fasta";

# Run AddConflicts.pm to include the conflict peptides
# This package also makes an uniprot including decoy fasta entries, where the sequence is reversed and the uniprot id changed to REV
# The first is called $outfile
# The latter is called decoy_$outfile

my $final = AddConflicts->new( -infile => "../tmp/sap_tr_sp_tmpfasta.fasta", -outfile => $outfile);
$final->AddConflicts();
unlink "../tmp/sap_tr_sp_tmpfasta.fasta" if -e "../tmp/sap_tr_sp_tmpfasta.fasta";

# Then we get the the SAP heterogeneity sequences on demand transmutator

#-------------------------------------------------------------------------------

sub help {
  print STDERR <<EOF;

SAPshodt.pl: generate the SAP heterogeneity sequences on demand transmutator

Usage: perl SAPshodt.pl -fasta <fasta_file> -sapfile <SAP file> -outfile <output filename>

Additional options:

	-h             : show this help
	-inc <file>    : specify whether only to handle some proteins in the fasta_file (default: all),
                         otherwise please use this option to input a file containing protein accessions.
                         These two kinds of accessions are both allowed: "sp|P31946|1433B_HUMAN" or "P31946"
	-targetp <Y/N> : specify whether to use TargetP (SignalP) to predict transit or signal peptides,
                         and then to append these peptides to the original sequence (default: N) 
	-spfile	<FILE> : specify whether to get conflicting sequences from Swiss-Prot conflict annotations,
                         and then to append these peptides to the original sequence (default: N)

EOF
  exit;

}

#-------------------------------------------------------------------------------
