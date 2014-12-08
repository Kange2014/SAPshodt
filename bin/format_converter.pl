#!/usr/bin/perl -w

#### perl format_converter.pl
#### to convert the ensembl IDs to uniprot IDs,
#### and filter the missense variations that don't result in SAPs, that means, amino acid indels will not be analyzied.

#### the input SNP annotation file format should be like this (tab separated):
#### Ensembl Protein ID	Protein location (aa)	Protein Allele	Reference ID	Variant Alleles
#### ENSP00000442112	4	L/Q	rs143466490	T/A
#### ENSP00000442112	7	M/T	TMP_ESP_15_63889611	T/C
#### ENSP00000442112	32	L/F	TMP_ESP_15_63889685	C/T
#### ENSP00000442112	42	V/M	TMP_ESP_15_63889715	G/A

##################################################
### - Lulu Zheng - March 2013 #####
##################################################

use strict;

use Bio::SeqIO;
use Getopt::Long;

my ( $help, $ensembl, $uniprot, $idmapping, $sapfile, $outfile );
GetOptions( 'help'         => \$help,
            'ensembl=s'      => \$ensembl,
			'sapfile=s'		=> \$sapfile,
			'idmapping=s'      => \$idmapping,
            'uniprot=s'        => \$uniprot,
            'outfile=s'    => \$outfile,
            'h'            => \$help,
);

help() if $help;
help() unless ( $ensembl or $uniprot or $idmapping or $sapfile or $outfile ); # required options

if ( $outfile and -s $outfile ){
	print "FATAL: output file $outfile already exists, please delete the original one or input a new filename\n";
	exit;
}

my %ensembl = ();
my $str = Bio::SeqIO->new(-file=> $ensembl,-format => "fasta");
while(my $seq = $str->next_seq()){
	$ensembl{$seq->display_id()} = $seq;
}

my %uniprot = ();
$str = Bio::SeqIO->new(-file=> $uniprot, format => "fasta");
while(my $seq = $str->next_seq()){
	my @ids = split /\|/,$seq->display_id();
	$uniprot{$ids[1]} = $seq;
}

my %converter = ();
open FILE, "$idmapping";
while(<FILE>){
	#next if $. == 1;
	chomp;
	my @columns = split /\t/;
	my @ensembls = ();
	next unless $columns[21];
	if($columns[21] =~ /;\s/){
		my @tmp = split /;\s/,$columns[21];
		push @ensembls,@tmp;
	}
	else{ push @ensembls,$columns[21]; }
	foreach my $id(@ensembls){
		#print "Ensembl\t$id\n" unless exists $ensembl{$id};
		next unless exists $ensembl{$id};
		foreach my $isform(2..100){
			if(exists $uniprot{$columns[0]."-".$isform}){
				if($ensembl{$id}->seq() eq $uniprot{$columns[0]."-".$isform}->seq() ){
					if(exists $converter{$id}){ $converter{$id} .= "|".$columns[0]."-".$isform; }
					else{ $converter{$id} = $columns[0]."-".$isform; }
				}
			}
			else{ last;}
		}
		if($ensembl{$id}->seq() eq $uniprot{$columns[0]}->seq()){
			if(exists $converter{$id}){ $converter{$id} .= "|".$columns[0]; }
			else{ $converter{$id} = $columns[0]; }
		}
	}
}
close FILE;

#my @keys = keys %converter;
#my $len = @keys;
#print "The number of mappings is: ".$len."\n";

open SAP,"$sapfile";
open OUT,">$outfile";
while(<SAP>){
	next if $. == 1;
	chomp;
	my @columns = split /\t/;
	next unless $columns[3];											#### every SAP should have a SNP reference ID to show its source
	next if $columns[4] =~ /\-/;										#### no indel will be analyzied
	my @snps = split /\//,$columns[4];
	my $flag = 0;
	for(my $i = 1;$i < @snps; $i++){
		if(length($snps[$i]) != length($snps[0])){ $flag = 1; last;}	#### Also, check whether there are indels and no indel is allowed 
	}
	next if $flag == 1;
	
	if(exists $converter{$columns[0]}){
			my @ids = ();
			if($converter{$columns[0]} =~ /\|/){
				my @tmp = split /\|/,$converter{$columns[0]};
				push @ids, @tmp;
			}
			else{ push @ids, $converter{$columns[0]}; }
			foreach my $id(@ids){
				#### to make large variations into some SAP/s if possible, e.g.,
				#### ENSP00000371119	15	RRYVQ/RRYVE	rs71263997	GTACATAGCGC/CTACGTAGCGT
				#### will be like this:
				#### ENSP00000371119	19	Q/E	rs71263997	GTACATAGCGC/CTACGTAGCGT
				if(length($columns[2]) > 3){
						my @mutations = split /\//,$columns[2];
						next if length($mutations[0]) != length($mutations[1]);	#### no indel will be analyzied
						for(my $k = 1; $k <= length($mutations[0]); $k++){
							if(substr($mutations[0],0,$k) eq substr($mutations[1],0,$k)){ next; }
							else{
								my $pos = $columns[1] + $k - 1;
								my $aa1 = substr($mutations[0],$k - 1,1);
								my $aa2 = substr($mutations[1],$k - 1,1);
								print OUT join("\t",($id,$pos,"$aa1/$aa2",@columns[3..4],$columns[0]))."\n";
							}
						}
				}else{
					print OUT join("\t",($id,@columns[1..4],$columns[0]))."\n";
				}
			}
	}
}
close SAP;
close OUT;

sub help {
  print STDERR <<EOF;

format_converter.pl: to convert ensembl-ID SAP dataset into uniprot-ID SAP dataset

Usage: perl format_converter.pl -ensembl <fasta_file from ensembl> -uniprot <fasta_file from uniprot> -idmapping <id mapping file from uniprot> -sapfile <SAP file from ensembl> -outfile <output filename>

EOF
  exit;

}
