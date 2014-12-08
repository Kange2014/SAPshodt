SAPshodt
========

Mass spectrometry (MS)-based proteomics analysis usually relies on database search. However, standard protein data-bases, such as NCBI-nr and Uniprot databases, don’t contain protein variants arising from SNPs that result in amino acid changes. Therefore, SAPshodt has been developed to introduce single amino acid polymorphism (SAP) information into an existing protein reference database on demand transmutator, enabling variant peptide detection in proteomics. 

1.	Preliminaries
To help understand this guide, the following keywords in the context mean:

SNP: A single-nucleotide polymorphism (SNP) is a DNA sequence variation occurring when a single nucleotide — A, T, C or G — in the genome (or other shared sequence) differs between members of a biological speciesor paired chromosomes.

SAP: A single amino-acid polymorphism (SAP) is an amino acid substitution caused by a non-synonymous SNP or mutation.

Trypsin: A serine protease found in the digestive system of many vertebrates, where it hydrolyses proteins. Trypsin is produced in the pancreas as the inactive proenzyme trypsinogen. Trypsin cleaves peptide chains mainly at the carboxyl side of the amino acids lysine (K) or arginine (R), except when either is followed by proline (P).

Tryptic peptide: Trypsin catalyses the hydrolysis of peptide bonds so that proteins can be broken down into smaller peptides. These peptides are tryptic peptides. In SAPshodt, for each SAP and their permutation occuring in the same tryptic peptide, we appended this enclosing tryptic peptide with substituted residues as well as the two flanking tryptic peptides to the original sequence. 

2.	Installation
SAPshodt will run on the most common UNIX (e.g., Linux etc.) and windows platforms. Make sure that Perl is available on your system. Decide where you wish to keep the software. Uncompress and untar the package in that location when under UNIX platforms as follows (users can use  compression tools  like winrar to uncompress it under windows platform):
>tar xvzf sapshodt.tar.gz

This will produce a directory “SAPshodt”, and it contains three folders: bin/, DATABASES/ and tmp/. You will run some scirpts from “bin/” and put your input or related datasets in “DATABASES/”. The temporary files produced by the tool will be put in the directory “tmp/”. There are two folders (“Add/” and “Bio/”) and two scripts in the “bin/” (“SAPshodt.pl” and “format_converter.pl”). “SAPshodt.pl” is a script for generating SAP sequence heterogeneity sequences on demand transmutator. It requires the local modules reposited in the “Add/” and “Bio/”. For each tryptic peptide cut from proteins by trypsin, this script will make a permutation for all SAPs in this peptide, and then adding the new generated peptides into the primary sequence. For “format_converter.pl”, it’s a script for converting Ensmbl-ID SAP dataset into Uniport-ID SAP dataset. You’ll see how to use these two scripts in the following part. Note: no raw dataset is provided with this package.

If you want to add novel N-terminal peptides, you should also install third party predictors: TargetP 1.1 and SignalP 3.0 (http://www.cbs.dtu.dk/services/software.php):

1).	Change the full path in the top of the signalp script 
a.	Do the change in this script: ./SAPshodt/bin/signalp-3.0/signalp
b.	See details in ./SAPshodt /bin/signalp-3.0/signalp-3.0.readme
2).	Change the paths in the top of the targetp script
a.	Do the change in this script: ./SAPshodt /bin/targetp-1.1/targetp
b.	See details in ./SAPshodt /bin/ targetp-1.1/targetp-1.1.readme
c.	It is important to set the tmp directory for TargetP to the tmp directory for SAPshodt ./SAPshodt/tmp

However, due to the TargetP 1.1 and SignalP 3.0’s limitation to UNIX platforms, for SAPshodt only UNIX platforms are feasible in this option.

3.	To prepare input files
1)	Human FASTA format file from Uniprot

The primary protein sequences used to generate a SAP heterogeneity database should be provided with FASTA format. It’s suggested that users get the primary protein sequences from Uniprot database. Users can either get the whole human proteome dataset from the websit (http://www.uniprot.org/) by advanced search or extracted from the UniProtKB/Swiss-Prot and UniProtKB/TrEMBL files reposited in ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/.Sometimes, you may only want to deal with a subset of this dataset, not the whole human proteome, you can list this subset’s all protein accessions in a new file, and one line, one accession. SAPshodt can automatically extract the subset sequences from the whole proteome. 

2)	SAP file from Ensembl

Single Amino acid Polymorphism (SAP) dataset can be downloaded from Ensembl by BioMart (http://www.ensembl.org/index.html). At present (to 2012-12-14), we can choose database “Ensembl Genes 69”, and then choose dataset “Homo sapiens genes (GRCh37.p8)”. In the “Filters”, we restrict the “Gene type” to “protein_coding” in the “GENE”, and “Variation type” to “missense_variant,stop_gained” in the “VARIATION”. In the “Attributes”, we choose “Variation” type and then choose detailed attributes to be included in the output. It should contain “Ensembl Protein ID”, “Protein location (aa)”, “Protein Allele”, “Reference ID” and “Variant Alleles” from “GERMLINE VARIATION INFORMATION”. The format is like this (tab delimited):

Ensembl Protein ID	Protein location (aa)	Protein Allele	Reference ID	Variant Alleles
ENSP00000442112	4	L/Q	rs143466490	T/A
ENSP00000442112	7	M/T	TMP_ESP_15_63889611	T/C
ENSP00000442112	32	L/F	TMP_ESP_15_63889685	C/T
ENSP00000442112	42	V/M	TMP_ESP_15_63889715	G/A
ENSP00000442112	44	E/D	rs142301663	G/C
ENSP00000442112	55	R/C	rs148140490	C/T
ENSP00000442112	78	W/R	rs201654150	T/C
ENSP00000442112	82	R/H	rs200328419	G/A
ENSP00000442112	88	I/N	rs150267919	T/A/C
ENSP00000442112	88	I/T	rs150267919	T/A/C

Note: the order of columns is important. Finally, you check “Unique results only” and then can get the compressed web file by email. 

Since the UniProt Knowledgebase (UniProtKB) is the central hub for the collection of functional information on proteins, with accurate, consistent and rich annotation, it’s suggested to use uniprot IDs instead of ensembl IDs in the following analysis. The ID mapping file can be accessed from   ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz. Considering this simple mapping doesn’t always make sure that the corresponding proteins have the identical sequences, we need two datasets more: human proteome sequence dataset in Ensembl and Uniprot, respectively. The former can be accessed from ftp://ftp.ensembl.org/pub/release-69/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.69.pep.all.fa.gz, and for the latter, if users have downloaded it in the before part as the input file, users can use this file; otherwise, users have to get it either from the websit (http://www.uniprot.org/) by advanced search or extracted from the UniProtKB/Swiss-Prot and UniProtKB/TrEMBL fasta files reposited in the ftp. Note: it’s a better choice to include isoform sequence data. When users get these files, it has to be uncompressed first:
>gunzip HUMAN_9606_idmapping_selected.tab.gz
>gunzip Homo_sapiens.GRCh37.69.pep.all.fa.gz

We have developed a perl script to convert the Ensembl-ID dataset into the Uniprot-ID dataset. You can run it like this:
> perl format_converter.pl -ensembl <fasta_file from ensembl> -uniprot <fasta_file from uniprot> -idmapping <id mapping file from uniprot> -sapfile <SAP file from ensembl> -outfile <output filename>

The process involves:
1) all the SNPs should have a SNP database reference ID (e.g., dbSNP) to show their sources;
2) all indels (amino acids inserts or deletions) are discarded;
3) Most of ensembl protein ID can be mapped to uniprot ID/s; if the corresponding sequences are the same, then we say a ensembl protein is converted into a uniprot protein. For instance, supposing ENSP00000461108 can be mapped to PQ0001, PQ0002, if the sequence of ENSP00000461108 is the same to the PQ0001's,  not PQ0002's, we say ENSP00000461108 is identical to PQ0001. Sometimes, we can find a uniprot protein may have more than one identical ensembl proteins. Then, the SNPs information for each ensembl protein will be merged.
The new output format is similar to the previous’ (tab delimited):

H0YG80	4	L/Q	rs143466490	T/A	ENSP00000442112
H0YG80	7	M/T	TMP_ESP_15_63889611	T/C	ENSP00000442112
H0YG80	32	L/F	TMP_ESP_15_63889685	C/T	ENSP00000442112
H0YG80	42	V/M	TMP_ESP_15_63889715	G/A	ENSP00000442112
H0YG80	44	E/D	rs142301663	G/C	ENSP00000442112
H0YG80	55	R/C	rs148140490	C/T	ENSP00000442112

The first three columns are required, which means protein accession, SAP location and SAP. And the last three columns here are only to indicate the source of SAPs. They’re not required any more in the following steps. Note: no header columns are allowed.
If users want to use their own SAP dataset, users should prepare the SAP dataset in the above’s format. In addition, it’s better to use protein IDs from Uniprot. These two kinds are both allowed:“sp|P31946|1433B_HUMAN” or “P31946”, and this “1433B_HUMAN” is not suggested.

3)	Swiss-Prot annotation file
Users can also choose whether to add information about conflicting sequence entries contained in the Swiss-Prot database. For humans, the Swiss-Prot annotation file can be accessed in ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_human.dat.gz . Then, you can uncompress it:
	>gunzip uniprot_sprot_human.dat.gz
This will produce a file “uniprot_sprot_human.dat”.

4.	To Run
Run it from ./SAPshodt/bin/. You can run it this way, e.g.: 

> nohup perl SAPshodt.pl –fasta <fasta_file> –sapfile <SAP file> –outfile <output file> &

-fasta <fasta_file>     : FASTA-format protein sequence file from Uniprot to be handled
-sapfile <SAP file>    : Uniprot-ID SAP file including “Protein Accession”, “Location” and “SAP” each line
-outfile <output file>  : output file

Additional options:

-h			          : show this help
-inc <file>	      : specify whether only to handle some proteins in the fasta_file (default: all), otherwise please use this                    option to input a file containing protein accessions. These two kinds of accessions are both allowed:                        “sp|P31946|1433B_HUMAN” or “P31946”
-targetp <Y/N>    : specify whether to use TargetP (SignalP) to predict transit or signal peptides, and then to append these                     peptides to the original sequence (default: N) 
-spfile <file>    : specify whether to get conflicting sequences from Swiss-Prot conflict annotations, and then to append                        these peptides to the original sequence (default: N)

5.	Output
1)	SAP heterogeneity database (also with a corresponding decoy database in the same directory, with a “decoy_” beginning filename):

a.	Example of modified uniprot entry:
>O60911 Cathepsin L2 OS=Homo sapiens GN=CTSL2 PE=1 SV=2 lng=334 # SP[336,V,18,V]SNP[371,A,81,G]SNP[411,R,81,G]CON[444,P,81,G] #
MNLSLVLAAFCLGIASAVPKFDQNLDTKWYQWKATHRRLYGANEEGWRRAVWEKNMKMIELHNGEYSQGKHGFTMAMNAFGDMTNEEFRQMMGCFRNQKFRKGKVFREPLFLDLPKSVDWRKKGYVTPVKNQKQCGSCWAFSATGALEGQMFRKTGKLVSLSEQNLVDCSRPQGNQGCNGGFMARAFQYVKENGGLDSEESYPYVAVDEICKYRPENSVANDTGFTVVAPGKEKALMKAVATVGPISVAMDAGHSSFQFYKSGIYFEPDCSSKNLDHGVLVVGYGFEGANSNNSKYWLVKNSWGPEWGSNGYVKIAKDKNNHCGIATAASYPNVJVPKFDQNLDTKJMIELHNGEYSQGKHGFTMAMNAFADMTNEEFRQMMGCFRJMIELHNGEYSQGKHGFTMAMNAFRDMTNEEFRJMIELHNGEYSQGKHGFTMAMNAFPDMTNEEFRQMMGCFR

To this protein four peptides were added: (added information in bold)
- One alternative N-terminal peptide corresponding to the removal of a signal peptide 17 residues long. Position of the appended peptide starts at 336 (residue=V). The original start position was at residue 18.
- Two alternative SNP peptides. Both for which the original position of the actual SNP is 81 (residue=G), and now at position 371 (residue A) and 411 (residue R), respectively.
- One alternative conflict peptide – also from original position 81 and now at position 444.
Note, that while signal and transit peptides positions specify the start of the peptide, the SNP and conflict positions specify the actual position of the SNP/conflict.
The occurrence of a “#” indicates that the sequence has been modified.
The occurrence of the second “#” may be used to detect if the header information has been truncated by e.g. the search engine.

b.	Another example of modified uniprot entry (stop codon gain):
>P12814-3_1 truncated isform 1: Isoform 3 of Alpha-actinin-1 OS=Homo sapiens GN=ACTN1 lng=914 # Variation: 232R/* #
MDHYDSQQTNDYMQPEEDWDRDLLLDPAWEKQQRKTFTAWCNSHLRKAGTQIENIEEDFRDGLKLMLLLEVISGERLAKPERGKMRVHKISNVNKALDFIASKGVKLVSIGAEEIVDGNVKMTLGMIWTIILRFAIQDISVEETSAKEGLLLWCQRKTAPYKNVNIQNFHISWKDGLGFCALIHRHRPELIDYGKLRKDDPLTNLNTAFDVAEKYLDIPKMLDAEDIVGTA

For each stop codon gain variation, it will produce a new entry for the primary entry protein. The ID is created with the primary entry plus a underline “_” and a number. And the sequence is the truncated sequence. “232R/*” means in the position 232, the amino acid R mutates into a stop codon.

2)	Problematic sequences
For some proteins, there are too many SAPs in their certain tryptic peptide. For example, the protein P68871 (http://www.uniprot.org/uniprot/P68871) is not a very long protein. It only contains 147 amino acids. If we cut it into trypsin peptides according to K|R rule that not followed by P, we can get the following peptides:
“9 18 31 41 60 62 66 67 83 96 105 121 133 145 147”
the first number 9 means 1-9 of P68871 is a trypsin peptide, and 18 means 10-18 is the second peptide...
Meanwhile, we're able to get its corresponding SAPs information from ensembl. So, if we permute each SAP in some certain peptide, the number of permutation will be very large. For example, for the peptide 42-60 of P68871, it contains following SAPs:
42:F/Y/C/S (it means at the 42nd position, a F can mutate into a Y or C or S)
43:F/S/V/L
44:E/A/K/Q
45:S/C
46:F/Y/C/S
47:G/E/R
48:D/V/G/A/Y/N/H
49:L/R/P
50:S/F/C
51:T/N/S
52:P/H/R/S/L
53:D/V/G/A/N/H
54:A/T
55:V/D
56:M/K
57:G/D/C/R
58:N/K/S/T/D/H
59:P/R
60:K/N/T/E

We find SAP can occur in each position. And we can compute the number of permutation easily: 4×4×4×2×4×3×7×3×3×3×4×6×2×2×2×4×6×2×4. This number is > 10 billion. Obviously, this is a very large number. And considering each peptide has at least 19 amino acids,  the new generated peptides will contain > 190 billion amino acids only for this peptide 42-60 of P68871. This can not be accepted, whatever for the hardware or software. So, it seems like it's not a good attempt to permute SNPs for any sequence. As we have seen, this will lead to the crash of softwares. We have applied a simplified strategy to handle this kind of peptides. For each SAP, we will produce a new peptide, but no any permutation between different SAPs. That is, one SAP, one new peptide. For the example above, it will only produce 3+3+3+1+3+2+6+2+2+2+3+5+1+1+1+3+5+1+3 = 50 new peptides. If there are such kind of proteins in the input file, we list all this kind of proteins in the “../tmp/problematic_proteins.txt”. Users can check it.

6.	Copyright & problems
Copyright (c) 2012: Proteome Center Rostock.

Author: Lu-Lu Zheng (mingkanghust@gmail.com).

This is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Please contact mingkanghust@gmail.com in case of problems.

