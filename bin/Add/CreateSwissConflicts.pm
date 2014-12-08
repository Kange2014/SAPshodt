
package CreateSwissConflicts;

# This package controls the parsing of SwissProt data
# It is part of the MS Friendly database scripting framework
# It parses the Swissprot part of UniProt for:
# Conflict annotation - being used to create variant tryptic peptides

##################################################
### - Lulu Zheng - March 2013 #####
##################################################

use warnings;
use strict;

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
	},$class;
	
	$class->_incr_count();
	return $self;
}

sub CreateSwissConflicts{
	my ($self) = @_;
	#Open files for Output
	open MYFILE1, '>../tmp/swiss_conflicts.txt';	# Feature annotation for conflicts

	#declare variables
	my $ac="";
	my $feat="";
	my $start="";
	my $end="";
	my $desc="";
	my $euk="";

	open REF, $self->{_file} or die;
	while (my $line=<REF>){
		chomp($line);
		# Pick the accession number
        if ($line =~ /^AC\s+(\w+);/x) {
        	$ac=$1;                                                   
		} 	

		#Pick an indicater that the current entry is human
		if ($line =~ /^OS\s+(Homo\ssapiens)/x) {
			$euk=$1;
		}

		# Pick the Conflict annotation FEATURES
	  	if ($line =~ /^FT\s+(CONFLICT)\s+(\d+)\s+(\d+)\s+(.*)$/x) {
			$feat=$1;
			$start=$2;
			$end=$3;
			$desc=$4;
			while($desc !~ /\)\.$/){ 
				$line = <REF>;
				chomp($line);
				my @tmp = split /\s+/,$line;
				last if $tmp[0] ne "FT";
				my $index = @tmp - 1;
				$desc .= " ".join(" ",@tmp[1..$index]) if $tmp[1] =~ /^\(/;
				$desc .= join(" ",@tmp[1..$index]) if $tmp[1] !~ /^\(/;
			}
            #Only print out conflict annotations if a human record - These are printed on the fly
	        if($euk) {
	                       print MYFILE1 "$ac $feat\t$start\t$end\t$desc\n";
	                       $feat="";
	                       $start="";
	                       $end="";
	                       $desc="";
            }
        }
		
 		# End of this record, print all general information and go on in the while loop
		if ($line =~ /^ID\s+/x) {
            #Reset variables
            $euk="";
			$ac="";
		}
	}


}

sub DESTROY {
    my($self) = @_;
    $self->_decr_count();
}

1;

__END__   
