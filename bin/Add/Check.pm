
package Check;

# This package check whether there are some problematic sequences that 
# containing too many SAPs in the trypsic region

##################################################
### - Lulu Zheng - March 2013 #####
##################################################

use strict;

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
		_snp => $args{-snp},
		_pep => $args{-pep},
	},$class;
	
	$class->_incr_count();
	return $self;
}

sub Check{
	my ($self) = @_;
	my ($snp,$pep_hash) = ($self->{_snp},$self->{_pep});
	my %problematic = ();
	foreach my $protein(keys %{$snp}){
		next unless exists $pep_hash->{$protein};
		
		my @locs = sort {$a <=> $b} keys %{$snp->{$protein}};
		my $string = $pep_hash->{$protein}->seq();
		my $pri_string = $string;
		$pri_string =~ s/J.*//;
		
		my @positions = ();
		
		while($pri_string =~ /([K|R])/g){
			next if substr($string,$-[1]+1,1) eq "P";
			push @positions, $-[1] + 1;
		}
		push @positions,length($pri_string);
		
		for(my $i = 0; $i < @positions; $i++){
			my @hits = ();
			for(my $j = 0; $j < @locs; $j++){
				my @mutations = split /\//,$snp->{$protein}->{$locs[$j]};
				if($locs[$j] <= $positions[$i] && $locs[$j] + length($mutations[0]) - 1 > $positions[$i]){
					print $protein.":".$locs[$j].":".$snp->{$protein}->{$locs[$j]}."\t cross the cleavage site\n";
					next;
				}
				elsif($i == 0){	#### for the first trypsin peptide in the primary sequence
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
				$num += @mutations - 1;
			}
			if($num >= 21){
				$problematic{$protein} = $pep_hash->{$protein}->desc();
				last;
			}
		}
	}
	return \%problematic;
}

# When an object is no longer being used, this will be automatically called
# and will adjust the count of existing objects
sub DESTROY {
    my($self) = @_;
    $self->_decr_count();
}

1;
