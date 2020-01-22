#!/usr/bin/perl -w
use strict;
=use
To find a peak for each identified cluster ...
...
=cut

sub findPeak {
	my ($region,$reads,$output) = @_;
	my %Locus;
	open IN,"$reads" or die $!;
	while(<IN>){
		chomp;
		my @Columns = split;
		if($Columns[5] =~ /\-/){
			if(!(exists $Locus{$Columns[0]}{$Columns[5]}{$Columns[2]})){
        		$Locus{$Columns[0]}{$Columns[5]}{$Columns[2]} = 1;
			}else{
				my $c = $Locus{$Columns[0]}{$Columns[5]}{$Columns[2]};
				$c ++;
				$Locus{$Columns[0]}{$Columns[5]}{$Columns[2]} = $c;
			}
		}elsif($Columns[5] =~ /\+/){
			if(!(exists $Locus{$Columns[0]}{$Columns[5]}{$Columns[1]})){
					$Locus{$Columns[0]}{$Columns[5]}{$Columns[1]} = 1;
			}else{
					my $c = $Locus{$Columns[0]}{$Columns[5]}{$Columns[1]};
					$c ++;
					$Locus{$Columns[0]}{$Columns[5]}{$Columns[1]} = $c;
			}       
		 }else{
		 	print "error: strand $Columns[5] is not recoglized\n";
		 }
	}
	close IN;
	open IN,"$region" or die $!;
	open OUT,">$output" or die $!;
	while(<IN>){
		chomp;
		my @Columns = split;
		my $c = 0;
		my $p = 0;
		my $sum = 0;
		for($Columns[1]..$Columns[2]){
			my $l = $_;
			if(exists $Locus{$Columns[0]}{$Columns[5]}{$l}){
				$sum = $sum + $Locus{$Columns[0]}{$Columns[5]}{$l};
				if($Columns[5] eq "-"){
					if(!($Locus{$Columns[0]}{$Columns[5]}{$l} < $c)){
						$p = $l;
						$c = $Locus{$Columns[0]}{$Columns[5]}{$l};
					}
				}elsif($Columns[5] eq "+"){
					if($Locus{$Columns[0]}{$Columns[5]}{$l} > $c){
						$p = $l;
						$c = $Locus{$Columns[0]}{$Columns[5]}{$l};
					}
				}else{
					print "strand information: $Columns[5] is not defined\n";
				}
			}
		}
		my $id = $Columns[0]."|".$Columns[5]."|NA|".$p."|".$sum.";".$Columns[3];
		if($Columns[5] eq "-"){
			my $p2 = $p - 1;
			print OUT "$Columns[0]\t$Columns[1]\t$Columns[2]\t$id\t$sum\t$Columns[5]\t$p2\t$p\n";
		}elsif($Columns[5] eq "+"){
			my $p2 = $p + 1;
			print OUT "$Columns[0]\t$Columns[1]\t$Columns[2]\t$id\t$sum\t$Columns[5]\t$p\t$p2\n";
		}else{
			print "error:strand $Columns[5] information can not recognized\n";
		}
	}
	close OUT;
	close IN;
}

if(@ARGV < 3){
	print "usage:idr_out\treads_bed\toutput\n";
	exit(0);
}else{
	&findPeak($ARGV[0],$ARGV[1],$ARGV[2]);
}
