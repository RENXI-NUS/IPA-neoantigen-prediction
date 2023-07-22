#!/usr/bin/perl
use strict;
use warnings;
my $sam_file=shift;
my $name=shift;
my $out_PAS_file=$name."_potential_polyA_sites_from_softclipped_reads.txt";
my @quality_list;
my $char;

sub average {
my @array = @_; # save the array passed to this function
my $sum; # create a variable to hold the sum of the array's values
foreach (@array) { $sum += $_; } # add each element of the array 
# to the sum
return $sum/@array; # divide sum by the number of elements in the
# array to find the mean
}

open(my $IN, "-|", "samtools view -F 4 -f 0x2 $sam_file") || die $!;
open(OUT, ">>$out_PAS_file")or die $!;
my $line1 = <$IN>;
my $line2 = <$IN>;

while(($line1)&&($line2)){
chomp $line1;
chomp $line2;

my @field1=split /\s+/,$line1;
my @field2=split /\s+/,$line2;
if($field1[0] =~ /^@/){
$line1 = $line2;
$line2 = <$IN>;
}
else{
if($field1[0] eq $field2[0]){
	if ($field1[5] =~ /^(([3-9]|[1-9]\d)S)/){
                my $index = $1;
                $index =~ s/\D//g;
                #$index = $index -1;
		if ($field1[9] =~ /^(T)\1{3,40}/){
                if ($field1[1] & 256) {
                $line1=<$IN>;
                $line2=<$IN>;
                }
                else{
                if ($field2[1] & 512) {
                $line1=<$IN>;
                $line2=<$IN>;
                }
                else{
                if ($field2[1] & 1024) {
                $line1=<$IN>;
                $line2=<$IN>;
                }
                else{
		my $seq = substr($field1[9],$index);
                print OUT $field1[2];
                print OUT "\t";
                print OUT ($field1[3]);
                print OUT "\t";
                print OUT "-";
                print OUT "\t";
		print OUT $field1[5]."\t";
                print OUT $field1[9]."\t";
                print OUT $field1[10]."\t";
                print OUT $field2[9]."\t";
                print OUT $field2[10]."\t";
                print OUT (length($seq))."\t";
		print OUT $seq;
		print OUT "\t";
                print OUT $field1[0]."\t";
		print OUT $field1[1];
		print OUT "\n";
                $line1=<$IN>;
                $line2=<$IN>;
                }
		}
		}
		}
		else{
                $line1=<$IN>;
                $line2=<$IN>;
                }
        }
	else{
	$line1=<$IN>;
	$line2=<$IN>;
	}
}
else {
$line1 = $line2;
$line2 = <$IN>;
	if($field1[0] eq $field2[0]){
        if ($field1[5] =~ /^(([3-9]|[1-9]\d)S)/){
                my $index = $1;
                $index =~ s/\D//g;
                #$index = $index -1;
		if ($field1[9] =~ /^(T)\1{3,40}/){
                if ($field1[1] & 256) {
                $line1=<$IN>;
                $line2=<$IN>;
                }
                else{
                if ($field2[1] & 512) {
                $line1=<$IN>;
                $line2=<$IN>;
                }
                else{
                if ($field2[1] & 1024) {
                $line1=<$IN>;
                $line2=<$IN>;
                }
                else{
		my $seq = substr($field1[9],$index);
                print OUT $field1[2];
                print OUT "\t";
                print OUT ($field1[3]);
                print OUT "\t";
                print OUT "-";
                print OUT "\t";
		print OUT $field1[5]."\t";
                print OUT $field1[9]."\t";
                print OUT $field1[10]."\t";
                print OUT $field2[9]."\t";
                print OUT $field2[10]."\t";
                print OUT (length($seq))."\t";
		print OUT $seq;
		print OUT "\t";
                print OUT $field1[0]."\t";
		print OUT $field1[1];
		print OUT "\n";
                $line1=<$IN>;
                $line2=<$IN>;
                }
		}
		}
		}
		else{
                $line1=<$IN>;
                $line2=<$IN>;
                }
        }
	else{
	$line1=<$IN>;
	$line2=<$IN>;
	}
	}
	else {
	last;
	print "error~~~";
	}
}
}
}
close($IN);
close(OUT);
