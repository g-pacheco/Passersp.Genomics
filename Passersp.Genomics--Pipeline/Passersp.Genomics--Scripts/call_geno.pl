#!/usr/local/bin/perl -w
#
# Cared for by Filipe G. Vieira <>
#
# Copyright Filipe G. Vieira
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

    call_geno.pl v0.0.3

=head1 SYNOPSIS

    perl call_geno.pl [--skip 0] [--ifs \t] [--ofs \t] [--min_prob -Inf] [--miss_data -1]

    --skip       INT    = Number of columns to initially skip [3]
    --ifs        STRING = Input field separator [\t]
    --ofs        STRING = Output field separator [\t]
    --min_prob   DOUBLE = Minimum probability to call a genotype; if lower set as missing data [-Inf]
    --miss_data  STRING = Missing data character [-1]

=head1 DESCRIPTION

    This script will call genotypes from either Genotype Likelihoods (GLs) 
    or Genotype Posterior Probabilities (PP). It reads either from a file 
    or STDIN and accepts several options to skip first N columns (e.g. if 
    labels), set input and output field separators, define minimum probability 
    to call a genotype (if lower will set genotype to missing), and define 
    missing data character.

    It can output several GENO formats:
                          HOM1   HET    HOM2
        angsd:             0      1      2
        ped/tped (plink): 2 2    1 2    1 1

=head1 AUTHOR

    Filipe G. Vieira - fgarrettvieira _at_ gmail _dot_ com

=head1 CONTRIBUTORS

    Additional contributors names and emails here

=cut


# Let the code begin...

use strict;
use Getopt::Long;
use List::MoreUtils qw(minmax firstidx);

my ($skip, $ifs, $ofs, $out_format, $min_prob, $miss_data);
my ($cnt, $geno, @line, @gl, @mm);

$skip = 3;
$ifs = "\t";
$ofs = "\t";
$out_format = 'angsd';
$min_prob = -Inf;

GetOptions('h|help'          => sub { exec('perldoc',$0); exit; },
           's|skip:s'        => \$skip,
	   'ifs:s'           => \$ifs,
	   'ofs:s'           => \$ofs,
	   'of|out_format:s' => \$out_format,
	   'p|min_prob:s'    => \$min_prob,
	   'm|miss_data:s'   => \$miss_data,
    );

if($out_format eq 'angsd') {
    $miss_data = -1 unless($miss_data);
} elsif($out_format eq 'tped') {
    $miss_data = "0 0" unless($miss_data);
} else {
    die("ERROR: invalid output format.\n");
}

while(<>) {
    $cnt=0;
    @line = split(/$ifs/,$_);

    # Skip and print first $skip columns
    print(join($ofs,splice(@line,0,$skip))."\t") if($skip > 0);

    while( @gl=splice(@line,0,3) ) {
	$cnt++;
	@mm = minmax(@gl);
	$geno = ($mm[0]==$mm[1] ? $miss_data : firstidx {$_==$mm[1]} @gl);
	$geno = $miss_data if($mm[1] < $min_prob);
	if($geno ne $miss_data){
	    if($out_format eq 'tped'){
		if($geno == 0) {
		    $geno = '2 2';
		} elsif($geno == 1) {
		    $geno = '1 2';
		} elsif($geno == 2) {
		    $geno = '1 1';
		}
	    }
	}
	print( ($cnt>1 ? "\t" : "").$geno );
    }
    print "\n";
}

exit(0);
