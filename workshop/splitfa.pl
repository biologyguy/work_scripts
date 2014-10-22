#!/usr/bin/perl
#use strict;   #Perl pragma to restrict unsafe constructs
#use warnings; #warnings about dubious constructs
use Getopt::Long; #Implements GetOptions()
use Bio::SeqIO;

#warning message
my $usage = <<_EOUSAGE_;

#########################################################################################
# splitfa.pl INPUT --split <NUMBER_OF_SEQUENCES_PER_FILE> --output <Output_Prefix>
# example command: splitfa.pl big_file.fa --split 100 -output small_file 
#
# Required:
# > BioPerl
#
##########################################################################################

_EOUSAGE_
	;
#calling the variables	
our $split;        
our $out;

#hash that puts user imput values into variables
&GetOptions( 'split=s' => \$split,
            'output=s' => \$out,);
unless ($split && $out) {
die $usage;
}
print $ARGV[0]."\n";
$input=$ARGV[0];

$num_Of_Fasta = $split;
my $x=1;
my $n=1;
$output= $out. $n . ".fa";

#open(my $fh, '>', $output);
$seq_out = Bio::SeqIO->new(-file => ">$output", -format => "fasta" );
$seqio_obj = Bio::SeqIO->new(-file => "<$input", -format => "fasta" );

while ($seq_obj = $seqio_obj->next_seq){
    if ($x>$num_Of_Fasta) {
        #close $fh;
        $x=1;
        $n=$n+1;
        $output= $out. $n . ".fa";
        $seq_out = Bio::SeqIO->new(-file => ">$output", -format => "fasta" );
#        open($fh, '>', $output);
    }
    $seq_out->write_seq($seq_obj);
#    print $fh ">".$seq_obj->id,"\n";
#    print $fh $seq_obj->seq,"\n";
    $x=$x+1;
}
#close $fh;
