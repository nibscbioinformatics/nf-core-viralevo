use strict;
use warnings;
use Getopt::Long;
use Cwd;
use IO::File;


##########################################
# Script options management
##########################################

#Pick up script name automatically for usage message
my $script=substr($0, 1+rindex($0,'/'));

#Set usage message
my $usage="Usage: $script -ivar original_ivarfile.tsv -vcf output_vcf_file
Please try again.\n\n\n";

#Declare all variables needed by GetOpt
my ($ivar, $vcf);

#Get command-line parameters with GetOptions, and check that all needed are there, otherwise die with usage message
die $usage unless
	&GetOptions(	'ivar:s' => \$ivar,
					'vcf:s' =>\$vcf
				)
	&& $ivar;

$vcf ||= $ivar.".vcf";

my $vcf_format = "GT:PVAL:AQ:DP:AF";


#############################
# file handling
#############################

my $total_lines=0;
open (FILE,$ivar) or die "Can't open $ivar: $!";
print STDERR "counting $ivar\n";
while (<FILE>){
	if ($_=~/#/) {next;}
	else{$total_lines++;}
}
close FILE;

open(my $ivarHandler, $ivar) or die("can't open original iVar file\n");
open(my $vcfHandler, ">$vcf") or die("can't write on the output vcf file\n");

my $headerLines = <<HEAD;
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="Passed iVar Filters">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PVAL,Number=1,Type=Float,Description="PValue of the variant call">
##FORMAT=<ID=AQ,Number=2,Type=Integer,Description="Allele Quality">
##FORMAT=<ID=DP,Number=2,Type=Integer,Description="Allelic Read Depth">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Alternative Allele Fraction">
##INFO=<ID=IVAR_DP,Number=1,Type=Integer,Description="Total Read Depth at site emitted by iVar">
##INFO=<ID=IVAR_GFF,Number=1,Type=String,Description="Annotation in iVar from GFF file">
##INFO=<ID=IVAR_REFAA,Number=1,Type=String,Description="iVar annotated reference aminoacid">
##INFO=<ID=IVAR_ALTAA,Number=1,Type=String,Description="iVar annotated alternative aminoacid">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$ivar
HEAD

## print the header in file
print $vcfHandler $headerLines."\n";

my $isInfo = 0;
my $fileLine=0;

my $previousID = "";

while(<$ivarHandler>){

	chomp();

	$fileLine++;
    my $percent = ($fileLine / $total_lines) *100;
    my $rounded = sprintf "%.3f", $percent;
    print STDERR "parsing and annotating results file, progress: ---  $rounded %       \r";

  if ($_=~/REGION/) {next;}

  my @data = split("\t", $_);

  ## handle the peculiar way of writing insertions and deletions
  ## assumption is that deletion or insertion is positioned downstream
  ## to reference allele  - this has been verified on UCSC browsert

  my $ref = "";
  my $alt = "";
  my $geno = "";

  if ($data[3]=~/\+/){
    $data[3]=~s/\+//;
    $ref = $data[2];
    $alt = $ref.$data[3];
  }
  elsif ($data[3]=~/-/){
    $data[3] =~s/-//;
    $ref = $data[2].$data[3];
    $alt = $data[2];
  }
  else {
    $ref = $data[2];
    $alt = $data[3];
  }
  $geno = "$ref/$alt";

	my $identifier = $data[0]."_".$data[1]."_".$ref."_".$alt;


  my $variantLine = "$data[0]\t$data[1]\t.\t$ref\t$alt\t.\t$data[13]\tIVAR_DP=$data[11];IVAR_GFF=$data[14];IVAR_REFAA=$data[16];IVAR_ALTAA=$data[18];";
  $variantLine .= "\tGT:PVAL:AQ:DP:AF\t$geno:$data[12]:$data[6],$data[9]:$data[4],$data[7]:$data[10]";

  if ($previousID ne $identifier){
    print $vcfHandler $variantLine."\n";
    $previousID = $identifier;
  }
  else {
    next;
  }


}

print STDERR "\n";
print STDERR "completed\n";
