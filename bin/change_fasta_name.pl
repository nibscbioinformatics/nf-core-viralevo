use warnings;
use Getopt::Long;
use IO::File;

#Pick up script name automatically for usage message
my $script=substr($0, 1+rindex($0,'/'));

#Set usage message
my $usage="Usage: $script -fasta sequence.fasta -name new_header -out output.fasta\nPlease try again.\n\n\n";

#Declare all variables needed by GetOpt
my ($fasta, $name, $out);

#Get command-line parameters with GetOptions, and check that all needed are there, otherwise die with usage message
die $usage unless
	&GetOptions(
					'fasta:s' => \$fasta,
					'name:s' => \$name,
          'out:s' => \$out
				)
	&& $fasta && $name;

$out ||= $fasta."_renamed.fasta";

open(FILE,$fasta) or die "Can't open $fasta: $!";
open(OUT, ">$out") or die "Can't write on $out: $!";

my $count=0;
my $append="";

while(<FILE>){
  chomp($_);
  if ($_=~/>/){
    if ($count >0){
      my $append = $count;
    }
    my $newline = ">".$name.$append;
    print OUT $newline."\n";
		$count++;
  }
  else {
    print OUT $_."\n";
  }
}
close FILE;
close OUT;
