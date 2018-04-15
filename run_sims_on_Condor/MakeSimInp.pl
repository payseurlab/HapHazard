
#!usr/bin/perl

$jobName = @ARGV[0];

#print "Type the name of the file that contains your input for the experiment.\n";
$inputFile = "HH_" .@ARGV[0] . ".inp"; 
chomp $inputFile;

open(INPUT, $inputFile) or die "Could not open file $inputFile!\n";
@inputFile = <INPUT>; foreach (@inputFile) { chomp $_; } close INPUT;

@input;
@lines;

foreach $line (@inputFile)
{
	@line = split( ' ' , $line);
	
	if( @line[0] eq '>' )
	{
		shift @line;
		push(@lines, $line);
		foreach $entry (@line)
		{
			#print $entry . "\n";
			push(@input, $entry);
		}
	}
}


open TEMPINP, ">$jobName/shared/$jobName.inp";

foreach (@input) { print TEMPINP "$_\n"; }

close TEMPINP;

print "Temp File \"" . @input[0] . ".inp\" Made \n";


