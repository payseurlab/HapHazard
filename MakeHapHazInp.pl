
#!usr/bin/perl

# A driver script to run the junction generator


#print "Type the name of the file that contains your input for the experiment.\n";
$inputFile = @ARGV[0]; 
chomp $inputFile;

$numJobs = @ARGV[1];
chomp $numJobs;

$memoryNeeded = @ARGV[2];
chomp $memoryNeeded;

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


open TEMPINP, ">@input[0].inp";

foreach (@input) { print TEMPINP "$_\n"; }

close TEMPINP;

print "Temp File \"" . @input[0] . ".inp\" Made \n";


open EXECUTABLE, ">run_hap_hazard";

print EXECUTABLE "#!/bin/sh
echo \"Running Hap_Hazard\"
./JunGen4 \$1 \$2 \$3
perl haphazard_makemarkers.pl
tar -cvzf @input[0].tar.gz *
echo \"Done\"";

close EXECUTABLE;

print "Executable, \"run_hap_hazard\" made.\n";

system('chmod +x run_hap_hazard');

$submitFile = @input[0] . "-chtc.sub";
open CHTCSUBMIT, ">$submitFile";

$seed = int(rand(9007199254740992));

print CHTCSUBMIT "#
# haphazard-chtc.sub
#
# The name of the job
#

job = haphazard

#
# Specify the HTCondor universe, executable, and log file
#

universe = vanilla
executable = run_hap_hazard
log = \$(job)_\$(Cluster).log

#
# Specify that HTCondor should transfer files to and from
# the remote execution hosts for us. 
#

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = JunGen4, " . @input[0] . ".inp, haphazard_makemarkers.pl
arguments = " . @input[0] . ".inp \$(Process) " . $seed ."
output = \$(job)_\$(Cluster)_\$(Process).out
transfer_output_files = " . @input[0] . "_\$(Process).tar.gz

#
# Request cpus, memory, and storage space
#

request_cpus = 1
request_memory = " . $memoryNeeded . "
request_disk = 14000

#
# Tell HTCondor to run 10 instance of the job
#

queue " . $numJobs . "
";
close CHTCSUBMIT;

print "CONDOR submit file, $submitFile, made.\n";

exit;


