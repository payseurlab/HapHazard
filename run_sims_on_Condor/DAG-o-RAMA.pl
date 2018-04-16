
#!usr/bin/perl

#DAG-o-RAMA.pl

$jobListFile = "DAG_jobs.txt";
#$runs_per_job = 1000; 

open DAGJOBS, $jobListFile or die "Could not open DAG jobs file $jobListFile!\n";
@dagjobs = <DAGJOBS>; foreach (@dagjobs) { chomp $_; } close DAGJOBS;

open MASTERDAG, ">master.dag";

foreach $job (@dagjobs)
{
	@job = split( ',' , $job);
	$job = @job[0];
	$memoryNeeded = "4.5GB";
	$runs_per_job = @job[1];
	print MASTERDAG "SPLICE $job $job/$job.dag\n";
	mkdir $job, 0755;
	mkdir "$job/shared", 0755;
	system("perl MakeSimInp.pl $job");
	

	$job_dagfilename = $job . "/" . $job . ".dag";
	open DAGFILE, ">$job_dagfilename";

	for( $i = 0 ; $i < $runs_per_job ; $i++ )
	{
		print DAGFILE "JOB $job" . "_" . $i . " $job/$job" . "_$i/$job" . "_$i.submit\n"; # DIR $job/$job" . "_$i/\n";
		$jobdir = $job . "/" . $job . "_" . $i;
		#print $jobdir . "\n";
		mkdir $jobdir, 0755;

		$Rversion = 'sl5-R.3.1.0';

		# make the submit file for the run
		$inpFileName = $jobdir . "/shared/" . $job . '.inp';
		$seed = int(rand(9007199254740992));

		$runFile = $jobdir . "/run_hap_hazard.R";
		open RUNFILE, ">$runFile";
		print RUNFILE "print(\"Running Hap_Hazard\")
seed <- " . $seed . "
input <- \"$job.inp\"
id <- " . $i . "
print(\"Running Simulation\")
command <- paste(\"./HapHazard\",input,id,seed,sep=\" \")
system(command)
print(\"Loading simulation parameters to R-workspace...\")
source(\"Rvars.R\")
print(\"Running ancestry block and junction analyses...\")
#source(\"HH_BJ_1s.R\")
print(\"Computing Hybrid Indices...\")
#source(\"HH_hybridIndex.R\")
#print(\"Computing junction ages and frequency spectra...\")
#source(\"HH_AgeFS.R\")
#print(\"Computing geographic clines...\")
#source(\"ClineMCMC2.R\")
print(\"Moving files...\")
system(\"mv Rvars.R " . $job . "_" . $i . "/Rvars.R\")
print(\"Tarring results folder...\")
system(\"tar -czf " . $job . "_" . $i . ".tar.gz " . $job ."*\" )
print(\"Done\")
";
		system("chmod +x $runFile");
		$subFileName = $jobdir . "/" . $job . "_" . $i . '.submit';
		open SUBMITFILE, ">$subFileName";

		#### several of the next lines violate the format scheme for readability and so that they appear
		### as they will in the file being printed to...
		print SUBMITFILE"universe = vanilla
executable = HH_executable.sh
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = ../../shared/, ../shared/, run_hap_hazard.R, http://proxy.chtc.wisc.edu/SQUID/mfrayer/SLIBS.tar.gz
#arguments = --type=R --version=sl5-R-3.1.0 --cmdtorun=run_hap_hazard.R --unique=\$(Cluster)\.\$(Process) --
requirements = (OpSysMajorVer =?= 6)
#+R = \"sl5-R-3.1.0\"
output = \$(Cluster)\.\$(Process).out
error = \$(Cluster)\.\$(Process).err
log = \$(Cluster)\.\$(Process).log
request_cpus = 1
request_memory = " . $memoryNeeded . "
request_disk = 500000
Initialdir = " . $job . "/" . $job . "_$i/\n
+wantGlidein = true
+wantFlocking = true\n";
#Initialdir = " . $job . "_" . $i . "/" . $job . "_$i\n";

	#	if($memoryNeeded <= 2000)
	#	{
	#		print SUBMITFILE "+wantGlideIn = true\n";
	#	}

		print SUBMITFILE "queue 1";

		close SUBMITFILE;

	}
}


