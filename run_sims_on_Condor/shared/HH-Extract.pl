#!usr/bin/perl


$expName = @ARGV[0];

$numOfExps = @ARGV[1];

for( $i = 0 ; $i < $numOfExps ; $i++ )
{
	print $expName . " " . $i . "\n";
	$tarCom = "tar -xvzf " . $expName . "_" . $i . ".tar.gz";
	system($tarCom);
	
}

exit;
