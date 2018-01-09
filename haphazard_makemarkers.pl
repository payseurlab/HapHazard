
#!usr/bin/perl

# This script will make marker files from the block data

# open and read in the list of files from the output parsing file
$fileListName = "block_file_list.txt";
open FILELIST, $fileListName or die "Could not open file list $fileListname!\n";
@fileList = <FILELIST>; foreach (@fileList) { chomp $_; } close FILELIST;

# process each
foreach $file (@fileList)
{
	@fileOp = split( ',' , $file );
	&MakeMarkerFile(@fileOp[0], @fileOp[1], @fileOp[2]);
}


exit;


sub MakeMarkerFile()
{
	
	my $blockFileName = @_[0];
	my $markerSpacing = @_[1];
	my $chromosomeLength = @_[2];

	#print $blockFileName . " " . $markerSpacing . " " . $chromosomeLength . "\n";


	# open the block file, or exit the subroutine
	unless (open BLOCKFILE, "$blockFileName" )
	{
		print "Could not open block file $blockFileName!\n";
		return ;
	}

	# read in the data from the block file and make a header row from the fields in the block file 
	my @blockFile = <BLOCKFILE>; foreach (@blockFile) { chomp $_; } close BLOCKFILE;

	my @header = split( ',' , @blockFile[0] );
	splice @header, -6;

	my $position = 0;

	# then append the positions of the markers to the header
	until ($position > $chromosomeLength)
	{
		push @header, sprintf("%.6f", $position);
		$position += $markerSpacing;
	}

	my $header = join( ',' , @header);

	# change the appendix of the block file name to .mkr, and create a markers file with that name
	@markerFileName = split( '' , $blockFileName );
	splice @markerFileName, -3;
	$markerFileName = join( '' , @markerFileName) . "mkr";
	
	open MARKERFILE, ">$markerFileName";
	print MARKERFILE $header . "\n";


	#scan the block file and add markers along blocks according to ancestry
	$position = 0;
	for( $i = 1 ; $i <= scalar @blockFile ; $i++ )
	{
		
		# get the current and previous block so we can check to see if they are on the same chromosome
		my @thisBlock = split( ',' , @blockFile[$i] );
		my @prevBlock = split( ',' , @blockFile[$i-1] );
		
		# if they are on different chromosomes
		if( $i == 1 || @thisBlock[2] ne @prevBlock[2] )
		{
			# and the current block is not the first element of the file
			unless ( $i == 1 )
			{
				# print the ancestry of the block for all of the markers whose positions
				# lie in that block
				while( $position <= @prevBlock[9] )
				{
					print MARKERFILE "," . @prevBlock[5] ;
					$position += $markerSpacing;
				}

				print MARKERFILE "\n";		
			}
	
			$position = 0;
			
			# make a tag carrying information about the haplotype for the new haplotype
			my $hapTag = join( ',' , @thisBlock[0..4]);
			print MARKERFILE "$hapTag" ; 
		}
		else
		{
			# print the ancestry of the block for all of the markers whose positions
			# lie in that block
			while($position <= @prevBlock[9] )
			{
				#print $i . " " . $position . " " . @prevBlock[9] . "\n";
				print MARKERFILE "," . @prevBlock[5] ;
				$position += $markerSpacing;
			}
		}
	}

	close MARKERFILE;

}
