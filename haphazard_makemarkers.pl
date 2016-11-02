
#!usr/bin/perl

# This script will make marker files from the block data

# open and read in the list of files from the output parsing file
$fileListName = "block_file_list.txt";
open FILELIST, $fileListName or die "Could not open file list $fileListname!\n";
@fileList = <FILELIST>; foreach (@fileList) { chomp $_; } close FILELIST;

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

	#print "Making file $blockFileName ($markerSpacing, $chromosomeLength)...\n";

	unless (open BLOCKFILE, "$blockFileName" )
	{
		print "Could not open block file $blockFileName!\n";
		return ;
	}
	my @blockFile = <BLOCKFILE>; foreach (@blockFile) { chomp $_; } close BLOCKFILE;

	my @header = split( ',' , @blockFile[0] );
	splice @header, -4;

	my $position = 0;

	until ($position > $chromosomeLength)
	{
		push @header, sprintf("%.6f", $position);
		$position += $markerSpacing;
	}

	my $header = join( ',' , @header);

	@markerFileName = split( '' , $blockFileName );
	splice @markerFileName, -3;
	$markerFileName = join( '' , @markerFileName) . "mkr";
	
	open MARKERFILE, ">$markerFileName";
	print MARKERFILE $header . "\n";

	$position = 0;
	for( $i = 1 ; $i <= scalar @blockFile ; $i++ )
	{
		my @thisBlock = split( ',' , @blockFile[$i] );
		my @prevBlock = split( ',' , @blockFile[$i-1] );
		

		if( $i == 1 || @thisBlock[3] ne @prevBlock[3] )
		{

			unless ( $i == 1 )
			{
				while( $position <= @prevBlock[7] )
				{
					print MARKERFILE "," . @prevBlock[5] ;
					$position += $markerSpacing;
				}

				print MARKERFILE "\n";		
			}
	
			$position = 0;

			my $hapTag = join( ',' , @thisBlock[0..4]);
			print MARKERFILE "$hapTag" ; 
		}
		else
		{
			while($position <= @prevBlock[7] )
			{
				print MARKERFILE "," . @prevBlock[5] ;
				$position += $markerSpacing;
			}
		}
	}

	close MARKERFILE;

}
