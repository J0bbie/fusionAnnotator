#!/usr/bin/perl
#Author: 			Job van Riet
#Date of creation:		10/9/12
#Date of modification:		14/1/13
#Known bugs:			None
#Function:			This module will house all of the functions (subroutines) of the fusion application.
#				It is therefore needed to import this module at every fusion script.

#Define package declaration and tags.
package fusionModule;

use strict;
use Exporter;
#Use DBI for database connectivity
use DBI;
#Diagnostics for debugging
use diagnostics;
#Use Bioperl and Ensembl API 67 for function
use lib '/usr/local/ensembl67/ensembl/modules';
use lib '/usr/local/ensembl67/ensembl-compara/modules';
use lib '/usr/local/ensembl67/ensembl-external/modules';
use lib '/usr/local/ensembl67/ensembl-functgenomics/modules';
use lib '/usr/local/bioperl1.2.3/bioperl-live';

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
#Used for file handling and temp folder generations
use File::Basename;
use File::Temp;	
use HTML::TokeParser;
use LWP::Simple;
#Use CGI
use CGI;
use CGI::Carp qw(fatalsToBrowser);    # Remove for production use

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
#Export nothing by default to prevent possible method overriding.
@EXPORT      = ();
#Define the functions that are exported on demand
@EXPORT_OK   = qw(saveTempFile 
	readFusionFile 
	openFusionConnection 
	openEnsemblConnection
	getBreakpointDataWhereSampleID 
	getFeaturesFromEnsembl 
	closeFusionConnectionWithCommit
	getSampleInfoWhereFusionID
	createFusionFolder
	getEnsemblSlices
	getSampleNames
	);
#Define functions that are exported using fusionModule qw( :Tag)
%EXPORT_TAGS = ();

##Global variables##
#If there is a break in a gene, define the max. position in which a promotor break is detected
our $promotorBreakMinBP = 1000;





######################
###Fusion Functions###
######################

########################################
###Make derivatives/Jbrowse functions###
########################################
#This function will create a new folder for multi-user use in which the deravitive data is stored, also add all the data by using calling the methods and creating the deravitive and its features
sub createFusionFolder{
	my ($fusID, $searchUp, $searchDown) = @_;
	
	#Get data based on fusID
	my $dbh = &openFusionConnection;
	my $sampleInfo = &getSampleInfoWhereFusionID($dbh, $fusID);
	my $fusionInfo = &getFusionPointsWhereFusionID($dbh, $fusID);
	my $registry = &openEnsemblConnection();
	
	&closeFusionConnection($dbh);
	
	#Simplify the variable names
	my ($chr1, $chr2, $chr1_Start, $chr1_End, $chr2_Start, $chr2_End, $species, $breakOrientation, $supportFilesSample, $supportFilesFusion, $sampleID) =
	($fusionInfo->{$fusID}->{'break_chr1'}, $fusionInfo->{$fusID}->{'break_chr2'}, $fusionInfo->{$fusID}->{'break_chr1_start'}, $fusionInfo->{$fusID}->{'break_chr1_end'}, $fusionInfo->{$fusID}->{'break_chr2_start'}, $fusionInfo->{$fusID}->{'break_chr2_end'},
		$sampleInfo->{$fusID}->{'sample_species'}, $fusionInfo->{$fusID}->{'break_orientation'}, $sampleInfo->{$fusID}->{'supportFilesSample'}, $sampleInfo->{$fusID}->{'supportFilesFusion'}, $sampleInfo->{$fusID}->{'sample_ID'});
	
	#Creates a new folder folder to store the data
	my $tmpFolder = File::Temp->newdir(
        DIR      => "../data/$species/fusionpoints",
        UNLINK => 0,
        CLEANUP => 0
        );
    	my $fusionFolderName = $tmpFolder->dirname;
    	
    	#Folder where the deravitive tracks are stored
    	system("mkdir $fusionFolderName/der");
    	system("mkdir $fusionFolderName/der/tracks");
    	system("mkdir $fusionFolderName/der/seq");

    	#Folder where the reference tracks are stored
    	system("mkdir $fusionFolderName/ref");
    	system("mkdir $fusionFolderName/ref/tracks");
    	system("mkdir $fusionFolderName/ref/seq");
    	
    	#Folder where the raw/unformatted datafiles are stores before being formatted.
    	system("mkdir $fusionFolderName/raw");
    	
    	#Symlink all the 'real' chromosomes to the fusionpoint instead of simply copying them to save space
    	opendir(D, "../data/$species/ref/seq") || die "Can't open dir: $!\n";
	while (my $f = readdir(D)) {
		system("ln -s /var/www/html/FusionAnnotator/data/$species/ref/seq/$f $fusionFolderName/ref/seq") if !($f =~ /json/);
	}
	closedir(D);
	
	#Symlink all the 'real' tracks to the fusionpoint instead of simply copying them to save space
    	opendir(D, "../data/$species/ref/tracks") || die "Can't open dir: $!\n";
	while (my $f = readdir(D)) {
		system("ln -s /var/www/html/FusionAnnotator/data/$species/ref/tracks/$f $fusionFolderName/ref/tracks");
	}
	closedir(D);
	
	#Make an copy of refseqs.json for this fusionpoint folder
    	system("cp ../data/$species/ref/seq/refSeqs.json $fusionFolderName/ref/seq");
	
	#Make true derivative chromosome, get the name of this derivative chromosome
	my $derName = &makeDerChromosomeTrue($registry, $chr1, $chr2, $chr1_Start, $chr1_End, $chr2_Start, $chr2_End, $species, $breakOrientation, $fusionFolderName);
	
	#Make derivate features of genes, transcripts, domains, and user defined support files.
	#Based on search distances defined by the user and break orientation, also add the features as tracks to Jbrowse
	my $breakPos1DerCoord = &makeDerFeatures($registry, $chr1, $chr2, $chr1_Start, $chr1_End, $chr2_Start, $chr2_End, $derName, $searchUp, $searchDown, $species, $breakOrientation, $supportFilesSample, $supportFilesFusion, $sampleID, $fusID, $fusionFolderName);
	
	#Delete the raw folder to save space since the data has already served its purpose.
	#system("rm -r $fusionFolderName/raw/");
	
	#Remove the .htaccess file generated by Jbrowse
	system("rm $fusionFolderName/.htaccess");
	system("rm $fusionFolderName/der/.htaccess");
	system("rm $fusionFolderName/ref/.htaccess");
	
	#Delete fusionpoints older than 1 day
	system("find /var/www/html/FusionAnnotator/data/$species/fusionpoints/* -mtime +1 -exec rm -rf {} \;");
	
	#Return the correct breakposition, based on orientation
	my ($o1, $o2);
	if ($breakOrientation =~ /^(H|T)(H|T)/i) {
		  $o1 = uc($1);
		  $o2 = uc($2);
	}
	
	#If orientation is TH
	my ($breakChrA,$breakChrB) = ($chr1_End, $chr2_Start);
	$breakChrA = $chr1_Start if $o1 eq "H";
	$breakChrB = $chr2_End if $o2 eq "T";
	
	my $confirmedBreakpoint = "Not-Confirmed";
	$confirmedBreakpoint = "Confirmed" if $chr1_Start eq $chr1_End;
	#Return fusion info
	return ($chr1, $chr2, $breakChrA, $breakChrB, $breakPos1DerCoord, $derName, $species, $breakOrientation, $fusionFolderName, $confirmedBreakpoint);
}

#This functions will add all the support tracks for a sample to the JBrowser
sub addSupportTracks{
	my ($derName, $species, $chr1, $chr2, $breakPos1, $breakPos2, $breakPos1DerCoord, $breakOrientation, $searchDown, $searchUp, $supportFilesSample, $supportFilesFusion, $fusionFolderName) = @_;
		
	#Make a derivative of the support files (.WIG/.GFF/.BED)
	foreach my $supportFile (@{ $supportFilesFusion }){
		makeBAMTrack($supportFile->[0],$derName, $species, $chr1, $chr2, $breakPos1, $breakPos2, $breakPos1DerCoord , $breakOrientation, $searchDown, $searchUp, $fusionFolderName) if $supportFile->[0] =~ /bam$/i;
		makeWIG_BED_GFF_Track($supportFile->[0],$derName, $species, $chr1, $chr2, $breakPos1, $breakPos2, $breakPos1DerCoord, $breakOrientation, $searchDown, $searchUp, $fusionFolderName) if $supportFile->[0] =~ /(wig$|bed$|gff$)/i;
	}
	
	foreach my $supportFile (@{ $supportFilesSample }){
		makeBAMTrack($supportFile->[0],$derName, $species, $chr1, $chr2, $breakPos1, $breakPos2, $breakPos1DerCoord , $breakOrientation, $searchDown, $searchUp, $fusionFolderName) if $supportFile->[0] =~ /bam$/i;
		makeWIG_BED_GFF_Track($supportFile->[0],$derName, $species, $chr1, $chr2, $breakPos1, $breakPos2, $breakPos1DerCoord, $breakOrientation, $searchDown, $searchUp, $fusionFolderName) if $supportFile->[0] =~ /(wig$|bed$|gff$)/i;
	}
}

#Makes a derivative of a WIG/BED or GFF support file and then adds both the original as the derivative as a track to JBrowser
sub makeWIG_BED_GFF_Track{
	my ($fileLink, $derName, $species, $chr1, $chr2, $breakPos1, $breakPos2, $breakPos1DerCoord, $breakOrientation, $searchDown, $searchUp, $fusionFolderName) = @_;
	($searchDown,$searchUp) = (100e3, 100e3);
	my ($fileName,$dir,$ext) = fileparse($fileLink, qr/\.[^.]*/);
	my $DerFileHandler;
	#Get the Break orientation (Determines the position the breakpoint is on)
	my ($o1, $o2);
	if ($breakOrientation =~ /^(H|T)(H|T)/i) {
	  $o1 = uc($1);
	  $o2 = uc($2);
	}
	
	#Open the correct file for editing
	my $fileText = get($fileLink) || die("Could'nt get the text of $fileLink");
	open(my $fileHandler, ">", "$fusionFolderName/raw/refsupport.$ext") || die "Could not make support File: $fileLink";
	print $fileHandler $fileText;
	close ($fileHandler);
	
	open(my $fileHandler, "<", "$fusionFolderName/raw/refsupport.$ext") || die "Could not read reference support File: $fileLink";


	open($DerFileHandler, ">", "$fusionFolderName/raw/Der_$fileName.wig") || die "Could not make derivative WIG file of $fileLink" if $ext =~ /wig$/i;
	open($DerFileHandler, ">", "$fusionFolderName/raw/Der_$fileName.bed") || die "Could not make derivative BED file of $fileLink" if $ext =~ /bed$/i;
	open($DerFileHandler, ">", "$fusionFolderName/raw/Der_$fileName.gff") || die "Could not make derivative GFF file of $fileLink" if $ext =~ /gff$/i;
	if($ext =~ /bed$/i or $ext =~ /wig$/){
		while(<$fileHandler>){
			my @line = split("\t", $_);
			#Get rid of chr-prefix
			$line[0] =~ s/chr//;
			#Print the line without the chr prefix for the reference tracks
			print $DerFileHandler (join("\t", @line)); 
			
			#Make deravitive tracks
			#If chromosome matches first breakchromosome
			if ($line[0] == $chr1){
				#Get only the features below or or above the breakpoint on chrA
				#If start position lower than break 1 and within search distance (T orientation)
				if($line[1] <= $breakPos1 and $line[1] >= ($breakPos1-$searchDown) and $o1 eq "T"){
					$line[0] = $derName;
					print $DerFileHandler (join("\t", @line)); 
				}
	
				if($line[1] >= $breakPos1 and $line[1] <= ($breakPos1 + $searchUp) and $o1 eq "H"){
					$line[0] = $derName;
					#Transpose the begin-end positions and flip the strand if GFF or BED
					($line[1], $line[2]) = &getTransPosition(1, $line[1], $line[2], $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
					#Flip the strand
					if($ext =~ /.bed/i){
						$line[5] = "-" if $line[5] eq "+";
						$line[5] = "+" if $line[5] eq "-";
					}
					print $DerFileHandler (join("\t", @line)); 	
				}
			}
			#If chromosome matches second breakchromosome
			if ($line[0] == $chr2){
				#Get only the features below or or above the breakpoint on chrB
				#If start position higer than break 2 and within search distance (H orientation)
				if($line[1] >= $breakPos2 and $line[1] <= ($breakPos2+$searchUp) and $o2 eq "H"){
					$line[0] = $derName;
					#Transpose the begin-end positions and flip the strand if GFF or BED
					($line[1], $line[2]) = &getTransPosition(2, $line[1], $line[2], $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
					print $DerFileHandler (join("\t", @line)); 
				}
	
				if($line[1] <= $breakPos2 and $line[1] >= ($breakPos2 - $searchDown) and $o2 eq "T"){
					$line[0] = $derName;
					#Transpose the begin-end positions and flip the strand if GFF or BED
					($line[1], $line[2]) = &getTransPosition(2, $line[1], $line[2], $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
					#Flip the strand
					$line[5] = "-" if $line[5] eq "+";
					$line[5] = "+" if $line[5] eq "-";
					print $DerFileHandler (join("\t", @line)); 
				}
			}
		}
	}#End of .BED and .WIG processing
	if($ext =~ /gff$/i){
		while(<$fileHandler>){
			my @line = split("\t", $_);
			#Get rid of chr-prefix
			$line[0] =~ s/chr//;
			#Print the line without the chr prefix for the reference tracks
			print $DerFileHandler (join("\t", @line)); 
			
			#Make deravitive tracks
			#If chromosome matches first breakchromosome
			if ($line[0] == $chr1){
				#Get only the features below or or above the breakpoint on chrA
				#If start position lower than break 1 and within search distance (T orientation)
				if($line[3] <= $breakPos1 and $line[3] >= ($breakPos1-$searchDown) and $o1 eq "T"){
					$line[0] = $derName;
					print $DerFileHandler (join("\t", @line)); 
				}
	
				if($line[3] >= $breakPos1 and $line[3] <= ($breakPos1 + $searchUp) and $o1 eq "H"){
					$line[0] = $derName;
					#Transpose the begin-end positions and flip the strand if GFF or BED
					($line[3], $line[4]) = &getTransPosition(1, $line[3], $line[4], $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
					#Flip the strand
					$line[6] = "-" if $line[6] eq "+";
					$line[6] = "+" if $line[6] eq "-";
					print $DerFileHandler (join("\t", @line)); 	
				}
			}
			#If chromosome matches second breakchromosome
			if ($line[0] == $chr2){
				#Get only the features below or or above the breakpoint on chrB
				#If start position higer than break 2 and within search distance (H orientation)
				if($line[3] >= $breakPos2 and $line[3] <= ($breakPos2+$searchUp) and $o2 eq "H"){
					$line[0] = $derName;
					#Transpose the begin-end positions and flip the strand if GFF or BED
					($line[3], $line[4]) = &getTransPosition(2, $line[1], $line[2], $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
					print $DerFileHandler (join("\t", @line)); 
				}
	
				if($line[3] <= $breakPos2 and $line[3] >= ($breakPos2 - $searchDown) and $o2 eq "T"){
					$line[0] = $derName;
					#Transpose the begin-end positions and flip the strand if GFF or BED
					($line[3], $line[4]) = &getTransPosition(2, $line[3], $line[4], $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
					#Flip the strand
					$line[6] = "-" if $line[6] eq "+";
					$line[6] = "+" if $line[6] eq "-";
					print $DerFileHandler (join("\t", @line)); 
				}
			}
		}
	}#End of .GFF processing	
	
	close($fileHandler);
	close($DerFileHandler);
	
	#Add the tracks to the JBrowser
	#WIG
	system("perl ../genomeBrowser/bin/wig-to-json.pl --wig $fusionFolderName/raw/Der_$fileName.wig --tracklabel '$fileName' --key '$fileName' --out $fusionFolderName/der/") if $ext =~ /wig$/i;
	system("perl ../genomeBrowser/bin/wig-to-json.pl --wig $fusionFolderName/raw/Der_$fileName.wig --tracklabel '$fileName' --key '$fileName' --out $fusionFolderName/ref/") if $ext =~ /wig$/i;

	#GFF
	system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff $fusionFolderName/raw/Der_$fileName.gff --tracklabel '$fileName' --key '$fileName' --out $fusionFolderName/der/") if $ext =~ /gff$/i;
	system("perl ../genomeBrowser/bin/wig-to-json.pl --gff $fusionFolderName/raw/Der_$fileName.gff --tracklabel '$fileName' --key '$fileName' --out $fusionFolderName/ref/") if $ext =~ /gff$/i;

	#BED
	system("perl ../genomeBrowser/bin/flatfile-to-json.pl --bed $fusionFolderName/raw/Der_$fileName.bed --tracklabel '$fileName' --key '$fileName' --out $fusionFolderName/der/") if $ext =~ /bed$/i;
	system("perl ../genomeBrowser/bin/wig-to-json.pl --bed $fusionFolderName/raw/Der_$fileName.bed --tracklabel '$fileName' --key '$fileName' --out $fusionFolderName/ref/") if $ext =~ /bed$/i;

}

#Make a deravative of the BAM by converting regions of the BAM to SAM, changing the coordinates and chromosome and back to BAM.
sub makeBAMTrack{
	my ($BAMLink, $derName, $species, $chr1, $chr2, $breakPos1, $breakPos2, $breakPos1DerCoord,$breakOrientation, $searchDown, $searchUp, $fusionFolderName) = @_;
	my ($o1, $o2);
	if ($breakOrientation =~ /^(H|T)(H|T)/i) {
		$o1 = uc($1);
		$o2 = uc($2);
	}
	
	my ($fileName,$dir,$ext) = fileparse($BAMLink, qr/\.[^.]*/);	

	#Replace the link with the wgs02 /data/ path
	my $BAMLinkLocal = $BAMLink;
	$BAMLinkLocal =~ s/http:\/\/genetics.genomicscenter.nl\/BAMS\//\/data\//;
	
	#Get the BAM-reads for the breaks, put them in a derSAM file (Overwrite older file with the same derName).
	system( "/usr/local/samtools/samtools view -h $BAMLinkLocal $chr1:".($breakPos1-$searchDown)."-$breakPos1 $chr2:$breakPos2-".($breakPos2+$searchDown)." > $fusionFolderName/raw/temp.sam") if uc($breakOrientation) eq "TH";
	system( "/usr/local/samtools/samtools view -h $BAMLinkLocal $chr1:$breakPos1-".($breakPos1+$searchDown)." $chr2:$breakPos2-".($breakPos2+$searchDown)." > $fusionFolderName/raw/temp.sam") if uc($breakOrientation) eq "HH";
	system( "/usr/local/samtools/samtools view -h $BAMLinkLocal $chr1:$breakPos1-".($breakPos1+$searchDown)." $chr2:".($breakPos2-$searchDown)."-$breakPos2 > $fusionFolderName/raw/temp.sam") if uc($breakOrientation) eq "HT";
	system( "/usr/local/samtools/samtools view -h $BAMLinkLocal $chr1:".($breakPos1-$searchDown)."-$breakPos1 $chr2:".($breakPos2-$searchDown)."-$breakPos2 > $fusionFolderName/raw/temp.sam") if uc($breakOrientation) eq "TT";
		
	#Open the SAM file for editing
	open(my $DerTempFileHandler, "<", "$fusionFolderName/raw/temp.sam") || die $!;
	open(my $DerSAMFileHandler, ">", "$fusionFolderName/raw/$fileName.sam") || die "Could not write SAM File: $fusionFolderName/raw/$fileName.sam";

	#Read through SAM File while making a new SAM file
	print $DerSAMFileHandler "\@SQ\tSN:$derName\tLN:9999999999999999\n";
	
	#Transpose the BAM reads
	while(<$DerTempFileHandler>){
		if(!($_ =~ /@/)){
			my ($tpos1, $tpos2);
			#Split SAM per tab
			my @lineSplit = split("\t");
			#Delete chr-prefix
			$lineSplit[2] =~ s/^chr//;
			#Transpose the element correctly
			#If start position lower than break 1 and == $chr1(T Orientation)
			if($lineSplit[2] eq $chr1 and $lineSplit[3] <= $breakPos1 and $o1 eq "T"){
				($tpos1, $tpos1) = &getTransPosition(1, $lineSplit[3], $lineSplit[3], $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
			}
			#If start position higher than break 1 and == $chr1(H Orientation)
			if($lineSplit[2] eq $chr1 and $lineSplit[3] >= $breakPos1 and $o1 eq "H"){
				($tpos1, $tpos1) = &getTransPosition(1, $lineSplit[3], $lineSplit[3], $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
			}
			
			#If start position higher than break 2 and == $chr2(T Orientation)
			if($lineSplit[2] eq $chr2 and $lineSplit[3] <= $breakPos2 and $o2 eq "T"){
				($tpos2, $tpos2) = &getTransPosition(2, $lineSplit[3], $lineSplit[3], $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
			}
			
			#If start position lower than break 2 and == $chr2(H Orientation)
			if($lineSplit[2] eq $chr2 and $lineSplit[3] >= $breakPos2 and $o2 eq "H"){
				($tpos2, $tpos2) = &getTransPosition(2, $lineSplit[3], $lineSplit[3], $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
			}
			
			#Change chromosome to derivative chromosome
			$lineSplit[2] = $derName;
			#Print the new feature to the SAM File
			if($tpos1){
				$lineSplit[3] = $tpos1;
				print $DerSAMFileHandler join("\t",@lineSplit);
			}
			if($tpos2){
				$lineSplit[3] = $tpos2;
				print $DerSAMFileHandler join("\t",@lineSplit);
			}
		}
	}
	#Close the modified SAM file.
	close($DerSAMFileHandler);
	close($DerTempFileHandler);
	
	#Convert to .SAM to .BAM
	system("/usr/local/samtools/samtools view -Sb  $fusionFolderName/raw/$fileName.sam  >  $fusionFolderName/raw/unsorted_$fileName.bam");
	#Sort the BAM
	system("/usr/local/samtools/samtools sort $fusionFolderName/raw/unsorted_$fileName.bam $fusionFolderName/der/".$fileName."_Der");
	#Make BAI index.
	system("/usr/local/samtools/samtools index $fusionFolderName/der/".$fileName."_Der.bam");

	
	#Add the BAM files to the Jbrowse as tracks
	&addBAMTracks($BAMLink, $fileName, $derName, $species, $fusionFolderName);
}

#Add the BAM files as tracks to the JBrowse by altering the tracklist.json and directly pointing the tracks to the BAM files
#Do this by substituting the existing link to the BAM file, not needed for the derivative of the BAM
sub addBAMTracks{
	my ($BAMLink, $fileName, $derName, $species, $fusionFolderName) = @_;
	
	#Open trackList.json for editing
	open(my $trackListHandlerDer, "<", "$fusionFolderName/der/trackList.json") || die "Could not read deravitive tracklist.json";

	open(my $tempHandler, ">", "$fusionFolderName/tempList.json") || die "Could not write to tempList.json";
	
	#Append the new .BAM tracks to the temp file for the deravitive view
	while (<$trackListHandlerDer>){
		print $tempHandler $_ if !($_ =~ /]/i);
		print $tempHandler '
		,
		      {
		    "storeClass" : "JBrowse/Store/SeqFeature/BAM",
		    "urlTemplate" : "'.$fileName.'_Der.bam",
		    "max_score" : 40,
		    "label" : "der_'.$fileName.' Coverage",
		    "type" : "JBrowse/View/Track/FeatureCoverage",
		    "min_score" : 0,
		    "key" : "der_'.$fileName.' Coverage"
		      },
		      {
		    "storeClass" : "JBrowse/Store/SeqFeature/BAM",
		    "urlTemplate" : "'.$fileName.'_Der.bam",
		    "style" : {
			"className" : "alignment",
			"arrowheadClass" : "arrowhead",
			"labelScale" : 100,
			"histScale" : 100
		    },
		    "label" : "der_'.$fileName.'",
		    "type" : "JBrowse/View/Track/Alignments",
		    "key" : "der_'.$fileName.'"
		      }
		      ],
		      ' if $_ =~ /]/;
	}
	
	close($trackListHandlerDer);
	close($tempHandler);
	
	#Write the temp content to the actual trackList,json
	open($tempHandler, "<", "$fusionFolderName/tempList.json") || die "Could not read tempList.json";
	open($trackListHandlerDer, ">", "$fusionFolderName/der/trackList.json") || die "Could not write to tracklist.json";
	
	while (<$tempHandler>){
		print $trackListHandlerDer $_;
	}
	
	close($trackListHandlerDer);
	close($tempHandler);
	
	#Append the new .BAM tracks to the ref file for the deravitive view
	open(my $tempHandler, ">", "$fusionFolderName/tempList.json") || die "Could not write to tempList.json";
	open(my $trackListHandlerRef, "<", "$fusionFolderName/ref/trackList.json") || die "Could not read reference tracklist.json";

	while (<$trackListHandlerRef>){
		print $tempHandler $_ if !($_ =~ /]/i);
		print $tempHandler '
		,
		      {
		    "storeClass" : "JBrowse/Store/SeqFeature/BAM",
		    "urlTemplate" : "'.$BAMLink.'",
		    "max_score" : 40,
		    "label" : "'.$fileName.' Coverage",
		    "type" : "JBrowse/View/Track/FeatureCoverage",
		    "min_score" : 0,
		    "key" : "'.$fileName.' Coverage"
		      },
		      {
		    "storeClass" : "JBrowse/Store/SeqFeature/BAM",
		    "urlTemplate" : "'.$BAMLink.'",
		    "style" : {
			"className" : "alignment",
			"arrowheadClass" : "arrowhead",
			"labelScale" : 100,
			"histScale" : 80
		    },
		    "label" : "'.$fileName.'",
		    "type" : "JBrowse/View/Track/Alignments",
		    "key" : "'.$fileName.'"
		      }
		      ],
		      ' if $_ =~ /]/;
	}
	
	close($tempHandler);
	close($trackListHandlerRef);
	
	#Write the temp content to the actual trackList,json
	open($tempHandler, "<", "$fusionFolderName/tempList.json") || die "Could not read tempList.json";
	open(my $trackListHandlerRef, ">", "$fusionFolderName/ref/trackList.json") || die "Could not read reference tracklist.json";

	
	while (<$tempHandler>){
		print $trackListHandlerRef $_;
	}
	
	close($trackListHandlerRef);
	close($tempHandler);
}

#Makes an gff track that can be used as a visualization aid for the current breakpoint and the surrounding breakpoints
#Colour each breaktype differently to differentiate between them more easily, clicking this feature opens a new genome browser with that breakpoint
sub addBreakpointVisualization{
	my ($derName, $species, $chr1, $chr2, $breakPos1, $breakPos2, $breakPos1DerCoord, $breakOrientation, $sampleID, $fusID, $fusionFolderName) = @_;
	my ($breakPosition, $breakPosition2, $breakPositionRef, $break2PositionRef);

	#Get the breakpoints in the vicinity of the current breakpoint
	my $dbh = openFusionConnection();
	#Open the file where the breakpoint will be written to
	open( my $breakpointGFF, ">", "$fusionFolderName/raw/breakpointsInVicinity.gff") || die "Could not write/create to breakpointsInVicinity.bed in $fusionFolderName/raw/";

	#Do not show the other breakpoints if the breakpoint is originating from the StrucVarDB
	my $breakpointsRef = $dbh->selectall_hashref("SELECT DISTINCT fusion_ID, break_chr1, break_chr1_start, break_chr1_end, break_chr2, break_chr2_start, break_chr2_end, break_orientation, T_breaktype_breaktype_ID
		FROM T_fusionpoint
		 INNER JOIN T_samples ON T_fusionpoint.T_samples_sample_ID=T_samples.sample_ID
		 WHERE T_samples_sample_ID = $sampleID", "fusion_ID");
	foreach my $fusKey (keys %$breakpointsRef){
		my $breaktype;

		#Get breakpoint type
		$breaktype = "Insertion" if $breakpointsRef->{$fusKey}->{'T_breaktype_breaktype_ID'} == 1;
		$breaktype = "Deletion" if $breakpointsRef->{$fusKey}->{'T_breaktype_breaktype_ID'} == 2;
		$breaktype = "Anti" if $breakpointsRef->{$fusKey}->{'T_breaktype_breaktype_ID'} == 3;
		$breaktype = "Inversion" if $breakpointsRef->{$fusKey}->{'T_breaktype_breaktype_ID'} == 4;
		$breaktype = "Remote" if $breakpointsRef->{$fusKey}->{'T_breaktype_breaktype_ID'} == 5;
		$breaktype = "Translocation" if $breakpointsRef->{$fusKey}->{'T_breaktype_breaktype_ID'} == 6;
		$breaktype = "Duplication" if $breakpointsRef->{$fusKey}->{'T_breaktype_breaktype_ID'} == 7;
                $breaktype = "Evertion" if $breakpointsRef->{$fusKey}->{'T_breaktype_breaktype_ID'} == 8;
		$breaktype = "Other" if $breakpointsRef->{$fusKey}->{'T_breaktype_breaktype_ID'} == 9;
		#Mark the selected breakpoint as selected for color purposes
		$breaktype = "Selected" if $fusID == $breakpointsRef->{$fusKey}->{'fusion_ID'};
		
		#Get the Break orientation (Determines the position the breakpoint is on)
		my ($o1, $o2);
		if ($breakpointsRef->{$fusKey}->{'break_orientation'} =~ /^(H|T)(H|T)/i) {
		  $o1 = uc($1);
		  $o2 = uc($2);
		}
		$breakPositionRef = $breakpointsRef->{$fusKey}->{'break_chr1_end'} if $o1 eq "T";
		$breakPositionRef = $breakpointsRef->{$fusKey}->{'break_chr1_start'} if $o1 eq "H";
		$break2PositionRef = $breakpointsRef->{$fusKey}->{'break_chr2_end'} if $o2 eq "T";
		$break2PositionRef = $breakpointsRef->{$fusKey}->{'break_chr2_start'} if $o2 eq "H";
		
		#Print the breakpoint to the GFF file for the reference tracks
		print $breakpointGFF $breakpointsRef->{$fusKey}->{'break_chr1'}."\t";
		print $breakpointGFF "FusionDB\t$breaktype\t";
		print $breakpointGFF $breakPositionRef."\t";
		print $breakpointGFF $breakPositionRef."\t.\t.\t+\t";
		print $breakpointGFF "Name=Break-$breaktype;URL=<a href = 'http://localhost/FusionAnnotator/cgi-bin/preVisual.pl?fusID=$fusKey' target='blank'>Visualize breakpoint</a>;\n";
		
		print $breakpointGFF $breakpointsRef->{$fusKey}->{'break_chr2'}."\t";
		print $breakpointGFF "FusionDB\t$breaktype\t";
		print $breakpointGFF $break2PositionRef."\t";
		print $breakpointGFF $break2PositionRef."\t.\t.\t+\t";
		print $breakpointGFF "Name=Break-$breaktype;URL=<a href = 'http://localhost/FusionAnnotator/cgi-bin/preVisual.pl?fusID=$fusKey' target='blank'>Visualize breakpoint</a>;\n";
		
		#Print the selected breakpoint
		if($breaktype eq "Selected"){
			print $breakpointGFF $derName."\t";
			print $breakpointGFF "FusionDB\t$breaktype\t";
			print $breakpointGFF $breakPosition."\t";
			print $breakpointGFF $breakPosition."\t.\t.\t+\t";
			print $breakpointGFF "Name=Break-Selected;URL=<a href = 'http://localhost/FusionAnnotator/cgi-bin/preVisual.pl?fusID=$fusKey' target='blank'>Visualize breakpoint</a>;\n";
		}
		#Print the deravative breakpoint to the GFF file
		#Transpose the positions accordingly
	
		#Are breakpoints from same chromosome as selected fusionpoint
		if($breakpointsRef->{$fusKey}->{'break_chr1'} eq $chr1){
			($breakPosition, $breakPosition2) = &getTransPosition(1, $breakPositionRef, $breakPositionRef, $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
			#Only print if it is not over the breakpoint
			if( ($breakPosition < $breakPos1DerCoord and $breaktype ne "Selected") or ($breakPosition == $breakPos1DerCoord and $breaktype eq "Selected")){
				print $breakpointGFF $derName."\t";
				print $breakpointGFF "FusionDB\t$breaktype\t";
				print $breakpointGFF $breakPosition."\t";
				print $breakpointGFF $breakPosition."\t.\t.\t+\t";
				print $breakpointGFF "Name=Break-$breaktype;URL=<a href = 'http://localhost/FusionAnnotator/cgi-bin/preVisual.pl?fusID=$fusKey' target='blank'>Visualize breakpoint</a>;\n";
			}
		}
			
		if($breakpointsRef->{$fusKey}->{'break_chr2'} eq $chr1){
			($breakPosition, $breakPosition2) = &getTransPosition(1, $break2PositionRef, $break2PositionRef, $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
			#Only print if it is not over the breakpoint
			if( ($breakPosition < $breakPos1DerCoord and $breaktype ne "Selected") or ($breakPosition == $breakPos1DerCoord and $breaktype eq "Selected")){
				print $breakpointGFF $derName."\t";
				print $breakpointGFF "FusionDB\t$breaktype\t";
				print $breakpointGFF $breakPosition."\t";
				print $breakpointGFF $breakPosition."\t.\t.\t+\t";
				print $breakpointGFF "Name=Break-$breaktype;URL=<a href = 'http://localhost/FusionAnnotator/cgi-bin/preVisual.pl?fusID=$fusKey' target='blank'>Visualize breakpoint</a>;\n";
			}
		}
			
		#Are breakpoints from same chromosome as selected fusionpoint
		if($breakpointsRef->{$fusKey}->{'break_chr1'} eq $chr2){
			($breakPosition, $breakPosition2) = &getTransPosition(2, $breakPositionRef, $breakPositionRef, $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
			#Only print if it is not over the breakpoint
			if($breakPosition > ($breakPos1DerCoord+1) and $breaktype ne "Selected"){
				print $breakpointGFF $derName."\t";
				print $breakpointGFF "FusionDB\t$breaktype\t";
				print $breakpointGFF $breakPosition."\t";
				print $breakpointGFF $breakPosition."\t.\t.\t+\t";
				print $breakpointGFF "Name=BreakB-$breaktype;URL=<a href = 'http://localhost/FusionAnnotator/cgi-bin/preVisual.pl?fusID=$fusKey' target='blank'>Visualize breakpoint</a>;\n";
			}
		}
			
		if($breakpointsRef->{$fusKey}->{'break_chr2'} eq $chr2){
			($breakPosition, $breakPosition2) = &getTransPosition(2, $break2PositionRef, $break2PositionRef, $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
			#Only print if it is not over the breakpoint
			if($breakPosition > ($breakPos1DerCoord+1) and $breaktype ne "Selected"){
				print $breakpointGFF $derName."\t";
				print $breakpointGFF "FusionDB\t$breaktype\t";
				print $breakpointGFF $breakPosition."\t";
				print $breakpointGFF $breakPosition."\t.\t.\t+\t";
				print $breakpointGFF "Name=BreakB-$breaktype;URL=<a href = 'http://localhost/FusionAnnotator/cgi-bin/preVisual.pl?fusID=$fusKey' target='blank'>Visualize breakpoint</a>;\n";
			}
		}
	}
	
	#Close the GFF file where the breakpoints have been written to.
	close($breakpointGFF);
	#Close connection to FusionDB
	&closeFusionConnection($dbh);
	
	#Add the GFF as a track to Jbrowse, use subfeatureClasses to give the exact breakposition a different color to "easily" recognize this site
	system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff $fusionFolderName/raw/breakpointsInVicinity.gff --tracklabel 'BreakpointsInVicinity' --key 'BreakpointsInVicinity' --subfeatureClasses '{\"Breakpoint\" : \"feature2\"}' --out $fusionFolderName/der/");
	system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff $fusionFolderName/raw/breakpointsInVicinity.gff --tracklabel 'BreakpointsInVicinity' --key 'BreakpointsInVicinity' --subfeatureClasses '{\"Breakpoint\" : \"feature2\"}' --out $fusionFolderName/ref/");
}


#Get the features from the Ensembl DB and put them in a file in the correct format (Faster than reading out of a big predefined feature file)
#Only show features close to the breakpoint (User Defined search distances also used to search the ensembl DB for features)
sub makeDerFeatures{
	my ($registry, $chr1, $chr2, $chr1_Start, $chr1_End, $chr2_Start, $chr2_End, $derName, $searchUp, $searchDown, $species, $breakOrientation, $supportFilesSample, $supportFilesFusion, $sampleID, $fusID, $fusionFolderName) = @_;
	
	($searchUp, $searchDown) = (200e3, 200e3);
	
	#Get slices
	my $sa = $registry->get_adaptor($species, 'Core', 'Slice');
	my $transcriptAdaptor = $registry->get_adaptor( $species, 'Core', 'Transcript' );


	#Get the Break orientation (Depends what strand to write a feature to)
	my ($o1, $o2);
	if ($breakOrientation =~ /^(H|T)(H|T)/i) {
		  $o1 = uc($1);
		  $o2 = uc($2);
	}	
	
	#Get the length of the chromosome if the orientation needs this (H for first and T for second break)
	my ($lengthChrA, $lengthChrB, $slice) = (0,0, undef);
	$slice = $sa->fetch_by_region( 'chromosome', $chr1 ) if $o1 eq "H";
	$slice = $sa->fetch_by_region( 'chromosome', $chr2 ) if $o1 eq "T";
	$lengthChrA = $slice->length() if $o1 eq "H";
	$lengthChrB = $slice->length() if $o2 eq "T";
	
	#Define the values for searching the features in ensembl and the exact breakpoint, based on orientation
	#breakPosxDerCoord is used to get absolute 5' -> 3' position on the deravative track of the breakpoints
	#If orientation is first orientation is T, second orientation is H
	my ($search1Start, $search1Stop, $breakPos1, $breakPos1DerCoord, $search2Start, $search2Stop, $breakPos2, $breakPos2DerCoord) = (($chr1_Start - $searchDown), $chr1_End, $chr1_End, $chr1_End, $chr2_Start, ($chr2_Start + $searchUp), $chr2_Start, $chr2_Start);
	
	#If first orientation is H
	($search1Start, $search1Stop, $breakPos1, $breakPos1DerCoord) = ($chr1_Start , ($chr1_Start + $searchUp), $chr1_Start, (($lengthChrA - $chr1_Start)+1)) if $o1 eq "H";
		
	#If second orientation is T
	($search2Start, $search2Stop, $breakPos2, $breakPos2DerCoord) = (($chr2_End - $searchUp), $chr2_End, $chr2_End, (($lengthChrB-$chr2_End)+1)) if $o2 eq "T";

	#Get the slices with the ensembl features for the "full deravitives"
	my ($derSlice1, $derSlice2) = ($sa->fetch_by_region( 'chromosome', $chr1, $search1Start, $search1Stop), $sa->fetch_by_region( 'chromosome', $chr2, $search2Start, $search2Stop));
	
	#Get slices only to check fusiongenes
	my ($fusionGeneSlice1, $fusionGeneSlice2) = ($sa->fetch_by_region( 'chromosome', $chr1, $breakPos1, $breakPos1), $sa->fetch_by_region( 'chromosome', $chr2, $breakPos2, $breakPos2));
	
	#Keep track of the features and slice to prevent duplicates
	my ($count, $geneStrand, %features);

	
	#Open feature files for writing
	open(my $DerGeneFileHandler, '>', "$fusionFolderName/raw/derGeneFeatures.gff") || die "Could not write feature file for the deravitive gene features: $fusionFolderName/raw/derGeneFeatures.gff\n";
	open(my $DerTranscriptFileHandler, '>', "$fusionFolderName/raw/derTranscriptFeatures.gff") || die "Could not write feature file for the deravitive transcript features: $fusionFolderName/raw/derTranscriptFeatures.gff\n";
	open(my $DerDomainFileHandler, '>', "$fusionFolderName/raw/derDomainFeatures.gff") || die "Could not write feature file for the deravitive domain features: $fusionFolderName/raw/derDomainFeatures.gff\n";
	
	#Open der reg. ele file
	open(my $DerRegEleFileHandler, '>', "$fusionFolderName/raw/derRegEleFeatures.gff") || die "Could not write feature file for the deravitive regEle features: $fusionFolderName/raw/derRegEleFeatures.gff\n";

	#Write the genes to DeravitiveFeatures.gff, transpose the features to their correct position and/or strand based on the break-orientation

	#First check for fusiongenes, if there are fusiongenes, fuse the genes together but cutt-off the gene @ the according breakpoint
	#Since a gene is possible on the - strand and on the + strand at the same position, check all genes for possible fusion genes
	foreach my $geneOn1(@{ $fusionGeneSlice1->get_all_Genes()} ){
		foreach my $geneOn2(@{ $fusionGeneSlice2->get_all_Genes() } ){			
			#Make true fusion-gene if orientation and strand permits this
			if(($o1 ne $o2 and $geneOn1->strand() eq $geneOn2->strand()) or ($o1 eq $o2 and $geneOn1->strand() ne $geneOn2->strand())){
				#Flip the strand if needed, based on orientation
				my $strand = "+";
				$strand = "-" if ($geneOn1->strand() == "-1" and $o1 eq "T") or ($geneOn1->strand() == "1" and $o1 eq "H");
				
				#Create fusion-name by simply joining the gene names
				my $fusionGeneHolder = join("::", $geneOn1->external_name, $geneOn2->external_name);
				my $fusionGeneName = $fusionGeneHolder;
				
				#Define the correct transposition coordinates based on orientation
				#If orientation of first break is T, second break is H
				my ($fusionGeneStart, $fusionGeneEnd) = ($geneOn1->seq_region_start(), ($geneOn2->seq_region_end() - $breakPos2) + ($breakPos1DerCoord + 1));
				
				#If first orientation is H
				$fusionGeneStart = ($breakPos1DerCoord - (($breakPos1 - $geneOn1->seq_region_end))*-1) 	if $o1 eq "H";
					
				#If second orientation is T
				$fusionGeneEnd = (($breakPos1DerCoord + 1) + ($breakPos2 - $geneOn2->seq_region_start))	  if $o2 eq "T";
					
				#Write the fusion-gene to the derGeneFilehandler
				print $DerGeneFileHandler join("\t", $derName, $species."_funcgen_69_37","FusionGene", $fusionGeneStart, $fusionGeneEnd, "." ,$strand, ".", "ID=$fusionGeneName;Name=$fusionGeneName;Note=FusionGene")."\n";
				#Make subfeature of gene on breakA, start from Gene-start to Breakpos
				print $DerGeneFileHandler join("\t", $derName, $species."_funcgen_69_37","GeneBreakA", $fusionGeneStart, $breakPos1DerCoord, "." ,$strand, ".", "Parent=$fusionGeneName;Name=".$geneOn1->external_name.";URL=<a href = http://www.ensembl.org/$species/Gene/Summary?g=".$geneOn1->stable_id." target='blank'>Go to Ensembl</a>")."\n";
				#Make subfeature of gene on breakB, start from Breakpos +1 to Gene-Stop (Minus the <--break bp)
				print $DerGeneFileHandler join("\t", $derName, $species."_funcgen_69_37","GeneBreakB", ($breakPos1DerCoord+1), $fusionGeneEnd, "." ,$strand, ".", "Parent=$fusionGeneName;Name=".$geneOn2->external_name.";URL=<a href = http://www.ensembl.org/$species/Gene/Summary?g=".$geneOn1->stable_id." target='blank'>Go to Ensembl</a>")."\n";	
				
				#Get exons, transcripts and domains for both GeneBreakA and GeneBreakB, cutt-of features at breakpoint position if overlapping
				#Get exons, transcripts and domains for (fusion)Gene on BreakA
				foreach my $transcript ( @{ $geneOn1->get_all_Transcripts() }){
					my ($transcriptStart, $transcriptEnd) = &getTransPosition(1, $transcript->seq_region_start(), $transcript->seq_region_end(), $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
					#Only print transcripts with the start < breaklocation
					if($transcriptStart <= ($breakPos1DerCoord -1)){
						#Print Transcripts
						print $DerTranscriptFileHandler join("\t", $derName, $species."_funcgen_69_37", "TranscriptBreakA", $transcriptStart, $transcriptEnd, ".", $strand , ".", "Name=".$transcript->stable_id().";Gene:=".$geneOn1->external_name().";URL=<a href = http://www.ensembl.org/$species/Transcript/Summary?t=".$transcript->stable_id()." target='blank'>Go to Ensembl</a>")."\n";

						#Get all exons
						foreach my $exon ( @{ $transcript->get_all_Exons() } ){
							my ($exonStart, $exonEnd) = &getTransPosition(1, $exon->seq_region_start(), $exon->seq_region_end(), $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
							#Only print exons with the start < breaklocation
							if($exonStart < ($breakPos1DerCoord)){
								print $DerGeneFileHandler join("\t", $derName, $species."_funcgen_69_37", "Exon", $exonStart, $exonEnd, ".", $strand , ".", "Parent=$fusionGeneName")."\n";
							}
						}
						
						#Get all domains
						my $domainfeat = $transcriptAdaptor->fetch_by_stable_id($transcript->stable_id())->translation()->get_all_DomainFeatures if defined($transcript->biotype()) and $transcript->biotype() eq "protein_coding";
						while ( my $pfeature = shift @{$domainfeat} ) {
							my ($domainStart, $domainEnd) = &getTransPosition(1, $pfeature->seq_region_start(), $pfeature->seq_region_end(), $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
							if($domainEnd <= ($breakPos1DerCoord -1)){
								print $DerDomainFileHandler join("\t", $derName, $species."_funcgen_69_37", "DomainBreakA", $domainStart, $domainEnd, ".", $strand , ".", "Name=".$pfeature->idesc().";Parent:=".$transcript->stable_id().";URL=<a href = http://pfam.sanger.ac.uk/family/".$pfeature->display_id()." target='blank'>Go to PFAM</a>")."\n" if $pfeature->idesc() ne "";
							}
						}	
					}
				}
				
				#Get exons for geneBreakB
				foreach my $transcript ( @{ $geneOn2->get_all_Transcripts() }){
					my ($transcriptStart, $transcriptEnd) = &getTransPosition(2, $transcript->seq_region_start(), $transcript->seq_region_end(), $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
					#Only print transcripts with the start < breaklocation
					if($transcriptEnd >= ($breakPos1DerCoord + 2)){
						#Print Transcripts
						print $DerTranscriptFileHandler join("\t", $derName, $species."_funcgen_69_37", "TranscriptBreakB", $transcriptStart, $transcriptEnd, ".", $strand , ".", "Name=".$transcript->stable_id().";Gene:=".$geneOn2->external_name().";URL=<a href = http://www.ensembl.org/$species/Transcript/Summary?t=".$transcript->stable_id()." target='blank'>Go to Ensembl</a>")."\n";

						#Get all exons
						foreach my $exon ( @{ $transcript->get_all_Exons() } ){
							my ($exonStart, $exonEnd) = &getTransPosition(2, $exon->seq_region_start(), $exon->seq_region_end(), $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
							#Only print exons with the start > breaklocation
							if($exonEnd > ($breakPos1DerCoord + 1)){
								print $DerGeneFileHandler join("\t", $derName, $species."_funcgen_69_37", "Exon2", $exonStart, $exonEnd, ".", $strand , ".", "Parent=$fusionGeneName")."\n";
							}
						}
						#Get all domains
						my $domainfeat = $transcriptAdaptor->fetch_by_stable_id($transcript->stable_id())->translation()->get_all_DomainFeatures if defined($transcript->biotype()) and $transcript->biotype() eq "protein_coding";
						while ( my $pfeature = shift @{$domainfeat} ) {
							my ($domainStart, $domainEnd) = &getTransPosition(2, $pfeature->seq_region_start(), $pfeature->seq_region_end(), $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
							if($domainEnd >= ($breakPos1DerCoord + 2)){
								print $DerDomainFileHandler join("\t", $derName, $species."_funcgen_69_37", "DomainBreakB", $domainStart, $domainEnd, ".", $strand , ".", "Name=".$pfeature->idesc().";Parent:=".$transcript->stable_id().";URL=<a href = http://pfam.sanger.ac.uk/family/".$pfeature->display_id()." target='blank'>Go to PFAM</a>")."\n" if $pfeature->idesc() ne "";
							}
						}	
					}
				}
				#Push the genes to the feature table so they wont be made again
				$features{$geneOn1->stable_id} = undef;
				$features{$geneOn2->stable_id} = undef;
			}
		}
	}
	###End of "true" fusion-gene making###
	
	###Make the non-fusion deravitives features
	
	#Get non-fusion genes, transcripts and domains and put them in the feature file
	#Do this per deravitive side
	foreach my $slice ($derSlice1, $derSlice2){
		$count++;
		#Genes, transcripts and domains
		while ( my $gene = shift @{$slice->get_all_Genes()} ){
				#Check if it has already been made (existing in a breakpoint for example)
				if(!(exists($features{$gene->stable_id()}))){
					my $strand = "+";
					#Flip the strand if orientation supports this
					my ($geneStart, $geneEnd, $transcriptStart, $transcriptEnd, $exonStart, $exonEnd, $domainStart, $domainEnd);
					if($count == 1){
						$strand = "-" if ($gene->strand() == "-1" and $o1 eq "T") or ($gene->strand() == "1" and $o1 eq "H");
						($geneStart, $geneEnd) = &getTransPosition(1, $gene->seq_region_start(), $gene->seq_region_end(), $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
					}
					if($count == 2){
						$strand = "-" if ($gene->strand() == "-1" and $o2 eq "H") or ($gene->strand() == "1" and $o2 eq "T");
						($geneStart, $geneEnd) = &getTransPosition(2, $gene->seq_region_start(), $gene->seq_region_end(), $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
					}
				
					#Get genes
					print $DerGeneFileHandler join("\t", $derName, $species."_funcgen_69_37", "Gene", $geneStart, $geneEnd , "." ,$strand, ".", "ID=".$gene->stable_id().";Name=".$gene->external_name().";Note=".$gene->stable_id().";URL=<a href = http://www.ensembl.org/$species/Gene/Summary?g=".$gene->stable_id()." target='blank'>Go to Ensembl</a>")."\n";
					
					print $DerGeneFileHandler join("\t", $derName, $species."_funcgen_69_37", "GeneBreakA", $geneStart, $geneEnd , "." ,$strand, ".", "Name=".$gene->external_name().";Parent=".$gene->stable_id().";URL=<a href = http://www.ensembl.org/$species/Gene/Summary?g=".$gene->stable_id()." target='blank'>Go to Ensembl</a>")."\n" if $count == 1;
					print $DerGeneFileHandler join("\t", $derName, $species."_funcgen_69_37", "GeneBreakB", $geneStart, $geneEnd , ".", $strand ,".", "Name=".$gene->external_name().";Parent=".$gene->stable_id().";URL=<a href = http://www.ensembl.org/$species/Gene/Summary?g=".$gene->stable_id()." target='blank'>Go to Ensembl</a>")."\n" if $count == 2;
					
					#Get domains, transcripts and exons
					while( my $transcript = shift @{ $gene->get_all_Transcripts() }){
						($transcriptStart, $transcriptEnd) = &getTransPosition(1, $transcript->seq_region_start(), $transcript->seq_region_end(), $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord) if $count == 1;
						($transcriptStart, $transcriptEnd) = &getTransPosition(2, $transcript->seq_region_start(), $transcript->seq_region_end(), $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord) if $count == 2;

						print $DerTranscriptFileHandler join("\t", $derName, $species."_funcgen_69_37", "TranscriptBreakA", $transcriptStart, $transcriptEnd, ".", $strand , ".", "Name=".$transcript->stable_id().";Gene=".$gene->external_name().";URL=<a href = http://www.ensembl.org/$species/Transcipt/Summary?t=".$transcript->stable_id()." target='blank'>Go to Ensembl</a>")."\n" if $count == 1;
						print $DerTranscriptFileHandler join("\t", $derName, $species."_funcgen_69_37", "TranscriptBreakB", $transcriptStart, $transcriptEnd, ".", $strand , ".", "Name=".$transcript->stable_id().";Gene=".$gene->external_name().";URL=<a href = http://www.ensembl.org/$species/Transcipt/Summary?t=".$transcript->stable_id()." target='blank'>Go to Ensembl</a>")."\n" if $count == 2;
						
						#Get exons
						foreach my $exon ( @{ $transcript->get_all_Exons() } ) {
							
							($exonStart, $exonEnd) = &getTransPosition(1, $exon->seq_region_start(), $exon->seq_region_end(), $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord) if $count == 1;
							($exonStart, $exonEnd) = &getTransPosition(2, $exon->seq_region_start(), $exon->seq_region_end(), $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord) if $count == 2;
								
							print $DerGeneFileHandler join("\t", $derName, $species."_funcgen_69_37", "Exon", $exonStart, $exonEnd, ".", $strand , ".", "Parent=".$gene->stable_id())."\n" if $count == 1 and $exonStart < $breakPos1DerCoord;
							print $DerGeneFileHandler join("\t", $derName, $species."_funcgen_69_37", "Exon2", $exonStart, $exonEnd, ".", $strand , ".", "Parent=".$gene->stable_id())."\n" if $count == 2 and $exonEnd > ($breakPos1DerCoord+1);
						}
						
						my $domainfeat = $transcriptAdaptor->fetch_by_stable_id($transcript->stable_id())->translation()->get_all_DomainFeatures if $transcript->biotype() eq "protein_coding";
						while ( my $pfeature = shift @{$domainfeat} ) {
								($domainStart, $domainEnd) = &getTransPosition(1, $pfeature->seq_region_start(), $pfeature->seq_region_end(), $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord) if $count == 1;
								($domainStart, $domainEnd) = &getTransPosition(2, $pfeature->seq_region_start(), $pfeature->seq_region_end(), $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord) if $count == 2;

								print $DerDomainFileHandler join("\t", $derName, $species."_funcgen_69_37", "DomainBreakA", $domainStart, $domainEnd, ".", $strand , ".", "Name=".$pfeature->idesc().";Transcript=".$transcript->stable_id().";URL=<a href = http://pfam.sanger.ac.uk/family/".$pfeature->display_id()." target='blank'>Go to PFAM</a>")."\n" if $pfeature->idesc() ne "" and $count == 1;
								print $DerDomainFileHandler join("\t", $derName, $species."_funcgen_69_37", "DomainBreakB", $domainStart, $domainEnd, ".", $strand , ".", "Name=".$pfeature->idesc().";Transcript=".$transcript->stable_id().";URL=<a href = http://pfam.sanger.ac.uk/family/".$pfeature->display_id()." target='blank'>Go to PFAM</a>")."\n" if $pfeature->idesc() ne "" and $count == 2;
						}
					}
					#Add the feature to the featurehash so it doesn't come up again.
					#$features{$gene->stable_id()} = 1;
				}
		}
				
		#Make a deravitive of the reg. elements for the homo sapiens and mouse
		if($species =~ /homo_sapiens|human|mouse|mus/i){
			makeDerRegElements($registry, $DerRegEleFileHandler, $slice , $derName, $species, $breakPos1, $breakPos2, $breakPos1DerCoord, $breakOrientation, $count);
		}
	}
	
	#Close the files
	close($DerGeneFileHandler);
	close($DerTranscriptFileHandler);
	close($DerDomainFileHandler);
	close($DerRegEleFileHandler);

	#Clean up the trackList.json to delete old tracks
	&makeTrackList($species, $fusionFolderName);
	
	#Add the derivative tracks to JBrowse	
	#Add the tracks for the genes/exons, transcripts and domains
	system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff $fusionFolderName/raw/derGeneFeatures.gff --tracklabel Genes_Der --key Genes_Der --type GeneBreakA,GeneBreakB,Gene,FusionGene,Exon,Exon2 --cssclass feature --subfeatureClasses '{\"Exon\" : \"exon\", \"Exon2\" : \"exon2\" , \"GeneBreakA\" : \"feature5\", \"GeneBreakB\" : \"feature3\"}' --out $fusionFolderName/der/");
	system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff $fusionFolderName/raw/derTranscriptFeatures.gff --tracklabel Transcripts_Der --key Transcripts_Der --type TranscriptBreakA,TranscriptBreakB --cssclass feature --subfeatureClasses '{\"TranscriptBreakA\" : \"feature5\", \"TranscriptBreakB\" : \"feature3\"}' --out $fusionFolderName/der/");
	system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff $fusionFolderName/raw/derDomainFeatures.gff --tracklabel Domains_Der --key Domains_Der --type DomainBreakA,DomainBreakB --cssclass feature --subfeatureClasses '{\"DomainBreakA\" : \"feature5\", \"DomainBreakB\" : \"feature3\"}' --out $fusionFolderName/der/");

	#Add the reg. element tracks
	#Methyl and acyl tracks
	if($species =~ /homo_sapiens|human|mouse|mus/i){
		system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff $fusionFolderName/raw/derRegEleFeatures.gff --tracklabel H3K27_Me_Der --key H3K27_Me_Der --type RegEleMethyl,H3K27me3,H3K27me2,H3K27me1 --out $fusionFolderName/der/ --cssClass feature --clientConfig '{ \"featureCss\": \"background-color:#FE2E2E;border-color:#FE2E2E\"}'");
		system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff $fusionFolderName/raw/derRegEleFeatures.gff --tracklabel H3K27_Ac_Der --key H3K27_Ac_Der --type RegEleAcetyl,H3K27ac3,H3K27ac2,H3K27ac1 --out $fusionFolderName/der/ --cssClass feature --clientConfig '{ \"featureCss\": \"background-color:#F5A9A9;border-color:#F5A9A9\"}'");
		system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff $fusionFolderName/raw/derRegEleFeatures.gff --tracklabel H3K4_Me_Der --key H3K4_Me_Der --type RegEleMethyl,H3K4me3,H3K4me2,H3K4me1 --out $fusionFolderName/der/ --cssClass feature --clientConfig '{ \"featureCss\": \"background-color:#2EFE2E;border-color:#2EFE2E\"}'");
		system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff $fusionFolderName/raw/derRegEleFeatures.gff --tracklabel H3K4_Ac_Der --key H3K4_Ac_Der --type RegEleAcetyl,H3K4ac3,H3K4ac2,H3K4ac1 --out $fusionFolderName/der/ --cssClass feature --clientConfig '{ \"featureCss\": \"background-color:#BCF5A9;border-color:#BCF5A9\"}'");
		system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff $fusionFolderName/raw/derRegEleFeatures.gff --tracklabel H3K9_Me_Der --key H3K9_Me_Der --type RegEleMethyl,H3K9me3,H3K9me2,H3K9me1 --out $fusionFolderName/der/ --cssClass feature --clientConfig '{ \"featureCss\": \"background-color:#0101DF;border-color:#0101DF\"}'");
		system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff $fusionFolderName/raw/derRegEleFeatures.gff --tracklabel H3K9_Ac_Der --key H3K9_Ac_Der --type RegEleAcetyl,H3K9me3,H3K9me2,H3K9me1 --out $fusionFolderName/der/ --cssClass feature --clientConfig '{ \"featureCss\": \"background-color:#819FF7;border-color:#819FF7\"}'");
		system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff $fusionFolderName/raw/derRegEleFeatures.gff --tracklabel H3K79_Me_Der --key H3K79_Me_Der --type RegEleMethyl,H3K79me3,H3K79me2,H3K79me1 --out $fusionFolderName/der/ --cssClass feature --clientConfig '{ \"featureCss\": \"background-color:#DF01D7;border-color:#DF01D7\"}'");
		system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff $fusionFolderName/raw/derRegEleFeatures.gff --tracklabel H3K79_Ac_Der --key H3K79_Ac_Der --type RegEleAcetyl,H3K79me3,H3K79me2,H3K79me1 --out $fusionFolderName/der/ --cssClass feature --clientConfig '{ \"featureCss\": \"background-color:#F5A9F2;border-color:#F5A9F2\"}'");
		
		#DNaseI track
		system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff $fusionFolderName/raw/derRegEleFeatures.gff --tracklabel DNaseI_Der --key DNaseI_Der --type RegEle,DNase1 --out $fusionFolderName/der/");
	}	
	#Make derivatives of supporting files (Link uploaded by user)
	&addSupportTracks($derName, $species, $chr1, $chr2, $breakPos1, $breakPos2, $breakPos1DerCoord, $breakOrientation, $searchDown, $searchUp, $supportFilesSample, $supportFilesFusion, $fusionFolderName);
	
	#Add a track that displays breakpoints around the selected breakpoint in the same sample, also colors the current breakpoint differtly as visualization aid
	&addBreakpointVisualization($derName, $species, $chr1, $chr2, $breakPos1, $breakPos2, $breakPos1DerCoord, $breakOrientation, $sampleID, $fusID, $fusionFolderName);
	
	#Return the position that should be zoomed into (breakpos1)
	return $breakPos1DerCoord;
}

#Returns the location of a feature for the deravative based on the orientation
#Needs the "Normal" 5->3 position of a feature and the der 5' -> 3' coordinate of the first breakpoint, returns transposition, also cuttsoff if overlapping
sub getTransPosition{
	my ($breakCount, $featStart, $featEnd, $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord) = @_;
	my ($o1, $o2, $derFeatStart, $derFeatEnd);
	if ($breakOrientation =~ /^(H|T)(H|T)/i) {
		  $o1 = uc($1);
		  $o2 = uc($2);
	}
	
	#On first break
	$derFeatStart = $breakPos1DerCoord - ($featEnd - $breakPos1) if $breakCount == 1 and $o1 eq "H";
	$derFeatEnd = $breakPos1DerCoord - ($featStart - $breakPos1) if $breakCount == 1 and $o1 eq "H";
	
	($derFeatStart,$derFeatEnd) = ($featStart,$featEnd) if $breakCount == 1 and $o1 eq "T";
	
	#On second break
	$derFeatStart = (($breakPos1DerCoord +1) + ($featStart - $breakPos2)) if $breakCount == 2 and $o2 eq "H";
	$derFeatEnd = (($breakPos1DerCoord +1)  + ($featEnd - $breakPos2)) if $breakCount == 2 and $o2 eq "H";
	
	$derFeatStart = (($breakPos1DerCoord +1) + ($breakPos2 - $featEnd)) if $breakCount == 2 and $o2 eq "T";
	$derFeatEnd = (($breakPos1DerCoord +1)  + ($breakPos2 - $featStart)) if $breakCount == 2 and $o2 eq "T";
	
	#Cut-off feature on breakpoint if overlapping
	$derFeatStart = $breakPos1DerCoord if $derFeatStart > $breakPos1DerCoord and $breakCount == 1;
	$derFeatEnd = $breakPos1DerCoord if $derFeatEnd > $breakPos1DerCoord and $breakCount == 1;
	
	$derFeatStart = ($breakPos1DerCoord+1) if $derFeatStart < $breakPos1DerCoord and $breakCount == 2;
	$derFeatEnd = ($breakPos1DerCoord+1) if $derFeatEnd < $breakPos1DerCoord and $breakCount == 2;
	
	return ($derFeatStart, $derFeatEnd);
}

#This functions makes a deravative track of the reg. elements based on a given slice
sub makeDerRegElements{
	my ($registry, $featFileHandler, $slice, $derName, $species, $breakPos1, $breakPos2, $breakPos1DerCoord, $breakOrientation, $count) = @_;	
	my ($regEleStart, $regEleEnd, $regAttributeStart, $regAttributeEnd, @split);

	#Make funcGen adaptor
	my $regfeat_adaptor = $registry->get_adaptor($species, 'funcgen', 'regulatoryfeature') || die "Couldn't get FuncGen Adaptor for species: $species";
	my $regType;
	#Get all reg. elements from slice
	my @reg_feats = @{$regfeat_adaptor->fetch_all_by_Slice($slice)};
	foreach my $cell_rf (@{@reg_feats}){
		my $makeParentMethyl = 0;
		my $makeParentAcetyl = 0;

	
		#Open the regulatory feature to get all the attributes like histone modification and pol. activity
		my $rfs = $regfeat_adaptor->fetch_all_by_stable_ID($cell_rf->stable_id); 
		#Only store the DNAseI, H3K4, H3K9, H3K27 en H3K79 (mono, di & tri variants and also acetyl)
		foreach my $regEle(@{$rfs}){
			#Get all the attributes of the reg. ele
			foreach my $attr_feat (@{$regEle->regulatory_attributes()}){
				if($attr_feat->display_label =~ /H3K4|H3K9|H3K27|H3K79/i){
					@split = split(" -", $attr_feat->display_label);
					$regType = $split[0];
					
					#Make the parent reg. ele. if not already
					if($makeParentMethyl == 0 and !($regType =~ /ac/i)){
						#Transpose the position
						($regEleStart, $regEleEnd) = &getTransPosition($count, $cell_rf->seq_region_start(), $cell_rf->seq_region_end(), $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
						print $featFileHandler join("\t", $derName, $species, "RegEleMethyl", $regEleStart, $regEleEnd, ".", ".", ".", "ID=".$cell_rf->stable_id."M");
						print $featFileHandler "\n";
						$makeParentMethyl++;
					}
					if($makeParentAcetyl == 0 and $regType =~ /ac/i){
						#Transpose the position
						($regEleStart, $regEleEnd) = &getTransPosition($count, $cell_rf->seq_region_start(), $cell_rf->seq_region_end(), $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
						print $featFileHandler join("\t", $derName, $species, "RegEleAcetyl", $regEleStart, $regEleEnd, ".", ".", ".", "ID=".$cell_rf->stable_id."A");
						print $featFileHandler "\n";
						$makeParentAcetyl++;
					}
					#Make the subfeature
					($regAttributeStart, $regAttributeEnd) = &getTransPosition($count, $attr_feat->seq_region_start(), $attr_feat->seq_region_end(), $breakOrientation, $breakPos1, $breakPos2, $breakPos1DerCoord);
					

					#Print the ID of the reg. ele region, what kind of ele and the start-end of the region
					print $featFileHandler join("\t", $derName, $species, $regType, $regAttributeStart, $regAttributeEnd, ".",".",".", "Name=".$attr_feat->display_label.";Parent=".$cell_rf->stable_id."M") if !($regType =~ /ac/i);
					print $featFileHandler join("\t", $derName, $species, $regType, $regAttributeStart, $regAttributeEnd, ".",".",".", "Name=".$attr_feat->display_label.";Parent=".$cell_rf->stable_id."A") if $regType =~ /ac/i;
					print $featFileHandler "\n";
				}
			}
		}
	}	
}

#Adds the real sequence to the Jbrowse, so must convert to .json (Takes a while)
sub addRefSeqTrue{
	my ($species, $fusionFolderName) = @_;
	system("perl ../genomeBrowser/bin/prepare-refseqs.pl --fasta $fusionFolderName/raw/derivative.fa --out $fusionFolderName/der/");
	system("rm $fusionFolderName/raw/derivative.fa");
}

#Make a real derivative chromosome by joining the nucleotide sequence of 2 chromosomes together into a derivative
sub makeDerChromosomeTrue{
	my ($registry, $chr1, $chr2, $chr1_Start, $chr1_End, $chr2_Start, $chr2_End, $species, $breakOrientation, $fusionFolderName) = @_;
	my $slice_adaptor = $registry->get_adaptor( $species, 'Core', 'Slice' ) || die "Couldn't get Slice adaptor for species: $species";
	my ($chr1NucSeq, $chr2NucSeq, $chr1DerNucSeq, $chr2DerNucSeq, $o1, $o2);
	if ($breakOrientation =~ /^(H|T)(H|T)/i) {
		 $o1 = uc($1);
		 $o2 = uc($2);
	}
	
	#Get the nucleotide sequence from the ensembl DB depending on the orientation
	$chr1NucSeq = $slice_adaptor->fetch_by_region( 'chromosome', $chr1, 1, $chr1_End)->seq() || die "Couldn't get chromosome of $chr1 of species: $species" if $o1 eq "T";
	$chr1NucSeq = $slice_adaptor->fetch_by_region( 'chromosome', $chr1, $chr1_Start)->seq() || die "Couldn't get chromosome of $chr1 of species: $species" if $o1 eq "H";

	$chr2NucSeq = $slice_adaptor->fetch_by_region( 'chromosome', $chr2, $chr2_Start)->seq() || die "Couldn't get chromosome of $chr2 of species: $species" if $o2 eq "H";
	$chr2NucSeq = $slice_adaptor->fetch_by_region( 'chromosome', $chr2, 1, $chr2_End)->seq() || die "Couldn't get chromosome of $chr2 of species: $species" if $o2 eq "T";
	
	#Path to the deravitive chromosome
	my $derFile =  "$fusionFolderName/raw/derivative.fa";
	open(my $newFile, '>', $derFile) || die "Could not write/create to derivate nuc.seq file: $derFile\n";

	#Open a dbAdaptor to get the karyotype bands
	my $dbAdaptor = &openDBAdaptor;
	#Get the correct nomenclature for this derivative based on the location of the centromere.
	my $derName = &getDerName($dbAdaptor, $chr1, $chr2, $chr1_Start, $chr1_End);
	
	if($chr1_Start == $chr1_End and $chr2_Start == $chr2_End){
		#DNA is complemented because of the flip of the strand at this orientation
		$chr1DerNucSeq = reverse($chr1NucSeq) if $o1 eq "H";
		$chr1DerNucSeq =~ tr/ACGTacgt/TGCAtgca/ if $o1 eq "H";
		$chr2DerNucSeq = reverse($chr2NucSeq) if $o2 eq "T";
		$chr2DerNucSeq =~ tr/ACGTacgt/TGCAtgca/ if $o2 eq "T";
		
		$chr1NucSeq = $chr1DerNucSeq if $o1 eq "H";
		$chr2NucSeq = $chr2DerNucSeq if $o2 eq "T";
	}
        
        #Print the derivative nuc. seq to the $newFile, with >dername \n (fasta)
	print $newFile ">$derName \n";
	print $newFile $chr1NucSeq;
	print $newFile $chr2NucSeq;
	
	#Close files
	close($newFile);
	
	#Add the deravitive chromosome to the refseqs, only add the nuc. seq to JBrowse if s1/e1 or s2/e2 are equal (Confirmed)
	if($chr1_Start == $chr1_End and $chr2_Start == $chr2_End){
		&addRefSeqTrue($species, $fusionFolderName);
	}else{
		&addRefSeqFake($species, $fusionFolderName);
	}
	#Return the filename of the derivative chromosome
	return $derName;
}

#Adds the noseq deravative chromosome as a refSeq into the JBrowse, overwrite older derivates.
sub addRefSeqFake{
	my ($species, $fusionFolderName) = @_;
	system("perl ../genomeBrowser/bin/prepare-refseqs.pl --noseq asd --fasta $fusionFolderName/raw/derivative.fa --out $fusionFolderName/der/");
	system("rm $fusionFolderName/raw/derivative.fa");

}

sub getDerName{
	my ($dbAdaptor, $chr1, $chr2, $chr1_Start, $chr1_End) = @_;
	#Variable to store the centrosome location
	my $centroLocChr1 = 0;
	#Get karyotype Adaptor
	my $karyAdaptor = $dbAdaptor->get_KaryotypeBandAdaptor();
	#Fetch the bands for one chromosome, one is sufficient since if the break if before the centromere, the deravative uses the centromere of the other chromosome.
	my @bandsChr1 = @{ $karyAdaptor->fetch_all_by_chr_name($chr1) };
	#Find the centromere, search for "acen" stain.
	foreach(@bandsChr1){
		$centroLocChr1 = $_->seq_region_end() if $_->stain() eq "acen" and $_->seq_region_end() > $centroLocChr1;
	}
	#Return the derivative name based on centromere position
	return "der($chr1)($chr1;$chr2)" if $centroLocChr1 < $chr1_End; #If the centromere is on the 'first' break
	return "der($chr2)($chr2;$chr1)" if $centroLocChr1 > $chr1_End; #If the centromere is on the 'second break'
}

#This subroutine will cleanse the tracklist of any old tracks with the exception of:
# -Genes, Transcripts, Domains (Add them twice, once with blue colors and once with green colors to differentiate between fusionpoint)
sub makeTrackList {
	my ($species, $fusionFolderName) = @_;
	open(my $trackFileHandlerRef, ">", "$fusionFolderName/ref/trackList.json") || die("Could'nt open $fusionFolderName/ref/trackList.json for editing");
	open(my $trackFileHandlerDer, ">", "$fusionFolderName/der/trackList.json") || die("Could'nt open $fusionFolderName/der/trackList.json for editing");

	#Make a tracklist for the reference tracks
	print $trackFileHandlerRef '{
   "tracks":[
      {
         "chunkSize":20000,
         "urlTemplate":"seq/{refseq_dirpath}/{refseq}-",
         "type":"SequenceTrack",
         "label":"DNA",
         "key":"DNA"
      },
      {
         "style":{
            "className":"feature3",
            "subfeatureClasses" : {
            	"Exon" : "exon2"
            }
         },
         "key":"GeneBreakB",
         "urlTemplate":"tracks/Genes/{refseq}/trackData.json",
         "phase":null,
         "compress":0,
         "type":"FeatureTrack",
         "label":"GeneBreakB",
         "subfeatures":null
      },
      {
         "style":{
            "className":"feature3"
         },
         "key":"TranscriptsBreakB",
         "urlTemplate":"tracks/Transcripts/{refseq}/trackData.json",
         "phase":null,
         "compress":0,
         "type":"FeatureTrack",
         "label":"TranscriptsBreakB",
         "subfeatures":null
      },
      {
         "style":{
            "className":"feature3"
         },
         "key":"DomainsBreakB",
         "urlTemplate":"tracks/Domains/{refseq}/trackData.json",
         "phase":null,
         "compress":0,
         "type":"FeatureTrack",
         "label":"DomainsBreakB",
         "subfeatures":null
      },
      {
         "style":{
            "className":"feature5",
            "subfeatureClasses" : {
            	"Exon" : "exon"
            }
         },
         "key":"GeneBreakA",
         "urlTemplate":"tracks/Genes/{refseq}/trackData.json",
         "phase":null,
         "compress":0,
         "type":"FeatureTrack",
         "label":"GeneBreakA",
         "subfeatures":null
      },
      {
         "style":{
            "className":"feature5"
         },
         "key":"TranscriptsBreakA",
         "urlTemplate":"tracks/Transcripts/{refseq}/trackData.json",
         "phase":null,
         "compress":0,
         "type":"FeatureTrack",
         "label":"TranscriptsBreakA",
         "subfeatures":null
      },
      {
         "style":{
            "className":"feature5"
         },
         "key":"DomainsBreakA",
         "urlTemplate":"tracks/Domains/{refseq}/trackData.json",
         "phase":null,
         "compress":0,
         "type":"FeatureTrack",
         "label":"DomainsBreakA",
         "subfeatures":null
      },
      {
         "style" : {
            "className" : "feature"
         },
         "key" : "Methyl_AcylElements",
         "urlTemplate" : "tracks/Methyl_AcylElements/{refseq}/trackData.json",
         "phase" : null,
         "compress" : 0,
         "type" : "FeatureTrack",
         "label" : "Methyl_AcylElements",
         "subfeatures" : null
      }

   ],
   "formatVersion":1
}';
	close($trackFileHandlerRef);
	
	#Make der tracklist
	print $trackFileHandlerDer '{
   "tracks":[
      {
         "chunkSize":20000,
         "urlTemplate":"seq/{refseq_dirpath}/{refseq}-",
         "type":"SequenceTrack",
         "label":"DNA",
         "key":"DNA"
      }   ],
   "formatVersion":1
}';
	close($trackFileHandlerDer);

	
}

##########################################
###Fill fusion DB with Ensembl features###
#####Mainly used by fusionPipe01.pl#######
##########################################
#These functions are used to get genomic features ((fusion-)genes) from the Ensembl DB
#and insert them into the fusion DB.

#Manages the feature search for the uploading of the breakpoints, calls the appropiate functions to get the genes(and check if fusion-gene) from the Ensembl DB	
sub getFeaturesFromEnsembl{
	my ($dbh, $registry, $refFusion, $searchUp, $searchDown, $species) = @_;
	#Get the genes with corresponding sub-features
	my ($breakSlice1, $breakSlice2) = getEnsemblSlices($registry, $refFusion->{"break_chr1"}, $refFusion->{"break_chr1_start"} ,$refFusion->{"break_chr1_end"}, $refFusion->{"break_chr2"}, $refFusion->{"break_chr2_start"},$refFusion->{"break_chr2_end"},$refFusion->{"break_orientation"}, $searchUp, $searchDown, $species);
	#Check for fusion genes, promotor break, geneslk
	findFusionGenesAndGeneBreaks($dbh, $registry, $species, $refFusion->{"break_chr1"}, $refFusion->{"break_chr1_start"}, $refFusion->{"break_chr1_end"}, $refFusion->{"break_chr2"}, $refFusion->{"break_chr2_start"}, $refFusion->{"break_chr2_end"} ,$refFusion->{"break_orientation"}, $refFusion->{"fusion_ID"});
}

#Gets a slice of the core or Funcgen DB based on the breakpoint coordinate, orientation and user-defined search distance
#Return 2 slices, one for the first breakpoint and one for the second breakpoint.
#The Core DB is used to get the genes, exons, transcripts, domains
#The Funcgen DB is used to get the regulatory elements
#$type is used as flag to determine whether to use Core DB or funcgen DB based on the data needed
sub getEnsemblSlices{
	my($registry, $chr1, $chr1Start, $chr1End, $chr2, $chr2Start, $chr2End , $breakOrientation, $searchUp, $searchDown, $species) = (@_);
	my $sa = $registry->get_adaptor($species, 'Core', 'Slice') || die("Could not get slice adaptor for species: $species");
	my ($breakSliceA, $breakSliceB, $o1, $o2);
	if ($breakOrientation =~ /^(H|T)(H|T)/i) {
		$o1 = uc($1);
		$o2 = uc($2);
	}

	#Get slices based on breakpoint coordinates and user defined search distance for gene search purposes
	#This retrieves slices that can be used for gene-searching and regulatory element searching
	$breakSliceA = $sa->fetch_by_region( 'chromosome', $chr1, ($chr1End-$searchDown), ($chr1End+$searchUp)) if $o1 eq "T";
	$breakSliceA = $sa->fetch_by_region( 'chromosome', $chr1, ($chr1Start-$searchDown), ($chr1Start+$searchUp)) if $o1 eq "H";

	$breakSliceB = $sa->fetch_by_region( 'chromosome', $chr2, ($chr2Start-$searchDown), ($chr2Start+$searchUp)) if $o2 eq "T";
	$breakSliceB = $sa->fetch_by_region( 'chromosome', $chr2, ($chr2End-$searchDown), ($chr2End+$searchUp)) if $o2 eq "H";
	
	return $breakSliceA, $breakSliceB;
}

#Checks the ensembl DB to see if there is a gene on the breakpoint position, then checks ORF and intronic/exonic breakposition
sub findFusionGenesAndGeneBreaks{
	my ($dbh, $registry, $species, $breakChr1, $breakChr1Start, $breakChr1End, $breakChr2, $breakChr2Start, $breakChr2End, $breakOrientation, $fusID) = @_;
	my ($breakSlice1, $breakSlice2, $o1, $o2, %foundGenes);
	#Query statement to alter the featuretype of a gene if needed (To Fusion-overlapping gene with out-frame exons or Fusion-overlapping gene with in-frame exons)
	my $pdh = $dbh->prepare("INSERT INTO T_features (T_fusionpoint_fusion_ID, feat_ENS_ID, T_featuretype_featuretype_ID, featName) VALUES(?,?,?,?)");

	#See if there is a gene on the position of the break
	my $sa = $registry->get_adaptor($species, 'Core', 'Slice');
	($breakSlice1, $breakSlice2) = ($sa->fetch_by_region( 'chromosome', $breakChr1, $breakChr1End, $breakChr1End), $sa->fetch_by_region( 'chromosome', $breakChr2, $breakChr2Start, $breakChr2Start)) if uc($breakOrientation) eq "TH";
	($breakSlice1, $breakSlice2) = ($sa->fetch_by_region( 'chromosome', $breakChr1, $breakChr1Start, $breakChr1Start), $sa->fetch_by_region( 'chromosome', $breakChr2, $breakChr2End, $breakChr2End)) if uc($breakOrientation) eq "HT";
	($breakSlice1, $breakSlice2) = ($sa->fetch_by_region( 'chromosome', $breakChr1, $breakChr1Start, $breakChr1Start), $sa->fetch_by_region( 'chromosome', $breakChr2, $breakChr2Start, $breakChr2Start)) if uc($breakOrientation) eq "HH";
	($breakSlice1, $breakSlice2) = ($sa->fetch_by_region( 'chromosome', $breakChr1, $breakChr1End, $breakChr1End), $sa->fetch_by_region( 'chromosome', $breakChr2, $breakChr2End, $breakChr2End)) if uc($breakOrientation) eq "TT";
	#Get genes on the breakpoints
	my ($genesBreak1, $genesBreak2) = (&getGenesFromSlice($breakSlice1), &getGenesFromSlice($breakSlice2));
	
	#Get orientation  
	if ($breakOrientation =~ /^(H|T)(H|T)/i) {
		  $o1 = uc($1);
		  $o2 = uc($2);
	}
	#Get the correct breakposition
	#For TH
	my ($breakPos1, $breakPos2) = ($breakChr1End, $breakChr2Start);
	#If first is H
	$breakPos1 = $breakChr1Start;
	#If second is T
	$breakPos2 = $breakChr2End;
	
	#Check in/out frame of fusiongenes if genes were found
	if(scalar keys(%$genesBreak1) != 0 and scalar keys(%$genesBreak2) != 0){
		#Since a gene is possible on the - strand and on the + strand at the same position, check all genes for possible fusion genes
		foreach my $geneOn1(keys %$genesBreak1){
			foreach my $geneOn2(keys %$genesBreak2){
				#If gene already found in fusiongene
				$foundGenes{$geneOn1}{'bp1'} = 0 unless $foundGenes{$geneOn1}{'bp1'};
				$foundGenes{$geneOn2}{'bp2'} = 0 unless $foundGenes{$geneOn2}{'bp2'};
				
				my $fusion = "Out-frame";
				my $geneName = join("::", $geneOn1, $geneOn2);

				#If the strand allow fusiongenes, mark it so
				$fusion = "Possible-In-frame" if $o1 ne $o2 and $genesBreak1->{$geneOn1}->{"strand"} eq $genesBreak2->{$geneOn2}->{"strand"};
				$fusion = "Possible-In-frame" if $o1 eq $o2 and $genesBreak1->{$geneOn1}->{"strand"} ne $genesBreak2->{$geneOn2}->{"strand"};
				
				#Check if the fusiongenes actually form an in-frame product if the breakpoint is confirmed(BreakChr1Start == BreakChr1End)
				#Only applies if one or both of the cuts are exonic and the frame stays the same or both cuts are intronic
				if($breakChr1Start == $breakChr1End and $fusion eq "Possible-In-frame"){
					my ($frameExonBreak1, $frameExonBreak2, $exonPos);
					#If the breakpoint on break1 is exonic
					$frameExonBreak1 = ($genesBreak1->{$geneOn1}->{exonEnd}-$breakPos1+$genesBreak1->{$geneOn1}->{frame})%3 if $genesBreak1->{$geneOn1}->{"type"} eq "Exonic" and $o1 eq "T";
					$frameExonBreak1 = ($genesBreak1->{$geneOn1}->{exonStart}-$breakPos1+$genesBreak1->{$geneOn1}->{frame})%3 if $genesBreak1->{$geneOn1}->{"type"} eq "Exonic" and $o1 eq "H";
					#If the breakpoint on break2 is exonic
					$frameExonBreak2 = ($genesBreak2->{$geneOn2}->{exonEnd}-$breakPos2+$genesBreak2->{$geneOn2}->{frame})%3 if $genesBreak2->{$geneOn2}->{"type"} eq "Exonic" and $o2 eq "T";
					$frameExonBreak2 = ($genesBreak2->{$geneOn2}->{exonStart}-$breakPos2+$genesBreak2->{$geneOn2}->{frame})%3 if $genesBreak2->{$geneOn2}->{"type"} eq "Exonic" and $o2 eq "H";
										
					#Get the exon closest to the first breakpoint if there wasn't a break in an exon
					if($genesBreak1->{$geneOn1}->{"type"} ne "Exonic"){
						my $geneSlice = $sa->fetch_by_gene_stable_id($genesBreak1->{$geneOn1}->{"ensID"});
						#Get the frame of the exon closest to the breakpoint (Based on orientation)
						foreach my $exon (@{$geneSlice->get_all_Exons()}){
							
							#Get last exon before the breakpoint
							last if $exon->seq_region_start > $breakPos1 and $o1 eq "T";
							$frameExonBreak1 = ($exon->seq_region_end-$breakPos1+$exon->phase())%3 if $o1 eq "T";
							
							#Get first exon after the breakpoint
							last if defined($frameExonBreak1) and $o1 eq "H";
							$frameExonBreak1 = ($exon->seq_region_start-$breakPos1+$exon->phase())%3 if $o1 eq "H" and $exon->seq_region_end > $breakPos1;
						}
					}
					
					#Get the exon closest to the second breakpoint if there wasn't a break in an exon
					if($genesBreak2->{$geneOn2}->{"type"} ne "Exonic"){
						my $geneSlice = $sa->fetch_by_gene_stable_id($genesBreak2->{$geneOn2}->{"ensID"});
						#Get the frame of the exon closest to the breakpoint (Based on orientation)
						foreach my $exon (@{$geneSlice->get_all_Exons()}){
							
							#Get last exon before the breakpoint
							last if $exon->seq_region_start > $breakPos2 and $o2 eq "T";
							$frameExonBreak2 = ($exon->seq_region_end-$breakPos2+$exon->phase())%3 if $o2 eq "T";
							
							#Get first exon after the breakpoint
							last if defined($frameExonBreak2);
							$frameExonBreak2 = ($exon->seq_region_start-$breakPos2+$exon->phase())%3 if $o2 eq "H" and $exon->seq_region_start < $breakPos2;
						}
					}
					#This fusion results in an in-frame protein product or not
					$fusion = "Real-In-frame" if $frameExonBreak2 == $frameExonBreak1;
					$fusion = "Real-Out-frame" if $frameExonBreak2 != $frameExonBreak1;
				}
				
				#Real in-frame fusiongene
				if( $fusion eq "Real-In-frame"){				
					foreach my $gene(%$genesBreak1->{$geneOn1}, %$genesBreak2->{$geneOn2}){
						#Update the database accordingly based on exonic/intronic and in/out frame
						$pdh->execute($fusID, $gene->{"ensID"}, 1, $geneName);
					}
					$foundGenes{$geneOn2}{'bp2'} = 1;
					$foundGenes{$geneOn1}{'bp1'} = 1;
				}
				
				#Real out-frame fusiongene
				if( $fusion eq "Real-Out-frame"){				
					foreach my $gene(%$genesBreak1->{$geneOn1}, %$genesBreak2->{$geneOn2}){
						#Update the database accordingly based on exonic/intronic and in/out frame
						$pdh->execute($fusID, $gene->{"ensID"}, 2, $geneName);
					}
					$foundGenes{$geneOn2}{'bp2'} = 1;
					$foundGenes{$geneOn1}{'bp1'} = 1;

				}
				#Possible fusion-gene (No exact breakposition)
				if( $fusion eq "Possible-In-frame"){				
					foreach my $gene(%$genesBreak1->{$geneOn1}, %$genesBreak2->{$geneOn2}){
						#Update the database accordingly based on exonic/intronic and in/out frame
						$pdh->execute($fusID, $gene->{"ensID"}, 3, $geneName) if $gene->{"type"} eq "Intronic";
						$pdh->execute($fusID, $gene->{"ensID"}, 4, $geneName) if $gene->{"type"} eq "Exonic";
					}
					$foundGenes{$geneOn2}{'bp2'} = 1;
					$foundGenes{$geneOn1}{'bp1'} = 1;
				}	
			}
		}
	}
	
	#Only report gene-breaks if it isn't a fusiongene with something else
	foreach my $gene(keys %foundGenes){
		#Check if a gene exist in a break1
		if(exists($foundGenes{$gene}{'bp1'})){
			#Check if gene was already found in fusiongene
			#If on break 1
			if($foundGenes{$gene}{'bp1'} == 0){
			   	#Update the database accordingly based on exonic/intronic and in/out frame
			   	$pdh->execute($fusID, $genesBreak1->{$gene}->{"ensID"}, 5, $genesBreak1->{$gene}->{"name"}) if $genesBreak1->{$gene}->{"type"} eq "Intronic";
			   	$pdh->execute($fusID, $genesBreak1->{$gene}->{"ensID"}, 6,  $genesBreak1->{$gene}->{"name"}) if $genesBreak1->{$gene}->{"type"} eq "Exonic";
			}
		}
		if(exists($foundGenes{$gene}{'bp2'})){
			#Check if gene was already found in fusiongene
			#If on break 1
			if($foundGenes{$gene}{'bp2'} == 0){
			   	#Update the database accordingly based on exonic/intronic and in/out frame
			   	$pdh->execute($fusID, $genesBreak2->{$gene}->{"ensID"}, 5, $genesBreak2->{$gene}->{"name"}) if $genesBreak2->{$gene}->{"type"} eq "Intronic";
			   	$pdh->execute($fusID, $genesBreak2->{$gene}->{"ensID"}, 6,  $genesBreak2->{$gene}->{"name"}) if $genesBreak2->{$gene}->{"type"} eq "Exonic";
			}
		}
		
	}
	
	
	#Check if there is a break in a gene without creating a fusiongene and in the first 1000 bp(Possibly leaving an "empty" promotor behind for another gene)
	#Search promotor breaks on first breakpoint
	if(scalar keys(%$genesBreak1) != 0){
		foreach my $geneOn1(keys %$genesBreak1){
			#Check promotor break on gene on +/- strand
			$pdh->execute($fusID, $genesBreak1->{$geneOn1}->{"ensID"}, 7 , $genesBreak1->{$geneOn1}->{"name"}) if ($breakPos1 - $genesBreak1->{$geneOn1}->{start}) <= $promotorBreakMinBP and $genesBreak1->{$geneOn1}->{strand} eq "+";
			$pdh->execute($fusID, $genesBreak1->{$geneOn1}->{"ensID"}, 7, $genesBreak1->{$geneOn1}->{"name"}) if ($genesBreak1->{$geneOn1}->{end} - $breakPos1) <= $promotorBreakMinBP and $genesBreak1->{$geneOn1}->{strand} eq "-";
		}
	}

	#Search promotor breaks on first breakpoint
	if(scalar keys(%$genesBreak2) != 0){
		foreach my $geneOn2(keys %$genesBreak2){
			#Check promotor break on gene on +/- strand
			$pdh->execute($fusID, $genesBreak2->{$geneOn2}->{"ensID"}, 7, $genesBreak2->{$geneOn2}->{"name"}) if ($breakPos2 - $genesBreak2->{$geneOn2}->{start}) <= $promotorBreakMinBP and $genesBreak2->{$geneOn2}->{strand} eq "+";
			$pdh->execute($fusID, $genesBreak2->{$geneOn2}->{"ensID"}, 7, $genesBreak2->{$geneOn2}->{"name"}) if ($genesBreak2->{$geneOn2}->{end} - $breakPos2) <= $promotorBreakMinBP and $genesBreak2->{$geneOn2}->{strand} eq "-";
		}	
	}
}

#Subroutine to get all the genes from a slice and return them in an array
sub getGenesFromSlice{
	my $slice = shift;
	#Make hash that keeps track of genes
	my %genesFound;
	my $genesOnSlice = $slice->get_all_Genes();
	#Check if there are genes, if yes -> concenate on string
	foreach my $gene (@{$genesOnSlice}){
		my $strand = "+";
		$strand = "-" if $gene->strand() == -1;
		#Check if it is an intronic or exonic break
		my $type = "Intronic";
		my ($exonStart, $exonEnd, $phase) = ("None", "None", "None");
		#Check if the exon found is from the actual gene (strand-wise)
		foreach my $exon (@{$slice->get_all_Exons()}){
			if($exon->strand() == $gene->strand()){
				$type = "Exonic";
				$phase = $exon->phase;
				$exonStart = $exon->seq_region_start;
				$exonEnd = $exon->seq_region_end;
			}
		}
		#Fill hash with gene info
		$genesFound{$gene->external_name} = {
			ensID => $gene->stable_id,
			name => $gene->external_name,
			strand => $strand,
			type => $type,
			frame => $phase,
			exonStart => $exonStart,
			exonEnd => $exonEnd,
			start => $gene->seq_region_start,
			end => $gene->seq_region_end
		};
	}
	#Return a string with the genes, strand of genes and breaktype
	return \%genesFound;
}

############################################
###Upload and handling of breakpoint file###
############################################
#These functions are used to facilitate the uploading and insertion of the mate-pair coordinates coordinates
#into the local fusion DB

#This opens and reads the submitted breakpoint data file
#This file should in the following format: (| == \t)
#Chr1|Start|End|Chr2|Start|End|Orientation|Breaktype or the new 13-tabbed 123SV output
sub readFusionFile{
		my ($q, $filename,  $species, $supportFiles) = @_;
		my $pdh;
		open(my $filehandler, '<', "$filename") || die("Could not read the breakpoint file\n");
		print "<p>Reading file and storing breakpoints into fusion DB</p>";
		#Open connection to DB
		my($dbh, $sample_ID) = &openFusionConnectionGetSampleID($filename,  $species, $supportFiles);
		#Per line (Also pass the $dbh and $sample_ID)
		while (<$filehandler>) {
			next if $_ =~ /^#/;
			#Split the line and check for format (8 or 13 tabs)
			my @fusData = (splitLine(chomp($_),$dbh));
			#Save data in fusionDB
			$pdh = uploadBreakpointData(\@fusData, $dbh, $sample_ID);
		}
		#Close and delete the temporary file
		close($filehandler);
		unlink($filename);
		#If nothing went wrong, commit the data and close the connection
		&closeFusionConnectionWithFinishAndCommit($dbh, $pdh);
		return $sample_ID;
}

#Splits (and checks) the line to get the relevant breakpoint data.
sub splitLine{
	#Split on Tab
	my @fusDataHolder = split(/\t/,$_);
	my @fusData;
	
	foreach(@fusDataHolder){
		chomp;
		push(@fusData, $_);
	}
	
	#Check if there are 8 or 11 (123SV) elements so that the correct format is used or else die
	if(scalar(@fusData) != 8 and scalar(@fusData) != 13){
		die("The submitted file was not in the correct format, either the columns did not count to 11 or no tab delimiter was used or no newline used between lines\n");
	}
	#Check if "chr" is used in front of the chromosomenumber and delete in the fields:Chr1(0) and Chr2(3)
	$fusData[0] =~ s/chr// if $fusData[0] =~ m/chr[x|y|\d]/;
	$fusData[3] =~ s/chr// if $fusData[3] =~ m/chr[\w|\d]/;
	#Return array with data
	return @fusData;
}

#Temporarily save the breakpoint file on the local server for reading purposes
sub saveTempFile{
	my ($q,$filename) = @_;
	my ($bytesread, $buffer, $totalbytes);
	my $num_bytes = 1024;
	my $untainted_filename;	
	#Check if the breakpoint file has been uploaded
	if (!$filename) {
		print $q->p("Upload the file with the breakpoint positions!");
		return;
	}
	
	# Untaint $filename
	if ($filename =~ /^([-\@:\/\\\w.]+)$/) {
		$untainted_filename = $1;
	} else {
		die <<"EOT";
		Unsupported characters in the filename "$filename". 
		Your filename may only contain alphabetic characters and numbers, 
		and the characters '_', '-', '\@', '/', '\\' and '.'
EOT
	}
	
	if ($untainted_filename =~ m/\.\./) {
		die <<"EOT";
		Your upload filename may not contain the sequence '..' 
		Rename your file so that it does not include the sequence '..', and try again.
EOT
	}
	
	#Store the file temporarily on the local server
	my $file = "/tmp/$untainted_filename";
	
	open (OUTFILE, ">", "$file") or die "Couldn't write temp file: $!";
	
	while ($bytesread = read($filename, $buffer, $num_bytes)) {
		$totalbytes += $bytesread;
		print OUTFILE $buffer;
	}
	die "Read failure" unless defined($bytesread);
	unless (defined($totalbytes)) {
		print "<p>Error: Could not read file ${untainted_filename}, ";
		print "or the file was zero length.";
	} else {
		print "<p>File $filename succesfully uploaded! ($totalbytes bytes)";
	}
	close OUTFILE or die "Couldn't close $file: $!";
	return $file;
}

#Opens connection to the fusionpoint DB to facilitate the upload of data (Non-auto committing)
#Also gets the new SampleID for use in &uploadData
sub openFusionConnectionGetSampleID{
	#Get filename, remove the extension of the file.
	my ($samplename, $species, $supportFiles) = @_;
	my $pdh;
	
	#Get the name of the sample from the file
	$samplename = basename($samplename);
	$samplename =~ s{\.[^.]+$}{};
	
	#Connect to local fusion DB, autocommit off to make use of a single transaction (Data will only be stored on DB with COMMIT)
	my $dbh = DBI->connect("dbi:mysql:FusionAnnotator:localhost:3306", "*", "*", { AutoCommit => 0 }) || die "Could not connect to fusionpoint database: $DBI::errstr \;";
	
	#Check if the sample has already been uploaded to the fusionDB through the StrucVarDB
	my @array = $dbh->selectrow_array("SELECT sample_ID FROM T_samples WHERE sample_name = '$samplename'");
	my $sample_ID = $array[0];
	
	#If the sample ID does not already exist
	unless($sample_ID){
		#Insert Sample
		$pdh = $dbh->prepare('INSERT INTO T_samples(sample_name, sample_species) VALUES(?,?)') || die "Could not prepare sample statement";
		$pdh->execute($samplename, $species) || die "Could not upload new sample";
		$dbh->commit() || die "Could not commit data";
		#Get the sample_ID
		$sample_ID = $pdh->{mysql_insertid};
	}
	
	#Insert the support files into the DB
	$pdh = $dbh->prepare('INSERT INTO T_support(support_link, T_samples_sample_ID, T_supType_supType_ID) VALUES(?,?,?)') || die "Could not prepare support File insertion";
	
	#Insert the other support files
	foreach my $supportFile (@{$supportFiles}){
			$pdh->execute($supportFile, $sample_ID, 1) || die "Could not upload BED File: $supportFile" if $supportFile =~ /(.bam$)/i;
			$pdh->execute($supportFile, $sample_ID, 2) || die "Could not upload BED File: $supportFile" if $supportFile =~ /(.bed$)/i;
			$pdh->execute($supportFile, $sample_ID, 3) || die "Could not upload WIG File: $supportFile" if $supportFile =~ /(.wig$)/i;
			$pdh->execute($supportFile, $sample_ID, 4) || die "Could not upload GFF File: $supportFile" if $supportFile =~ /(.gff$)/i;
	}
	print "<p>Given sample ID: $sample_ID </p>";

	#Return DBhandler, PREPARE statement for insertion and sample_ID
	return $dbh, $sample_ID;
}

#Uploads a selected breakpoint from the StrucVarDB to the fusionDB with sample_ID 1 which holds all the samples submitted from the StrucVarDB
sub uploadBreakpointDataFromStrucVarDB{
	my ($dbh, $sampleName, $species, $chr1, $s1, $e1, $chr2, $s2, $e2, $or, $breakType) = (@_);
	my (@array, $fusID);
	#Check if the breakpoint already exist in the fusionDB, if so return that fusion ID
	@array = $dbh->selectrow_array("SELECT fusion_ID FROM T_fusionpoint 
		INNER JOIN T_samples ON T_fusionpoint.T_samples_sample_ID=T_samples.sample_ID
		WHERE break_chr1 = $chr1 AND break_chr1_start = $s1 AND break_chr1_end = $e1 AND break_chr2 = $chr2 AND break_chr2_start = $s2 AND break_chr2_end = $e2 AND sample_name = '$sampleName'");
	$fusID = $array[0];
	return $fusID if $fusID != "";
	
	#Get the correct species name
	$species = "homo_sapiens" if $species =~ /human|homo/i;
	$species = "danio_rerio" if $species =~ /danio|rerio|zebra/i;
	$species = "c_elegans" if $species =~ /elegans/i;
	$species = "rattus_rattus" if $species =~ /rat|rattus/i;
	
	#If it does'nt exist, insert into fusionDB
	my $pdh = $dbh->prepare('INSERT INTO T_fusionpoint (break_chr1, break_chr1_start, break_chr1_end, break_chr2, break_chr2_start, break_chr2_end, break_orientation, strucVarOrigin, T_breaktype_breaktype_ID, T_samples_sample_ID)  
				VALUES(?,?,?,?,?,?,?,?,?,?)');
	
	#Get the breaktype (Insertion/Deletion/Remote/etc.), corresponds with database entry
	#This is done to prevent data redundancy, all types can be clearly categorized
	my $breakTypeID = 9;
	$breakTypeID = 1 if $breakType =~ /Ins/i; 		#Insertion synonyms
	$breakTypeID = 2 if $breakType =~ /Del/i; 		#Deletion synonyms
	$breakTypeID = 3 if $breakType =~ /Anti/i;		#Anti synonyms
	$breakTypeID = 4 if $breakType =~ /(Inversion|Invrs)/i;	#Inversion synonyms
	$breakTypeID = 5 if $breakType =~ /(Remote|Rem)/i;	#Remote synonyms
	$breakTypeID = 6 if $breakType =~ /Trans/i; 		#Translocation
	$breakTypeID = 7 if $breakType =~ /Dup/i;		#Duplication
	$breakTypeID = 8 if $breakType =~ /Evertion/i; 		#Evertion
	
	my $orientation;
	$orientation = "TH" if $or =~ /th/i;
	$orientation = "TT" if $or =~ /tt/i;
	$orientation = "HH" if $or =~ /hh/i;
	$orientation = "HT" if $or =~ /ht/i;
	
	#Find if the sample has already been added once, if so get that sample ID. Else make a new sample.
	@array = $dbh->selectrow_array("SELECT sample_ID FROM T_samples WHERE sample_name = '$sampleName'");
	my $sampleID = $array[0];
		
	unless($sampleID){
		my $pdhSample = $dbh->prepare("INSERT INTO T_samples (sample_species, sample_name) VALUES(?,?)");
		$pdhSample->execute($species, $sampleName);
		#Get the sampleID back
		$sampleID = $pdhSample->{mysql_insertid};
	}
	
	#Upload data to fusionpoint DB using the prepare statement with placeholders to handle the quoting and prevent "bad" data
	$pdh->execute( $chr1, $s1, $e1, $chr2, $s2, $e2, $orientation, "1", $breakTypeID, $sampleID) || die "Could not upload the selected StrucVarDB breakpoint";
	#Get the fus_ID back
	$fusID = $pdh->{mysql_insertid};
	$pdh->finish();
	$dbh->commit();
	
	return $fusID;
	
}


#Uploads the data from the splitLine subroutine into the database (Without committing the data)
#Also keeps track of the given fusionpointIDs for use in the Analyze_Data.pl script.
#Gets the data in Array (split on tab), the prepareHandler and the sample_ID.
sub uploadBreakpointData{
	my($fusdata ,$dbh ,$sample_ID) = (@_);
	
	#Prepare with placeholders (Less load on DB and prevents SQL Injection etc.)
	my $pdh = $dbh->prepare('INSERT INTO T_fusionpoint (break_chr1, break_chr1_start, break_chr1_end, break_chr2, break_chr2_start, break_chr2_end, break_orientation, T_breaktype_breaktype_ID, T_samples_sample_ID)  
				VALUES(?,?,?,?,?,?,?,?,?)');
	
	#Skip the fusionpoint if it has already been uploaded to the FusionDB through the StrucVarDB
	my @array = $dbh->selectrow_array("SELECT fusion_ID FROM T_fusionpoint 
		INNER JOIN T_samples ON T_fusionpoint.T_samples_sample_ID=T_samples.sample_ID
		WHERE break_chr1 = $fusdata->[0] AND break_chr1_start = $fusdata->[1] AND break_chr1_end = $fusdata->[2] AND break_chr2 = $fusdata->[3] AND break_chr2_start = $fusdata->[4] AND break_chr2_end = $fusdata->[5] AND T_samples.sample_ID = '$sample_ID'");
	
	my $fusID = $array[0];
		
	unless($fusID){
		#Get the breaktype (Insertion/Deletion/Remote/etc.), corresponds with database entry
		#This is done to prevent data redundancy, all types can be clearly categorized
		my $breaktype = 9;
		$breaktype = 1 if $fusdata->[-1] =~ /Ins/i; 			#Insertion synonyms
		$breaktype = 2 if $fusdata->[-1] =~ /Del/i; 			#Deletion synonyms
		$breaktype = 3 if $fusdata->[-1] =~ /Anti/i;			#Anti synonyms
		$breaktype = 4 if $fusdata->[-1] =~ /(Inversion|Invrs)/i;	#Inversion synonyms
		$breaktype = 5 if $fusdata->[-1] =~ /(Remote|Rem)/i;		#Remote synonyms
		$breaktype = 6 if $fusdata->[-1] =~ /Trans/i; 			#Translocation
		$breaktype = 7 if $fusdata->[-1] =~ /Dup/i; 			#Duplication
		$breaktype = 8 if $fusdata->[-1] =~ /Evertion/i; 		#Evertion
		
		
		my ($orHolder,$orientation);
		$orHolder = $fusdata->[6] if scalar(@{$fusdata}) == 8;
		$orHolder = $fusdata->[7] if scalar(@{$fusdata}) == 13;
	
		$orientation = "TH" if $orHolder =~ /th/i;
		$orientation = "TT" if $orHolder =~ /tt/i;
		$orientation = "HH" if $orHolder =~ /hh/i;
		$orientation = "HT" if $orHolder =~ /ht/i;
		
		#Upload data to fusionpoint DB using the prepare statement with placeholders to handle the quoting and prevent "bad" data
		#Upload 8 column file
		if(scalar(@{$fusdata}) == 8){
			$pdh->execute($fusdata->[0], $fusdata->[1],$fusdata->[2],$fusdata->[3],$fusdata->[4], $fusdata->[5], $orientation, $breaktype, $sample_ID) || die "Could not upload the breakpoint data from the breakpoint file";	
		}
		#Upload 11 column file from 123SV
		if(scalar(@{$fusdata}) == 13){
			$pdh->execute($fusdata->[0], $fusdata->[1],$fusdata->[2],$fusdata->[3],$fusdata->[4], $fusdata->[5], $orientation, $breaktype, $sample_ID) || die "Could not upload the breakpoint data from the breakpoint file";	
		}
	}
	
	return $pdh;
}

#########################
###Check User Interest\Features###
##Mainly used by the View pages###
#########################
#These functions will check with the Ensembl DB if the provided fusionpoints have features of interest to the user

#This function will sort the fusionpoint on importance, it will place the fusionpoints with in-frame/out-frame fusion-overlapping genes and found genes and domains of user interest higher than others.
sub sortImportance{
	my ($refFeatureHash, $pointsInFusionGene, $pointsOutFusionGene, $pointsPosInFusionGene, $pointsUserGene, $pointsUserDomain, $pointsBiotype, $pointsUserRegEle, $pointsPromotorBreak, $pointsGeneBreak)= @_;

	#Array in which the sorted fusionID's are stored
	my %scoreHash;

	#Loop through the featurehash
	foreach my $fusID (keys %$refFeatureHash){
		#Add the found user-defined gens to the score
		if($refFeatureHash->{$fusID}->{"geneString"} ne "No"){
			$scoreHash{$fusID} += (($refFeatureHash->{$fusID}->{"geneString"} =~ tr/,//)*$pointsUserGene);
			#If there is only one user-gene found, there is no comma to separate but it is not empty.
			if($refFeatureHash->{$fusID}->{"geneString"} =~ tr/,// == 0){
				$scoreHash{$fusID} += $pointsUserGene;
			}
		}
		#Add the found user-defined domains to the score
		if($refFeatureHash->{$fusID}->{"domainString"} ne "No"){
			$scoreHash{$fusID} += (($refFeatureHash->{$fusID}->{"domainString"} =~ tr/,//)*$pointsUserDomain);
			if($refFeatureHash->{$fusID}->{"domainString"} =~ tr/,// == 0){
				$scoreHash{$fusID} += $pointsUserDomain;
			}
		}
		#Add the user-specified reg. ele to the score
		if($refFeatureHash->{$fusID}->{"regEleString"} ne "No"){
			$scoreHash{$fusID} += (($refFeatureHash->{$fusID}->{"regEleString"} =~ tr/,//)*$pointsUserRegEle);
			if($refFeatureHash->{$fusID}->{"regEleString"} =~ tr/,// == 0){
				$scoreHash{$fusID} += $pointsUserRegEle;
			}
		}
		#Add the fusionGenes/promotorbreaks/genebreaks and biotype to the score
		$scoreHash{$fusID} += $pointsInFusionGene*$refFeatureHash->{$fusID}->{"countFusInFrameGenes"};
		$scoreHash{$fusID} += $pointsOutFusionGene*$refFeatureHash->{$fusID}->{"countFusOutFrameGenes"};
		
		$scoreHash{$fusID} += $pointsPosInFusionGene*$refFeatureHash->{$fusID}->{"countFusPosInFrameGenes"};

		#Gene breaks
		$scoreHash{$fusID} += $pointsGeneBreak*$refFeatureHash->{$fusID}->{"countIntronicBreak"};
		$scoreHash{$fusID} += $pointsGeneBreak*$refFeatureHash->{$fusID}->{"countExonicBreak"};
		#Promotorbreaks
		$scoreHash{$fusID} += $pointsPromotorBreak*$refFeatureHash->{$fusID}->{"countPromotorBreaks"};
		#Biotype
		$scoreHash{$fusID} += $pointsBiotype*$refFeatureHash->{$fusID}->{"countBiotype"} if $refFeatureHash->{$fusID}->{"countBiotype"} ne "No";		
	}
	return %scoreHash;
}
#This function will get the features of the fusionpoints
sub AnalyzeFeatures{
	my ($dbh, $registry, $fusionpoints , $sampleInfo, $searchUp, $searchDown, $refInterestGenes, $refInterestDomains, $interestBiotype, $refInterestRegEles) = (@_);

	my ($foundRegEles, %featureHash, @fusionGeneNames);

	#Foreach fusionpoint
	foreach my $fusID (keys %$fusionpoints) {
		my $species = $sampleInfo->{$fusID}->{'sample_species'};
		#Get Slice adaptor
		my $sa = $registry->get_adaptor($species, 'Core', 'Slice');
		#Get slices for gene search purposes
		my ($breakSlice1, $breakSlice2) = getEnsemblSlices($registry, $fusionpoints->{$fusID}->{"break_chr1"}, $fusionpoints->{$fusID}->{"break_chr1_start"} ,$fusionpoints->{$fusID}->{"break_chr1_end"}, $fusionpoints->{$fusID}->{"break_chr2"}, $fusionpoints->{$fusID}->{"break_chr2_start"},$fusionpoints->{$fusID}->{"break_chr2_end"},$fusionpoints->{$fusID}->{"break_orientation"}, $searchUp, $searchDown, $species);
		#Gets a hash with the count of fusionGenes, gene-breaks without fusionproducts and promotor breaks
		%featureHash->{$fusID} = &getFusionBreakInfoWhereFusionID($dbh, $fusID);
		
		#Get an array of all the genes and domains that fit the criteria of the user.
		my ($foundGenes, $foundDomains, $countBiotype) = &checkUserInterestGenes($registry, $dbh, $breakSlice1, $breakSlice2, $refInterestGenes, $refInterestDomains, $interestBiotype);
		
		#Get an array of all the reg. eles that fit the criteria of the user if species is Human or mouse (Only species with funcgen DB)
		if($species =~ /homo_sapiens|human|mouse|mus/i){
			$foundRegEles = &checkUserInterestRegEles($registry, $dbh, $species, $breakSlice1, $breakSlice2, $refInterestRegEles);
		}
		
		#Make a text string of the found features
		my ($geneString, $domainString, $regEleString) = ("No", "No", "No");
		$geneString = join(', ', @$foundGenes) if !scalar(@$foundGenes) == 0;
		$domainString = join(', ', @$foundDomains) if !scalar(@$foundDomains) == 0;
		if($species eq "homo_sapiens"){
			$regEleString = join(', ', @$foundRegEles) if !scalar(@$foundRegEles) == 0;
		}else{
			$regEleString = "No";
		}

		$countBiotype = "No" if $countBiotype == 0;
		$countBiotype = "Yes" if $countBiotype != 0;
		
		#Save these features in the hash.		
		%featureHash->{$fusID}->{"geneString"} = $geneString;
		%featureHash->{$fusID}->{"domainString"} = $domainString;
		%featureHash->{$fusID}->{"regEleString"} = $regEleString;
		%featureHash->{$fusID}->{"countBiotype"} = $countBiotype;
		
		
	}
	return \%featureHash;
}

#Checks the vicinity of the breakpoints for used-specified reg. elements
sub checkUserInterestRegEles{
	my ($registry, $dbh, $species, $breakSlice1, $breakSlice2, $refInterestRegEles)= @_;
	#Keep track of reg. ele, so they don't come up twice
	my @foundRegEle;
	#Only if user specified reg. Elements
	if(!(scalar(@$refInterestRegEles)) == 0){
		#Make funcgen Adaptor
		my $regfeat_adaptor = $registry->get_adaptor($species, 'funcgen', 'regulatoryfeature');
		foreach my $slice($breakSlice1, $breakSlice2){
			my $regType;
			#Get all reg. elements from slice
			my @reg_feats = @{$regfeat_adaptor->fetch_all_by_Slice($slice)};
			foreach my $cell_rf (@{@reg_feats}){
				#Open the regulatory feature to get all the attributes like histone modification and pol. activity
				my $rfs = $regfeat_adaptor->fetch_all_by_stable_ID($cell_rf->stable_id); 
				foreach my $regEle(@{$rfs}){
					my $regex = "(".join("|", @$refInterestRegEles).")";

					#Get all the attributes of the reg. ele
					foreach my $attr_feat (@{$regEle->regulatory_attributes()}){
						if($attr_feat->display_label =~ /$regex/i){
							my @split = split(" -", $attr_feat->display_label);
							#Add the reg. ele if it hasn't already been added before
							push(@foundRegEle, $split[0]) if (grep { $_ eq $split[0] } @foundRegEle);
						}
					}
				}
			}
		}	
	}
	return \@foundRegEle;
}

#Checks the found genes and domains against the user genes and domains of interest
sub checkUserInterestGenes{
	my ($registry, $dbh, $breakSlice1, $breakSlice2, $refInterestGenes, $refInterestDomains, $interestBiotype) = @_;
	#Arrays to keep track of found genes and domains of interest, which are used to later show the user what was found
	my (@foundGenes, @foundDomains,$regex, $countBiotype);

	if(!(scalar(@$refInterestGenes)) == 0 or !(scalar(@$refInterestDomains)) == 0 or $interestBiotype ne "none"){
		foreach my $slice ($breakSlice1, $breakSlice2){
			my $genesOnSlice = $slice->get_all_Genes();
			#Check if there are genes, if yes -> concenate on string
			foreach my $gene (@{$genesOnSlice}){
				my $geneAdaptor = $registry->get_adaptor('Human', 'Core', 'Gene');
				#Check if a transcript adaptor is needed, only when the user specified domains of interest
				my $transcriptAdaptor = $registry->get_adaptor( 'Human', 'Core', 'Transcript' ) if scalar(@$refInterestDomains) != 0;
				#Add the gene to an array to prevent data redundancy of the same gene on both slices
				push(@foundGenes, $gene->external_name()) if grep { $_ eq $gene->external_name() } @$refInterestGenes and !(grep { $_ eq $gene->external_name() } @foundGenes);
	
				#Biotype filtering, checks if the gene is of the biotype the user has specified
				if($interestBiotype ne "None"){
					$countBiotype++ if $gene->biotype() eq $interestBiotype;				
				}
					
				#Check user-defined genes of interest
				if(!(scalar(@$refInterestGenes)) == 0){
					$regex = "(".join("|", @$refInterestGenes).")";
					#Search through all synonyms of the gene
					my $db_entry_adaptor = $registry->get_adaptor( 'Human', 'Core', 'DBEntry' );
					my @db_entries = @{ $db_entry_adaptor->fetch_all_by_Gene($gene) };
					foreach my $dbentry( @db_entries){
						foreach my $syn (@{$dbentry->get_all_synonyms}) {
							#Check synonym, only save it once			
							push(@foundGenes, $syn) if $syn =~ /$regex/i and !(grep { $_ eq $syn } @foundGenes);
						}  
					}
				}
				#Domain search
				if(!(scalar(@$refInterestDomains)) == 0){
					$regex = "(".join("|", @$refInterestDomains).")";
					while( my $transcript = shift @{ $gene->get_all_Transcripts() }){
						my $domainfeat = $transcriptAdaptor->fetch_by_stable_id($transcript->stable_id())->translation()->get_all_DomainFeatures if $transcript->biotype eq "protein_coding";
						while ( my $pfeature = shift @{$domainfeat} ) {
							my $domain = $pfeature->idesc();
							push(@foundDomains, $domain) if $domain  =~ /$regex/i and !(grep { $_ eq $domain } @foundDomains);
						}
					}
				}
			}
		}
	}
	
	#Return all the found genes and domains of interest
	return \@foundGenes, \@foundDomains, $countBiotype;
}

########################################
###General functions for FusionAnnotator#######
########################################

#Inserts or deletes support files per fusionpoint or sample as given as parameter by detailOverview.pl
sub changeSupportFiles{
	my ($dbh, $fusID, $sampleID, $supportFiles) = @_;	
	my (@splitLine, $specific, $typeName, $typeID, $link, $pdh);
	#Delete ALL the supportfiles for the given fusionpoint and ALL the sample-specific supportfiles for that sample (Except for the StrucVarDB samples ( <= 4)
	if($sampleID > 4){
		$dbh->do("DELETE FROM T_support WHERE T_samples_sample_ID = $sampleID OR T_fusionpoint_fusion_ID = $fusID");
	}
	if($sampleID < 4){
		$dbh->do("DELETE FROM T_support WHERE T_fusionpoint_fusion_ID = $fusID");
	}
		
	#Insert the supporfiles (Inserting them to the specific fusionpoint or the whole sample)
	foreach(@{$supportFiles}){
		@splitLine = split(/\|\|/, $_);

		$specific = @splitLine[0];
		$typeName = @splitLine[1];
		$link = @splitLine[2];
		
		$typeID = 1 if $typeName =~ /(bam$)/i;
		$typeID = 2 if $typeName =~ /(bed$)/i;
		$typeID = 3 if $typeName =~ /(wig$)/i;
		$typeID = 4 if $typeName =~ /(gff$)/i;
		
		#Add support for all fusionpoint from the same sample
		if($specific eq "sample"){
			$pdh = $dbh->prepare("INSERT INTO T_support (support_link, T_samples_sample_ID, T_supType_supType_ID) VALUES(?,?,?)");
			$pdh->execute($link, $sampleID, $typeID);
			$pdh->finish;
		}
		#Add support file for the selected breakpoint only
		else{
			$pdh = $dbh->prepare("INSERT INTO T_support (support_link, T_fusionpoint_fusion_ID, T_supType_supType_ID) VALUES(?,?,?)") || die "Cant upload Support file";
			$pdh->execute($link, $fusID, $typeID);
			$pdh->finish;
		}
	}
	$dbh->commit();
}


#Gets the genes stored in the fusion DB for a given fusionpoint
sub getGenesFromFusionDB{
	return $_[0]->selectall_hashref("SELECT feat_ENS_ID
		FROM T_features WHERE  T_fusionpoint_fusion_ID  = $_[1] AND  T_featuretype_featuretype_ID < '6'", "feat_ENS_ID") || die "Could not get gene features for fusionpoint $_[1] from DB";
}

#Returns a hash of information about the products of the fusionpoint
sub getFusionBreakInfoWhereFusionID{
	my ($dbh, $fusID) = @_;
	my (%infoHash, @fusionGeneNames, @fusionGeneNamesHolder, @geneBreakNamesHolder, @geneBreakNames);
	
	%infoHash->{$fusID}->{'countFusInFrameGenes'} = ($dbh->selectrow_array("SELECT COUNT(feat_ID) FROM T_features WHERE T_fusionpoint_fusion_ID = '$fusID' AND T_featuretype_featuretype_ID = 1"))/2;
	%infoHash->{$fusID}->{'countFusOutFrameGenes'} = ($dbh->selectrow_array("SELECT COUNT(feat_ID) FROM T_features WHERE T_fusionpoint_fusion_ID = '$fusID' AND T_featuretype_featuretype_ID = 2"))/2;
	%infoHash->{$fusID}->{'countFusPosInFrameGenes'} = ($dbh->selectrow_array("SELECT COUNT(feat_ID) FROM T_features WHERE T_fusionpoint_fusion_ID = '$fusID' AND (T_featuretype_featuretype_ID = 3 OR T_featuretype_featuretype_ID = 4)"))/2;
	%infoHash->{$fusID}->{'countIntronicBreak'} = $dbh->selectrow_array("SELECT COUNT(feat_ID) FROM T_features WHERE T_fusionpoint_fusion_ID = '$fusID' AND T_featuretype_featuretype_ID = 5");
	%infoHash->{$fusID}->{'countExonicBreak'} = $dbh->selectrow_array("SELECT COUNT(feat_ID) FROM T_features WHERE T_fusionpoint_fusion_ID = '$fusID' AND T_featuretype_featuretype_ID = 6");
	%infoHash->{$fusID}->{'countPromotorBreaks'} = $dbh->selectrow_array("SELECT COUNT(feat_ID) FROM T_features WHERE T_fusionpoint_fusion_ID = '$fusID' AND T_featuretype_featuretype_ID = 7");
	
	#Get all the fusionnames for the specific breakpoint from the database
	my $fusionGeneNames = $dbh->selectall_arrayref("SELECT DISTINCT featName FROM T_features WHERE T_fusionpoint_fusion_ID = '$fusID' AND T_featuretype_featuretype_ID <5");
	
	foreach(@{$fusionGeneNames}){
		push(@fusionGeneNamesHolder, $_->[0]);
	}
	
	%infoHash->{$fusID}->{'fusionGeneNames'} = join(", ", @fusionGeneNamesHolder);
	if(scalar(@fusionGeneNamesHolder) == 0){
		%infoHash->{$fusID}->{'fusionGeneNames'} = "None";
	}
	
	#Get all the geneNames for geneBreaks for the specific breakpoint from the database
	my $geneBreakNames = $dbh->selectall_arrayref("SELECT DISTINCT featName FROM T_features WHERE T_fusionpoint_fusion_ID = '$fusID' AND (T_featuretype_featuretype_ID = 5 OR T_featuretype_featuretype_ID = 6)");
	
	foreach(@{$geneBreakNames}){
		push(@geneBreakNamesHolder, $_->[0]);
	}
	
	%infoHash->{$fusID}->{'geneBreakNames'} = join(", ", @geneBreakNamesHolder);
	if(scalar(@geneBreakNamesHolder) == 0){
		%infoHash->{$fusID}->{'geneBreakNames'} = "None";
	}
	
	return %infoHash->{$fusID};
}

#Gets the fusion genes stored in the fusion DB for a given fusionpoint based on the fusionpoint_ID
sub getFusionGenesWhereFusionID{
	return $_[0]->selectall_hashref("SELECT feat_ENS_ID,  T_featuretype_featuretype_ID 
		FROM T_features WHERE  T_fusionpoint_fusion_ID  = $_[1] AND T_featuretype_featuretype_ID >=2 AND T_featuretype_featuretype_ID <=3", "feat_ENS_ID") || die "Could not get gene features for fusionpoint $_[1] from DB";
}

#Gets the breakpoint coordinates from the fusion DB based on the given sample_ID
sub getBreakpointDataWhereSampleID{
	return $_[0]->selectall_hashref("SELECT fusion_ID, break_chr1, break_chr1_start, break_chr1_end, break_chr2, break_chr2_start, break_chr2_end, break_orientation, T_breaktype_breaktype_ID, sample_species FROM T_fusionpoint
		INNER JOIN T_samples ON T_fusionpoint.T_samples_sample_ID=T_samples.sample_ID
		WHERE T_samples_sample_ID = $_[1]","fusion_ID") || die "Could not get fusionpoint data from DB";
}

#Gets the breakpoint coordinates from the fusion DB based on the given sample_ID
sub getBreakpointDataWhereFusionID{
	return $_[0]->selectall_hashref("SELECT fusion_ID, break_chr1, break_chr1_start, break_chr1_end, break_chr2, break_chr2_start, break_chr2_end, break_orientation, T_breaktype_breaktype_ID, sample_species FROM T_fusionpoint
		INNER JOIN T_samples ON T_fusionpoint.T_samples_sample_ID=T_samples.sample_ID
		WHERE fusion_ID = $_[1]","fusion_ID") || die "Could not get fusionpoint data from DB";
}

#Gets the breakpoint coordinates from the fusion DB based on the given fusionpoint
sub getFusionPointsWhereFusionID{
	return $_[0]->selectall_hashref("SELECT fusion_ID, break_chr1,  break_chr1_start, break_chr1_end, break_chr2, break_chr2_start, break_chr2_end, break_orientation, T_breaktype_breaktype_ID FROM T_fusionpoint WHERE fusion_ID = $_[1]", "fusion_ID") || die "Could not get fusionpoint from DB";
}

#Gets sample Info based on sample_ID
sub getSampleInfoWhereSampleID{
	return $_[0]->selectall_hashref("SELECT DISTINCT sample_ID, sample_species, sample_name FROM T_samples 
          WHERE sample_ID = $_[1]", "sample_ID") || die "Could not get sample info from DB";
}

#Gets sample Info based on sample_ID
sub getSampleNames{
	return $_[0]->selectall_hashref("SELECT DISTINCT sample_name FROM T_samples", "sample_name") || die "Could not get sample info from DB";
}
#Get the sample information using a fusion ID
sub getSampleInfoWhereFusionID{
	#Get the sample ID and species
	my $sampleInfo = $_[0]->selectall_hashref("SELECT DISTINCT sample_ID, sample_species, sample_name, fusion_ID FROM T_fusionpoint  
        INNER JOIN T_samples ON T_fusionpoint.T_samples_sample_ID=T_samples.sample_ID
          WHERE fusion_ID = $_[1]", "fusion_ID") || die "Could not get sample info from DB";
        
        #Get the support files
        #Sample-Specific
        my $supportFiles = $_[0]->selectall_arrayref("SELECT support_link FROM T_support WHERE T_samples_sample_ID = ".$sampleInfo->{$_[1]}->{'sample_ID'});
        $sampleInfo->{$_[1]}->{"supportFilesSample"} = $supportFiles;
        #Fusionpoint-specific
        my $supportFiles = $_[0]->selectall_arrayref("SELECT support_link FROM T_support WHERE T_fusionpoint_fusion_ID = $_[1]");
        $sampleInfo->{$_[1]}->{"supportFilesFusion"} = $supportFiles;

        
        #Return a hashref of the sample information
        return $sampleInfo;
}

#Opens connection to the fusionpoint DB to facilitate the upload of data
sub openFusionConnection{
	#Local fusion DB
	return DBI->connect("dbi:mysql:FusionAnnotator:localhost:3306", "*", "*", { AutoCommit => 0 }) || die "Could not connect to fusionpoint database: $DBI::errstr \;";
}

#Close succesfull transaction and close fusion DB connection
sub closeFusionConnection{
	$_[0]->disconnect() || die "Could not close connection";
}

#Close succesfull transaction and close fusion DB connection
sub closeFusionConnectionWithCommit{
	$_[0]->commit() || die "Could not commit data";
	$_[0]->disconnect() || die "Could not close connection";
	print "Data committed & database connection closed\n";
}

#Closes the connection to the fusionpoint DB, finish the prepare statement and commits the data.
sub closeFusionConnectionWithFinishAndCommit{
	my ($dbh, $pdh) = (@_);
	$pdh->finish()|| die "Could not close PREPARE statement";
	$dbh->commit() || die "Could not commit breakfile data";
	$dbh->disconnect() || die "Could not close connection to FusionDB";
}

#This function will get the fusion data from the local fusion DB depending on the user specified search parameters on the ViewExisting page
#It return an feature hash of the fusionpoints with the score
sub getFusionpointsBasedOnSearchParameters{
	my ($dbh, $registry, $species, $searchUp, $searchDown, $fusID, $sampleID, $sampleName, $fusionGenesCheck, $confirmedCheck, $geneBreakCheck, $promotorCheck, $refInterestGenes, $refInterestDomains, $refInterestRegEles, $interestBiotype) = @_;
	my ($fusionpoints,@searchParam,$featureHash);
	#Add the search parameters upon the standard select all query
	my $sqlQuery = "SELECT DISTINCT sample_species, sample_name,  fusion_ID, break_chr1, break_chr1_start, break_chr1_end, break_chr2, break_chr2_start, break_chr2_end, break_orientation, T_breaktype_breaktype_ID FROM T_fusionpoint  
        INNER JOIN T_samples ON T_fusionpoint.T_samples_sample_ID=T_samples.sample_ID";
     	if($fusID ne ''){
     		push(@searchParam, "fusion_ID = '$fusID'") if (split(";",$fusID)) == 1;
     		push(@searchParam, "fusion_ID BETWEEN ".(split(';',$fusID))[0]. " AND ".(split(';',$fusID))[1]) if (split(";",$fusID) != 1);		
     	}
     	if($sampleID ne ''){
     		push(@searchParam, "T_samples_sample_ID = '$sampleID'") if (split(";",$sampleID)) == 1;
     		push(@searchParam, "T_samples_sample_ID BETWEEN ".(split(';',$sampleID))[0]. " AND ".(split(';',$sampleID))[1]) if (split(";",$sampleID) != 1);		
     	}
	push(@searchParam, "T_samples.sample_name LIKE '$sampleName\%'") if $sampleName ne '';
	push(@searchParam, "break_chr1_start = break_chr1_end") if $confirmedCheck == 1;
	push(@searchParam, "T_featuretype_featuretype_ID BETWEEN 1 AND 4") if $fusionGenesCheck == 1;
	push(@searchParam, "T_featuretype_featuretype_ID = 7") if $promotorCheck == 1;
	push(@searchParam, "T_featuretype_featuretype_ID BETWEEN 5 AND 6") if $geneBreakCheck == 1;
	
	#Add the T_Feature to the select query through INNER JOIN and thus eliminate all non-feature fusionpoints
	if($fusionGenesCheck == 1 or $promotorCheck == 1 or $geneBreakCheck == 1){
		$sqlQuery = "SELECT DISTINCT sample_species, fusion_ID, break_chr1, break_chr1_start, break_chr1_end, break_chr2, break_chr2_start, break_chr2_end, break_orientation, T_breaktype_breaktype_ID, T_featuretype_featuretype_ID FROM T_fusionpoint  
		INNER JOIN T_samples ON T_fusionpoint.T_samples_sample_ID=T_samples.sample_ID
		INNER JOIN T_features ON T_fusionpoint.fusion_ID =  T_fusionpoint_fusion_ID";
	}
	
	push(@searchParam, "sample_species = \"$species\"") if $species ne "none";
	
	$sqlQuery = join("",$sqlQuery," WHERE ", join(" AND ",@searchParam)) if (scalar(@searchParam)) != 0;
	
	#Print the sqlQuery for debugging
	#print $sqlQuery;
	
	#Get the fusionpoints based on the search parameters
	$fusionpoints = $dbh->selectall_hashref($sqlQuery, "fusion_ID");
	
	#Get the features
	$featureHash = &AnalyzeFeatures($dbh, $registry, $fusionpoints , $fusionpoints, $searchUp, $searchDown, $refInterestGenes, $refInterestDomains, $interestBiotype, $refInterestRegEles);
	#Return 
	return $fusionpoints, $featureHash;
}

##############################
###General functions for the Ensembl DB###
##############################
#Open connection (Registry) to Ensembl DB
sub openEnsemblConnection{
	my $registry = 'Bio::EnsEMBL::Registry';
	#Local Ensembl DB
	$registry->load_registry_from_db(
    -host => 'wgs10.op.umcutrecht.nl',
    -user => '*') || die "Could not open ensembl registry connection to wgs10 using ensembl user";
    	return $registry;
}

#Open DB Adaptor to Ensembl DB to get the karyotype bands for the centromere
sub openDBAdaptor{
	return Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -user   => '*',
    -host   => 'wgs10.op.umcutrecht.nl',
    -dbname => 'homo_sapiens_core_67_37',
    -driver => 'mysql');
}
	##########
###End of Functions###
	##########
	
#Need to close module with a true value for validation, thus 1;
1;
