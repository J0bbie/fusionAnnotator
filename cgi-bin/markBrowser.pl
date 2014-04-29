#!/usr/bin/perl
#Author: 		Job van Riet
#Date of creation:	22/11/12
#Data of modification:	22/11/12
#Known bugs:		None
#Function:		Uploads multiple .BED files to a special JBrowse instance for use in Mark's script

use strict;
#Diagnostics for debugging
use diagnostics;
#Use CGI and other HTML modules
use CGI;
use CGI::Carp qw(fatalsToBrowser);    # Remove for production use
use File::Basename;
use HTML::TokeParser;
use LWP::Simple;
use File::Temp;	

#Get the the path of the folder where the .BED files are that are to be displayed
my $q = new CGI;
my ($folderPath, $species) = ($q->param('folderPath'),$q->param('species'));
#Set content type and print HTML
print  "Content-type:text/html\r\n\r\n";
#Test variables
#$folderPath = "http://genetics.genomicscenter.nl/HomozygosityMapping/results/MarkTest_p2_r3_h20_m0_lNo_dNo_cNo_aNo_s1000000/";
#$species = "homo_sapiens";

#Makes the data-folder where the data is stored
my $homoZygFolder = &makeHomoZygfolder($species);

#Make tracks of all the .BED files in the specified folder, get a string of all the .BED tracks
my $BEDtracks = &makeBEDTracks($homoZygFolder, $folderPath);

#Delete .htacces file and tempBED
system("rm $homoZygFolder/.htaccess");
system("rm $homoZygFolder/tempBED.BED");


print '
<HTML>
<HEAD>
	<title>Visualization of .BED Files</title>
	<LINK href="../VisualBrowserLayout.css" rel="stylesheet" type="text/css">
</HEAD>
	<BODY>
		<!--Genome with the .BED data displayed-->
		<p>Visualization of .BED files</p>
		<i>Patients are <FONT COLOR="#FF0000">red</FONT>, references are <FONT COLOR="#F7FE2E">yellow</FONT> and overlap is <FONT COLOR="#8904B1">purple</FONT>!</i>
		<iframe align="center" style="border: 1px solid black" src="../genomeBrowser/fusionBrowser.html?&data='.$homoZygFolder.'&tracklist=1&nav=1&overview=1&tracks=DNA%2CGenes%2C'.$BEDtracks.'" width="90%" height="90%"></iframe>
		<br>
	</body>
</html>
';

#Gets all the .bed files in a folder and then makes them a track in the data/mark/data folder for viewing in JBrowse
sub makeBEDTracks{
	my ($homoZygFolder, $folderPath) = @_;
	my (@trackNames, $BEDFileHandler, $tempFileHandler);
	#Check if file is the overlap file (Isn't mentioned specifically)
	
	#Clean up the trackList to remove old tracks
	&cleanTrackList($homoZygFolder);
	
	#Reads the settings.txt file to see what files correspond to what type of patients/references
	my ($patients, $references) = &getPatientsAndRefs($folderPath);
	
	#Order the tracks, first overlap then patients and then reference
	foreach(@{$references}){
		push(@trackNames, $_);
	}
	foreach(@{$patients}){
		push(@trackNames, $_);
	}
	
	#Read the HTML of the page/directory supplied
	my $page = get($folderPath);
	my $p = HTML::TokeParser->new(\$page);
	
	#Iterate through files in directory through the use of the <a> tag, only act upon .bed files
	while (my $file = $p->get_tag("a")) {
		my $text = $p->get_trimmed_text("/a");
		
		#Split the link on filename/dir and extension
		my ($fileName,$dir,$ext) = fileparse($text, qr/\.[^.]*/);
		#Only act upon .BED files
		if ($ext =~ /(.bed)/i){
			my $overlapCheck = 0;

			#Get the text of the .BED file
			my $BEDFileText = get($folderPath.$text) || die("Could'nt get $folderPath.$text");
			#Store the text in a temporary local file
			open($tempFileHandler, ">", "$homoZygFolder/tempBED.BED") || die("Could not open temp BED file in $homoZygFolder");
			print $tempFileHandler $BEDFileText;
			close($tempFileHandler);			
			
			#Add the BED file as a track to JBrowse, color patients and references differently
			if(grep { $_ eq $fileName } @{$patients}){
				system("perl ../genomeBrowser/bin/flatfile-to-json.pl --bed $homoZygFolder/tempBED.BED --tracklabel '$fileName' --key '$fileName' --out $homoZygFolder --cssClass basic --clientConfig '{ \"featureCss\": \"background-color:#FF0000;border-color:#FF0000\"}'");
				$overlapCheck = 1;
			}
			if(grep { $_ eq $fileName } @{$references}){
				system("perl ../genomeBrowser/bin/flatfile-to-json.pl --bed $homoZygFolder/tempBED.BED --tracklabel '$fileName' --key '$fileName' --out $homoZygFolder --cssClass basic --clientConfig '{ \"featureCss\": \"background-color:#F7FE2E;border-color:#F7FE2E\"}'");
				$overlapCheck = 1;
			}
			#If not patient nor reference then overlap
			system("perl ../genomeBrowser/bin/flatfile-to-json.pl --bed $homoZygFolder/tempBED.BED --tracklabel '$fileName' --key '$fileName' --out $homoZygFolder --cssClass basic --clientConfig '{ \"featureCss\": \"background-color:#8904B1;border-color:#8904B1\"}'") if $overlapCheck == 0;
			
			#Add the overlap tracks to the tracklist
			push(@trackNames, $fileName) if $overlapCheck == 0;
		}
	}
	#Return a string of the tracks that will be displayed
	return join("\%2C", @trackNames);
}

#Reads the settings.txt file to see what files belong to the patients or reference, this is later used to color the tracks accordingly
sub getPatientsAndRefs{
	my $folderPath = shift;
	my (@patients, @references);
	#Get the settings text and split on newline
	my @splitText = split("\n", get($folderPath."/settings.txt"));
	#Find the lines with Patients (If not a patient than a reference)
	foreach(@splitText){
		#Split on /t
		my $pat = (split("\t", $_))[1] if $_ =~ /Patient/;
		my $ref = (split("\t", $_))[1] if $_ =~ /Reference/;

		#Delete extension
		$pat =~ s/\..*// if $_ =~ /Patient/;
		$ref =~ s/\..*// if $_ =~ /Reference/;

		#Push in stack
		push(@patients, $pat)  if $_ =~ /Patient/;
		push(@references, $ref)  if $_ =~ /Reference/;

	}
	return \@patients, \@references;
}

#This function makes a new folder for the selected homozygosity view, also symlinks the tracks to that folder from the homo_sapiens/ref folder where the tracks for the  reference chromosomes are stored
sub makeHomoZygfolder{
	my $species = shift;
	
	#Get the  correct species
	$species = "homo_sapiens" if $species =~ /homo|human/i;         
	#Sample ID for StrucVarDB samples(Danio Rerio)
	$species = "danio_rerio" if $species =~ /danio|zebra/i;
	#Sample ID for StrucVarDB samples(Danio Rerio)
	$species = "rattus_norvegicus" if $species =~ /rattus|rat/i;
	#Sample ID for StrucVarDB samples(Danio Rerio)
	$species = "c_elegans" if $species =~ /elegans/i;
		
	#Creates a new folder folder to store the data
	my $tmpFolder = File::Temp->newdir(
        DIR      => "../data/mark/",
        UNLINK => 0,
        CLEANUP => 0
        );
    	my $homoZygFolder = $tmpFolder->dirname;
    	
    	#Folder where the formatted nuc. seq are stored
    	system("mkdir $homoZygFolder/seq");
    	#Folder where the formatted tracks are stores
    	system("mkdir $homoZygFolder/tracks");
    	#Folder where the raw/unformatted datafiles are stores before being formatted.
    	system("mkdir $homoZygFolder/raw");
    	
    	#Symlink the reference nuc. seqs to this folder, do the same for the genes, transcript and domain track
    	opendir(D, "../data/$species/ref/seq") || die "Can't open refSeq $species dir: $!\n";
	while (my $f = readdir(D)) {
		system("ln -s /var/www/html/FusionAnnotator/data/$species/ref/seq/$f $homoZygFolder/seq");
	}
	closedir(D);
	
	#Symlink all the tracks
    	opendir(D, "../data/$species/ref/tracks") || die "Can't open refTrack dir: $!\n";
	while (my $f = readdir(D)) {
		system("ln -s /var/www/html/FusionAnnotator/data/$species/ref/tracks/$f $homoZygFolder/tracks");
	}
	closedir(D);
	
	#Delete folders older than 1 day
	system("find /var/www/html/FusionAnnotator/data/mark/* -mtime +1 -exec rm {} \;");
	
    	
    	#Return the name of the folder
    	return $homoZygFolder;
}


#Remake tracklist.json to clean up old .BED tracks
sub cleanTrackList{
	my $homoZygFolder = shift;
	open(my $trackFileHandler, ">", "$homoZygFolder/trackList.json") || die("Could'nt open $homoZygFolder/trackList.json for editing");
	print $trackFileHandler '{
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
            "className":"feature",
            "subfeatureClasses" : {
            	"Exon" : "exon"
            }
         },
         "key":"Genes",
         "urlTemplate":"tracks/Genes/{refseq}/trackData.json",
         "phase":null,
         "compress":0,
         "type":"FeatureTrack",
         "label":"Genes",
         "subfeatures":null
      },
      {
         "style":{
            "className":"feature"
         },
         "key":"Transcripts",
         "urlTemplate":"tracks/Transcripts/{refseq}/trackData.json",
         "phase":null,
         "compress":0,
         "type":"FeatureTrack",
         "label":"Transcripts",
         "subfeatures":null
      }
   ],
   "formatVersion":1
}';
	close($trackFileHandler);
}
