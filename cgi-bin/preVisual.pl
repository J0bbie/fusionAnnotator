#!/usr/bin/perl
#Author: 				Job van Riet
#Date of creation:		13/10/12
#Data of modification:	13/10/12
#Known bugs:			None
#Function:			This page will make the derivative chromosomes by splitting the nuc. seq fasta files and derivative genomic features by querying the EnsEMBL DB
#					It will then give the user a button to be redirected to the visualization browser with the correct view.

use lib '../modules/';

use strict;
#Diagnostics for debugging
use diagnostics;
#Use fusionModule for functions
use fusionModule;
#Use CGI
use CGI;
use CGI::Carp qw(fatalsToBrowser);    # Remove for production use

#Get the fusionID from the form
my $q = new CGI;
my $fusID = $q->param('fusID');
my $searchUp = $q->param('searchUp');
my $searchDown = $q->param('searchDown');

print  "Content-type:text/html\r\n\r\n";

#Make the derivative chromosome with all its features
my ($chr1, $chr2, $breakChrA, $breakChrB, $breakPos1DerCoord, $derName, $species, $breakOrientation, $fusionFolderName, $confirmedBreakpoint) = fusionModule::createFusionFolder($fusID, $searchUp, $searchDown);


#Give the user a button to go to the correct view
print '
<HTML>
<HEAD>
<title>Pre-Visualization of Fusionpoint</title>
<LINK href="VisualBrowserLayout.css" rel="stylesheet" type="text/css">
</HEAD>
<BODY>
<div class="options" align="center" style="color: #ff7f00">
	<p>Derivative chromosome, BAM linkage and genomic features are made, click the button to go to the genome browser</p>
	<form action="fusBrowser.pl" method="post" enctype="multipart/form-data">
			<input type="hidden" name="fusID" value='.$fusID.' />
			<input type="hidden" name="chr1" value='.$chr1.' />
			<input type="hidden" name="chr2" value='.$chr2.' />
			<input type="hidden" name="breakChrA" value='.$breakChrA.' />
			<input type="hidden" name="breakChrB" value='.$breakChrB.' />
			<input type="hidden" name="breakPos1DerCoord" value='.$breakPos1DerCoord.' />
			<input type="hidden" name="derName" value='.$derName.' />
			<input type="hidden" name="species" value='.$species.' />
			<input type="hidden" name="breakOrientation" value='.$breakOrientation.' />
			<input type="hidden" name="confirmedBreakpoint" value='.$confirmedBreakpoint.' />
			<input type="hidden" name="fusionFolderName" value='.$fusionFolderName.' />
			<br><input type="submit" name="submit" value="Go to genome browser"/><br>
	
	</form>
</div>
</BODY>
</HTML>';
