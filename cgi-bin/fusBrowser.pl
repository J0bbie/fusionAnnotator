#!/usr/bin/perl
#Author: 		Job van Riet
#Date of creation:	2/10/12
#Data of modification:	2/10/12
#Known bugs:		None
#Function:		This page features 3 Jbrowse browsers, one for breakpoint 1 and another for breakpoint 2 and the bottom module depicts
#			the fusionpoints (The end of break 1 is "welded" onto the end of break 2). The user is able to add/remove or edit the tracks shown
#			in the modules. The fusionpoint presented is chosen by the user by clicking the
#			"View" button.

use lib '../modules/';
use strict;
#Diagnostics for debugging
use diagnostics;
#Use fusionModule for functions
use fusionModule;
#Use CGI
use CGI;
use CGI::Carp qw(fatalsToBrowser);    # Remove for production use

#Get the fusion parameters from the form
my $q = new CGI;
my ($fusID, $chr1, $chr2, $breakChrA, $breakChrB, $breakPos1DerCoord, $derName, $species, $breakOrientation, $fusionFolderName, $confirmedBreakpoint) = ($q->param('fusID'),$q->param('chr1'),$q->param('chr2'), $q->param('breakChrA'), $q->param('breakChrB'), $q->param('breakPos1DerCoord'), $q->param('derName'), $q->param('species'), $q->param('breakOrientation'), $q->param('fusionFolderName'), $q->param('confirmedBreakpoint'));
my ($o1, $o2);
if ($breakOrientation =~ /^(H|T)(H|T)/i) {
	  $o1 = uc($1);
	  $o2 = uc($2);
}
#Set content type and print HTML
print  "Content-type:text/html\r\n\r\n";
print '
<HTML>
<HEAD>
<title>Visualization of Fusionpoint</title>
<LINK href="../VisualBrowserLayout.css" rel="stylesheet" type="text/css">
</HEAD>
<BODY>
	<!--Origal genomes-->
	<table width="100%"><tr><td align="left"><b>Original Reference 1</b><br><i>'.$chr1.';'.$breakChrA.'</i> ('.$o1.') </td><td align="center"><i>Visualization of fusion ID: '.$fusID.'</i></td><td align="right"><b>Original Reference 2</b><br><i>'.$chr2.';'.$breakChrB.'</i> ('.$o2.') </td></tr></table>
	<br>
	<iframe align="left" style="border: 1px solid black" src="../genomeBrowser/fusionBrowser.html?loc='.$chr1.':'.($breakChrA-10000).'..'.($breakChrA +10000).'&data='.$fusionFolderName.'/ref&tracklist=1&nav=1&overview=1&tracks=DNA%2CBreakpointsInVicinity%2CGeneBreakA" width="49%" height="25%"></iframe>
	<iframe align="right" style="border: 1px solid black" src="../genomeBrowser/fusionBrowser.html?loc='.$chr2.':'.($breakChrB-10000).'..'.($breakChrB +10000).'&data='.$fusionFolderName.'/ref&tracklist=1&nav=1&overview=1&tracks=DNA%2CBreakpointsInVicinity%2CGeneBreakB" width="49%" height="25%"></iframe>
	<!--Deritative genome-->
	<p align ="center"><b>Derivitave chromosome </b><br><i>'.$derName.'</i> ('.$breakOrientation.') <br> Breakpoint status: '.$confirmedBreakpoint.'</p>
	<iframe style="border: 1px solid black" src="../genomeBrowser/fusionBrowser.html?loc='.$derName.':'.($breakPos1DerCoord - 10000).'..'.($breakPos1DerCoord +10000).'&data='.$fusionFolderName.'/der/&tracklist=1&nav=1&overview=1&tracks=DNA%2CBreakpointsInVicinity%2CGenes_Der" width="100%" height="55%"></iframe>
	<br>
	<form align="center" action="fusBrowser.pl" method="post">
	<input type="hidden" name="fusID" value='.$fusID.' />
	<input type="hidden" name="chr1" value='.$chr1.' />
	<input type="hidden" name="chr2" value='.$chr2.' />
	<input type="hidden" name="breakChrA" value='.$breakChrA.' />
	<input type="hidden" name="breakChrB" value='.$breakChrB.' />
	<input type="hidden" name="breakPos1DerCoord" value='.$breakPos1DerCoord.' />
	<input type="hidden" name="derName" value='.$derName.' />
	<input type="hidden" name="species" value='.$species.' />
	<input type="hidden" name="breakOrientation" value='.$breakOrientation.' />
	<input type="hidden" name="fusionFolderName" value='.$fusionFolderName.' />
	<input type="hidden" name="confirmedBreakpoint" value='.$confirmedBreakpoint.' />
	<input type= "submit" name="submit" value="Reset View"/></form>
	</body>
</html>
';

