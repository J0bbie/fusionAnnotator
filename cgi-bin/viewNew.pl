#!/usr/bin/perl
#Author:					Job van Riet
#Date of creation:			22/10/12
#Last date of modification:	22/10/12
#Modifications:

#Known bugs:			None
#Function:			This script will get the newly-uploaded fusionpoints from the fusion DB based on a sample ID.
#					It will display the features found in a scrollable table, the fusionpoints are ordered on user defined importance
#					The domains, genes and regulatory elements the user has defined will be looked up in these fusionpoints and are also listed.
#					This scripts requires the sampleID, and the list of user defined features of importance.
#					After every fusionpoint, there is a possibility to 'view' the fusionpoint in the genome browser or to open a new screen/tab where additional 'info' can be found about the selected fusionpoint.

use lib '../modules/';

#Use strict programming
use strict;
#Use DBI for database connectivity
use fusionModule;
#Diagnostics for debugging
use diagnostics;
#Use CGI
use CGI;
use CGI::Carp qw(fatalsToBrowser);    # Remove for production use

### Main script ###

#Get these variables from the fusion01.pl
my $q = new CGI;

my (@interestGenes, @interestDomains, @interestRegEles);
if($q->param('interestGenes') ne "/"){
	#If only 1 gene
	if(split(/;/,$q->param('interestGenes')) == 0){
		push (@interestGenes, $q->param('interestGenes'));
	}
	else{
		@interestGenes = split(/;/,$q->param('interestGenes'));
	}
}

if($q->param('interestDomains') ne "/"){
	#If only 1 domains
	if(split(/;/,$q->param('interestDomains')) == 0){
		push (@interestDomains, $q->param('interestDomains'));
	}
	else{
		@interestDomains = split(/;/,$q->param('interestDomains'));
	}
}

if($q->param('interestRegEles') ne "/"){
	#If only 1 regEles
	if(split(/;/,$q->param('interestRegEles')) == 0){
		push (@interestRegEles, $q->param('interestRegEles'));
	}
	else{
		@interestRegEles = split(/;/,$q->param('interestRegEles'));
	}
}
	
my $pointsInFusionGene = $q->param('pointsInFusionGene');
my $pointsOutFusionGene = $q->param('pointsOutFusionGene');
my $pointsPosInFusionGene = $q->param('pointsPosInFusionGene');

my $pointsUserGene = $q->param('pointsUserGene');
my $pointsPromotorBreak = $q->param('pointsPromotorBreak');
my $pointsGeneBreak = $q->param('pointsGeneBreak');
my $pointsUserDomain = $q->param('pointsUserDomain');
my $pointsUserRegEle = $q->param('pointsUserRegEle');
my $interestBiotype = $q->param('interestBiotype');
my $pointsBiotype = $q->param('pointsBiotype');
my $sampleID = $q->param('sampleID');
my $species = $q->param('species');
my $searchUp = $q->param('searchUp');
my $searchDown = $q->param('searchDown');


#Open the connections to the fusion DB & Ensembl DB
my ($dbh,$registry) = (fusionModule::openFusionConnection,fusionModule::openEnsemblConnection);
#Gets the breakpoint coordinates (fusion_ID, chr coordinates, orientation and type (Ins/Del etc.)
my $fusionpoints = fusionModule::getBreakpointDataWhereSampleID($dbh, $sampleID);
#Get sample information
my $sampleInfo = fusionModule::getSampleInfoWhereSampleID($dbh, $sampleID);
#Get a hash with the fusionpoint features
my $featureHash = fusionModule::AnalyzeFeatures($dbh, $registry, $fusionpoints , $fusionpoints, $searchUp, $searchDown, \@interestGenes, \@interestDomains, $interestBiotype, \@interestRegEles);
#Get a hash with the sorting scores
my %scoreHash = fusionModule::sortImportance($featureHash, $pointsInFusionGene, $pointsOutFusionGene, $pointsPosInFusionGene, $pointsUserGene, $pointsUserDomain, $pointsBiotype, $pointsUserRegEle, $pointsPromotorBreak, $pointsGeneBreak);

#Set content type and print HTML
print  "Content-type:text/html\r\n\r\n";

print '<HTML>
<HEAD>
<title>Overview Uploaded Fusionpoints</title>
<LINK href="../VisualBrowserLayout.css" rel="stylesheet" type="text/css">
</HEAD>
<BODY style="background-color:#EEEEEE;">
<div class="horizontalMenu">
	<table >
		<tr> 
			<td><a href="../uploadScreen.html">Upload new breakpoints</a></td>
			<td><a href="viewExisting.pl">View existing breakpoints</a></td>
		</tr>
	</table>
</div>

<div class="content">
	<p align="center"><i>Uploaded fusionpoints and number of features, ordered on possible biological relevance:</i	></p>
	Sample ID: <b>'.$sampleID.'</b> Samplename: <b>'.$sampleInfo->{$sampleID}->{'sample_name'}.'</b> <i>('.$species.')</i><br>
	<i>Total records found: '.scalar(keys %scoreHash).' </i><br>
	<!--Form to make a .CSV file of the selected fusionpoints-->
	<form action="makeCSV.pl" method="post" target="_blank">
		<input type="hidden" name="sampleID" value="'.$sampleID.'"/>
		<input type= "submit" name="submit" value="Save as .CSV"/>

	</form>

<div class="fusionTable" >
	<table >	
		<tr> 
			<td>Fusion #</td>
			<td>Chr1+Start&amp;End</td>
			<td>Chr2+Start&amp;End</td>
			<td>Orientation</td>
			<td>SV Type</td>
			<td>Fusion-genes? (In|Out) frame</td>
			<td>Fusion-Names</td>
			<td>Gene breaks?(Intronic|Exonic)</td>
			<td>Genebreak name(s)</td>
			<td>Promotor breaks?</td>
			<td>User Genes found:</td>
			<td>User Domains found:</td>
			<td>User Reg. Eles. found:</td>
			<td>User Biotypes found:</td>
			<td>Score:</td>
		</tr>';
#Fill the table based on the score
foreach my $fusID (sort {$scoreHash{$b} <=> $scoreHash{$a}} keys %scoreHash) {
#Fill the table
	print "<tr>
			<td>".$fusID."</td>
			<td>".$fusionpoints->{$fusID}->{'break_chr1'}.":".$fusionpoints->{$fusID}->{'break_chr1_start'}.":".$fusionpoints->{$fusID}->{'break_chr1_end'}."</td>
			<td>".$fusionpoints->{$fusID}->{'break_chr2'}.":".$fusionpoints->{$fusID}->{'break_chr2_start'}.":".$fusionpoints->{$fusID}->{'break_chr2_end'}."</td>
			<td>".$fusionpoints->{$fusID}->{'break_orientation'}."</td>";
			print "<td>Insertion</td>" if $fusionpoints->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 1;
			print "<td>Deletion</td>" if $fusionpoints->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 2;
			print "<td>Anti</td>" if $fusionpoints->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 3;
			print "<td>Inversion</td>" if $fusionpoints->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 4;
			print "<td>Remote</td>"if $fusionpoints->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 5;
			print "<td>Translocation</td>" if $fusionpoints->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 6;
			print "<td>Duplication</td>" if $fusionpoints->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 7;
			print "<td>Evertion</td>" if $fusionpoints->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 8;
			print "<td>Other</td>" if $fusionpoints->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 9;
			if($fusionpoints->{$fusID}->{'break_chr1_start'} == $fusionpoints->{$fusID}->{'break_chr1_end'}){
			print "<td>".$featureHash->{$fusID}->{'countFusInFrameGenes'}."|".$featureHash->{$fusID}->{"countFusOutFrameGenes"}."</td>"
			}else{
			print "<td>".$featureHash->{$fusID}->{'countFusPosInFrameGenes'}." (Possible)</td>"
			}
			print"
			<td>".$featureHash->{$fusID}->{"fusionGeneNames"}."</td>
			<td>".$featureHash->{$fusID}->{'countIntronicBreak'}."|".$featureHash->{$fusID}->{"countExonicBreak"}."</td>
			<td>".$featureHash->{$fusID}->{"geneBreakNames"}."</td>
			<td>".$featureHash->{$fusID}->{"countPromotorBreaks"}."</td>
			<td>".$featureHash->{$fusID}->{"geneString"}."</td>
			<td>".$featureHash->{$fusID}->{"domainString"}."</td>
			<td>".$featureHash->{$fusID}->{"regEleString"}."</td>
			<td>".$featureHash->{$fusID}->{"countBiotype"}."</td>
			<td>".$scoreHash{$fusID}."</td>
			<td><form action='preVisual.pl' method='post' target='_blank'><input type='hidden' name='fusID' value=".$fusID." /><input type='hidden' name='searchUp' value=".$searchUp." /><input type='hidden' name='searchDown' value=".$searchDown." /><input type= 'submit' name='submit' value='View'/></form></td>
			<td><form action='detailOverview.pl' method='post' target='_blank'><input type='hidden' name='fusID' value=".$fusID." /><input type= 'submit' name='submit' value='Info'/></form></td>
			";
}

#Last bit of HTML (Closing the table, div etc.)
print '
	</table>
</div>
</BODY>
</HTML>
';

#Close the connection to the fusion DB
fusionModule::closeFusionConnection($dbh);

### End of main script ###
1;
