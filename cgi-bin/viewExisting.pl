#!/usr/bin/perl
#Author:			Job van Riet
#Date of creation:		29/10/12
#Last date of modification:	23/1/13
#Modifications:

#Known bugs:			None
#Function:			This script will get the newly-uploaded fusionpoints from the fusion DB based on a sample ID.
#				It will display the features found in a scrollable table, the fusionpoints are ordered on user defined importance
#				The domains, genes and regulatory elements the user has defined will be looked up in these fusionpoints and are also listed.
#				This scripts requires the sampleID, and the list of user defined features of importance.
#				After every fusionpoint, there is a possibility to 'view' the fusionpoint in the genome browser or to open a new screen/tab where additional 'info' can be found about the selected fusionpoint.

#Get the fusionModule
use lib '../modules/';

#Use strict programming
use strict;
#Use the fusion module
use fusionModule;
#Use CGI
use CGI;
use CGI::Carp qw (fatalsToBrowser);    # Remove for production use

### Main script ###

#Get the search parameters
my $q = new CGI;
my ($fusID, $sampleID, $sampleName)  = ($q->param('fusID'), $q->param('sampleID'), $q->param('sampleName'));

my $searchUp = $q->param('searchUp')*10e3;
my $searchDown =  $q->param('searchDown')*10e3;
my $species =  $q->param('species');

#Get user interest
my $fusionGenesCheck = $q->param('fusionGeneCheck');
my $confirmedCheck = $q->param('confirmedCheck');
my $promotorCheck = $q->param('promotorCheck');
my $geneBreakCheck = $q->param('geneBreakCheck');

my $interestGenes = $q->param('interestGenes');
my $interestDomains = $q->param('interestDomains');
my $interestRegEles = $q->param('interestRegEles');
my $interestBiotype = $q->param('interestBiotype');

#Get points for user-interests
my $pointsInFusionGene = $q->param('pointsInFusionGene');
my $pointsOutFusionGene = $q->param('pointsOutFusionGene');
my $pointsPosInFusionGene = $q->param('pointsPosInFusionGene');
my $pointsPromotorBreak = $q->param('pointsPromotorBreak');
my $pointsGeneBreak = $q->param('pointsGeneBreak');
my $pointsUserGene = $q->param('pointsUserGene');
my $pointsUserDomain = $q->param('pointsUserDomain');
my $pointsUserRegEle = $q->param('pointsUserRegEle');
my $pointsBiotype = $q->param('pointsBiotype');

#Set fusID 0 if the page is requested for the first time, to prevent all the fusionpoints from loading
$fusID = 0 if ! defined($q->param('fusID'));

$species = "none" if ! defined($q->param('species'));
$pointsInFusionGene = "1500" if ! defined($q->param('pointsInFusionGene'));
$pointsOutFusionGene = "1250" if ! defined($q->param('pointsOutFusionGene'));
$pointsPosInFusionGene = "1000" if ! defined($q->param('pointsPosInFusionGene'));
$pointsPromotorBreak = "100" if ! defined($q->param('pointsPromotorBreak'));
$pointsGeneBreak = "10" if ! defined($q->param('pointsGeneBreak'));
$pointsUserGene = "5" if ! defined($q->param('pointsUserGene'));
$pointsUserDomain = "3" if ! defined($q->param('pointsUserDomain'));
$pointsUserRegEle = "5" if ! defined($q->param('pointsUserRegEle'));
$pointsBiotype = "5" if ! defined($q->param('pointsBiotype'));

#Make an array of the user-specified features
my (@interestGenes, @interestDomains, @interestRegEles);

if($q->param('interestGenes') ne ""){
	#If only 1 gene
	if(split(/;/,$q->param('interestGenes')) == 0){
		push (@interestGenes, $q->param('interestGenes'));
	}
	else{
		@interestGenes = split(/;/,$q->param('interestGenes'));
	}
}
if($q->param('interestDomains') ne ""){
	#If only 1 domains
	if(split(/;/,$q->param('interestDomains')) == 0){
		push (@interestDomains, $q->param('interestDomains'));
	}
	else{
		@interestDomains = split(/;/,$q->param('interestDomains'));
	}
}

if($q->param('interestRegEles') ne ""){
	#If only 1 reg Ele
	if(split(/;/,$q->param('interestRegEles')) == 0){
		push (@interestRegEles, $q->param('interestRegEles'));
	}
	else{
		@interestRegEles = split(/;/,$q->param('interestRegEles'));
	}
}

#Open connection to fusion and ensembl DB
my $dbh =fusionModule::openFusionConnection;
my $registry = fusionModule::openEnsemblConnection;

#Get the fusionpoints
print  "Content-type:text/html\r\n\r\n";

print '<HTML>
<HEAD>
<title>View Existing Fusionpoints</title>
<LINK href="../VisualBrowserLayout.css" rel="stylesheet" type="text/css">
</HEAD>';
my ($fusionpoints, $featureHash, %scoreHash);

#Get all the breakpoints from the fusionDB and search for user-specified features in the ensembl DB
my ($fusionpoints, $featureHash) = fusionModule::getFusionpointsBasedOnSearchParameters($dbh, $registry, $species, $searchUp, $searchDown, $fusID, $sampleID, $sampleName, $fusionGenesCheck, $confirmedCheck, $geneBreakCheck, $promotorCheck, \@interestGenes, \@interestDomains, \@interestRegEles, $interestBiotype);
#Sort the fusionpoints according to the user-specified scores
my %scoreHash = fusionModule::sortImportance($featureHash, $pointsInFusionGene, $pointsOutFusionGene, $pointsPosInFusionGene, $pointsUserGene, $pointsUserDomain, $pointsBiotype, $pointsUserRegEle, $pointsPromotorBreak, $pointsGeneBreak);

#Get the names of all the uploaded samples to show in the dropdown menu
my $sampleNames = fusionModule::getSampleNames($dbh);


print '
<BODY style="background-color:#eeeeee;">
<!--Menu-->
<div class="horizontalMenu">
	<table >
		<tr> 
			<td><a href="../uploadScreen.html">Upload new breakpoints</a></td>
			<td>View existing breakpoints</td>
		</tr>
	</table>
</div>
<!--Filter options-->
<div class="content">
	<form action="#View_Existing.pl" method="post" enctype="multipart/form-data">
		<div class="options" align="center">
			<i>Query the DB with the following commands:</i><br>
			<b>Get by sample or fusion ID</b><br>
			<i>Use ; for "between" like 1;4, leave all empty for all fusionpoints.</i><br>';
			print '<label for="fusID">Fusion ID</label><input type="text" name="fusID" value="'.$fusID.'"/><br>' if $fusID != 0;
			print q(<label for="fusID">Fusion ID</label><input type="text" name="fusID" value=""/><br>) if $fusID == 0;


			print '
			<label for="sampleID">Sample ID</label><input type="text" name="sampleID" value="'.$sampleID.'"/><br>
			<label for="sampleName">Sample Name</label><select name="sampleName">
			<option value="">All Samples</option>';
			foreach(keys %{$sampleNames}){
				if($_ eq $sampleName){
					print '<option value="'.$_.' " selected >'.$_.'</option>';
				}
				else{
					print '<option value="'.$_.'">'.$_.'</option>';
				}
			}			
			print '
			</select><br>
			
			<b>Get by features.</b><br>
			<i>Use ; as a delimiter in the defining of the domains and/or genes that should be searched for.</i><br>';
			print q(<label for="confirmedCheck">Only show confirmed breakpoints</label><input type="checkbox" name="confirmedCheck" value="1"/><br>) if $confirmedCheck != 1;
			print q(<label for="confirmedCheck">Only show confirmed breakpoints</label><input type="checkbox" name="confirmedCheck" value="1" checked="yes"/><br>) if $confirmedCheck == 1;

			print q(<label for="fusionGeneCheck">Only show fusion genes</label><input type="checkbox" name="fusionGeneCheck" value="1"/><br>) if $fusionGenesCheck != 1;
			print q(<label for="fusionGeneCheck">Only show fusion genes</label><input type="checkbox" name="fusionGeneCheck" value="1" checked="yes" /><br>) if $fusionGenesCheck == 1;
			
			print q(<label for="promotorCheck">Only show promotor breaks</label><input type="checkbox" name="promotorCheck" value="1"/><br>) if $promotorCheck != 1;
			print q(<label for="promotorCheck">Only show promotor breaks</label><input type="checkbox" name="promotorCheck" value="1" checked="yes" /><br>) if $promotorCheck == 1;
			
			
			print q(<label for="geneBreakCheck">Only show gene-breaks</label><input type="checkbox" name="geneBreakCheck" value="1"/><br>) if $geneBreakCheck != 1;
			print q(<label for="geneBreakCheck">Only show gene-breaks</label><input type="checkbox" name="geneBreakCheck" value="1"/><br>) if $geneBreakCheck == 1;

			print '<label for="species">Show only organism:</label>
			<select name="species">';
			print q(<option value="none">No species filtering</option>) if $species ne "none";
			print q(<option value="none" selected>No species filtering</option>) if $species eq "none";
			print q(<option value="homo_sapiens">Homo sapiens</option>) if $species ne "homo_sapiens";
			print q(<option value="homo_sapiens" selected>Homo sapiens</option>) if $species eq "homo_sapiens";
			print q(<option value="danio_rerio">Danio rerio</option>) if $species ne "danio_rerio";
			print q(<option value="danio_rerio" selected>Danio rerio</option>) if $species eq "danio_rerio";
			print '</select>
			<br>
			<label for="interestGenes">Search for genes</label><input type="text" name="interestGenes" value="'.$interestGenes.'"/><br>
			<label for="interestDomains">Search for domains</label><input type="text" name="interestDomains" value="'.$interestDomains.'"/><br>
			<label for="interestRegEles">Search for reg Eles</label><input type="text" name="interestRegEles" value="'.$interestRegEles.'"/><br>
			<label for="interestBiotype">Only show features with biotype</label>
			<select name="interestBiotype">';
			print q(<option value="none">No biotype filtering</option>) if $interestBiotype ne "none";
			print q(<option value="none" selected>No biotype filtering</option>) if $interestBiotype eq "none";
			print q(<option value="protein_coding">Protein Coding</option>) if $interestBiotype ne "protein_coding";
			print q(<option value="protein_coding" selected>Protein Coding</option>) if $interestBiotype eq "protein_coding";
			print q(<option value="non_coding">Non-Coding</option>) if $interestBiotype ne "non_coding";
			print q(<option value="non_coding" selected>Non-Coding</option>) if $interestBiotype eq "non_coding";			
			print q(<option value="pseudogene">Pseudogene</option>) if $interestBiotype ne "pseudogene";
			print q(<option value="pseudogene" selected>Pseudogene</option>) if $interestBiotype eq "pseudogene";

			print'</select><br><br>
			<i>Search distances from fusionpoint in kb for the features (Default = 20 kb)</i><br>
			<label for="searchUp"><b>Upstream</b></label><input type="text" name="searchUp" value="'.($searchUp/10e3).'"/><br>
			<label for="searchDown"><b>Downstream</b></label><input type="text" name="searchDown" value="'.($searchDown/10e3).'"/><br>
			<br>
			<!-- Extra options for fusionpoint sorting -->
			<i>Adjust the points given to genomic features and found user-definitions, this will adjust how the fusionpoints are sorted in the display.</i><br>
			<label for="pointsInFusionGene"><b>"Real" In-frame Fusion gene:</b></label> <input type="text" name="pointsInFusionGene" value="'.$pointsInFusionGene.'"/><br>
			<label for="pointsOutFusionGene"><b>"Real" Out-frame Fusion gene:</b></label> <input type="text" name="pointsOutFusionGene" value="'.$pointsOutFusionGene.'"/><br>
			<label for="pointsPosInFusionGene"><b>Possible fusion gene:</b></label> <input type="text" name="pointsPosInFusionGene" value="'.$pointsPosInFusionGene.'"/><br>
			<label for="pointsPromotorBreak"><b>Promotor break:</b></label> <input type="text" name="pointsPromotorBreak" value="'.$pointsPromotorBreak.'"/><br>
			<label for="pointsGeneBreak"><b>Gene-break (without fusiongene):</b></label> <input type="text" name="pointsGeneBreak" value="'.$pointsGeneBreak.'"/><br>
			<label for="pointsUserGene"><b>User-Genes:</b></label> <input type="text" name="pointsUserGene" value="'.$pointsUserGene.'"/><br>
			<label for="pointsUserDomain"><b>User-Domains:</b></label> <input type="text" name="pointsUserDomain" value="'.$pointsUserDomain.'"/><br>
			<label for="pointsUserRegEle"><b>User-Reg. Ele.:</b></label> <input type="text" name="pointsUserRegEle" value="'.$pointsUserRegEle.'"/><br>
			<label for="pointsBiotype"><b>Points biotype:</b></label> <input type="text" name="pointsBiotype" value="'.$pointsBiotype.'"/><br>

			<br>
			<input type="submit" name="submit" value="Get fusionpoints"/><br>
			<br>
		</div>
	</form>
<div class="fusionTable" >
	<i>Total records found: '.scalar(keys %scoreHash).' </i><br>
	<!--Form to make a .CSV file of the selected fusionpoints-->
	<form action="makeCSV.pl" method="post" target="_blank">
		<input type="hidden" name="fusID" value="'.$fusID.'"/>
		<input type="hidden" name="sampleID" value="'.$sampleID.'"/>
		<input type="hidden" name="sampleName" value="'.$sampleName.'"/>
		<input type="hidden" name="fusionGeneCheck" value="'.$fusionGenesCheck.'"/>
		<input type="hidden" name="promotorCheck" value="'.$promotorCheck.'"/>
		<input type="hidden" name="geneBreakCheck" value="'.$geneBreakCheck.'"/>
		<input type="hidden" name="confirmedCheck" value="'.$confirmedCheck.'"/><br>
		<input type="hidden" name="species" value="'.$species.'"/>
		<input type="hidden" name="interestGenes" value="'.$interestGenes.'"/>
		<input type="hidden" name="interestDomains" value="'.$interestDomains.'"/>
		<input type="hidden" name="interestRegEles" value="'.$interestRegEles.'"/>
		<input type="hidden" name="interestBiotype" value="'.$interestBiotype.'"/>
		<input type="hidden" name="searchUp" value="'.$searchUp.'"/>
		<input type="hidden" name="searchDown" value="'.$searchDown.'"/>
		<input type="hidden" name="pointsInFusionGene" value="'.$pointsInFusionGene.'"/>
		<input type="hidden" name="pointsOutFusionGene" value="'.$pointsOutFusionGene.'"/>
		<input type="hidden" name="pointsPosInFusionGene" value="'.$pointsPosInFusionGene.'"/>
		<input type="hidden" name="pointsPromotorBreak" value="'.$pointsPromotorBreak.'"/>
		<input type="hidden" name="pointsGeneBreak" value="'.$pointsGeneBreak.'"/>
		<input type="hidden" name="pointsUserGene" value="'.$pointsUserGene.'"/>
		<input type="hidden" name="pointsUserDomain" value="'.$pointsUserDomain.'"/>
		<input type="hidden" name="pointsUserRegEle" value="'.$pointsUserRegEle.'"/>
		<input type="hidden" name="pointsBiotype" value="'.$pointsBiotype.'"/>
		<input type= "submit" name="submit" value="Save as .CSV"/>

	</form>
	<table >	
		<tr> 
			<td>Fusion #</td>
			<td>Sample Name</td>
			<td>Chr1+Start&amp;End</td>
			<td>Chr2+Start&amp;End</td>
			<td>Orientation</td>
			<td>SV Type</td>
			<td>Fusion-genes? (In|Out) frame</td>
			<td>Fusion-name(s)</td>
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
	#Only show the fusionpoints that have the features the user has requested
	if(((scalar(@interestGenes) != 0 and $featureHash->{$fusID}->{"geneString"} ne "No") or scalar(@interestGenes) == 0)){
		if(((scalar(@interestDomains) != 0 and $featureHash->{$fusID}->{"domainString"} ne "No") or scalar(@interestDomains) == 0)){
			if(((scalar(@interestRegEles) != 0 and $featureHash->{$fusID}->{"regEleString"} ne "No") or scalar(@interestRegEles) == 0)){
				if((($interestBiotype != "" and $featureHash->{$fusID}->{"countBiotype"} ne "No") or $interestBiotype == "")){
					my $sampleInfo = fusionModule::getSampleInfoWhereFusionID($dbh, $fusID);
					#Fill the table
					#Fill the table
					print "<tr>
					<td>".$fusID."</td>
					<td>".$sampleInfo->{$fusID}->{'sample_name'}."</td>
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
			}
		}
	}
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

#End of main Script
