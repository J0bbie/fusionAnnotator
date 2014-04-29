#!/usr/bin/perl
#Author:			Job van Riet
#Date of creation:		7/1/13
#Last date of modification:	7/1/13
#
#Function:			This script will accept one breakpoint from the StrucVarDB and store it in the fusionDB under the correct sample (depending on species), it then creates a table in which the features
#				of the fusionpoint can be seen. (Derivative of viewNew.pl and FusionPipe01.pl)

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

#Get the breakpoint data through POST parameters
my $q = new CGI;

my $chr1 = $q->param('chr1');
my $chr2 = $q->param('chr2');
my $s1 = $q->param('s1');
my $e1 = $q->param('e1');
my $s2 = $q->param('s2');
my $e2 = $q->param('e2');
my $species = $q->param('species');
my $orientation = $q->param('orientation');
my $breakType = $q->param('breakType');
my $sampleName = $q->param('sampleName');



#($species, $chr1, $s1, $e1, $chr2, $s2, $e2, $orientation, $breakType) = ("homo_sapiens", 1, 12212,12212, 4, 12345,12345, "TH", "Evertion");
#Other parameters needed for functions
my $searchUp = 20e3;
my $searchDown =  20e3;
my @interestGenes = undef;
my @interestDomains = undef;
my @interestRegEles = undef;
my $pointsInFusionGene = 1500;
my $pointsOutFusionGene = 1250;
my $pointsPosInFusionGene = 1000;
my $pointsPromotorBreak = 100;
my $pointsGeneBreak = 10;
my $pointsUserGene = 5;
my $pointsUserDomain = 3;
my $pointsUserRegEle = 5;
my $interestBiotype = "protein_coding";
my $pointsBiotype = 3;


#Set content type and print HTML
print  "Content-type:text/html\r\n\r\n";


#Store the breakpoint in the fusionAnnotator DB in sample ID 1 (StrucVarBreakpoints)
#Open the connections to the fusion DB & Ensembl DB
my ($dbh,$registry) = (fusionModule::openFusionConnection,fusionModule::openEnsemblConnection);
my $fusID = fusionModule::uploadBreakpointDataFromStrucVarDB($dbh, $sampleName, $species, $chr1, $s1, $e1, $chr2, $s2, $e2, $orientation, $breakType);

#Gets the breakpoint coordinates (fusion_ID, chr1, chr1_end, chr2, chr2_start)
my $refFusionData = fusionModule::getBreakpointDataWhereFusionID($dbh, $fusID);
#For each fusionpoint
foreach my $id (keys %$refFusionData) {
	fusionModule::getFeaturesFromEnsembl($dbh, $registry, $refFusionData->{$id}, $searchDown, $searchUp, $species);
}

##Create fusionTable
#Gets the breakpoint coordinates (fusion_ID, chr coordinates, orientation and type (Ins/Del etc.)
my $fusionpoints = fusionModule::getBreakpointDataWhereFusionID($dbh, $fusID);
#Get a hash with the fusionpoint features
my $featureHash = fusionModule::AnalyzeFeatures($dbh, $registry, $fusionpoints , $fusionpoints, $searchUp, $searchDown, \@interestGenes, \@interestDomains, $interestBiotype, \@interestRegEles);
#Get a hash with the sorting scores
my %scoreHash = fusionModule::sortImportance($featureHash, $pointsInFusionGene, $pointsOutFusionGene, $pointsPosInFusionGene, $pointsUserGene, $pointsUserDomain, $pointsBiotype, $pointsUserRegEle, $pointsPromotorBreak, $pointsGeneBreak);


print '<HTML>
<HEAD>
<title>Overview StrucVarDB Breakpoint</title>
<LINK href="../VisualBrowserLayout.css" rel="stylesheet" type="text/css">
</HEAD>
<BODY style="background-color:#EEEEEE;">

<div class="content">
	<p align="center"><i>Selected breakpoint and its features:</i></p>
	Sample name: <b>StrucVarDB-Breakpoint</b> <i>('.$species.')</i><br>

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
			<td>User Biotypes found:</td>
			<td>Score</td>
		</tr>';
#Fill the table based on the score
foreach my $fusID (reverse sort {$scoreHash{$a} <=> $scoreHash{$b}} keys %scoreHash) {
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
