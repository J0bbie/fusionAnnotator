#!/usr/bin/perl -T
#Author:					Job van Riet
#Date of creation:			17/10/12
#Last date of modification:	17/10/12
#Modifications:
# - Changed the fusion DB to localhost DB
#Known bugs:			Cant use compacted code to get the HTML parameters because it messes up to correct declaration to the specific variables, solved by using individual declarations
#Function:			This script will query the (local) Ensembl DB for genes, exon, domains, and regulatory element
#					It will also check if the exons of the genes overlapping the fusionpoint are in-frame and store this accordingly

use lib '../modules/';

#Use strict programming
use strict;
#Get the fusion application functions
use fusionModule;
#Use CGI
use CGI;
use CGI::Carp qw(fatalsToBrowser);    # Remove for production use

#Main script
$CGI::POST_MAX = 1024 * 100;  # maximum upload filesize is 100mb

my $q = new CGI;
#Start HTML
print $q->header;
print $q->start_html(
	-title => "Upload Fusionpoints",
	-style=>"../VisualBrowserLayout.css",
	-bgcolor => "#FFFFFF"
	);

print '<div class="options" align="center" style="color: #ff7f00">';
#Get the passed form parameters from the upload screen
#Do it individually instead of compacted or else it messes up somehow...
my $breakFile = $q->upload('file');
my @supportFiles = split(/\,/, $q->param('supportListHolder'));
my $searchUp = $q->param('searchUp')*10e3;
my $searchDown =  $q->param('searchDown')*10e3;
my $species =  $q->param('species');
my $interestGenes = $q->param('interestGenes');
my $interestDomains = $q->param('interestDomains');
my $interestRegEles = $q->param('interestRegEles');
my $pointsInFusionGene = $q->param('pointsInFusionGene');
my $pointsOutFusionGene = $q->param('pointsOutFusionGene');
my $pointsPosInFusionGene = $q->param('pointsPosInFusionGene');
my $pointsPromotorBreak = $q->param('pointsPromotorBreak');
my $pointsGeneBreak = $q->param('pointsGeneBreak');
my $pointsUserGene = $q->param('pointsUserGene');
my $pointsUserDomain = $q->param('pointsUserDomain');
my $pointsUserRegEle = $q->param('pointsUserRegEle');
my $interestBiotype = $q->param('interestBiotype');
my $pointBiotype = $q->param('pointsBiotype');

#Call the methods
#Save the uploaded file temporarily for reading purposes
my $tempFile = fusionModule::saveTempFile($q, $breakFile);
#Read the breakpoint file and store the breakpoints in the DB
my $sample_ID = fusionModule::readFusionFile($q, $tempFile, $species, \@supportFiles);
#Get the genes, regulatory elements and fusion-overlapping (in-frame/out-frame) exons and store them in the DB.
my $dbh = fusionModule::openFusionConnection;
my $registry = fusionModule::openEnsemblConnection;

#Gets the breakpoint coordinates (fusion_ID, chr1, chr1_end, chr2, chr2_start)
my $refFusionData = fusionModule::getBreakpointDataWhereSampleID($dbh, $sample_ID);
#For each fusionpoint
foreach my $id (keys %$refFusionData) {
	fusionModule::getFeaturesFromEnsembl($dbh, $registry, $refFusionData->{$id}, $searchDown, $searchUp, $species);
}

#Save the data if everything went right
fusionModule::closeFusionConnectionWithCommit($dbh);

print '<form action="viewNew.pl" method="post" enctype="multipart/form-data">
		<input type="hidden" name="interestGenes" value='.$interestGenes.' />
		<input type="hidden" name="searchUp" value='.$searchUp.' />
		<input type="hidden" name="searchDown" value='.$searchDown.' />
		<input type="hidden" name="interestDomains" value='.$interestDomains.' />
		<input type="hidden" name="interestRegEles" value='.$interestRegEles.' />
		<input type="hidden" name="interestBiotype" value='.$interestBiotype.' />
		<input type="hidden" name="pointsInFusionGene" value='.$pointsInFusionGene.' />
		<input type="hidden" name="pointsOutFusionGene" value='.$pointsOutFusionGene.' />
		<input type="hidden" name="pointsPosInFusionGene" value='.$pointsPosInFusionGene.' />
		<input type="hidden" name="pointsGeneBreak" value='.$pointsGeneBreak.' />
		<input type="hidden" name="pointsUserGene" value='.$pointsUserGene.' />
		<input type="hidden" name="pointsPromotorBreak" value='.$pointsPromotorBreak.' />
		<input type="hidden" name="pointsUserDomain" value='.$pointsUserDomain.' />
		<input type="hidden" name="pointsUserRegEle" value='.$pointsUserRegEle.' />
		<input type="hidden" name="pointsBiotype" value='.$pointBiotype.' />
		<input type="hidden" name="sampleID" value='.$sample_ID.' />
		<input type="hidden" name="species" value='.$species.' />
		<br><input type="submit" name="submit" value="Go to overview"/><br>
		</div></form>';

#End of HTML
print $q->end_html;

#End of main script
1;
