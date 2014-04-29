#!/usr/bin/perl
#Author:			Job van Riet
#Date of creation:		4/1/13
#Last date of modification:	4/1/13
#Modifications:

#Known bugs:			None
#Function:			This script will make a .csv file from the selected breakpoints

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
my $fusID = $q->param('fusID');
my $sampleID = $q->param('sampleID');
my $sampleName = $q->param('sampleName');

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
$species = "none" if ! defined($q->param('species'));


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
my ($fusionpoints, $featureHash, %scoreHash);

#Get all the breakpoints from the fusionDB and search for user-specified features in the ensembl DB
my ($fusionpoints, $featureHash) = fusionModule::getFusionpointsBasedOnSearchParameters($dbh, $registry, $species, $searchUp, $searchDown, $fusID, $sampleID, $sampleName, $fusionGenesCheck, $confirmedCheck, $geneBreakCheck, $promotorCheck, \@interestGenes, \@interestDomains, \@interestRegEles, $interestBiotype);
#Sort the fusionpoints according to the user-specified scores
my %scoreHash = fusionModule::sortImportance($featureHash, $pointsInFusionGene, $pointsOutFusionGene, $pointsPosInFusionGene, $pointsUserGene, $pointsUserDomain, $pointsBiotype, $pointsUserRegEle, $pointsPromotorBreak, $pointsGeneBreak);

print "Content-type: application/octet-stream\n";
print "Content-Disposition: attachment;filename=fusionPoints.csv\n\n";

#Print header
print join("\t", "Fusion #", "Chr1+Start&;End", "Chr2+Start&;End","Orientation","SV Type","Fusion-genes? (In|Out) frame","Fusion-name(s)","Gene breaks?(Intronic|Exonic)","Genebreak name(s)","Promotor breaks?", "User Genes found:", "User Domains found:", "User Reg. Eles. found:","User Biotypes found:")."\n";
	
#Make CSV File
foreach my $fusID (reverse sort {$scoreHash{$a} <=> $scoreHash{$b}} keys %scoreHash) {
	#Only show the fusionpoints that have the features the user has requested
	
	if(((scalar(@interestGenes) != 0 and $featureHash->{$fusID}->{"geneString"} ne "No") or scalar(@interestGenes) == 0)){
		if(((scalar(@interestDomains) != 0 and $featureHash->{$fusID}->{"domainString"} ne "No") or scalar(@interestDomains) == 0)){
			if(((scalar(@interestRegEles) != 0 and $featureHash->{$fusID}->{"regEleString"} ne "No") or scalar(@interestRegEles) == 0)){
				if((($interestBiotype != "" and $featureHash->{$fusID}->{"countBiotype"} ne "No") or $interestBiotype == "")){
					my $sampleInfo = fusionModule::getSampleInfoWhereFusionID($dbh, $fusID);
					print $fusID ."\t";
					print $fusionpoints->{$fusID}->{'sample_name'}."\t";
					print $fusionpoints->{$fusID}->{'break_chr1'}.":".$fusionpoints->{$fusID}->{'break_chr1_start'}.":".$fusionpoints->{$fusID}->{'break_chr1_end'}."\t";
					print $fusionpoints->{$fusID}->{'break_chr2'}.":".$fusionpoints->{$fusID}->{'break_chr2_start'}.":".$fusionpoints->{$fusID}->{'break_chr2_end'}."\t";
					print $fusionpoints->{$fusID}->{'break_orientation'}."\t";
					print "Insertion\t" if $fusionpoints->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 1;
					print "Deletion\t" if $fusionpoints->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 2;
					print "Anti\t" if $fusionpoints->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 3;
					print "Inversion\t" if $fusionpoints->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 4;
					print "Remote\t"if $fusionpoints->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 5;
                       			print "Translocation\t" if $fusionpoints->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 6;
                        		print "Duplication\t" if $fusionpoints->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 7;
                        		print "Other\t" if $fusionpoints->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 8;

					if($fusionpoints->{$fusID}->{'break_chr1_start'} == $fusionpoints->{$fusID}->{'break_chr1_end'}){
						print $featureHash->{$fusID}->{'countFusInFrameGenes'}."|".$featureHash->{$fusID}->{"countFusOutFrameGenes"}."\t";
					}else{
						print $featureHash->{$fusID}->{'countFusPosInFrameGenes'}." (Possible)\t";
					}
					print $featureHash->{$fusID}->{"fusionGeneNames"}."\t";
					print $featureHash->{$fusID}->{'countIntronicBreak'}."|".$featureHash->{$fusID}->{"countExonicBreak"}."\t";
					print $featureHash->{$fusID}->{"geneBreakNames"}."\t";
					print $featureHash->{$fusID}->{"countPromotorBreaks"}."\t";
					print $featureHash->{$fusID}->{"geneString"}."\t";
					print $featureHash->{$fusID}->{"domainString"}."\t";
					print $featureHash->{$fusID}->{"regEleString"}."\t";
					print $featureHash->{$fusID}->{"countBiotype"}."\t";
					print $scoreHash{$fusID}."\n";
				}
			}
		}
	}
}



