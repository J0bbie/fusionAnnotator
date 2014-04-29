#!/usr/bin/perl
#Author: 			Job van Riet
#Date of creation:		19/9/12
#Data of modification:		19/9/12
#Known bugs:			None
#Function:			This page will show a more detailed overview of the features found in a certian fusionpoint (Declared by the user by clicking on the "Info" button)
#				The tables on this page are scrollable to allow vast amounts of features yet still hold some compact layout. The IDs/names of the features are clickable
#				and will direct users to the Ensembl page of the feature.

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

#Get the fusionID
my $q = new CGI;
my $fusID = $q->param('fusID');
my $searchUp = $q->param('searchUp');
my $sampleID = $q->param('sampleID');
my $searchDown = $q->param('searchDown');
my @supportFilesNew = split(/\,/, $q->param('supportListHolder'));


#Open connection to fusion and ensembl DB
my $dbh =fusionModule::openFusionConnection;
my $registry = fusionModule::openEnsemblConnection;

#Change the support files if the user has requested this (By using the change support form)
if(defined($q->param('supportListHolder'))){
	fusionModule::changeSupportFiles($dbh, $fusID, $sampleID, \@supportFilesNew);
}

#Get the features of the fusionpoint
my $fusionInfo = fusionModule::getFusionPointsWhereFusionID($dbh, $fusID);
#Get sample info and support files
my $sampleInfo = fusionModule::getSampleInfoWhereFusionID($dbh, $fusID);
#Get the fusionGenes/geneBreaks and promotor breaks
my %fusionGenes->{$fusID} = fusionModule::getFusionBreakInfoWhereFusionID($dbh, $fusID);
my $fusionGenes = \%fusionGenes;

#Get slices to gather all the genes surrounding the fusionpoints (Displayed is most-lowest table.
my $sa = $registry->get_adaptor($sampleInfo->{$fusID}->{'sample_species'}, 'Core', 'Slice');
my ($breakSlice1, $breakSlice2) = fusionModule::getEnsemblSlices($registry, $fusionInfo->{$fusID}->{"break_chr1"}, $fusionInfo->{$fusID}->{"break_chr1_start"} ,$fusionInfo->{$fusID}->{"break_chr1_end"}, $fusionInfo->{$fusID}->{"break_chr2"}, $fusionInfo->{$fusID}->{"break_chr2_start"},$fusionInfo->{$fusID}->{"break_chr2_end"},$fusionInfo->{$fusID}->{"break_orientation"}, $searchUp, $searchDown, $sampleInfo->{$fusID}->{'sample_species'});

#Set content type and print HTML
print  "Content-type:text/html\r\n\r\n";

print q(
<HTML>
<HEAD>
<title>Detailed Overview</title>
<LINK href="../VisualBrowserLayout.css" rel="stylesheet" type="text/css">

<!--Script to have a dynamic table -->
<SCRIPT language="javascript">
    //Array to save the link in to pass on
    var supportStack = new Array();
    var idCounter = 0;
    //String used to keep the id's
    var idString = "";

    //Function to add rows to a table
    function addSupportFile(tableID, type, link, chkSpecificity) {
        var table = document.getElementById(tableID);

        var rowCount = table.rows.length;
        var row = table.insertRow(rowCount);
        //Make new ID (Keep in mind the rows from the DB)
        idCounter++;

        var rowID = (document.getElementById('rowCountDB').value * 1);
        var ID = idCounter + rowID;

        var idS = ID.toString();
        var idStringHolder = idString.concat(idS+"||||");
        idString = idStringHolder;

        //Selected support file
        var cell1 = row.insertCell(0);
        var element1 = document.createElement("input");
        element1.type = "checkbox";
        element1.name = "chk";
        element1.value = ID;
        cell1.appendChild(element1);

        //Cell used to indicate whether or not the link is applied to all the fusionpoints in a sample or just to the selected fusionpoint
        var cell2 = row.insertCell(1);
        var element2 = document.createElement("input");
        element2.type = "checkbox";
        element2.id = ("spec-" + ID);
        if (chkSpecificity.match(/(Y|J)/i)) {
            element2.checked = true;
            element2.value = "fusion";
        } else {
            element2.value = "sample";
        }
        element2.disabled = "disabled";
        cell2.appendChild(element2);

        //Cell used to display the type of the link
        var cell3 = row.insertCell(2);
        var element3 = document.createElement("input");
        element3.type = "text";
        element3.id = ("type-" + ID);
        element3.value = type;
        element3.disabled = "disabled";
        cell3.appendChild(element3);

        //Cell to write link to
        var cell4 = row.insertCell(3);
        cell4.style.width = "100%";
        var element4 = document.createElement("input");
        element4.type = "text";
        element4.id = ("link-" + ID);
        element4.value = link;
        element4.disabled = "disabled";
        element4.style.width = "100%";
        cell4.appendChild(element4);
    }

    function addLink(tableID) {
        var link = prompt("Copy/Paste link of support file", "HTTP link here");
        var chkSpecificity = prompt("Fusion-point specific?(Yes/No)", "No");

        if (link.match(/(\.bed$|\.wig$|\.gff$|\.bam$)/i)) {
            if (link.match(/http:\/\//i)) {
                var type;
                if (link.match(/(\.bed$)/i)) {
                    var type = "BED";
                }
                if (link.match(/(\.bam$)/i)) {
                    var type = "BAM";
                }
                if (link.match(/(\.gff$)/i)) {
                    var type = "GFF";
                }
                if (link.match(/(\.wig$)/i)) {
                    var type = "WIG";
                }
                //Add the support file to the table
                addSupportFile(tableID, type, link, chkSpecificity);
            } else {
                alert("Not a supported file type! (Only .BED, .GFF , .WIG or .BAM in http links)");
            }
        } else {
            alert("Not a supported file type! (Only .BED, .GFF , .WIG or .BAM in http links)");
        }
    }

    //Delete rows from a table
    function deleteRow(tableID) {
        try {
            var table = document.getElementById(tableID);
            var rowCount = table.rows.length;

            for (var i = 0; i < rowCount; i++) {
                var row = table.rows[i];
                var chkbox = row.cells[0].childNodes[0];
                var rowIDString = document.getElementById('rowIDString').value;

                if (null != chkbox && true == chkbox.checked) {
                    //Delete ID from idString
                    var idStringHolder = idString.replace("||" + chkbox.value + "||", "");
                    idString = idStringHolder;

                    //Delete ID from original rowCount from DB if needed
                    idStringHolder = rowIDString.replace(chkbox.value + "||||", "");
                    document.getElementById('rowIDString').value = idStringHolder;
                    table.deleteRow(i);
                    rowCount--;
                    i--;
                }


            }
        } catch (e) {
            alert(e);
        }
    }

    //Get all the rows and push their values into a stack that is POST passable
    function getSupportFiles() {
        var table = document.getElementById('dataTable');
        //Make an array of all the IDs from the DB
        var rowIDString = document.getElementById('rowIDString').value;
        var rowIDHolder = rowIDString.split("||||");

        //Push all the left IDs from the DB to the stack
        if(rowIDHolder.length != 0){
		for (var i = 0; i < (rowIDHolder.length-1); i++) {
			var id = rowIDHolder[i].replace(/\|/g,"");
			//Get values
			var chkSpecificity = document.getElementById('spec-' + id).value;
			var type = document.getElementById('type-' + id).value;
			var link = document.getElementById('link-' + id).value;
			//Push support file to stack
			supportStack.push(chkSpecificity + "||" + type + "||" + link);
		}
        }
        //Push all the added support files to the stack
        var addedIDHolder = idString.split("||||");
        if(addedIDHolder.length != 1){
		for (var i = 0; i < (addedIDHolder.length-1); i++) {
			var id = addedIDHolder[i];
			//Get values from fields
			var chkSpecificity = document.getElementById('spec-' + id).value;
			var type = document.getElementById('type-' + id).value;
			var link = document.getElementById('link-' + id).value;
			//Push support file to stack
			supportStack.push(chkSpecificity + "||" + type + "||" + link);
		}
	}
        
        document.getElementById('supportListHolder').value = supportStack.join();
    }
</SCRIPT>
);

print '


</HEAD>	
<BODY style="background-color:#eeeeee;">

<div class="content">
	<p align="left">This screen will provide detailed information about the requested fusionpoint</p>
<!--Shows general information about the fusionpoint-->
	<div class="fusionInfo">
		<table>
		<i>General Fusionpoint info</i>
		<tr><td>Fusion ID:</td><td>'.$fusID.'</td></tr>
			<tr><td>Chr1+Pos:</td><td>'.$fusionInfo->{$fusID}->{'break_chr1'}.':'.$fusionInfo->{$fusID}->{'break_chr1_start'}.'|'.$fusionInfo->{$fusID}->{'break_chr1_end'}.'</td></tr>	
			<tr><td>Chr2+Pos:</td><td>'.$fusionInfo->{$fusID}->{'break_chr2'}.':'.$fusionInfo->{$fusID}->{'break_chr2_start'}.'|'.$fusionInfo->{$fusID}->{'break_chr2_end'}.'</td></tr>	
			<tr><td>Orientation:</td><td>'.$fusionInfo->{$fusID}->{'break_orientation'}.'</td></tr>';
			print "<td>Breaktype:</td><td>Insertion</td>" if $fusionInfo->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 1;
			print "<td>Breaktype:</td><td>Deletion</td>" if $fusionInfo->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 2;
			print "<td>Breaktype:</td><td>Anti</td>" if $fusionInfo->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 3;
			print "<td>Breaktype:</td><td>Inversion</td>" if $fusionInfo->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 4;
			print "<td>Breaktype:</td><td>Remote</td>"if $fusionInfo->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 5;
			print "<td>Breaktype:</td><td>Translocation</td>"if $fusionInfo->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 6;
			print "<td>Breaktype:</td><td>Duplication</td>"if $fusionInfo->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 7;
			print "<td>Breaktype:</td><td>Evertion</td>"if $fusionInfo->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 8;
			print "<td>Breaktype:</td><td>Other</td>"if $fusionInfo->{"$fusID"}->{"T_breaktype_breaktype_ID"} == 9;
	
			#If there is a fusionpoint-overlapping gene found.
			if(scalar(keys % {$fusionGenes}) != 0){
				my ($strand,$type);
				foreach my $fusID (keys %$fusionGenes){
					if($fusionInfo->{$fusID}->{'break_chr1_start'} == $fusionInfo->{$fusID}->{'break_chr1_end'}){
						print '<tr><td># of In-frame fusion-genes</td><td>'.$fusionGenes->{$fusID}->{'countFusInFrameGenes'}.'</td></tr>';
						print '<tr><td># of Out-frame fusion-genes</td><td>'.$fusionGenes->{$fusID}->{'countFusOutFrameGenes'}.'</td></tr>';
						print '<tr><td>FusionGenes:</td><td><tr><td></tr></td>';
						foreach (split(", ", $fusionGenes->{$fusID}->{'fusionGeneNames'})){
							print "<tr><td><i>Fusiongene: </i></td><td>$_</td></tr>";		
						}
					}else{
						print '<tr><td># of possible fusion-genes</td><td>'.$fusionGenes->{$fusID}->{'countFusPosInFrameGenes'}.'</td></tr>';
						foreach (split(", ", $fusionGenes->{$fusID}->{'fusionGeneNames'})){
							print "<tr><td><i>Fusiongene: </i></td><td>$_</td></tr>" ;		
						}
					}
					print '<tr><td># of intronic non-fusion gene breaks</td><td>'.$fusionGenes->{$fusID}->{'countIntronicBreak'}.'</td></tr>';
					print '<tr><td># of exonic non-fusion gene breaks</td><td>'.$fusionGenes->{$fusID}->{'countExonicBreak'}.'</td></tr>';
					foreach (split(", ", $fusionGenes->{$fusID}->{'geneBreakNames'})){
						print "<tr><td>Gene broken: </td><td>$_</td></tr>" ;		
					}
					print '<tr><td># of possible promotor breaks</td><td>'.$fusionGenes->{$fusID}->{'countPromotorBreaks'}.'</td></tr>';
				}
			}	
			#If no fusionpoint-overlapping gene found
			else{
				print '<tr><td>No fusion gene(s)</td><td></td></tr>';
			}
			print '
			<tr><td>Sample ID: </td><td>'.$sampleInfo->{$fusID}->{'sample_ID'}.'</td></tr>			
			<tr><td>Sample name: </td><td>'.$sampleInfo->{$fusID}->{'sample_name'}.'</td></tr>			
			<tr><td>Sample species: </td><td>'.$sampleInfo->{$fusID}->{'sample_species'}.'</td></tr>
			';

			print ' 
		</table>
	</div>';
			
			print q(	
	
	<!--Dynamic table of support Files-->
	<br>
	<i>Add/Remove support files, indicate whether the support-file is fusion-point specific or applies to all fusionpoints in sample.</i><br>
	<i>If breakpoint is from StrucVarDB, the support files can only affect the specific requested breakpoint</i><br>
	<b>Press "Submit Changes" to let changes take affect!!</b><br>
	<INPUT type="button" value="Add support file" onclick="addLink('dataTable')" />
	<INPUT type="button" value="Delete selected support file(s)" onclick="deleteRow('dataTable')" />
	<!--Form to update the support files-->
	<form action="detailOverview.pl" method='post' onsubmit='getSupportFiles();'>
	<INPUT type="submit" value="Submit changes"/>
	<input type='hidden' name='fusID' value=).$fusID.q( />
	<input type='hidden' name='sampleID' value=).$sampleInfo->{$fusID}->{'sample_ID'}.q( />
	<input type='hidden' name='searchDown' value=).$searchDown.q( />
	<input type='hidden' name='searchUp' value=).$searchUp.q( />		
	<input type="hidden"  name="supportListHolder" id="supportListHolder"/>


 	<div class="featureTable">

	<table id="dataTable" width="100%" border="1">

	<tr>
		<th width="25%">Selected</th>
		<th width="25%">Fusionpoint-specific?</th>
		<th width="25%">Type</th>
		<th width="25%">Link</th>
	</tr>
	);
	my $rowCountDB = 0;
	foreach my $supportFile (@{$sampleInfo->{$fusID}->{'supportFilesSample'}}){
		$rowCountDB++;
		print '<tr><TD><INPUT type="checkbox" name="chk" value='.$rowCountDB.' /></TD><TD><INPUT type="checkbox" name="chkSpecific" id="spec-'.$rowCountDB.'" value="sample" disabled/></TD><TD><INPUT type="txt" id="type-'.$rowCountDB.'" value="BAM" disabled></TD><TD width="100%"><INPUT size="100%" type="text" id="link-'.$rowCountDB.'" disabled value='.$supportFile->[0].' /></TD></tr>' if $supportFile->[0] =~ /(.bam$)/i;
		print '<tr><TD><INPUT type="checkbox" name="chk" value='.$rowCountDB.' /></TD><TD><INPUT type="checkbox" name="chkSpecific" id="spec-'.$rowCountDB.'" value="sample" disabled/></TD><TD><INPUT type="txt" id="type-'.$rowCountDB.'" value="GFF" disabled></TD><TD width="100%"><INPUT size="100%" type="text" id="link-'.$rowCountDB.'" disabled value='.$supportFile->[0].' /></TD></tr>' if $supportFile->[0] =~ /(.gff$)/i;
		print '<tr><TD><INPUT type="checkbox" name="chk" value='.$rowCountDB.' /></TD><TD><INPUT type="checkbox" name="chkSpecific" id="spec-'.$rowCountDB.'" value="sample" disabled/></TD><TD><INPUT type="txt" id="type-'.$rowCountDB.'" value="BED" disabled></TD><TD width="100%"><INPUT size="100%" type="text" id="link-'.$rowCountDB.'" disabled value='.$supportFile->[0].' /></TD></tr>' if $supportFile->[0] =~ /(.bed$)/i;
		print '<tr><TD><INPUT type="checkbox" name="chk" value='.$rowCountDB.' /></TD><TD><INPUT type="checkbox" name="chkSpecific" id="spec-'.$rowCountDB.'" value="sample" disabled/></TD><TD><INPUT type="txt" id="type-'.$rowCountDB.'" value="WIG" disabled></TD><TD width="100%"><INPUT size="100%" type="text" id="link-'.$rowCountDB.'" disabled value='.$supportFile->[0].' /></TD></tr>' if $supportFile->[0] =~ /(.wig$)/i;
	}
	foreach my $supportFile (@{$sampleInfo->{$fusID}->{'supportFilesFusion'}}){
		$rowCountDB++;
		print '<tr><TD><INPUT type="checkbox" name="chk" value='.$rowCountDB.' /></TD><TD><INPUT type="checkbox" checked name="chkSpecific" id="spec-'.$rowCountDB.'" value="fusion" disabled/></TD><TD><INPUT type="txt" id="type-'.$rowCountDB.'" value="BAM" disabled></TD><TD width="100%"><INPUT size="100%" type="text" id="link-'.$rowCountDB.'" disabled value='.$supportFile->[0].' /></TD></tr>' if $supportFile->[0] =~ /(.bam$)/i;
		print '<tr><TD><INPUT type="checkbox" name="chk" value='.$rowCountDB.' /></TD><TD><INPUT type="checkbox" checked name="chkSpecific" id="spec-'.$rowCountDB.'" value="fusion" disabled/></TD><TD><INPUT type="txt" id="type-'.$rowCountDB.'" value="GFF" disabled></TD><TD width="100%"><INPUT size="100%" type="text" id="link-'.$rowCountDB.'" disabled value='.$supportFile->[0].' /></TD></tr>' if $supportFile->[0] =~ /(.gff$)/i;
		print '<tr><TD><INPUT type="checkbox" name="chk" value='.$rowCountDB.' /></TD><TD><INPUT type="checkbox" checked name="chkSpecific" id="spec-'.$rowCountDB.'" value="fusion" disabled/></TD><TD><INPUT type="txt" id="type-'.$rowCountDB.'" value="BED" disabled></TD><TD width="100%"><INPUT size="100%" type="text" id="link-'.$rowCountDB.'" disabled value='.$supportFile->[0].' /></TD></tr>' if $supportFile->[0] =~ /(.bed$)/i;
		print '<tr><TD><INPUT type="checkbox" name="chk" value='.$rowCountDB.' /></TD><TD><INPUT type="checkbox" checked name="chkSpecific" id="spec-'.$rowCountDB.'" value="fusion" disabled/></TD><TD><INPUT type="txt" id="type-'.$rowCountDB.'" value="WIG" disabled></TD><TD width="100%"><INPUT size="100%" type="text" id="link-'.$rowCountDB.'" disabled value='.$supportFile->[0].' /></TD></tr>' if $supportFile->[0] =~ /(.wig$)/i;
	}
	#Make an string of all the IDs
	my (@arrRowCount,$rowIDString);
	for (my $count = 1; $count <= $rowCountDB; $count++) {
		push(@arrRowCount, $count);
 	}
 	$rowIDString=join("||||",@arrRowCount);
 	if($rowIDString ne ""){
 		$rowIDString=$rowIDString."||||";
 	}
	print'
	<input type="hidden"  name="rowCountDB" id="rowCountDB" value='.$rowCountDB.' />
	<input type="hidden"  name="rowIDString" id="rowIDString" value='.$rowIDString.' />


	</table>
	</div>

	<br>
	<i>Name of genes(Exons), strand and type:</i><br>
		<div class="featureTable">
			<table>
				<tr>
					<th>Ensembl Gene Name</th>
					<th>Description</th>
					<th>Strand</th>
					<th>Type</th>
				</tr>';
					#Print the genes
					foreach my $slice ($breakSlice1, $breakSlice2){
						my $genesOnSlice = $slice->get_all_Genes();
						#Check if there are genes, if yes -> concenate on string
						foreach my $gene (@{$genesOnSlice}){
							print '<tr><td STYLE="text-align: left;"><a href = http://www.ensembl.org/'.$sampleInfo->{$fusID}->{'sample_species'}.'/Gene/Summary?g='.$gene->stable_id.'>'.$gene->external_name().' </a></td>';
							print '<td>No Description</td>' if $gene->description() eq '';
							print '<td>'.$gene->description().'</td>' if $gene->description() ne '';
							print '<td>+</td>' if $gene->strand() == 1;
							print '<td>-</td>' if $gene->strand() == -1;
							print '<td>'.$gene->biotype().'</td></tr>';
						}
					}
print q(		</table>
	</div>
	</form>

</BODY>
</HTML>);

#Close connection to fusion DB
fusionModule::closeFusionConnection($dbh);
