#!/usr/bin/perl
#Author: 			Job van Riet
#Date of creation:		5/11/12
#Date of modification:		5/11/12
#Known bugs:			None
#Function:			Gets all the genes, transcripts, domains from the EnsEMBL DB
#				Puts them in a flat file for use of multiple tracks in Jbrowser, also load these tracks in



use strict;
#Use Bioperl and Ensembl API 67 for function
use lib '/usr/local/ensembl67/ensembl/modules';
use lib '/usr/local/ensembl67/ensembl-compara/modules';
use lib '/usr/local/ensembl67/ensembl-external/modules';
use lib '/usr/local/ensembl67/ensembl-functgenomics/modules';
use lib '/usr/local/bioperl1.2.3/bioperl-live';

#Diagnostics for debugging
use diagnostics;
#Use CGI
use CGI;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

#Main
#Define the species!
my $species = shift;

die("Species is empty") unless($species);

my $registry = &openEnsemblConnection;

#Make the folders for the new species
system("mkdir /var/www/html/FusionAnnotator/data/$species/ref/");
system("mkdir /var/www/html/FusionAnnotator/data/$species/ref/seq");
system("mkdir /var/www/html/FusionAnnotator/data/$species/ref/tracks");


#Get the chromosome seqs and format these into .JSON
&getRefSequences($registry, $species);

#Get all the genes, transcripts and domains from the ensembl DB
&getGenomeFeatures($registry, $species);
#Get the regulatory features from the ensembl DB if there is a funcgen DB for that species
if($species =~ /homo_sapiens|human|mouse|mus/i){
	getRegElements($registry, $species);
}

#Add the tracks to Jbrowse
&addTracks($species);

#Gets the sequences of the reference chromosomes and formats/add these to the correct species
sub getRefSequences{
	my ($registry, $species) = @_;
	my $slice_adaptor = $registry->get_adaptor( $species, 'Core', 'Slice' );

	#Get the chromosomes of the species
	my @slices = @{ $slice_adaptor->fetch_all('chromosome') };
	
	#Do this for each chromosome
	foreach my $chromosome (@slices){
		#Open file for writing
		open(my $seqFileHandler, '>', "/var/www/html/FusionAnnotator/data/$species/ref/seq.fasta") || die "Could not write or make feature File\n";
		#Print seq to file
		print $seqFileHandler ">".$chromosome->seq_region_name()."\n";
		print $seqFileHandler $chromosome->seq();
		#Close file
		close($seqFileHandler);
		#Add the reference seq to the species
		system("perl ../genomeBrowser/bin/prepare-refseqs.pl --fasta /var/www/html/FusionAnnotator/data/$species/ref/seq.fasta --out ../data/$species/ref/");
	}
	#Remove the file used for seq storage
	system("rm /var/www/html/FusionAnnotator/data/$species/ref/seq.fasta");

}


#Add the feature tracks to Jbrowse (To the correct species), also delete the allEnsemblGenomeFeatures.gff file after formatting to JSON
sub addTracks{
	$species = shift;
	system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff /var/www/html/FusionAnnotator/data/$species/ref/allEnsemblGenomeFeatures.gff --tracklabel Genes --type Gene,Exon --out ../data/$species/ref/");
	system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff /var/www/html/FusionAnnotator/data/$species/ref/allEnsemblGenomeFeatures.gff --tracklabel Transcripts --type Transcript --out ../data/$species/ref/");
	system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff /var/www/html/FusionAnnotator/data/$species/ref/allEnsemblGenomeFeatures.gff --tracklabel Domains --type Domain --out ../data/$species/ref/");
	system("rm /var/www/html/FusionAnnotator/data/$species/ref/allEnsemblGenomeFeatures.gff");
	#Only homo_sapiens has funcgen DB
	if($species =~ /homo_sapiens|human|mouse|mus/i){
		system("perl ../genomeBrowser/bin/flatfile-to-json.pl --gff /var/www/html/FusionAnnotator/data/$species/ref/EnsemblFuncgenFeatures.gff --tracklabel Methyl_AcylElements --key Methyl_AcylElements --type H3K27me3,H3K27me2,H3K27me1,H3K4me3,H3K4me2,H3K4me1,H3K9me3,H3K9me2,H3K9me1,H3K79me3,H3K79me2,H3K79me1 --out ../data/$species/ref/");
		system("rm /var/www/html/FusionAnnotator/data/$species/ref/EnsemblFuncgenFeatures.gff");
	}
}


#Gets specified enriched-sites with certain reg. elements from the Ensembl DB en stores them in a GFF format for use in Jbrowse
sub getRegElements{
	my ($registry, $species) = @_;
	
	#Open file for writing
	open(my $featFileHandler, '>', "/var/www/html/FusionAnnotator/data/$species/ref/allEnsemblGenomeFeatures.gff") || die "Could not write or make feature File\n";

	#Get slices of all chromosomes
	my $sa = $registry->get_adaptor($species, 'Core', 'Slice') || die "Could not get slice adaptor for species: $species";
	my @slices = @{ $sa->fetch_all('chromosome') } || die "Could not retrieve chromosomes of species: $species";
	my ($regType, @reg_feats);
	
	#Funcgen Adaptor
	my $regfeat_adaptor = $registry->get_adaptor($species, 'funcgen', 'regulatoryfeature') || "No funcgen for species: $species";

	#For each chromosome
	foreach my $slice(@slices){
		#Search reg features for first slice
		my @reg_feats = @{$regfeat_adaptor->fetch_all_by_Slice($slice)};
		foreach my $cell_rf (@{@reg_feats}){
			#Print the Parent reg. element
			print $featFileHandler join("\t", $cell_rf->seq_region_name, $species, "RegEle", $cell_rf->seq_region_start, $cell_rf->seq_region_end, ".", ".", ".", "ID=".$cell_rf->stable_id);
			print $featFileHandler "\n";
			#Open the regulatory feature to get all the attributes like histone modification and pol. activity
			my $rfs = $regfeat_adaptor->fetch_all_by_stable_ID($cell_rf->stable_id); 
			#Only store the DNAseI, H3K4, H3K9, H3K27 en H3K79 (mono, di & tri variants and also acetyl)
			foreach my $regEle(@{$rfs}){
				#Get all the attributes of the reg. ele
				foreach my $attr_feat (@{$regEle->regulatory_attributes()}){
					if($attr_feat->display_label =~ /H3K4|H3K9|H3K27|H3K79|DNase/i){
						my @split = split(" -", $attr_feat->display_label);
						$regType = $split[0];
						#Print the ID of the reg. ele region, what kind of ele and the start-end of the region
						print $featFileHandler join("\t", $attr_feat->seq_region_name(), $species, $regType, $attr_feat->seq_region_start, $attr_feat->seq_region_end, ".",".",".", "Name=".$attr_feat->display_label.";Parent=".$cell_rf->stable_id);
						print $featFileHandler "\n";
					}
				}
			}
		}
	}
	close($featFileHandler);
}


#Gets genes, transcripts and domains from the EnsEMBL DB and writes them into a file
sub getGenomeFeatures{
	my ($registry, $species )= @_;
	#File where the features will be stored
	open(my $featFileHandler, '>', "/var/www/html/FusionAnnotator/data/$species/ref/allEnsemblGenomeFeatures.gff") || die "Could not write or make feature File\n";
	
	#Get slices of every chromosome
	my $slice_adaptor = $registry->get_adaptor( $species, 'Core', 'Slice' ) || die "Could not get slice Adaptor for species: $species";
	my $transcriptAdaptor = $registry->get_adaptor( $species, 'Core', 'Transcript' ) || die "Could not get transcript Adaptor for species: $species";;
	
	# Retrieve slices of every chromosome in the database
	my @slices = @{ $slice_adaptor->fetch_all('chromosome') };
	
	my $geneStrand;
	
	#Get all genes, transcripts and domains
	foreach my $slice(@slices){
		for my $gene (@{ $slice->get_all_Genes() }){
			$geneStrand = "-" if $gene->strand() == -1;
			$geneStrand = "+" if $gene->strand() == 1;
			print $featFileHandler join("\t", $gene->seq_region_name(), $species,"Gene", $gene->seq_region_start(), $gene->seq_region_end(), "." ,$geneStrand, ".", "ID=".$gene->stable_id().";Name=".$gene->external_name().";URL=<a href = http://www.ensembl.org/$species/Gene/Summary?g=".$gene->stable_id()." target='blank'>Go to Ensembl</a>")."\n";
			#Get  transcripts
			while( my $transcript = shift @{ $gene->get_all_Transcripts() }){
				print $featFileHandler join("\t", $gene->seq_region_name(), $species, "Transcript", $transcript->seq_region_start(), $transcript->seq_region_end(), ".", $geneStrand , ".", "Name=".$transcript->stable_id().";Gene:=".$gene->external_name().";URL=<a href = http://www.ensembl.org/$species/Transcript/Summary?t=".$transcript->stable_id()." target='blank'>Go to Ensembl</a>")."\n";
				#Get Exons
				foreach my $exon ( @{ $transcript->get_all_Exons() } ) {
					print $featFileHandler join("\t", $gene->seq_region_name(), $species, "Exon", $exon->seq_region_start(), $exon->seq_region_end(), ".", $geneStrand , ".", "Parent=".$gene->stable_id())."\n";
				}
				#Get domains
				my $domainfeat = $transcriptAdaptor->fetch_by_stable_id($transcript->stable_id())->translation()->get_all_DomainFeatures if defined($transcript->biotype()) and $transcript->biotype() eq "protein_coding";
				while ( my $pfeature = shift @{$domainfeat} ) {
						print $featFileHandler join("\t", $gene->seq_region_name(), $species, "Domain", ($transcript->seq_region_start()+$pfeature->start()), ($transcript->seq_region_start()+$pfeature->end()), ".", $geneStrand , ".", "Name=".$pfeature->idesc().";Transcript:=".$transcript->stable_id().";URL=<a href = http://pfam.sanger.ac.uk/family/".$pfeature->display_id()." target='blank'>Go to PFAM</a>")."\n" if $pfeature->idesc() ne "";
				}
			}
		}
	}
	close($featFileHandler);
}

#Open connection (Registry) to Ensembl DB
sub openEnsemblConnection{
	my $registry = 'Bio::EnsEMBL::Registry';
	#Local Ensembl DB
		$registry->load_registry_from_db(
    -host => 'wgs10.op.umcutrecht.nl',
    -user => 'ensembl');
    	return $registry;
}
