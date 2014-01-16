#!/usr/bin/perl -w

# **********************************************************************
# **********************************************************************
# 
# Date: January 11, 2014
# The purpose of this program is to calculate the correlation between
# CCLE gene expression data and COLT RNAi knockdown scores across multiple
# cell lines.
# In processData a single gene's correlation is calculated against 
# every other gene for which there is data. Just set the 'geneToAnalyze'
# the 'sourceData' (i.e. what data to use for the geneToAnalyze), and
# the 'compareDataTo'. SourceData and CompareDataTo can be either 'CCLE' 
# or 'COLT' or 'GEO', but they cannot be the same. 'pValueCutoff' can be set to the 
# correlation p-value threshold you would like to outputted results to be limited to
# 'tumorTypesToConsider' is a hash that can be used to specify the types of 
# cell lines to consider (source tissue). If this parameter is omitted all 
# tissues will be considered.
# **********************************************************************
# **********************************************************************

use strict;
use Modules::UtilityFunctions;
use Modules::Maths; 

&processData( {'geneToAnalyze' => 'PLK1',
								'sourceData' => 'colt',
								'compareDataTo'=>'geo',
								'limitToYeastData'=>0,
								'pValueCutoff' => 1,
								# 'tumorTypesToConsider' => {'BREAST' => 1}
							});


sub processData{
	use Data::Dumper;
	my ($options) = @_;
	my ($geneToAnalyze, $compareDataTo, $sourceData, $limitToYeastData, $pValueCutoff, $yeastDataset);

	if($options->{'geneToAnalyze'}){	$geneToAnalyze = $options->{'geneToAnalyze'};	}
	else{	die "no geneToAnalyze";	}

	if($options->{'compareDataTo'}){	$compareDataTo = $options->{'compareDataTo'};	}
	else{	die "no compareDataTo";	}
	if($options->{'sourceData'}){	$sourceData = $options->{'sourceData'};	}
	else{	die "no sourceData";	}


	if($options->{'limitToYeastData'}){	
		$limitToYeastData = $options->{'limitToYeastData'};	
		if($options->{'yeastDataset'}){	$yeastDataset = $options->{'yeastDataset'};	}
		else{	die "If limitToYeastData is true, yeastDataset must be defined";	}
	}
	else{	$limitToYeastData = 0;}

	if($options->{'pValueCutoff'}){	$pValueCutoff = $options->{'pValueCutoff'};	}
	else{	$pValueCutoff = 0.05;}

	my %datasets = ("colt" => {'sub' => \&loadColtData, "param" => 'source_data/COLT_data/GARP-score.txt'}, 
									# "ccle" => {'sub' => \&load_CCLE_data, "param" => ['source_data/CCLE_data/rma - CCLE_Expression_2012-09-29.res']},
									"ccle" => {'sub' => \&load_CCLE_data, "param" => ['source_data/CCLE_data/rma_CCLE_Expression_Entrez_2012-10-18.res']}, 
									"geo" => {'sub' => \&load_GEO_data, "param" => ['source_data/GEO_data/rma_GSE12777_NUSE.txt']}
									);
	# what should the drug data be compared to? options are defined in %datasets
	
	if(!defined $datasets{$compareDataTo}){warn "Invalid value for compareDataTo ($compareDataTo), it should be one of the following:"; warn Dumper keys %datasets; exit(0);}
	if(!defined $datasets{$sourceData}){warn "Invalid value for sourceData ($sourceData), it should be one of the following:"; warn Dumper keys %datasets; exit(0);}
	if($compareDataTo eq $sourceData){
		warn "Source data ($sourceData) cannot be the same as compareDataTo ($compareDataTo)."; exit(0);
	}
		
	my $outfilename = "correlation of $options->{'geneToAnalyze'} $sourceData to all $compareDataTo data";

	my $data = &getRelevantSourceData($sourceData, $options->{'geneToAnalyze'}, \%datasets);
	my %cellLinesToConsider = ();
	foreach my $dataset(keys %{$data}){
		foreach my $cellLine(keys %{$data->{$dataset}}){
			$cellLinesToConsider{$cellLine}=1;
		}
	}
	# open dataset, store data in hash --> {dataset}->{'genes'}->{geneName}->{cellLine} = score
	# {'cellLinesToConsider'} = array of cell line names
	my $comparisonData = $datasets{$compareDataTo}->{'sub'}->(\%cellLinesToConsider, $datasets{$compareDataTo}->{'param'});
	# warn Dumper $comparisonData;
	my $count=0;
	my $start = 'Gene';
	if($limitToYeastData){	$start = "Yeast Gene\tHuman Gene";}
	$start = "$start\tSorting correlation\tSorting P-value\tActually Significant?\tSpearman Correlation\tSpearman P-value\tCorrelation\tP-Value\tCook Corrected Correlation\tP-value\tExplanation\t# cell lines\t";
	my $numTabs = 0;
	while ($start =~ /\t/g) { $numTabs++; }
	# warn "numTabs = $numTabs";

	foreach my $sourceDataset(keys %{$data}){
		foreach my $dataset(keys %{$comparisonData}){
			my @cellLinesToCompare=();
			if(defined $options->{'tumorTypesToConsider'}){
				foreach my $cl(sort keys %{$comparisonData->{$dataset}->{'cellLinesToConsider'} }){
					if(defined $options->{'tumorTypesToConsider'}->{$comparisonData->{$dataset}->{'cellLinesToConsider'}->{$cl}}){
						push(@cellLinesToCompare, $cl);
					}
				}
			}
			else{@cellLinesToCompare = sort keys %{$comparisonData->{$dataset}->{'cellLinesToConsider'}};}

			$start.=join("\t",@cellLinesToCompare)."\n";

			my $numberTumors = scalar(@cellLinesToCompare);
			if($numberTumors > 10){
				open (my $OUT, ">$outfilename - $count.tab") || die "Couldn't open $outfilename - $!\n";	
				print $OUT "Comparing\n$sourceDataset\nto\n$dataset\n\n";
				$count++;

				my @sourceValues=();
				foreach my $cl(@cellLinesToCompare){
					push(@sourceValues,$data->{$sourceDataset}->{$cl}); #log($scores->{$cl})/log(10));
				}

				print $OUT $start;
				print $OUT "$options->{'geneToAnalyze'} $sourceData Values:". ("\t" x $numTabs) .join("\t",@sourceValues)."\n\n";

				foreach my $gene(keys %{$comparisonData->{$dataset}->{'genes'}}){
					my @comparisonScores = ();
					foreach my $cellLine(@cellLinesToCompare){
						if(defined $comparisonData->{$dataset}->{'genes'}->{$gene}->{$cellLine}){
							push(@comparisonScores, $comparisonData->{$dataset}->{'genes'}->{$gene}->{$cellLine});
						}
						else{
							warn "dataset = $dataset, gene = $gene, cell line = $cellLine";
							warn $comparisonData->{$dataset}->{'genes'}->{$gene}->{$cellLine};
							die "major problemo";
						}
					} # end for each gene in comparison data

					my $outliers = &Maths::determineInfluentialPoints(\@sourceValues,\@comparisonScores);
					my ($outlierPearlmanCorrelation, $outlierPvalue, $outlierExplanation) = ('na',1,'na');
					if( keys %{$outliers}){
						my @temp1 = @comparisonScores;
						my @temp2 = @sourceValues;
						$outlierExplanation = '';
						#  sort in descending order so we do not run into any issues when we splice up the array
						my @indices = sort {$b<=>$a} keys %{$outliers};
						for(my $i=0;$i<@indices;$i++){
							splice(@temp1, $indices[$i],1);
							splice(@temp2, $indices[$i],1);
							$outlierExplanation .= "[Removed index $indices[$i], ($outliers->{$indices[$i]}->[0],$outliers->{$indices[$i]}->[1]) - cook dist = $outliers->{$indices[$i]}->[2]] ";
						}
						($outlierPearlmanCorrelation, $outlierPvalue) = &Maths::calculateCorrelation(\@temp1,\@temp2);
					}
					my ($correlation,$pValue) = &Maths::calculateCorrelation(\@comparisonScores,\@sourceValues);
					
					my ($sCorrelation,$sPvalue) = &Maths::calculateSpearmanCorrelation(\@comparisonScores,\@sourceValues);

					if($pValue < $pValueCutoff || $sPvalue < $pValueCutoff || $outlierPvalue < $pValueCutoff){
						my $actuallySignificant = 'yes';
						my $sortingCorrelation = "$correlation\t$pValue";
						if($outlierPearlmanCorrelation ne 'na'){
							if(abs($correlation-$outlierPearlmanCorrelation)>=0.2){
								if($outlierPvalue < $pValueCutoff){
									if(abs($outlierPearlmanCorrelation) > abs($correlation)){
										$actuallySignificant='yes - outliers in set (large)';
										$sortingCorrelation = "$outlierPearlmanCorrelation\t$outlierPvalue";
									}
								}
								else{
									$actuallySignificant='no';
									$sortingCorrelation="0.0\tna";
								}
							}
							elsif(abs($outlierPearlmanCorrelation) > abs($correlation)){
								$actuallySignificant='yes - outliers in set (small)';
								$sortingCorrelation = "$outlierPearlmanCorrelation\t$outlierPvalue";
							}
						}
						print $OUT "$gene\t$sortingCorrelation\t$actuallySignificant\t$sCorrelation\t$sPvalue\t$correlation\t$pValue\t$outlierPearlmanCorrelation\t$outlierPvalue\t$outlierExplanation\t$numberTumors\t".join("\t",@comparisonScores)."\n";
					}
				}
				close $OUT;
			}
			else{
				warn "Did not considered the dataset '$dataset', too few cell lines:";
				warn Dumper \@cellLinesToCompare;
			}
		}
	}
	warn "DONE!";
}

sub getRelevantSourceData {
	my ($sourceData, $gene, $datasets) = @_;
	my $data = $datasets->{$sourceData}->{'sub'}->(undef, $datasets->{$sourceData}->{'param'});
	
	# if this is defined then the datasource is colt, else it is CCLE and may contain several datasets
	my %data = ();
	foreach my $dataset(keys %{$data}){
		if(!defined $data->{$dataset}->{'genes'}->{uc($gene)}){
			warn "Could not find '$gene' in $sourceData dataset.";
			exit(0);
		}
		$data{$dataset}=$data->{$dataset}->{'genes'}->{uc($gene)};
	}
	return \%data;
}

sub loadColtData{
	# open COLT data, store data in hash --> {geneName}->{cellLine} = score
	my %coltData = ();
	my ($cellLinesToConsider, $filename) = @_;
	$filename = '../coltData.txt' if(!defined $filename);
	open (my $colt, "<$filename") || die "Couldn't open '$filename' $!\n";
	$/ = &UtilityFunctions::line_break_check( $colt );
	my $header = <$colt>;
	$header = uc($header);
	my @header = split(/\t/, $header);

	if(!($header[1] =~ /Gene/i && $header[3] eq 'DESCRIPTION' )){
		die "Could not parse colt data file header";
	}
	# warn "cell lines to consider = ";
	# warn Dumper $cellLinesToConsider;
	my (@indicesToConsider)=();
	for (my $i = 4; $i < @header; $i++) {	
		$header[$i] = &normalizeCellLineName($header[$i]);
		# if cellLinesToConsider is not defined, consider everything
		if(!$cellLinesToConsider){
			push(@indicesToConsider, $i);
			$coltData{'colt'}->{'cellLinesToConsider'}->{$header[$i]}=1;
		}
		else{
			# warn $header[$i];
			if(defined $cellLinesToConsider->{$header[$i]}){
				# warn "found it!";
				push(@indicesToConsider, $i);
				$coltData{'colt'}->{'cellLinesToConsider'}->{$header[$i]}=1;
			}
			else{$coltData{'colt'}->{'cellLinesNotConsidered'}->{$header[$i]}=1; }
		}
	}
	
	while(<$colt>){
		chomp;
		$_=~ s/"//g; # remove quotes
		$_=&UtilityFunctions::trimErroneousCharactersAtEnds($_);
		my @data = split(/\t/);
		# first index is gene, after that are the scores
		foreach my $i(@indicesToConsider){
			if(!defined $header[$i]){
				warn "colt problem!\ni = $i";
				warn scalar(@header);
				use Data::Dumper;
				warn Dumper \@header;
				warn Dumper \@data;
				exit(0);
			}
			#  if no gene name defined, default back to ref seq id
			if(!defined $data[1] || $data[1] eq ''){$data[1]=$data[0];}
			$coltData{'colt'}->{'genes'}->{uc($data[1])}->{$header[$i]} = $data[$i];
		}
	}
	close $colt;
	return \%coltData;
}

# load_GEO_data
# input: cellLinesToConsider == a hash reference that contains the cell lines we are interested in as keys
# 			 datasetsToAnalyze == an array reference containing the filenames of the datasets we are considering. Most of thie time the size of this array will be 1
# Processes GEO gene expression data (e.g. GSE12777) -- this subroutine is similar to load_CCLE_data
# data returned:
# {datasetName}->{'genes'}->{geneName}->{cellLine} = score
# {datasetName}->{'cellLinesToConsider'} = array of desired cell lines actually considered
# {datasetName}->{'desiredCellLinesNotConsidered'} = array of desired cell lines not considered (e.g. those in cellLinesToConsider that are not in a dataset)
# {datasetName}->{'cellLinesNotConsidered'} = array of other cell lines not considered (i.e. those in the dataset that are not in cellLinesToConsider)
# only pull data for cell lines in the $cellLinesToConsider array ref, if $cellLinesToConsidered is undef, consider everything
sub load_GEO_data {
	my ($cellLinesToConsider, $datasetsToAnalyze) = @_;
	# will load all files in the current directory preceding with 'RMA'
	# if an argument is specified it will search in the argument directory instead of the current dir
	if(!$datasetsToAnalyze){die "no GEO datasets to analyze";}
	
	my (%scores);
	foreach my $dataset(@{$datasetsToAnalyze}){
		open (my $datasetFile, "<$dataset") || die "in load_GEO_data → Couldn't open '$dataset' $!\n";	
		$/ = &UtilityFunctions::line_break_check( $datasetFile );
		&processGEOFile($dataset, $datasetFile, \%scores, $cellLinesToConsider);
	}
	return \%scores;
} 

# load_CCLE_data
# processes CCLE expression datafile
# data returned:
# {datasetName}->{'genes'}->{geneName}->{cellLine} = score
# {datasetName}->{'cellLinesToConsider'} = array of desired cell lines actually considered
# {datasetName}->{'desiredCellLinesNotConsidered'} = array of desired cell lines not considered (e.g. those in cellLinesToConsider that are not in a dataset)
# {datasetName}->{'cellLinesNotConsidered'} = array of other cell lines not considered (i.e. those in the dataset that are not in cellLinesToConsider)
# only pull data for cell lines in the $cellLinesToConsider array ref, if $cellLinesToConsidered is undef, consider everything
sub load_CCLE_data {
	my ($cellLinesToConsider, $datasetsToAnalyze) = @_;
	# will load all files in the current directory preceding with 'RMA'
	# if an argument is specified it will search in the argument directory instead of the current dir
	if(!$datasetsToAnalyze){$datasetsToAnalyze = &getCCLE_Datasets('source_data/CCLE_data');}
	
	# return hash with the following keys:
	# $cellline{'label'}->{ccleCellLineLabel}->{'primaryName'} = primary name of cell line
	# $cellline{'label'}->{ccleCellLineLabel}->{'site'} = site / tissue type of cell line
	# $cellline{'alias'}->{cellLineAlias}=ccleCellLineLabel
	my $ccleCellLines = &loadCCLE_CellLines('source_data/CCLE_data');
	my (%scores);
	foreach my $dataset(@{$datasetsToAnalyze}){
		open (my $datasetFile, "<$dataset") || die "in load_CCLE_data → Couldn't open '$dataset' $!\n";	
		$/ = &UtilityFunctions::line_break_check( $datasetFile );
		if($dataset =~ /\.res$/){	&processResFile($dataset, $datasetFile, \%scores, $ccleCellLines, $cellLinesToConsider);	}
		else{die "I don't know how to process $dataset";}
	}
	return \%scores;
}

# return hash with the following keys:
# $cellline{'label'}->{ccleCellLineLabel}->{'primaryName'} = primary name of cell line
# $cellline{'label'}->{ccleCellLineLabel}->{'site'} = site / tissue type of cell line
# $cellline{'alias'}->{cellLineAlias}=ccleCellLineLabel
sub loadCCLE_CellLines{
	my ($dir,$fileName) = @_;
	$dir = '.' if ! defined $dir; # use current dir if $dir is not defined.
	$fileName = 'source_data/CELL_LINE_annotations/all_ccle_Cell_Line_Annotations.txt' if ! defined $fileName;
	open(my $cellLines, "<$fileName") || die "Could not open '$fileName'$!\n";
	$/ = &UtilityFunctions::line_break_check(\*$cellLines);
	my $header = <$cellLines>;
	my %cellLines;
	foreach my $line(<$cellLines>){
		chomp($line);
		my @data = split("\t",uc($line));
		# $data[0] = ccle label for cell line
		# data[1] = primary cell name
		# data[2] == alias for cell line
		# data[3] == gender
		# data[4] == $site primate (i.e. tissue type)
		$data[1] = &normalizeCellLineName($data[1]);
		$cellLines{'label'}->{$data[0]}->{'primaryName'}=$data[1];
		$cellLines{'label'}->{$data[0]}->{'site'}=$data[4];
		$cellLines{'alias'}->{$data[2]} = $data[0];
		$cellLines{'siteTypes'}->{$data[4]}=1;
		$cellLines{'primaryName'}->{$data[1]} = $data[0];
	}
	return \%cellLines;
}

sub getCCLE_Datasets {
	my $dir = shift;
	$dir = '.' if ! defined $dir; # use current dir if $dir is not defined.
	opendir(TEMP_DIR, $dir) || die "Could not open dataset directory: $!";
	my @files=grep(/^rma/i, readdir TEMP_DIR); # exclude files starting with . and ..
	close TEMP_DIR;
	for (my $i = 0; $i < @files; $i++) {
		$files[$i] = $dir."/$files[$i]";
	}
	return \@files;
}

# processGEOFile
# inputs -- dataset = filename currently being processed
# 					datasetfile = file handle 
# 					$scores == hash reference that handles will store all the expression values
# 					$cellLinesToConsider == hash reference to look up the cell lines too consider data from
# function:
# 	process a tab delimited geo file that contains cell lines on the first line and expression values for the indicated probe on subsequent rows 
# 	(probe id = index 0 and each subsequent index is a score).
sub processGEOFile{
	my ($dataset, $datasetFile, $scores, $cellLinesToConsider) = @_;
	# for this type of file the first line is cell lines
	my $cellLineNames = '#';
	my $lineNum=0;
	# check for and skip over comment rows - there shouldn't be any, but just in case...
	while($cellLineNames =~ /^\#/){
		$cellLineNames = <$datasetFile>;
		$cellLineNames = &UtilityFunctions::trimErroneousCharactersAtEnds($cellLineNames);
		$lineNum++;
	}
	$cellLineNames = uc($cellLineNames);
	my @cellLineNames = split(/\t/, $cellLineNames);
	# first index of cellLineNames should be blank since this column of the file should contain probe ids (i.e. not scores)
	if($cellLineNames[0] && $cellLineNames[0] ne ''){	unshift(@cellLineNames, '');	}

	my %cellLinesNotConsidered = ();
	my %cellLineCancerTypes = ();
	my @columnsToConsider;
	
	my $considerAllCellLines = 0;
	if(!defined $cellLinesToConsider){	$considerAllCellLines=1;}

	HEAD:for (my $i = 1; $i < @cellLineNames; $i++) {
		# skip blanks...
		if(! defined $cellLineNames[$i]){next HEAD;}
		$cellLineNames[$i] = &normalizeCellLineName($cellLineNames[$i]);

		if($considerAllCellLines){
			$cellLinesToConsider->{$cellLineNames[$i]} = 1;
			push(@columnsToConsider,$i);	
			$cellLineCancerTypes{$cellLineNames[$i]}->{'all'}=1;
		}
		#  only consider those in hash
		else{
			if(!defined $cellLinesToConsider->{$cellLineNames[$i]}){
				$cellLinesNotConsidered{$cellLineNames[$i]}=1;
			}
			else{	
				push(@columnsToConsider,$i);
				$cellLineCancerTypes{$cellLineNames[$i]}->{'all'}=1;
			}
		}
	}
	my $probeTranslations = &loadProbeTranslations('source_data/GEO_data/hgu133plus2_probeToGene.txt');
	foreach my $probeScores(<$datasetFile>){
		$lineNum++;
		chomp($probeScores);
		$probeScores =~ s/"//g; # remove quotes
		$probeScores=uc($probeScores);
		my @data = split(/[\t|\s]/, $probeScores);
		my $gene = $probeTranslations->{$data[0]};
		foreach my $i(@columnsToConsider) {
			$scores->{$dataset}->{'genes'}->{$gene}->{$cellLineNames[$i]} = $data[$i];
		}
	}	
	$scores->{$dataset}->{'cellLinesToConsider'} = \%cellLineCancerTypes;
	$scores->{$dataset}->{'cellLinesNotConsidered'} = \%cellLinesNotConsidered;
	close $datasetFile;
	return 1;
}

# loadProbeTranslations
# simple subroutine that processes a probe to gene translation file and returns the results as a hash ref
# the file it processes should be a simple tab or space delimited file in which the 1st column contains a probeID and the 2nd
# contains the corresponding gene name.
# input: fileName
# returns: hash ref
sub loadProbeTranslations{
	my $fileName = shift;
	die "no probe translation fileName defined!" if (!defined $fileName);
	open (my $handle, "<$fileName") || die "in load_CCLE_data → Couldn't open '$fileName' $!\n";	
	$/ = &UtilityFunctions::line_break_check( $handle );
	my %probeLookup = ();
	foreach my $probeTranslation(<$handle>){
		chomp($probeTranslation);
		$probeTranslation =~ s/"//g; # remove quotes
		$probeTranslation=uc($probeTranslation);
		my @data = split(/[\t|\s]/, $probeTranslation);
		$probeLookup{$data[0]}=$data[1];
	}
	return \%probeLookup;
}

sub processResFile {
	my ($dataset, $datasetFile, $scores, $ccleCellLines, $cellLinesToConsider) = @_;
	# for these dataset files, the header line contains "name", "description", then cell line name. Each subsequent line has
	# the corresponding expression level for the indicated probe (probe id = index 0, description / gene name = index 1, each subsequent index is a score).
	my @header;
	my $headerFound = -100; # keep track of if the header has been found and how many rows we've iterated over (don't want to
	# get stuck in an infinate loop)
	my $lineNum = 0;
	while($headerFound < 0){
		$headerFound++;
		$lineNum++;
		my $header = <$datasetFile>; 
		$header = uc($header);
		$header = &UtilityFunctions::trimErroneousCharactersAtEnds($header);
		next if($header =~ /^\#/); # skip over comment rows
		@header = split(/\t/, $header);

		# header should contain the following column names
		if($header[1] eq 'ACCESSION' && ($header[0] eq 'DESCRIPTION' || $header[0] eq 'GENE')){
			$headerFound = 100;
			$lineNum++;$lineNum++;
			 # get rid of the 2 lines after the header;
			$header = <$datasetFile>; 
			$header = <$datasetFile>; 
		} #header found!
		# else, header not found
	}

	# if we exitted the while loop do to getting kicked out of an infinate loop, kill program
	if($headerFound < 10){die "Could not find $dataset header!";}


	my @allowedCellLineLabels = keys %{$ccleCellLines->{'label'}};

	# remove spaces from header, also remove anything after the 1st underscore
	my %columns;
	# my $oHeader = join("\t",@header);
	# open (my $OUT, ">cks1Data.tab") || die "Couldn't open cks1Data - $!\n";
	# print $OUT $oHeader."\n";
	my %cellLinesNotConsidered = ();
	my $correctedIndex = 4;
	# warn Dumper \@header;
	my $considerAllCellLines = 0;
	if(!$cellLinesToConsider){$considerAllCellLines=1;}
	HEAD:for (my $i = $correctedIndex; $i < @header; $i++) {	

		my @temp = split(/_/,$header[$i]);

		# skip blanks...
		if(! defined $temp[0]){next HEAD;}

		#  else we are good, proceed as we normally would
		my $cellLineName = &normalizeCellLineName($temp[0]);
		if(defined $ccleCellLines->{'primaryName'}->{$cellLineName}){
			# we're good, keep moving
		}
		elsif(defined $ccleCellLines->{'label'}->{$header[$i]}){
			#  this shouldn't be a situation
			$cellLineName = $ccleCellLines->{'label'}->{$header[$i]}->{'primaryName'};
		}
		else{	die "Error! could not find $header[$i] in ccle cell lines";	}
		if($considerAllCellLines){
			$cellLinesToConsider->{$cellLineName} = 1;
			$columns{$cellLineName}->{'col'}=$correctedIndex;	
		}
		#  only consider those in hash
		else{
			if(!defined $cellLinesToConsider->{$cellLineName}){
				$cellLinesNotConsidered{$cellLineName}=1;
			}
			else{	$columns{$cellLineName}->{'col'}=$correctedIndex;	}
		}
		
		$correctedIndex+=2;

		# make copy of header, remove the front portion, check to see of rest is an
		# actual site type → if not then this label is likely 2 cell line names stuck together,
		# so fix it.
		my $restOfName = $header[$i];
		$restOfName =~ s/^$temp[0]_//;
		if(!defined $ccleCellLines->{'siteTypes'}->{$restOfName}){
			# warn "Combo found! → '$restOfName'";
			# this is likely 2 cell lines smushed together
			my $foundIt = 0;
			foreach my $label(@allowedCellLineLabels){
				if($restOfName =~ /$label/){
					# found it!
					my $otherCellLineName = $ccleCellLines->{'label'}->{$label}->{'primaryName'};
					if(!defined $cellLinesToConsider->{$otherCellLineName}){
						$cellLinesNotConsidered{$otherCellLineName}=1;
					}
					else{	$columns{$otherCellLineName}->{'col'}=$correctedIndex;	}
					$correctedIndex+=2;
					$foundIt = 1;
				}
			}
			#  if the other cell line couldn't be found, kill the program
			if(! $foundIt){	die "Could not find a cell line for $restOfName in res file, $dataset";	}
		}
	}

	# iterate over $cellLinesToConsider, gather column numbers containing data for the cell lines we
	# wish to consider
	# %columnsToConsider is a hash whose keys are the column numbers we are considering and whose values define the 
	# cancer type of that column - actual cell line names are derived from @header
	my (%columnsToConsider,%desiredCellLinesNotConsidered) = ((),()); 
	foreach my $cellLine(keys %{$cellLinesToConsider}){
		if(defined $columns{$cellLine}){
			$columnsToConsider{$columns{$cellLine}->{'col'}} = $cellLine;
			delete $columns{$cellLine};
		}
		else{$desiredCellLinesNotConsidered{$cellLine}=1};
	}


	my %cellLineCancerTypes = ();
	foreach my $probeScores(<$datasetFile>){
		$lineNum++;
		chomp($probeScores);
		$probeScores =~ s/"//g; # remove quotes
		$probeScores=uc($probeScores);
		my @data = split(/[\t|\s]/, $probeScores);
		my $gene = $data[0];
		# if ($gene eq 'CKS1B'){
		# 	print $OUT "$probeScores\n";
		# }
		#warn Dumper \@data if $gene eq 'IRX3';
		# data starts at index #2
		foreach my $i(keys %columnsToConsider){
			# use header, store cancer cell type
			#$scores->{$dataset}->{'genes'}->{$gene}->{$columnsToConsider{$i}}->{$header[$i]}=$data[$i];
			# ignore cancer cell type
			# if($gene eq 'IRX3'){
			# 	warn "$gene\t$i\t$columnsToConsider{$i}\t$data[$i]";
			# }
			$scores->{$dataset}->{'genes'}->{$gene}->{$columnsToConsider{$i}}=$data[$i];
			$cellLineCancerTypes{$columnsToConsider{$i}}=$ccleCellLines->{'label'}->{$ccleCellLines->{'primaryName'}->{$columnsToConsider{$i}}}->{'site'};
		}
	}
	close $datasetFile;
	$scores->{$dataset}->{'cellLinesToConsider'} = \%cellLineCancerTypes;
	$scores->{$dataset}->{'desiredCellLinesNotConsidered'} = \%desiredCellLinesNotConsidered;
	$scores->{$dataset}->{'cellLinesNotConsidered'} = \%cellLinesNotConsidered;
	warn "done loading $dataset";
}

sub normalizeCellLineName{
	my $cellLine = shift;
	$cellLine =~ s/"+|\-+|\_+|\.+|\s+//g; # remove quotes and some other stuff
	return $cellLine;
}

