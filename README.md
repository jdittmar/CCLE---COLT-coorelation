CCLE---COLT-coorelation
=======================


 The purpose of this program is to calculate the correlation between
CCLE gene expression data and COLT RNAi knockdown scores across multiple cell lines.
In processData a single gene's correlation is calculated againstevery other gene for which there is data. Just set the 'geneToAnalyze' the 'sourceData' (i.e. what data to use for the geneToAnalyze), and the 'compareDataTo'. SourceData and CompareDataTo can be either 'CCLE' or colt, but they cannot be the same. 'pValueCutoff' can be set to the correlation p-value threshold you would like to outputted results to be limited to 'tumorTypesToConsider' is a hash that can be used to specify the types of cell lines to consider (source tissue). If this parameter is omitted all tissues will be considered.

Data sources:

CCLE
    http://www.broadinstitute.org/ccle/data/browseData?conversationPropagation=begin

COLT
    http://colt.ccbr.utoronto.ca/cancer/download.html
