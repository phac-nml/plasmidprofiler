#' Example Table of Blast Results
#'
#'
#' @docType data
#'
#' @usage data(blastdata)
#'
#' @format Dataframe.
#'
#' @keywords datasets
#'
#' @references None Yet
#' (\href{http://www.ncbi.nlm.nih.gov/}{PubMed})
#'
#' @source Strains graciously provided by the authors of the following papers:
#' Complete Genome and Plasmid Sequences of Three Canadian Isolates of
#' Salmonella enterica subsp. enterica Serovar Heidelberg from Human
#' and Food Sources. 2016 Labbe et al.
#' PMID: 26769926
#'
#' Complete Sequence of Four Multidrug-Resistant MOBQ1 Plasmids Harboring
#' blaGES-5 Isolated from Escherichia coli and Serratia marcescens
#' Persisting in a Hospital in Canada. 2015 Boyd et al.
#' PMID: 25545311
#'
#' Colistin-Nonsusceptible Pseudomonas aeruginosa Sequence Type 654 with
#' blaNDM-1 Arrives in North America. 2016 Mataseje et al.
#' PMID: 26824951
#'
#' @examples
#' data(blastdata)
"blastdata"

#' Example Table of SRST2 Results
#'
#'
#' @docType data
#'
#' @usage data(srst2data)
#'
#' @format Dataframe.
#'
#' @keywords datasets
#'
#' @references None Yet
#' (\href{http://www.ncbi.nlm.nih.gov/}{PubMed})
#'
#' @source Strains graciously provided by the authors of the following papers:
#' Complete Genome and Plasmid Sequences of Three Canadian Isolates of
#' Salmonella enterica subsp. enterica Serovar Heidelberg from Human
#' and Food Sources. 2016 Labbe et al.
#' PMID: 26769926
#'
#' Complete Sequence of Four Multidrug-Resistant MOBQ1 Plasmids Harboring
#' blaGES-5 Isolated from Escherichia coli and Serratia marcescens
#' Persisting in a Hospital in Canada. 2015 Boyd et al.
#' PMID: 25545311
#'
#' Colistin-Nonsusceptible Pseudomonas aeruginosa Sequence Type 654 with
#' blaNDM-1 Arrives in North America. 2016 Mataseje et al.
#' PMID: 26824951
#'
#' @examples
#' data(srst2data)
"srst2data"

#' Example Complete Report after the following steps.
#' Blast data from attached blastdata table
#' SRST2 data from attached srst2data table
#'
#' read_blast Import the blast file, add column names
#' blast_parser Parse imported file
#' amr_positives Detect AMR positive plasmids
#' read_srst2 Import SRST2 file
#' combine_results Combine SRST2 and Blast
#' zetner_score Add Sureness value
#' amr_presence Add detected AMR to report
#' order_report Arrange report
#'
#'
#' @docType data
#'
#' @usage data(report)
#'
#' @format Dataframe.
#'
#' @keywords datasets
#'
#' @references None Yet
#' (\href{http://www.ncbi.nlm.nih.gov/}{PubMed})
#'
#' @source Strains graciously provided by the authors of the following papers:
#' Complete Genome and Plasmid Sequences of Three Canadian Isolates of
#' Salmonella enterica subsp. enterica Serovar Heidelberg from Human
#' and Food Sources. 2016 Labbe et al.
#' PMID: 26769926
#'
#' Complete Sequence of Four Multidrug-Resistant MOBQ1 Plasmids Harboring
#' blaGES-5 Isolated from Escherichia coli and Serratia marcescens
#' Persisting in a Hospital in Canada. 2015 Boyd et al.
#' PMID: 25545311
#'
#' Colistin-Nonsusceptible Pseudomonas aeruginosa Sequence Type 654 with
#' blaNDM-1 Arrives in North America. 2016 Mataseje et al.
#' PMID: 26824951
#'
#' @examples
#' data(report)
"report"
