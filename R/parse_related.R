#  =============================================================================
#
# Copyright Government of Canada 2015-2016
#
# Written by: Adrian Zetner, Public Health Agency of Canada,
#     National Microbiology Laboratory
#
# Funded by the National Microbiology Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use
# this file except in compliance with the License. You may obtain a copy of the
# License at:
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed
# under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
# CONDITIONS OF ANY KIND, either express or implied. See the License for the
# specific language governing permissions and limitations under the License.
#
#  =============================================================================


#' Blast file import function
#'
#' This function imports the 25 column blast file and adds column headers
#'
#' @param br.file System location of the blast file, no default.
#' @return Dataframe of blast data with correct column headers.
#' @examples
#' \dontrun{
#' read_blast("/data/blast_results.tsv")
#' }
#' @importFrom stringr str_extract str_split_fixed
#' @importFrom utils read.table
#' @export

read_blast <- function(br.file){
  # Define column names, import blast output tabular file
  blast_cols <- c("qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
                "sallseqid",
                "score",
                "nident",
                "positive",
                "gaps",
                "ppos",
                "qframe",
                "sframe",
                "qseq",
                "sseq",
                "qlen",
                "slen",
                "salltitles")
  br <- read.table(br.file, sep = "\t", col.names = blast_cols)
  if (dim(br)[2] != 25){
    stop("Input blast file must be 25 columns")
  }else{
    br
  }
}


#' Blast Results Parser Function
#'
#' This function loads the imported blast results, extracts desired columns, Create new column of ratio between hit
#' length to query length - higher as denominator, adjusts pID by this ratio. Any AMR results are removed from the returned df.
#'
#' @param blast.results Blast results loaded from read_blast
#' @return Blast table with pID adjusted by ratio of hit length to query length (larger as denominator)
#' @examples
#' \dontrun{
#' blast_parser(blastdata)
#' }
#' @export

blast_parser <- function(blast.results){
  if (colnames(blast.results)[1] == "V1"){
    colnames(blast.results) <- c("qseqid",
                                 "sseqid",
                                 "pident",
                                 "length",
                                 "mismatch",
                                 "gapopen",
                                 "qstart",
                                 "qend",
                                 "sstart",
                                 "send",
                                 "evalue",
                                 "bitscore",
                                 "sallseqid",
                                 "score",
                                 "nident",
                                 "positive",
                                 "gaps",
                                 "ppos",
                                 "qframe",
                                 "sframe",
                                 "qseq",
                                 "sseq",
                                 "qlen",
                                 "slen",
                                 "salltitles")
  }
  #pull first 4 and 23/24 columns (qseqid, sseqid, pident, length, qlen, slen)
  adj.br <- blast.results[, c(1:4, 23, 24)]
  # Create new column of ratio between hit length to query length
  # (higher as denominator), use to adjust pID
  for (i in 1:nrow(adj.br)){
    if (adj.br$length[i] <= adj.br$qlen[i]){
      adj.br$ratio[i] <- adj.br$length[i] / adj.br$qlen[i]
    }else{
      adj.br$ratio[i] <- adj.br$qlen[i] / adj.br$length[i]
    }
  }

  #multiply percent ID by ratio and format properly
  adj.br$ADJpID <- adj.br$pident * adj.br$ratio
  adj.br$ADJpID <- format(adj.br$ADJpID, digits = 3)
  adj.br$ratio <- format(adj.br$ratio, digits = 3)
  adj.br$ADJpID <- as.numeric(adj.br$ADJpID)
  if (length(grep("AMR", adj.br$qseqid)) > 0){
    adj.br[-grep("(AMR)", adj.br$qseqid), ]
  }else{
    adj.br
  }

}

#' Identify Antimicrobial Resistance Positive Plasmids from Blast Results
#'
#' This function loads the imported blast results, identifies which plasmids
#' carry AMR genes at highest identity. May have issues with multiple genes per
#' plasmid, currently optimized for identifying one of two genes
#'
#' @param blast.results Blast results loaded from read_blast
#' @return Two column DF of plasmid names and genes present
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by select top_n bind_rows
#' @examples
#' \dontrun{
#' amr_positives(blastdata)
#' }
#' @export

amr_positives <- function(blast.results){

  # Find all AMR positives, select best matches only, append to pos.samples

  # Pull the best per plasmid
  blast.results <- (blast.results[grep("AMR", blast.results$qseqid), ] %>%
                      group_by(sseqid) %>%
                      select(sseqid, qseqid, pident) %>%
                      top_n(n = 1))
  blast.results$qseqid <- as.character(blast.results$qseqid)
  pos.samples <- data_frame()
  if (length(blast.results$qseqid) == 0){ # Return empty if none
    print("No match to AMR genes in DB")
    stop(pos.samples)
  }


  for (i in 1:length(blast.results$sseqid)){
    splt <- strsplit(blast.results$qseqid[i], split = ")")[[1]][2]
    splt <- strsplit(splt, split = "_")[[1]][1]
    tempplas <- as.character(blast.results$sseqid[i])
    tempgene <- paste(splt, blast.results$pident[i], sep = "_%")
    pos.samples <- bind_rows(pos.samples,
                             data_frame(Plasmid=tempplas, Gene=tempgene))
  }

  # This is dumb, replacing with better method
  # if (length(grep("AMR", blast.results$qseqid)) > 0){
  #   for (i in grep("AMR", blast.results$qseqid)){
  #     if (blast.results[i, 3] == 100){
  #       splt <- strsplit(blast.results[i, 1], split = ")")[[1]][2]
  #       splt <- strsplit(splt, split = "_")[[1]][1]
  #       pos.samples[ii, 1] <- as.character(blast.results[i, 2])
  #       pos.samples[ii, "Gene"] <- splt
  #       ii <- ii + 1
  #     }
  #   }
  #   pos.samples <- unique(pos.samples)
  #   pos.samples[, 1] <- as.factor(pos.samples[, 1])
  #   pos.samples[, 2] <- as.factor(pos.samples[, 2])
  #
  # }else {
  #   print("No match to AMR genes in DB")
  #   }

  pos.samples
}

#' SRST2 file import function
#'
#' This function imports the 14 column SRST2 file. Kind of superfluous
#'
#' @param srst2.file System location of the srst2 file, no default.
#' @return Dataframe of srst2 data with correct column headers.
#' @importFrom utils read.delim
#'
#' @examples
#' \dontrun{
#' read_srst2("/data/srst2_results.tsv")
#' }
#' @export
read_srst2 <- function(srst2.file){
  read.delim(srst2.file, sep = "\t", stringsAsFactors = FALSE) # New data
}

#' Combines SRST2 and Blast results into a single dataframe
#'
#' Cuts to desired columns, matches plasmids to BR and appends simplified INC names,
#' all future modifications are done to this dataframe
#'
#' @param br Blast results loaded from read_blast
#' @param sr SRST2 results loaded from read_srst2
#' @return Seven column dataframe of SRST2 results now including INC groups
#' @examples
#' \dontrun{
#' combine_results(example_srst2_results, example_blast_results)
#' }
#' @export

combine_results <- function(sr, br){
  # Create temp variable for reporting
  report <- sr[, c("Sample",
                  "gene",
                  "coverage",
                  "divergence",
                  "clusterid",
                  "length")]
  # Match plasmids to BR Inc Group and append to report
  report$Plasmid <- as.character(br$qseqid[match(report$gene, br$sseqid)])
  # Replace NAs (no Inc match) with Hyphen
  report$Plasmid[is.na(report$Plasmid)] <- "-"
  # Simplify names of Inc groups (remove the sequence data)
  report$Plasmid <- str_split_fixed(report$Plasmid, "_", 3)[, 1]
  colnames(report) <- c("Sample",
                        "Plasmid",
                        "Coverage",
                        "Divergence",
                        "Clusterid",
                        "Length",
                        "Inc_group")
  report <- report[, c("Sample",
                      "Plasmid",
                      "Inc_group",
                      "Coverage",
                      "Divergence",
                      "Length",
                      "Clusterid")]
  report
}
