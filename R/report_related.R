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


#' Adds the Zetner Score column to report
#'
#' Runs mimmax function on Coverage and Divergence, returns sum of normalized Coverage with negative normalized Divergence
#' a value which is then normalized from 0 to 1.
#'
#' @seealso \code{\link{subsampler}}, \code{\link{combine_results}}
#' @param report Dataframe of results produced by \code{\link{subsampler}} or \code{\link{combine_results}}
#' @return Report with zetner score added
#' @examples
#' \dontrun{
#' zetner_score(report)
#' }
#' @export
zetner_score <- function(report){
  # Combined maximum of coverage and minimum of divergence
  report <- report %>% minmax("Coverage", "Divergence") %>% group_by(Sample)
  report <- arrange(report, Sample)
  report$mmsum <- normalize(report$mmsum) # Improves final HM
  report <- ungroup(report)
  plyr::rename(report, c(mmsum = "Sureness"))
}


#' Adds the AMR_gene column to report
#'
#' Appends the results of amr_positives to the report
#' in column AMR_gene, missing have "-" instead
#'
#' @seealso \code{\link{subsampler}}, \code{\link{combine_results}}
#' @param report Dataframe of results produced
#' by \code{\link{subsampler}} or \code{\link{combine_results}}
#' @param pos.samples Two column DF of plasmid names and genes present produced
#' by \code{\link{amr_positives}}
#' @return Report with AMR_genes added
#' @examples
#' \dontrun{
#' amr_presence(report, pos.samples)
#' }
#' @export
amr_presence <- function(report, pos.samples){
  # Mark plasmids with AMR present and type
  report$AMR_gene <- as.character(pos.samples$Gene[match(report$Plasmid,
                                                         pos.samples$Plasmid)])
  report$AMR_gene[is.na(report$AMR_gene)] <- "-"
  report
}


#' Subsetting Results
#'
#' \preformatted{Several filters can be applied:
#'    Coverage: Filters results below percent read coverage specified
#'                eg. 95.9 cuts results where reads covered less than 95.9\% of the total length
#'    Sureness: Filters results below sureness specified
#'                eg. 0.9 cuts results where the sureness falls below 0.9
#'    Length:   Filters plasmid sequences shorter than length specified
#'                eg. 10000 cuts out results where the plasmid was less than 10kb
#'    Incompatibility groups can also be combined (eg. Fii(S) and Fii(K) are combined into Fii)}
#'
#' @seealso \code{\link{subsampler}}, \code{\link{combine_results}}
#' @param report Dataframe of results produced
#' by \code{\link{subsampler}} or \code{\link{combine_results}}
#' @param cov.filter Filters results below percent read coverage specified (eg. 80)
#' @param sure.filter Filters results below sureness specified (eg. 0.75)
#' @param len.filter Filters plasmid sequences shorter than length specified (eg. 10000)
#' @param inc.combine Flag to ombine incompatibility sub-groups into their main type (set to 1)
#' @return Report with filters applied
#' @importFrom gdata drop.levels
#' @importFrom stringr str_split_fixed
#' @examples
#' \dontrun{
#' subsampler(report, sureness.filter = 0.75, len.filter = 10000)
#' }
#' @export
subsampler <- function(report,
                       cov.filter = NA,
                       sure.filter = NA,
                       len.filter = NA,
                       inc.combine = NA){

  if (!is.na(cov.filter)){
    filename <<- paste(filename, "_cov", cov.filter, sep = "")
    report <- report[report$Coverage > cov.filter, ]
  }
  if (!is.na(sure.filter)){
    filename <<- paste(filename, "_sure", sure.filter, sep = "")
    report <- report[report$Sureness > sure.filter, ]
  }
  if (!is.na(len.filter)){
    filename <<- paste(filename, "_len", len.filter, sep = "")
    report <- report[report$Length > len.filter, ]
  }
  if (!is.na(inc.combine)){
    # Simplify names of Inc groups (remove the subsetting)
    report$Inc_group <- str_split_fixed(report$Inc_group, "\\(", 2)[, 1]
    # Replace all individual Col-type plasmids with just Col
    report$Inc_group[grep("Col", report$Inc_group)] <- "Col"
  }
  drop.levels(report)
}


#' Create Dendrogram Based on Plasmid Content
#'
#' Reads report, converts to matrix of Sample ~ Plasmid with Sureness as cell values.
#' Performs a hierarchical cluster analysis on a set of dissimilarities derived from the matrix.
#' Creates a dendrogram from this data. Returns either the HC data or the dendrogram plot
#'
#' @seealso \code{\link{subsampler}}, \code{\link{combine_results}}
#' @param report Dataframe of results produced
#' by \code{\link{subsampler}} or \code{\link{combine_results}}
#' @param hc.only Flag to return only hierarchical clustering
#' results instead of dendrogram plot (set to 1)
#' @return Dendrogram object or hierarchical clustering results
#' @examples
#' \dontrun{
#' tree_maker(report)
#' }
#' @importFrom reshape2 dcast
#' @importFrom ape as.phylo
#' @importFrom ggplot2 ggplot
#' @importFrom grid unit
#' @importFrom ggdendro dendro_data segment theme_dendro
#' @importFrom stats as.dendrogram dist hclust
#' @export
tree_maker <- function(report, hc.only = NA){
  reportable.wide <- dcast(report, Sample ~ Plasmid, value.var = "Sureness")
  reportable.wide[is.na(reportable.wide)] <- 0

  # Put results into matrix format for HM generation
  reportable.matrix <- data.matrix(reportable.wide[, 2:ncol(reportable.wide)])
  # Use sample names from col 1 of DCW for rownames in matrix
  rnames <- reportable.wide[, 1]
  rownames (reportable.matrix) <- rnames


  reportable.hc <- hclust(dist(reportable.matrix))

  tree.data <- dendro_data(as.dendrogram(reportable.hc), type = "rectangle")
  tree <- ggplot(segment(tree.data)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    coord_flip() +
    scale_y_reverse(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0.001, 0)) +
    theme_dendro() +
    if (nlevels(report$Sample) < 20){
      # Add some padding around the upper and lower edges of the tree.
      # 0.0027 per sample
      theme(plot.margin =
              unit(c(0.0027 * nlevels(report$Sample),
                     0,
                     0.0027 * nlevels(report$Sample),
                     0),
                   "null"))
    }else{
      # Add some padding around the upper and lower edges of the tree.
      # 0.0015 per sample
      theme(plot.margin = unit(c(0.00015 * nlevels(report$Sample),
                                 0,
                                 0.00015 * nlevels(report$Sample),
                                 0),
                               "null"))
    }
  if (!is.na(hc.only)){
    reportable.hc
  }else{
    tree
  }
}


#' Order the Report
#'
#' Order the report first by sample order (tree), then by incompatibility group, then by sureness on each plasmid
#'
#' @seealso \code{\link{subsampler}}, \code{\link{combine_results}}
#' @param report Dataframe of results produced by \code{\link{subsampler}} or \code{\link{combine_results}}
#' @param anonymize Flag to anything other than NA to replace plasmid and sample names with generic names
#' @return Ordered report
#' @examples
#' \dontrun{
#' order_report(report)
#' }
#' @importFrom plyr mapvalues
#' @importFrom dplyr arrange mutate group_by ungroup
#' @importFrom magrittr %>%
#' @export
order_report <- function(report, anonymize = NA){

  reportable.hc <- tree_maker(report, hc.only = 1)
  reportable.phylo <- as.phylo(reportable.hc)

  # First order Samples based on HC / Phylo
  report$Sample <- ordered(report$Sample,
                    levels = reportable.phylo$tip.label[reportable.hc$order])

  report <- report %>%
    group_by(Plasmid) %>%
    mutate(average = mean(Sureness))

  report <- ungroup(report)

  # Arrange by Inc_group then average Sureness
  report <- arrange(report, Inc_group, average)

  # Order the Plasmids based on order of appearance (ie. by inc group)
  report$Plasmid <- ordered(report$Plasmid,
                            levels = unique(report$Plasmid))

  report$Inc_group <- as.factor(report$Inc_group)
  report$AMR_gene <- as.factor(report$AMR_gene)

  if (!is.na(anonymize)){
    report$Sample <- mapvalues(report$Sample,
                               from = levels(report$Sample),
                               to = paste("Sample",
                                          1:length(levels(report$Sample))))
    report$Plasmid <- mapvalues(report$Plasmid,
                                from = levels(report$Plasmid),
                                to = paste("Plasmid",
                                           1:length(levels(report$Plasmid))))
    filename <- paste(filename, "_anon", sep = "")
  }

  subset(report, select = -c(average))
}
