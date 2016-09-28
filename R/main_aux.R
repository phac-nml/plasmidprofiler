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


# Create filecache environment to store filename outside of all functions
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome")
  file_cacher()
}

#' Main: Run everything
#'
#' Run all the interim functions to produce outputs.
#' \enumerate{
#'   \item \code{\link{read_blast}} Import the blast file, add column names
#'   \item \code{\link{blast_parser}} Parse imported file
#'   \item \code{\link{amr_positives}} Detect AMR positive plasmids
#'   \item \code{\link{read_srst2}} Import SRST2 file
#'   \item \code{\link{combine_results}} Combine SRST2 and Blast
#'   \item \code{\link{zetner_score}} Add Sureness value
#'   \item \code{\link{amr_presence}} Add detected AMR to report
#'   \item \code{\link{subsampler}} Apply filters to report
#'   \item \code{\link{order_report}} Arrange report
#'   \item \code{\link{save_files}} Save JPG and CSV
#'   \item \code{\link{create_plotly}} Creates plot
#'   \item \code{\link{save_files}} Save HTML plot
#' }
#'
#' @param blast.file Either system location of blast results (tsv) or dataframe
#' @param srst2.file Either system location of srst2 results (tsv) or dataframe
#' @param coverage.filter Filters results below percent read coverage specified (eg. 80)
#' @param sureness.filter Filters results below sureness specified (eg. 0.75)
#' @param length.filter Filters plasmid sequences shorter than length specified (eg. 10000)
#' @param combine.inc Flag to combine incompatibility sub-groups into their main type (set to 1)
#' @param plotly.user Enter your plotly info to upload to (\href{https://plot.ly/feed/}{Plotly})
#' @param plotly.api Enter your plotly info to upload to (\href{https://plot.ly/feed/}{Plotly})
#' @param post.plotly Flag to post to (\href{https://plot.ly/feed/}{Plotly})
#' @param anonymize Flag to post to anonymize plasmids and samples (set to 1)
#' @param main.title A title for the figure
#' @return Saves output files in working directory
#' @examples
#' main(blastdata,
#' srst2data,
#' coverage.filter=NA,
#' sureness.filter=0.75,
#' length.filter=10000,
#' main.title="Example Results")
#' @export
main <- function(blast.file,
                 srst2.file,
                 coverage.filter=NA,
                 sureness.filter = NA,
                 length.filter = NA,
                 combine.inc=NA,
                 plotly.user,
                 plotly.api,
                 post.plotly=NA,
                 anonymize=NA,
                 main.title="Plasmid Profiles") {

  file_cacher()

  if (typeof(blast.file) == "character"){
    blast.file <- read_blast(blast.file)
  }

  blast_results <- blast_parser(blast.file)
  pos.samples <- amr_positives(blast.file)

  if (typeof(srst2.file) == "character"){
    srst2.file <- read_srst2(srst2.file)
  }
  cr <- combine_results(srst2.file, blast_results)
  report <- zetner_score(cr)
  report <- amr_presence(report, pos.samples)
  report <- subsampler(report,
                       cov.filter = coverage.filter,
                       sure.filter = sureness.filter,
                       len.filter = length.filter,
                       inc.combine = combine.inc)
  report <- order_report(report, anonymize)
  save_files(report,
             plot.jpg = 1,
             report.csv = 1,
             title = main.title)
  create_plotly(report,
                user = plotly.user,
                api.key = plotly.api,
                post = post.plotly,
                title = main.title)
  save_files(report, webpage = 1, title = main.title)
}

#' Save Files Produced
#'
#' This function uses RColorBrewer to produce palettes based
#' on the factor levels of the identified column in a report.
#'
#' @param report Dataframe of results
#' @param plot.jpg Do you want to save a jpg? (Anything but NA)
#' @param report.csv Do you want to save a text report? (Anything but NA)
#' @param webpage Do you want to save an interactive heatmap as html? (Anything but NA)
#' @param title Enter a title for the plot
#' @return Named vector of colours, names are factor levels of column supplied
#' @import ggplot2
#' @import dplyr
#' @importFrom plotly ggplotly plotly_POST as.widget
#' @importFrom htmlwidgets saveWidget
#' @importFrom utils write.csv
#' @examples
#' \dontrun{
#'  save_files(report, plot.jpg=1, report.csv=1, webpage=NA)
#' }
#' @export
save_files <- function(report,
                       plot.jpg = NA,
                       report.csv = NA,
                       webpage = NA,
                       title = "Plasmid Profiles" ){

  filename <- get("name", envir = filecache)

  if (!is.na(plot.jpg)){
    g <- create_grob(report, grob.title = title)
    ggsave(paste(filename, ".jpg", sep = ""), g, device = "jpg", width = 12)
  }

  if (!is.na(report.csv)){
    report <- arrange(report, Sample, Inc_group, desc(Sureness))
    write.csv(report[, c(1:9)], paste(filename, ".csv", sep = ""))
  }

  # Write offline HTML object
  if (!is.na(webpage)){
    ppp <- create_plotly(report)
    htmlwidgets::saveWidget(as.widget(ppp), paste(filename, ".html", sep = ""))
  }
}

#' Normalize
#'
#' Normalizes a vector of values to a range of 0-1
#'
#' @param x Vector of values
#' @return Normalized vector of values
#' @examples
#' \dontrun{
#'  normalize(x)
#'  }
#' @export
normalize <- function(x){
  (x - min(x)) / (max(x) - min(x))
}

#' Minmax
#'
#' Takes two columns of numerical data,
#' normalizes it to ranges from 0 to 1 (0 to -1 for minimums),
#' sums them, arranges by sum, then returns the sorted dataframe
#'
#' @param df Dataframe
#' @param maxcol Column to normalize from 0 to 1
#' @param mincol Column to normalize from 0 to -1
#' @return Dataframe sorted by sum of maxcol and mincol
#' @importFrom dplyr arrange
#' @examples
#' \dontrun{
#'  minmax(df, maxcol, mincol)
#'  }
#' @export

# returns the sorted dataframe
minmax <- function(df, maxcol, mincol){
  mincolnorm <- - normalize(df[, mincol])
  maxcolnorm <- normalize(df[, maxcol])
  mmsum <- maxcolnorm + mincolnorm
  df$mmsum <- mmsum
  df <- arrange(df, desc(mmsum))
  df
}

#' Filecacher
#'
#' Creates filecache environment if needed
#'
file_cacher <- function(){
  if (!exists("filecache")){
    packageStartupMessage("No filecache found, creating...")
    filecache <<- new.env(parent = .GlobalEnv)
    filename <- "P2Run"
    assign("name", filename, envir = filecache)
  }
}
