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


#' Defining Colours Based on a Column of Data
#'
#' \preformatted{
#' This function uses RColorBrewer to produce palettes based on the factor levels of the identified column in a report.
#'    }
#'
#' @param report Dataframe of results produced by \code{\link{subsampler}} or \code{\link{combine_results}}
#' @param column Specify a column by name
#' @return Named vector of colours, names are factor levels of column supplied
#' @importFrom viridis inferno
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @examples
#' \dontrun{
#' define_colours(report, "AMR_gene")
#' }
#' @export
define_colours <- function(report, column){
  report[[column]] <- as.factor(report[[column]])
  levs <- levels(report[[column]])
  if (column == "AMR_gene"){
    colours <- colorRampPalette(inferno(20, end = 0.7))(length(levs))
  }else{
    colours <- colorRampPalette(brewer.pal(9, "Set1"))(length(levs))
  }
  names(colours) <- levs
  colours
}

#' Create GGPLOT Heatmap
#'
#' Using a ggplot2 tile geometry this function will create a heatmap of values in the report
#' coloured by incompatibility group, with alpha values from the sureness Zetner score. The order of
#' samples is determined by order_reporty and plasmids by incompatibility group and Zetner score.
#'
#' @param report Dataframe of results
#' @param len.highlight If anything but NA will highlight the largest plasmid hit per incompatibility group
#' @return GGPLOT plotted heatmap
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' plot_heatmap(report)
#' }
#' @export
plot_heatmap <- function(report, len.highlight=NA){

  colours.amr <- define_colours(report, "AMR_gene")
  colours.inc <- define_colours(report, "Inc_group")

  cnames <- colnames(report)
  cnames[8] <- "Sureness"
  colnames(report) <- cnames

  font.size <-  0.9 - 0.001 * length(levels(report$Plasmid))

  pp <- ggplot(report, aes(Plasmid,
                           Sample,
                           alpha = Sureness,
                           fill = Inc_group,
                           text = paste("AMR Gene: ",
                                        report$AMR_gene))
               ) +
    geom_tile(colour = "white") +
    scale_fill_manual(values = colours.inc, name = "Incompatibility Group") +
    guides(fill=guide_legend(ncol=2)) +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.45,
                                     size = rel(font.size)),
          axis.text.y = element_text(size = rel(0.8)),
          #axis.title.x=element_blank(),
          legend.title = element_text(),
          #axis.title.y=element_blank(),
          panel.background = element_rect(fill = "grey95")) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))

  #Add colours based on AMR for plasmid names

  report$Colours2 <- colours.amr[match(report$AMR_gene,
                                       levels(report$AMR_gene))]
  colours.amr2 <- report$Colours2[match(levels(report$Plasmid),
                                        report$Plasmid)]

  pp <- pp + theme(axis.text.x = element_text(colour = colours.amr2)) +
    geom_text(aes(label = "", colour = AMR_gene)) +
    scale_colour_manual(values = colours.amr, name = "Resistance Gene")

  if (!is.na(len.highlight)){
    seq.lengths <- report[report$Inc_group != "-", ] %>%
      group_by(Sample, Inc_group) %>%
      filter(Coverage >= 99) %>%
      arrange(Sample, Inc_group, desc(Length)) %>%
      slice(c(1))
    seq.lengths$maxed <- "Attention"
    report <- full_join(seq.lengths, report, by = colnames(report))
    pp <- pp + geom_point(aes(x = Plasmid, y = Sample),
                          seq.lengths[, 1:2],
                          inherit.aes = F,
                          alpha = 0.5,
                          size = 0.5
                          )
  }
  pp
}

#' Create Heatmap Graphical Object
#'
#' Combines the tree, heatmap, and titles to create final heatmap image.
#'
#' @param report Dataframe of results
#' @param grob.title Title of heatmap
#' @return Composite image
#' @import ggplot2
#' @import grid
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @import gtable
#'
#' @examples
#' \dontrun{
#' create_grob(report, grob.title="Plasmid Profiles")
#' }
#' @export
create_grob <- function(report, grob.title = "Plasmid Profiles"){
  pp <- plot_heatmap(report)
  tree <- tree_maker(report)
  # Build a graphical property table of the plot / tree
  pp.gt <- ggplot_gtable(ggplot_build(pp))
  tree.gt <- ggplot_gtable(ggplot_build(tree))

  # Get maximum widths and heights
  max.height <- unit.pmax(pp.gt$heights[4:5], tree.gt$heights[4:5])

  # Set the maximums in the gtables for tree and pp
  pp.gt$heights[4:5] <- as.list(max.height)
  tree.gt$heights[4:5] <- as.list(max.height)

  # Create a new gtable
  gt <- gtable(widths = unit(c(1, 7), "null"), heights = unit(c(7), "null"))

  title <- textGrob(grob.title, gp = gpar(fontsize = 24))

  mods <- get("mods", envir = filecache)

  footnote <- textGrob(mods,
                       x = 0,
                       hjust = 0,
                       gp = gpar(fontface = "italic"))

  padding <- unit(0.5, "line")

  # Insert tree and heatmap into the new gtable
  gt <- gtable_add_grob(gt, tree.gt, 1, 1) # Dend
  gt <- gtable_add_grob(gt,
                        pp.gt,
                        t = 1,
                        l = ncol(gt),
                        b = 1,
                        r = ncol(gt)) # Heatmap

  # Rows for the title / footnote
  gt <- gtable_add_rows(gt, heights = grobHeight(title) + padding, pos = 0)
  gt <- gtable_add_rows(gt, heights = grobHeight(footnote) + padding)
  gt <- gtable_add_grob(gt,
                        list(title, footnote),
                        t = c(1, nrow(gt)), l = c(1, 2), r = ncol(gt))

  # And render the plot
  grid.newpage()
  grid.draw(gt)

  g <- grid.arrange(gt)
  g
}


#' Create Plotly Object
#'
#' Combines the tree, heatmap, and titles to create final heatmap image.
#'
#' @param report Dataframe of results
#' @param user User ID for plotly web publishing
#' @param api.key API key for plotly web publishing
#' @param post Flag determines whether or not to post to plotly (default NA, no post)
#' @param title Title of heatmap
#' @param len.highlight If anything but NA will highlight the largest plasmid hit per incompatibility group
#' @return plotly object
#' @importFrom plotly ggplotly plotly_POST config
#' @importFrom graphics layout
#' @import ggplot2
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' create_plotly(report, title="Outbreak Plasmid Profiles")
#' }
#' @export
create_plotly <- function(report,
                          user,
                          api.key,
                          post=NA,
                          title="Plasmid Profiles",
                          len.highlight=NA){

  colours.amr <- define_colours(report, "AMR_gene")
  colours.inc <- define_colours(report, "Inc_group")
  report$Colours2 <- colours.amr[match(report$AMR_gene,
                                       levels(report$AMR_gene))]
  colours.amr2 <- report$Colours2[match(levels(report$Plasmid),
                                        report$Plasmid)]

  pp.noalpha <- ggplot(report,
                       aes(Plasmid,
                           Sample,
                           label = AMR_gene,
                           fill = Inc_group,
                           text = paste("Sureness: ",
                                        round(report$Sureness, 2)))) +
    geom_tile(colour = "white", alpha = 1) +
    scale_fill_manual(values = colours.inc, name = "Inc Group") +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.45,
                                     size = rel(0.7)),
          axis.text.y = element_text(size = rel(0.7)),
          axis.title.x = element_blank(),
          legend.title = element_text(),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "grey95")) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))

  # if (!is.na(len.highlight)){
  #   seq.lengths <- report[report$Inc_group != "-", ] %>%
  #     group_by(Sample, Inc_group) %>%
  #     filter(Coverage >= 99) %>%
  #     arrange(Sample, Inc_group, desc(Length)) %>%
  #     slice(c(1))
  #   seq.lengths$maxed <- "Attention"
  #   report <- full_join(seq.lengths, report, by = colnames(report))
  #   pp.noalpha <- pp.noalpha +
  #     geom_point(aes(x = Plasmid,
  #                    y = Sample),
  #                seq.lengths[, 1:2],
  #                inherit.aes = F,
  #                alpha = 0.5,
  #                size = 0.5)
  # }

  pp.noalpha <- pp.noalpha +
    theme(axis.text.x = element_text(colour = colours.amr2)) +
    guides(fill=guide_legend(ncol=2))

  pp.noalpha <- pp.noalpha +
    geom_tile(aes(x = Plasmid,
                   y = Sample,
                   label = AMR_gene,
                   fill = Inc_group,
                   text = paste("Sureness: ",
                                round(Sureness, 2))),
                    inherit.aes = F,
                    alpha = 0.001,
                    width = 0.99,
                    height = 0.99,
                    show.legend = FALSE,
                    colour = "white")

  # Originally included to show which have AMR gene present. Not working as desired
  # pp.noalpha <- pp.noalpha +
  #   geom_point(aes(x = Plasmid,
  #                  y = Sample,
  #                  label = AMR_gene,
  #                  text = paste("Sureness: ",
  #                               round(Sureness, 2))),
  #              data = subset(report, AMR_gene != "-"),
  #              inherit.aes = F,
  #              alpha = 0.25,
  #              size = 1)

  m <- list(l = 200,
           r = 100,
           b = 200,
           t = 100,
           pad = 10)

  f <- list(family = "Arial, Helvetica, sans-serif",
            size = 18,
            color = "#7f7f7f")

  hacktitle <- paste(".", paste(rep(" ", 75), collapse=""),
                     "Plasmid",
                     paste(rep(" ", 75), collapse=""),
                     ".",
                     collapse="")
  x <- list(title = hacktitle,
            titlefont = f)

  y <- list(title = "Sample",
            titlefont = f)

  # Make ggplotly object
  ppp <- ggplotly(pp.noalpha, tooltip = c("y",
                                          "x",
                                          "text",
                                          "label",
                                          "fill")) %>%
    plotly::layout(autosize = F,
                   width = 1600,
                   height = 1000,
                   font = f,
                   margin = m,
                   xaxis = x,
                   yaxis = y,
                   title = title) %>%
    plotly::config(., displaylogo = FALSE,
                   modeBarButtonsToRemove = list('sendDataToCloud',
                                                 'select2d',
                                                 'lasso2d',
                                                 'hoverClosestCartesian',
                                                 'hoverCompareCartesian'))
                                                 # 'toImage'))


  # Publish
  if (!is.na(post)){
    Sys.setenv("plotly_username" = user)
    Sys.setenv("plotly_api_key" = api.key)
    filename <- get("name", envir = filecache)
    plotly_POST(ppp, filename)
  }
  ppp
}
