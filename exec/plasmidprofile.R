#' RScript capable
#'

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Plasmidprofiler"))


cl_arguments <- function(){
  # CL arguments ####
  option_list = list(
    make_option(c("-b", "--blastfile"), type="character", default=NULL,
                help="BLAST TSV file name", metavar="character"),
    make_option(c("-s", "--srst2file"), type="character", default=NULL,
                help="SRST2 TSV file name", metavar="character"),
    make_option(c("-u", "--sureness"), type="numeric", default=0.75,
                help="Sureness cut off, defaults to 0.75", metavar="numeric"),
    make_option(c("-c", "--coverage"), type="numeric", default=NA,
                help="Percent coverage cut off", metavar="numeric"),
    make_option(c("-l", "--length"), type="numeric", default=NA,
                help="Plasmid length cut off", metavar="numeric"),
    make_option(c("-a", "--anonymize"), action="store_true", default=NA,
                help="Anonymize plasmid and sample names"),
    make_option(c("-t", "--title"), type="character", default="Plasmids",
                help="Title of image [default= %default]", metavar="character")
  );

  opt_parser <- OptionParser(option_list=option_list);
  opt <- parse_args(opt_parser);

  if (is.null(opt$blastfile) | is.null(opt$srst2file)){
    print_help(opt_parser)
    stop("SRST2 and BLAST files must be supplied.", call.=FALSE)
  }
  opt
}


opt <- cl_arguments()

Plasmidprofiler::main(blast.file = opt$blastfile,
     srst2.file = opt$srst2file,
     coverage.filter = opt$coverage,
     sureness.filter = opt$sureness,
     length.filter = opt$length,
     anonymize = opt$anonymize,
     main.title = opt$title)
