# Support functions


# Create filecache environment to store filename outside of all functions
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome")
  file_cacher()
}

# On load fix the notes:

.onLoad <- function(libname, pkgname){
  utils::globalVariables(c(".", "AMR_gene", "Coverage", "Inc_group", "Length",
                    "Plasmid", "Sample", "Sureness", "average", "filecache",
                    "pident", "qseqid", "sseqid", "x", "xend", "y", "yend"))
  file_cacher()
}



