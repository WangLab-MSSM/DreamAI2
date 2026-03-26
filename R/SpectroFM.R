# requireNamespace("Rcpp")

#' @title Imputation of Missing Protein Abundances using Matrix Factorization
#' @description The function impute.SpectroFM imputes a dataset with missing values or NA's using factorization
#' @param input_table A data frame with missing values or Nas
#' @param iter Number of iterations
#' @param verbose verbose in mcmc procedure
#'
#' @return Imputed data frame
#' @export
#'
#' @examples
#' \dontrun{
#' data<-data.DIA[1:100,1:50]
#' impute.SpectroFM(input_table=as.data.frame(data), iter = 40)
#' }
#' @importFrom Rcpp sourceCpp

impute.SpectroFM <- function(input_table, iter, verbose=FALSE){
  package.dir <- find.package(utils::packageName())
  print(package.dir)
  libfm_cpp_path = paste(sep="", package.dir, "/src/src/libfm/libfm.cpp")
  if(!file.exists(libfm_cpp_path)) {
    stop(paste("[Error] LibFM C++ file path'", libfm_cpp_path, "' is wrong./n",
               "/tCurrent Directory : ", getwd()))
  }

  verbose_int = "0"
  if(verbose) {verbose_int="1"}

  strvectors = utils::capture.output(utils::write.table(input_table, sep="\t"))

  command = paste("argv0 -task r -train dummy -test dummy -dim 1,1,40 -iter", iter, "-method mcmc -init_stdev 0.1", "-verbose", verbose_int, sep=" ")

  args = unlist(strsplit(command, " "), use.names=FALSE)

  #print("Compling LibFM...")

  sourceCpp(libfm_cpp_path)

  ret = do_libfm(strvectors, args)

  table_ret = utils::read.table(text=ret)

  return(table_ret)
}

