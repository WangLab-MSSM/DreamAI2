#' Imputation of Missing Protein Abundances with Iterative Prediction Model
#' @description The function DreamAI2 imputes a dataset with missing values or NA's using individual or ensemble methods.
#' Individual methods include: "KNN": k nearest neighbor, "MissForest": nonparametric Missing Value Imputation using Random Forest, "ADMIN": abundance dependent missing imputation, "Birnn": imputation using IRNN-SCAD algorithm, "SpectroFM": imputation using matrix factorization, "RegImpute": imputation using Glmnet ridge regression, and "MICE": Multivariate Imputation by Chained Equations.
#' Ensemble methods include: "Ensemble": average of the selected methods among the 7, "Ensemble.Fast": average of the selected methods without "MissForest".
#'
#' @param data dataset in the form of a matrix or dataframe with missing values or NA's. The function throws an error message and stops if any row or column in the dataset is missing all values.
#' @param k number of neighbors to be used in the imputation by KNN and ADMIN (default is 10).
#' @param maxiter_MF maximum number of iteration to be performed in the imputation by "MissForest" if the stopping criteria is not met beforehand.
#' @param ntree number of trees to grow in each forest in "MissForest".
#' @param maxnodes maximum number of terminal nodes for trees in the forest in "MissForest", has to equal at least the number of columns in the given data.
#' @param maxiter_ADMIN maximum number of iteration to be performed in the imputation by "ADMIN" if the stopping criteria is not met beforehand.
#' @param tol convergence threshold for "ADMIN".
#' @param gamma_ADMIN parameter for ADMIN to control abundance dependent missing. Set gamma_ADMIN=0 for log ratio intensity data. For abundance data put gamma_ADMIN=NA, and it will be estimated accordingly.
#' @param gamma parameter of the supergradients of popular nonconvex surrogate functions, e.g. SCAD and MCP of L0-norm for Birnn.
#' @param CV a logical value indicating whether to fit the best gamma with cross validation for "Birnn". If CV=FALSE, default gamma=50 is used, while if CV=TRUE gamma is calculated using cross-validation.
#' @param fillmethod a string identifying the method to be used to initially filling the missing values using simple imputation for "RegImpute". That could be "row_mean" or "zeros", with "row_mean" being the default. It throws an warning if "row_median" is used.
#' @param maxiter_RegImpute maximum number of iterations to reach convergence in the imputation by "RegImpute".
#' @param conv_nrmse convergence threshold for "RegImpute".
#' @param n_train number of predictor used for training for "RegImpute" , default is 50
#' @param nfolds number of folds in cross validation of glmnet fitting for "RegImpute", default is 10.
#' @param nlambda The number of lambda values for "RegImpute", default is 100.
#' @param iter_SpectroFM number of iterations for "SpectroFM".
#' @param m_mice Number of multiple imputations in MICE. The default is m=1.
#' @param method_mice Specifying the imputation method to be used for each column in MICE. The default is 'pmm'.
#' @param maxiter_mice A scalar giving the number of iterations in MICE. The default is 20.
#' @param seed_mice random seed used in MICE.
#' @param method a vector of imputation methods: ("KNN", "MissForest", "ADMIN", "Birnn", "SpectroFM, "RegImpute") based on which "Ensemble" imputed matrix will be obtained.
#' @param out a vector of imputation methods for which the function will output the imputed matrices. Default is "Ensemble".
#'
#'@note If all methods are specified for obtaining "Ensemble" imputed matrix, the approximate time required to output the imputed matrix for a dataset of dimension 26000 x 200 is ~50 hours.
#'
#' @return a list of imputed datasets by different methods as specified by the user. Always returns imputed data by "Ensemble"
#' @export
#' @examples
#' \dontrun{
#' data<-data.DIA[1:100,1:50]
#' impute<- DreamAI(data,k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,
#' maxiter_ADMIN=30,tol=10^(-2),gamma_ADMIN=NA,
#' gamma=50,CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40,
#' method = c("KNN", "MissForest", "ADMIN", "Birnn", "SpectroFM", "RegImpute", "MICE"),
#' out="Ensemble.Fast")
#' impute$Ensemble
#' }
DreamAI2<-function(data,k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,maxiter_ADMIN=30,tol=10^(-2),gamma_ADMIN=NA,gamma=50,
                   CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,n_train = 50,nfolds = 10,nlambda = 100,iter_SpectroFM=40,
                   m_mice = 1, method_mice = 'pmm', maxiter_mice = 20,seed_mice = 123,
                   method=c("KNN","MissForest","ADMIN","Birnn","SpectroFM","RegImpute","MICE"),out=c("Ensemble.Fast"))
{
  TimeStart<-proc.time()
  pkg.all = c("cluster"
            ,"survival"
            ,"randomForest"
            ,"missForest"
            ,"glmnet"
            ,"Rcpp"
            ,"foreach"
            ,"itertools"
            ,"iterators"
            ,"Matrix"
            ,"devtools"
            ,"impute"
            ,"mice")

  pkg.req = pkg.all[!sapply(pkg.all,base::requireNamespace, quietly = T)]
  if(length(pkg.req)>0)
  {
    base::stop(paste0("\nSome packages is not yet installed:\n  ",
                    paste(pkg.req,collapse = ', '),
                    ".\nPlease install them before running DreamAI. \n"),
             call. = FALSE)
  }

  # requireNamespace("cluster")
  # requireNamespace("survival")
  # requireNamespace("randomForest")
  # requireNamespace("missForest")
  # requireNamespace("glmnet")
  # requireNamespace("Rcpp")
  # requireNamespace("foreach")
  # requireNamespace("itertools")
  # requireNamespace("iterators")
  # requireNamespace("Matrix")
  # requireNamespace("devtools")
  # requireNamespace("impute")
  # requireNamespace("mice")

  missing_rows = (which(rowSums(is.na(data))==dim(data)[2]))
  if(length(missing_rows)>0){
    stop("Some Rows Are Missing All Values !!\n")
  }
  missing_cols = (which(colSums(is.na(data))==dim(data)[1]))
  if(length(missing_cols)>0){
    stop("Some Columns Are Missing All Values !!\n")
  }
  options(warn=-1)
  ### subsetting data without rows/cols with all missing ###

  # data = subset(data, (rowSums(is.na(data))!=dim(data)[2]))
  # data = subset(data, (colSums(is.na(data))!=dim(data)[1]))

  ### method of imputation ###
  method.all<-c("KNN","MissForest","ADMIN","Birnn","SpectroFM","RegImpute","MICE")

  if(length(setdiff(method,method.all))>0)
  {
    stop(paste0('\nNot identifiable methods selection: ',
              paste(setdiff(method,method.all),collapse = ', '),
              '.\nPlease select methods from the 7 components: ',
              paste(method.all,collapse = ', '),'.\n'))
  }

  # methods.match is the matched method names from input to the default 7 methods
  methods.match<- method.all[which(method.all %in% method)]

  n.method = length(methods.match)

  if(n.method==0)
  {
    sink()
    return(print("specify method"))
  }

  if(mean(out%in%c(methods.match,'Ensemble','Ensemble.Fast'))!=1)
  {
    stop(paste0('\nNot identifiable output method name: ',
              paste(setdiff(out,c(methods.match,'Ensemble','Ensemble.Fast')),collapse = ', '),
              '.\nPlease select among from the specified method: ',
              paste(methods.match,collapse = ', '),'.\n'))
  }

  message(paste('\n',n.method,'methods specified, ensemble imputation will be generated with those algorithms:\n',
              paste0(methods.match,collapse = ', '),'\n'))

  # out.match is the matched output names from out to the input default methods (include all individual methods to form the ensemble methods)
  out.match = intersect(out,methods.match)

  if("Ensemble.Fast" %in% out)
    out.match = union(out.match,setdiff(methods.match,"MissForest"))
  if("Ensemble" %in% out)
    out.match = union(out.match,methods.match)

  ensemble<-matrix(0,nrow(data),ncol(data),dimnames = dimnames(data))
  method.idx<-1
  imputed_matrix=list()
  out_matrix=list()

  ## KNN ##

  if("KNN" %in% out.match)
  {
    name = "KNN"
    sink("NULL")
    d.impute.KNN = impute.KNN(data = as.matrix(data),k = k)
    # ensemble<-ensemble+d.impute.knn
    sink()
    # print(paste("Method",method.idx,"complete"))
    # method.idx<-method.idx+1
    # imputed_matrix<-c(imputed_matrix,list("KNN"=as.matrix(d.impute.knn)))

    print(paste("Method",name,"complete"))
    imputed_matrix[[name]] = as.matrix(d.impute.KNN)
  }

  ## MF ##
  if("MissForest" %in% out.match)
  {
    name = "MissForest"
    if(is.null(maxnodes))
    {
      maxnodes<-ncol(as.matrix(data))
    }
    sink("NULL")
    d.impute.MF = impute.MF(data = as.matrix(data),maxiter_MF=maxiter_MF,ntree=ntree,maxnodes=maxnodes)
    # ensemble<-ensemble+d.impute.MF
    sink()
    # print(paste("Method",method.idx,"complete"))
    # method.idx<-method.idx+1
    # imputed_matrix<-c(imputed_matrix,list("MissForest"=as.matrix(d.impute.MF)))

    print(paste("Method",name,"complete"))
    imputed_matrix[[name]] = as.matrix(d.impute.MF)
  }

  ## ADMIN ##
  if("ADMIN" %in% out.match)
  {
    name = "ADMIN"
    sink("NULL")
    d.impute.ADMIN = impute.ADMIN(data = as.matrix(data),k = k, gamma=gamma_ADMIN,maxiter_ADMIN = maxiter_ADMIN, tol = tol)
    # ensemble<-ensemble+d.impute.ADMIN
    sink()
    # print(paste("Method",method.idx,"complete"))
    # method.idx<-method.idx+1
    # imputed_matrix<-c(imputed_matrix,list("ADMIN"=as.matrix(d.impute.ADMIN)))

    print(paste("Method",name,"complete"))
    imputed_matrix[[name]] = as.matrix(d.impute.ADMIN)
  }


  ## Birnn ##
  if("Birnn" %in% out.match)
  {
    name = "Birnn"
    sink("NULL")
    d.impute.Birnn=impute.Birnn(as.matrix(data), gamma = 50, CV = CV)
    # ensemble<-ensemble+d.impute.Birnn
    sink()
    # print(paste("Method",method.idx,"complete"))
    # method.idx<-method.idx+1
    # imputed_matrix<-c(imputed_matrix,list("Birnn"=as.matrix(d.impute.Birnn)))

    print(paste("Method",name,"complete"))
    imputed_matrix[[name]] = as.matrix(d.impute.Birnn)
  }

  ## SpectroFM ##
  if("SpectroFM" %in% out.match)
  {
    name = "SpectroFM"
    sink("NULL")
    d.impute.SpectroFM=impute.SpectroFM(input_table=as.data.frame(data),iter=iter_SpectroFM, verbose=FALSE)
    # ensemble<-ensemble+d.impute.SpectroFM
    sink()
    # print(paste("Method",method.idx,"complete"))
    # method.idx<-method.idx+1
    # imputed_matrix<-c(imputed_matrix,list("SpectroFM"=as.matrix(d.impute.SpectroFM)))

    print(paste("Method",name,"complete"))
    imputed_matrix[[name]] = as.matrix(d.impute.SpectroFM)
  }


  ## RegImpute ##
  if("RegImpute" %in% out.match)
  {
    name = "RegImpute"
    sink("NULL")
    d.impute.RegImpute=impute.RegImpute(data=as.matrix(data),fillmethod=fillmethod,
                                        maxiter_RegImpute = maxiter_RegImpute,conv_nrmse = conv_nrmse,
                                        n_train = n_train,nfolds = nfolds,nlambda = nlambda)
    # ensemble<-ensemble+d.impute.RegImpute
    sink()
    # print(paste("Method",method.idx,"complete"))
    # method.idx<-method.idx+1
    # imputed_matrix<-c(imputed_matrix,list("RegImpute"=as.matrix(d.impute.RegImpute)))

    print(paste("Method",name,"complete"))
    imputed_matrix[[name]] = as.matrix(d.impute.RegImpute)
  }

  ## MICE ##
  if("MICE" %in% out.match)
  {
    if ("Matrix" %in% .packages()) {
      detach("package:Matrix", unload=TRUE)
    }
    name = "MICE"
    sink("NULL")
    d.impute.MICE = impute.mice(data = as.matrix(data),m = m_mice ,method = method_mice ,maxit = maxiter_mice, seed = seed_mice)
    # ensemble<-ensemble+d.impute.knn
    sink()
    # print(paste("Method",method.idx,"complete"))
    # method.idx<-method.idx+1
    # imputed_matrix<-c(imputed_matrix,list("MICE"=as.matrix(d.impute.knn)))

    print(paste("Method",name,"complete"))
    imputed_matrix[[name]] = as.matrix(d.impute.MICE)
  }

  ## Ensemble ##

  df.impute.Ensemble = Reduce('+',imputed_matrix)/length(imputed_matrix)
  df.impute.Ensemble.Fast = Reduce('+',imputed_matrix[setdiff(names(imputed_matrix),"MissForest")])/length(setdiff(names(imputed_matrix),"MissForest"))

  imputed_matrix[['Ensemble']] = df.impute.Ensemble
  imputed_matrix[['Ensemble.Fast']] = df.impute.Ensemble.Fast

  out_matrix = imputed_matrix[out]

  sink()
  # print(paste("Method",method.idx,"complete"))


  # sink("NULL")
  # ensemble<- (d.impute.knn+d.impute.MF+d.impute.ADMIN+d.impute.Birnn+d.impute.SpectroFM+d.impute.RegImpute)/n.method
  # sink()
  #
  # imputed_matrix=list("KNN"=as.matrix(d.impute.knn),"MissForest"=as.matrix(d.impute.MF),"ADMIN"=as.matrix(d.impute.ADMIN),"Birnn"=as.matrix(d.impute.Birnn),"SpectroFM"=as.matrix(d.impute.SpectroFM),"RegImpute"=as.matrix(d.impute.RegImpute),"Ensemble"=as.matrix(ensemble))
  #
  # num<-which(method %in% out)
  #
  # output<-out_matrix

  TimeTaken<-round((proc.time()-TimeStart)[3]/60,3)
  cat("\n\n")
  print(paste("Imputation complete - time taken ",TimeTaken," min",sep=""))
  cat("\n\n")

  return(out_matrix)
}



