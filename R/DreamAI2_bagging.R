
#' @noRd
avg_batch = function(data,SamplesPerBatch)
{
  t(apply(data,1,function(x){apply(matrix(x,SamplesPerBatch),2,mean,na.rm=T)}));
}

#' @noRd
diff_mr = function(x,gm.pnnl,data.obs.b,N){
  p.temp = exp(x-gm.pnnl[2]*data.obs.b);
  p.temp[is.na(p.temp)] = 1;
  return(abs(mean(p.temp)-(sum(is.na(data.obs.b))+N)/sum(is.na(data.obs.b))*mean(is.na(data.obs.b))));
}



#' Bag Imputation of Missing Protein Abundances with Iterative Prediction Model
#' @description The function DreamAI_bagging imputes a dataset with missing values or NA's by bag imputaion with help of parallel processing. Pseudo datasets are generated having true missing (as in the original dataset) and pseudo missing and every such pseudo dataset is imputed by 7 different methods: KNN, MissForest, ADMIN, Birnn, SpectroFM, RegImpute and Ensemble (descriptions are included in the documentation of the function DreamAI).
#' @details This function can be run as parallel job in cluster. It generates and saves a .RData file containing the output from the current process in the location provided by the user, with the process number in the file name. If the user runs it in local computer multiple times, then changing the ProcessNumber everytime will generate and save .RData file with the given ProcessNumber.
#'
#' @param data dataset in the form of a matrix or dataframe with missing values or NA's. The function throws an error message and stops if any row or column in the dataset is missing all values
#' @param k number of neighbors to be used in the imputation by KNN and ADMIN (default is 10)
#' @param maxiter_MF maximum number of iteration to be performed in the imputation by "MissForest" if the stopping criteria is not met beforehand
#' @param ntree number of trees to grow in each forest in "MissForest"
#' @param maxnodes maximum number of terminal nodes for trees in the forest in "MissForest", has to equal at least the number of columns in the given data
#' @param maxiter_ADMIN maximum number of iteration to be performed in the imputation by "ADMIN" if the stopping criteria is not met beforehand
#' @param tol convergence threshold for "ADMIN"
#' @param gamma_bagging parameter to control abundance-dependence for pseudo missing values generation in bagging sets. Set gamma_bagging = NA, to learn the abundance-dependence from the observed data matrix. Set gamma_bagging = 0 to generate abundance-independent pseudo missing values.
#' @param gamma_ADMIN parameter for ADMIN to control abundance dependent missing. Set gamma_ADMIN=0 for log ratio intensity data. For abundance data put gamma_ADMIN=NA, and it will be estimated accordingly
#' @param gamma parameter of the supergradients of popular nonconvex surrogate functions, e.g. SCAD and MCP of L0-norm for Birnn
#' @param CV a logical value indicating whether to fit the best gamma with cross validation for "Birnn". If CV=FALSE, default gamma=50 is used, while if CV=TRUE gamma is calculated using cross-validation.
#' @param fillmethod a string identifying the method to be used to initially filling the missing values using simple imputation for "RegImpute". That could be "row_mean" or "zeros", with "row_mean" being the default. It throws an warning if "row_median" is used.
#' @param maxiter_RegImpute maximum number of iterations to reach convergence in the imputation by "RegImpute"
#' @param conv_nrmse convergence threshold for "RegImpute"
#' @param iter_SpectroFM number of iterations for "SpectroFM"
#' @param m_mice Number of multiple imputations in MICE. The default is m=1.
#' @param method_mice Specifying the imputation method to be used for each column in MICE. The default is 'pmm'.
#' @param maxiter_mice A scalar giving the number of iterations in MICE. The default is 20.
#' @param method a vector of imputation methods: ("KNN", "MissForest", "ADMIN", "Birnn", "SpectroFM, "RegImpute") based on which "Ensemble" imputed matrix will be obtained.
#' @param out a vector of imputation methods for which the function will output the imputed matrices. Default is "Ensemble".
#' @param seed.bags random seed used for generating missing values in bagging sets
#' @param SamplesPerBatch number of samples per batch (batch size in the original data)
#' @param n.bag number of pseudo datasets to generate and impute in the current process
#' @param path location to save the output file from the curent process.  Path only needs to be specified when save.out=TRUE
#' @param ProcessNum process number starting from 1 when run in cluster, e.g. 1 - 10, 1 - 100 etc. Needs to be specified only if the output is saved
#' @param save.out logical indicator whether or not to save the output. When TRUE output is saved, when FALSE output is only returned
#'
#' @return list of imputed dataset(average over all pseudo imputed data matrices) by different methods as specified by the user, n.bag and a matrix containing gene name, sample name, true and imputed values of every pseudo missing combined from n.bag datasets
#' @export
#' @examples
#' \dontrun{
#' data<-data.DIA[1:100,1:50]
#' impute<- DreamAI_Bagging(data=data,k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,
#' maxiter_ADMIN=30,tol=10^(-2),gamma_bagging=NA,gamma_ADMIN=NA,
#' gamma=50,CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40,
#' method = c("KNN", "MissForest", "ADMIN", "Birnn", "SpectroFM", "RegImpute","MICE"),
#' out=c("Ensemble.Fast"),
#' SamplesPerBatch=1,n.bag=2,save.out=TRUE,path="path_of_bagging_results",ProcessNum=1)
#' impute$Ensemble
#' }
DreamAI2_Bagging<-function(data,k=10,maxiter_MF = 10, ntree = 100,maxnodes = NULL,maxiter_ADMIN=30,tol=10^(-2),
                           gamma_bagging = NA, gamma_ADMIN=NA,gamma=50,
                           CV=FALSE,fillmethod="row_mean",maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40,
                           m_mice = 1, method_mice = 'pmm', maxiter_mice = 20,
                           method=c("KNN","MissForest","ADMIN","Birnn","SpectroFM","RegImpute","MICE"),out=c("Ensemble.Fast"),
                           seed.bags = NULL,SamplesPerBatch,n.bag,save.out=TRUE,path=NULL,ProcessNum=1)
{
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

  if(sum(rowMeans(is.na(data))==1)>0){cat(paste('removing',sum(rowMeans(is.na(data))==1),'features with all missing values'))}
  data.obs <- data[rowMeans(is.na(data))!=1,];

  data.obs.b <- avg_batch(data.obs,SamplesPerBatch=SamplesPerBatch)
  if(is.na(gamma_bagging)){gm.pnnl <- gamma_est(data.obs.b)}else{gm.pnnl <- c(NA,gamma_bagging)}

  N <- sum(is.na(data.obs.b[!(apply(is.na(data.obs.b), 1, mean)==1),]))


  gm1.new <- stats::optimize(f = diff_mr,interval = c(-10,10),gm.pnnl=gm.pnnl,data.obs.b=data.obs.b,N=N)$minimum

  p.pnnl <- exp(gm1.new-gm.pnnl[2]*data.obs.b)
  p.pnnl[is.na(p.pnnl)] <- 1
  p.pnnl[p.pnnl>1] <- 1
  # results<-NULL
  data.v<-c(as.matrix(data.obs))

  # true.miss <- which(is.na(data.v))
  # bag.knn.sum<-0
  # bag.MissForest.sum<-0
  # bag.ADMIN.sum<-0
  # bag.BruinGo.sum<-0
  # bag.DMIS.sum<-0
  # bag.Jeremy.sum<-0
  # bag.ensemble.sum<-0

  summary.all<-list()

  if(length(seed.bags)!=n.bag)
    seed.bags <- 1:n.bag+1000000

  # out.temp = union(method,out)
  # out.temp is the matched output names from out to the input default methods (include all individual methods to form the ensemble methods)
  out.temp = intersect(out,methods.match)

  if("Ensemble.Fast" %in% out)
    out.temp = union(out.temp,setdiff(methods.match,"MissForest"))
  if("Ensemble" %in% out)
    out.temp = union(out.temp,methods.match)

  out.temp = union(out.temp,out)

  impute.sum = 0

  for (i in 1:n.bag)
  {
    # TimeStart<-proc.time()

    set.seed(seed.bags[i]);
    m.missing.b <- matrix(stats::rbinom(length(p.pnnl), 1, c(p.pnnl)),dim(p.pnnl)[1]);
    m.missing <- t(apply(m.missing.b,1,rep,each=SamplesPerBatch));
    # sum(m.missing.b-is.na(data.obs.b))
    #print(sum(m.missing!=t(is.na(apply(data.obs.b,1,rep,each=SamplesPerBatch)))));

    data.obs.new <- data.obs;
    data.obs.new<-as.matrix(data.obs.new)
    data.obs.new[m.missing!=t(is.na(apply(data.obs.b,1,rep,each=SamplesPerBatch)))] <- NA;

    DreamAI.bagging.temp<-DreamAI2(data=data.obs.new,k=k,maxiter_MF = maxiter_MF, ntree = ntree,maxnodes = maxnodes,maxiter_ADMIN=maxiter_ADMIN,tol=tol,gamma_ADMIN=gamma_ADMIN,gamma=gamma,
                                    CV=CV,fillmethod=fillmethod,maxiter_RegImpute=maxiter_RegImpute,conv_nrmse = conv_nrmse,iter_SpectroFM=iter_SpectroFM,
                                    m_mice = m_mice, method_mice = method_mice, maxiter_mice = maxiter_mice,
                                    method=method,out=out.temp)

    ########## summary ###############
    ind.pseudo <- c(is.na(data.obs.new)-is.na(data.obs))==1;

    d.s <- data.frame(bag = i,
                      gene = rep(rownames(data.obs),ncol(data.obs))[ind.pseudo],
                      sample = rep(colnames(data.obs),each = nrow(data.obs))[ind.pseudo],
                      true = (as.matrix(data.obs))[ind.pseudo],
                      sapply(DreamAI.bagging.temp,function(x){x[ind.pseudo]}),
                      check.rows = F,check.names = F);

    summary.all[[i]]<-d.s

    impute.i = lapply(DreamAI.bagging.temp,function(x){x[is.na(data.obs)]})[out]

    impute.sum = impute.sum + do.call(cbind,impute.i)

    # methods<-c("KNN","MissForest","ADMIN","Birnn","SpectroFM","RegImpute")
    #
    # if("KNN" %in% out){
    #   bag.knn.sum <- bag.knn.sum+(c(ResultDreamImputation$KNN))[true.miss]
    # }
    #
    # if("MissForest" %in% out){
    #   bag.MissForest.sum <- bag.MissForest.sum+(c(ResultDreamImputation$MissForest))[true.miss]
    # }
    #
    # if("ADMIN" %in% out){
    #   bag.ADMIN.sum <- bag.ADMIN.sum+(c(ResultDreamImputation$ADMIN))[true.miss]
    # }
    #
    # if("Birnn" %in% out){
    #   bag.BruinGo.sum <- bag.BruinGo.sum+(c(ResultDreamImputation$Birnn))[true.miss]
    # }
    #
    # if("SpectroFM" %in% out){
    #   bag.DMIS.sum <- bag.DMIS.sum+(c(ResultDreamImputation$SpectroFM))[true.miss]
    # }
    #
    # if("RegImpute" %in% out){
    #   bag.Jeremy.sum <- bag.Jeremy.sum+(c(ResultDreamImputation$RegImpute))[true.miss]
    # }
    # if("Ensemble" %in% out){
    #   bag.ensemble.sum <- bag.ensemble.sum+(c(ResultDreamImputation$Ensemble))[true.miss]
    # }
    cat("\n\n")
    print(paste("imputation on bagging data ",i," completed",sep=""))
    cat("\n\n")
  }

  impute.bag.mean = impute.sum/n.bag

  summary.all.df <- do.call(rbind,summary.all);

  imputed_matrix<-list()

  for(method.out.i in out)
  {
    imputed_matrix[[method.out.i]] = data.obs
    imputed_matrix[[method.out.i]][is.na(data.obs)] = impute.bag.mean[,method.out.i]
  }

  # if("KNN" %in% out){
  #   bag.knn.v<- bag.knn.sum/n.bag
  #   data.temp <- data.v
  #   data.temp[true.miss]<-bag.knn.v
  #   bag.knn<-matrix(data.temp,nrow(data.obs),dimnames = dimnames(data.obs))
  #   imputed_matrix<-c(imputed_matrix,list("KNN"=bag.knn))
  # }
  #
  # if("MissForest" %in% out){
  #   bag.MissForest.v<- bag.MissForest.sum/n.bag
  #   data.temp <- data.v
  #   data.temp[true.miss]<-bag.MissForest.v
  #   bag.MissForest<-matrix(data.temp,nrow(data.obs),dimnames = dimnames(data.obs))
  #   imputed_matrix<-c(imputed_matrix,list("MissForest"=bag.MissForest))
  # }
  #
  # if("ADMIN" %in% out){
  #   bag.ADMIN.v<- bag.ADMIN.sum/n.bag
  #   data.temp <- data.v
  #   data.temp[true.miss]<-bag.ADMIN.v
  #   bag.ADMIN<-matrix(data.temp,nrow(data.obs),dimnames = dimnames(data.obs))
  #   imputed_matrix<-c(imputed_matrix,list("ADMIN"=bag.ADMIN))
  # }
  #
  # if("Birnn" %in% out){
  #   bag.BruinGo.v<- bag.BruinGo.sum/n.bag
  #   data.temp <- data.v
  #   data.temp[true.miss]<-bag.BruinGo.v
  #   bag.BruinGo<-matrix(data.temp,nrow(data.obs),dimnames = dimnames(data.obs))
  #   imputed_matrix<-c(imputed_matrix,list("Birnn"=bag.BruinGo))
  # }
  #
  # if("SpectroFM" %in% out){
  #   bag.DMIS.v<- bag.DMIS.sum/n.bag
  #   data.temp <- data.v
  #   data.temp[true.miss]<-bag.DMIS.v
  #   bag.DMIS<-matrix(data.temp,nrow(data.obs),dimnames = dimnames(data.obs))
  #   imputed_matrix<-c(imputed_matrix,list("SpectroFM"=bag.DMIS))
  # }
  #
  # if("RegImpute" %in% out){
  #   bag.Jeremy.v<- bag.Jeremy.sum/n.bag
  #   data.temp <- data.v
  #   data.temp[true.miss]<-bag.Jeremy.v
  #   bag.Jeremy<-matrix(data.temp,nrow(data.obs),dimnames = dimnames(data.obs))
  #   imputed_matrix<-c(imputed_matrix,list("RegImpute"=bag.Jeremy))
  # }
  # if("Ensemble" %in% out){
  #   bag.ensemble.v<- bag.ensemble.sum/n.bag
  #   data.temp <- data.v
  #   data.temp[true.miss]<-bag.ensemble.v
  #   bag.ensemble<-matrix(data.temp,nrow(data.obs),dimnames = dimnames(data.obs))
  #   imputed_matrix<-c(imputed_matrix,list("Ensemble"=bag.ensemble))
  # }

  bag.output<-list(impute=imputed_matrix,n.bag=n.bag,summary=summary.all.df,out.method=out)

  # num<-which(methods %in% method)
  # bag.imputed_matrix[[1]]<-bag.imputed_matrix[[1]][c(num,7)]
  # bag.output<-bag.imputed_matrix

  if(save.out){
    saveRDS(bag.output,file=paste(path,"bag_imputed_",sprintf("%03d",ProcessNum),".rds",sep=""))
  }else{
    return(bag.output)
  }
}


