#' @importFrom mice mice
#' @importFrom mice mice.impute.pmm
#' @importFrom mice complete
#' @importFrom impute impute.knn
#' @importFrom missForest missForest


#' @noRd
my.normalize = function(data)
{
  median.temp = apply(data,1,stats::median,na.rm=T);
  sd.temp = apply(data,1,stats::sd,na.rm=T);
  # apply(t(t(data)-median.temp),2,median,na.rm=T);
  # apply(t((t(data)-median.temp)/sd.temp),2,sd,na.rm=T);
  data.new = ((data-median.temp)/sd.temp);
  return(list(data = data.new, median.list = median.temp, sd.list=sd.temp));
}


#' @noRd
my.normalize.rev = function(data, median.old, sd.old)
{
  data.new = (data*sd.old+median.old);
  return(data.new)
}


#' Estimate abundance dependent missing parameter
#' @description The function estimate log linear dependency between missing rate and protein abundance.
#' @details Input of this function is the data matrix with rows as proteins and columns as samples. Out put is the parameter of abundance missing dependence.
#'
#' @param data dataset in the form of a matrix or dataframe with missing values or NA's.
#'
#' @return intercept and slope of the log linear trend between missing rate and protein mean abundance.
#' @export
#' @examples
#' \dontrun{
#' data<-data.DIA[1:100,1:50]
#' gm_temp = gamma_est(data)
#' }
gamma_est = function(data)
{
  ### data dim Lxn
  Y.mean = apply(data,1,mean,na.rm=T);
  q = apply(is.na(data),1,mean);
  par = stats::lsfit(Y.mean[q>0],log(1/q)[q>0])[[1]];
  names(par) = c('gamma0','gamma');
  return(par);
}


### 2 knn linear regression
#' @noRd
# knn.est.it2 = function(data, k, m.ind=T)
# {
#   L = dim(data)[1];
#   ccc = stats::cor(t(data));
#   order.k = apply(ccc,1,function(x){order(x,decreasing = T)[1:k+1]});
#
#   data.est = data;
#   data.est[m.ind,] = t(sapply((1:L)[m.ind],function(l){stats::predict(stats::lm(data[l,]~t(data[order.k[,l],])))}));
#   dimnames(data.est) = dimnames(data);
#   data.est[data.est<0] = 0;
#   return(data.est);
# }
knn.est.it2 = function(data, k, m.ind=T)
{
  L = dim(data)[1];
  ccc = stats::cor(t(data));
  # order.k = apply(ccc,1,function(x){order(x,decreasing = T)[1:k+1]});
  order.k = apply(ccc,1,order,decreasing = T)[1:k+1,]

  data.est = data;
  for(l in (1:L)[m.ind])
  {
    data.est[l,] = stats::predict(stats::lm(data[l,]~t(data[order.k[,l],])))
  }
  dimnames(data.est) = dimnames(data);
  data.est[data.est<0] = 0;
  return(data.est);
}

#########################################

#' @title Imputation of Missing Protein Abundances using K Nearest Neighbour
#' @description The function impute.KNN imputes a dataset with missing values or NA's using k nearest neighbour
#' @param data dataset in the form of a matrix or data frame with NAs as missings
#' @param k number of neighbors to be used in the imputation (default=10)
#'
#' @return the imputed version of the dataset
#' @export
#'
#' @examples
#' \dontrun{
#' data<-data.DIA[1:100,1:50]
#' impute.KNN(data=as.matrix(data), k=10)
#' }
impute.KNN = function(data,k)
{
  requireNamespace("impute")

  norm.temp = my.normalize(data);

  data.new = norm.temp$data;
  median.new = norm.temp$median.list;
  sd.new = norm.temp$sd.list;

  impu.temp = impute.knn(data.new,k,rowmax = 0.9,colmax = 0.9,maxp=dim(data)[1]);

  X.new = my.normalize.rev((impu.temp[[1]]), median.new, sd.new);
  # X.new[X.new<0] = 0;
  return((X.new));

}



#' @title Imputation of Missing Protein Abundances using MissForest
#' @description The function impute.MF imputes a dataset with missing values or NA's using MissForest
#' @param data dataset in the form of a matrix or data frame with NAs as missing values
#' @param maxiter_MF maximum number of iteration to be performed if the stopping criteria is not met beforehand
#' @param ntree number of trees to grow in each forest
#' @param maxnodes maximum number of terminal nodes for trees in the forest, has to equal at least the number of columns in the given data
#' @return the imputed version of the dataset
#' @export
#'
#' @examples
#' \dontrun{
#' data<-data.DIA[1:100,1:50]
#' impute.MF(data=as.matrix(data), maxiter_MF = 10, ntree = 100,maxnodes = NULL)
#' }
impute.MF = function(data = as.matrix(data),maxiter_MF, ntree, maxnodes)
{
  requireNamespace("missForest")

  if(is.null(maxnodes))
  {
    maxnodes<-ncol(as.matrix(data))
  }
  mF.temp = missForest(t(data), maxiter = maxiter_MF, ntree = ntree, variablewise = FALSE,
                       decreasing = FALSE, verbose = FALSE,
                       mtry = ceiling(sqrt(dim(data)[1])), replace = TRUE,
                       classwt = NULL, cutoff = NULL, strata = NULL,
                       sampsize = NULL, nodesize = NULL, maxnodes = maxnodes,
                       xtrue = NA, parallelize = c('no', 'variables', 'forests'));
  X.temp = t(mF.temp[[1]]);
  return(X.temp);
}




#' @title Imputation of Missing Protein Abundances using ADMIN (Abundance  Dependent  Missing Imputation Mechanism)
#'@description The function impute.ADMIN imputes a dataset with missing values or NA's using ADMIN
#' @param data dataset in the form of a matrix or data frame with NAs as missings
#' @param data.ini initial dataset, set to be NA as default
#' @param gamma parameter of the non-ignorable missing mechanism
#' @param k number of neighbors to be used in the imputation (default=10)
#' @param maxiter_ADMIN maximum number of iteration to be performed if the stopping criteria is not met beforehand
#' @param tol convergence threshold
#'
#' @return the imputed version of the dataset
#' @export
#'
#' @examples
#' \dontrun{
#' data<-data.DIA[1:100,1:50]
#' impute.ADMIN(data = as.matrix(data),k = 10, maxiter_ADMIN = 30, tol = 10^(-2))
#' }
# impute.ADMIN = function(data,data.ini=NA,gamma, k, maxiter_ADMIN,tol)
# {
#   L = dim(data)[1];
#   iter = 0;
#   diff = 999;
#   llh = Inf;
#   if(is.na(gamma))
#   {gamma = gamma_est(data)[2]}
#
#   M.t = data.ini;
#   if(sum(is.na(data.ini))>0){
#     set.seed(1234);
#     fit.0 = impute.knn(data,k=5,rowmax = 0.9,colmax = 0.9,maxp=dim(data)[1])[[1]];
#     M.t = fit.0;
#   }
#
#   d.t = 0;
#
#   while ((iter<maxiter_ADMIN) & (diff>tol))
#   {
#     gc();
#     iter = iter+1;
#
#     X.temp = (M.t-is.na(data)*gamma*d.t);
#
#     M.t = knn.est.it2(X.temp, k);
#
#     d.t = apply(M.t-X.temp,1,stats::var);
#
#     # X.temp = M.t;
#     M.diff = (M.t-data);
#     M.diff[is.na(M.diff)]=0;
#
#     llh = c(llh,-norm(M.diff,type = 'F')^2-sum((M.t*gamma*d.t)[is.na(data)]));
#
#     # diff = abs(diff(llh)[iter]/llh[iter+1]);
#     diff = abs(diff(llh)[iter]);
#
#     M.t[!is.na(data)] = data[!is.na(data)];
#
#     print(c(iter,llh[iter+1]));
#   }
#   dimnames(M.t) = dimnames(data);
#   # return(list(M.t,llh))
#   return(M.t);
# }
impute.ADMIN = function(data,data.ini=NA,gamma=NA, k=10, maxiter_ADMIN=30,tol=10e-3)
{
  requireNamespace("impute")

  L = dim(data)[1];
  iter = 0;
  diff = 999;
  llh = Inf;
  if(is.na(gamma))
  {gamma = gamma_est(data)[2]}

  M.t = data.ini;
  if(sum(is.na(data.ini))>0){
    set.seed(1234);
    M.t = impute.knn(data,k=5,rowmax = 0.9,colmax = 0.9,maxp=dim(data)[1])[[1]];
    # fit.0 = impute.knn(data,k=5,rowmax = 0.9,colmax = 0.9,maxp=dim(data)[1])[[1]];
    # M.t = fit.0;
  }

  d.t = 0;

  llh.new = 0

  while ((iter<maxiter_ADMIN) & (diff>tol))
  {
    iter = iter+1;

    X.temp = (M.t-is.na(data)*gamma*d.t);

    M.t = knn.est.it2(X.temp, k);

    d.t = apply(M.t-X.temp,1,stats::var);

    # X.temp = M.t;
    M.diff = (M.t-data);
    M.diff[is.na(M.diff)]=0;

    # llh = c(llh,-norm(M.diff,type = 'F')^2-sum((M.t*gamma*d.t)[is.na(data)]));
    #
    # # diff = abs(diff(llh)[iter]/llh[iter+1]);
    # diff = abs(diff(llh)[iter]);
    #
    # M.t[!is.na(data)] = data[!is.na(data)];
    #
    # print(c(iter,llh[iter+1]));

    llh.old = llh.new

    llh.new = -norm(M.diff,type = 'F')^2-sum((M.t*gamma*d.t)[is.na(data)])

    diff = abs(llh.new-llh.old)
    print(c(iter,llh.new))

    gc()
  }
  dimnames(M.t) = dimnames(data);
  # return(list(M.t,llh))
  return(M.t);
}





#' @title Imputation of Missing Protein Abundances using MICE
#' @description The function impute.mice imputes a dataset with missing values or NA's using mice.
#' @param data dataset in the form of a matrix or data frame with NAs as missings.
#' @param m Number of multiple imputations. The default is m=1.
#' @param method Specifying the imputation method to be used for each column in data. The default is 'pmm'.
#' @param maxit A scalar giving the number of iterations. The default is 20.
#' @param ... other parameters in mice() function
#' @return the imputed version of the dataset
#' @export
#'
#' @examples
#' \dontrun{
#' data<-data.DIA[1:100,1:50]
#' impute.mice(data=as.matrix(data))
#' }

impute.mice = function(data,m = 1,method = 'pmm',maxit = 20,...)
{
  requireNamespace("mice")
  requireNamespace("impute")
  # impute.temp = mice::mice(data, m = m, method = method, maxit = maxit,...)
  # impute.list <- mice::complete(impute.temp, 'all')
  # impute.out = Reduce('+',impute.list)/length(impute.list)
  # rownames(impute.out) = rownames(data)

  pred = 1-diag(1,ncol(data))
  colnames(pred) = colnames(data)
  rownames(pred) = colnames(data)

  if(ncol(data)>50)
  {
    q.cut = 0.9
    cor.s = cor(data,use = 'pairwise.complete.obs')
    cor.q.cut = as.numeric(quantile(cor.s[upper.tri(cor.s)],probs = q.cut))
    pred <- mice::quickpred(data, mincor = cor.q.cut)

    cor.top50 = as.matrix(t(apply(cor.s,1,function(x){rank(-abs(x))<52})))
    diag(cor.top50) = 0
    rownames(cor.top50) = colnames(data)
    colnames(cor.top50) = colnames(data)

    pred[rowSums(pred)<50,] = cor.top50[rowSums(pred)<50,]
  }


  impute.temp = mice::mice(data, m = m, method = method, maxit = maxit,predictorMatrix = pred,...)

  impute.list <- mice::complete(impute.temp, 'all')
  impute.out = Reduce('+',impute.list)/length(impute.list)
  rownames(impute.out) = rownames(data)
  colnames(impute.out) = colnames(data)

  if(sum(is.na(impute.out))>0)
  {impute.out = t(impute::impute.knn(t(impute.out),k = 10)$data)}

  return(impute.out);
}


