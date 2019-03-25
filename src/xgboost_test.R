require(xgboost)

X <- matrix(rnorm(1000*200),1000,200)
w <- matrix(rnorm(200*1),200,1)

Y <- X%*%w+matrix(rnorm(1000*1),1000,1)


mseregobj<- function(preds, dtrain) {
  labels <- getinfo(dtrain, "label")
  grad <- sum(preds - labels)
  hess <- 1
  return(list(grad = grad, hess = hess))
}

logregobj <- function(preds, dtrain) {
  labels <- getinfo(dtrain, "label")
  preds <- 1/(1 + exp(-preds))
  grad <- preds - labels
  hess <- preds * (1 - preds)
  return(list(grad = grad, hess = hess))
}

evalerror <- function(preds, dtrain) {
  labels <- getinfo(dtrain, "label")
  err <- as.numeric(sum(labels != (preds > 0)))/length(labels)
  return(list(metric = "error", value = err))
}

bst <- xgboost(data = X,label = Y, max.depth = 2, eta = 1, nthread = 2, nrounds = 2, objective = "reg:linear")
xgb.dump(bst,with_stats = T)
xgb.plot.tree(model = bst,trees = 0)

library("RGBM")
library(doParallel)
cl <- makeCluster(4)
aaa=matrix(rnorm(100*30),100,30)
bbb=matrix(0,100,30)
ccc=matrix(1,10,20)
colnames(aaa) <- c(paste0('tf',1:10),paste0('target',1:20))
colnames(bbb) <- colnames(aaa)
colnames(ccc) <- paste0('target',1:20)
rownames(ccc) <- paste0('tf',1:10)


# run network inference on a 100-by-100 dummy expression data.
A = RGBM(E = aaa,K = bbb,g_M = ccc,tfs = paste0('tf',1:10),targets = paste0('target',1:20))
stopCluster(cl)
