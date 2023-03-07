
#############################################################################
# Computational Model for Foot-and-Mouth Disease (FMD) Virus classification #
# and serotype prediction                                                   #
#############################################################################
# Feature matrix using k-mer technque
#' @title Feature generation from sequence data
#'
#' @description Returns k-mer counts data.
#'
#' @param x Nucleotide sequence data of virus isolates in fasta format.
#' @param k Integral number (>0) indicating the size of mers to generate features.
#'
#' @return Data frame containing feature counts of k-mers, where rows are virus isolates & columns are features.
#'
#' @author Samarendra Das
#'
#' @import kmer
#'
#' @export

#Feature generation
Feature.dna <- function(x, k) {
  res <- LETTERS[-c(1, 3, 7, 20)]
  if(missing(k)) k =3
  if(is.null(k)) {k=3}
  else(k = k)
  data <- kcount(x, k = k, residues = res, encode = FALSE)
  data <- as.data.frame(data)
  return(data)
}

#######
# Feature matrix using k-mer technque
#' @title Normalization of k-mer count data
#'
#' @description Normalizes the k-mer counts data to remove bias due to varying sequence lengths.
#'
#' @param dat Dataframe of k-mer counts of the virus isolates usually obtained from kmer techniques.
#'
#' @return Data frame containing normalized feature values of k-mers, where rows are virus isolates & columns are features.
#'
#' @author Samarendra Das
#'
#'
#' @export

#Feature generation
feat.norm <- function(dat) {
  id <- rownames(dat)
  minMax <- function(x) (x - min(x)) / (max(x) - min(x))
  dta.norm <- as.data.frame(lapply(dat, minMax))
  row.names(dta.norm) <- id
  return(dta.norm)
}

############
# Relevant Feature Selection
#' @title Releveant k-mer feature selection.
#'
#' @description Returns the list of relevant k-mer features required for model training.
#'
#' @param x Nucleotide sequence data of virus isolates in fasta format.
#' @param serotypes Serotype information of the FMD virus isolates, whose sequences are given in x
#' @param n Integral number (>0) indicating the number of relevant features to be included in the model.
#'
#' @return Data frame containing feature counts of k-mers, where rows are virus isolates & columns are features.
#' \itemize{
#'   \item1 feat.sel is the list of selected features.
#'   \item2 x.sel is the k-mer count of the virus isolates over the selected features.
#'  }
#'
#' @author Samarendra Das
#'
#' @import FSelector
#' @importFrom FSelector gain.ratio
#'
#' @export

#Feature selection
Feature.sel <- function (x, serotypes, n) {
  if(missing(n)) n=8
  if(is.null(n)) n=8
  dat <- Feature.dna(x, k=3)
  id <- rownames(dat)
  dat <- as.data.frame(dat)
  serotypes <- as.factor(serotypes)
  weights.gr <- gain.ratio (serotypes~., dat, unit = "log2")
  idd.gr <- sort(weights.gr[,1], decreasing = T, index.return=T)$ix
  x.sel <- dat[, idd.gr[1:n]]
  feat.sel <- colnames(dat)[idd.gr][1:n]
  out <- list(feat.sel=feat.sel, x.sel=x.sel)
  return(out)
}

#### Model training
#' @title Model training for serotype prediction and classification.
#'
#' @description Returns the trained model along with results on test data.
#'
#' @param x Nucleotide sequence data of virus isolates in fasta format.
#' @param serotyp Serotype information of the virus isolates whose sequence information was supplied in fasta format.
#' @param ind vector of two-length whose elements reprsent proportion of data used for model training and testing.
#' @param m Integral number (>0) indicating the number of runs performed on test data to evaluate the trained model.
#' @param k Integral number (>0) indicating the size of mers to generate features.
#' @param n number of features to be considered in model building.
#'
#'
#' @return Returns the prediction model along with model fitting measures including mean and standard error in accuracy and prediction error over the runs on test data.
#'
#' @author Samarendra Das
#'
#' @import randomForest
#' @importFrom randomForest randomForest
#' @importFrom stats predict
#' @importFrom stats sd
#' @import stats
#'
#' @export

SerotypPredModel <- function(x, serotyp, ind, m, k, n) {

  if(missing(k)) k =3
  if(is.null(k)) {k=3}
  else(k = k)

  if(missing(n)) n =8
  if(is.null(n)) {n=8}
  else(n = n)

  if(missing(m)) m = 100
  if(is.null(m)) {m = 100}
  else(m = m)

  if(missing(ind)) ind = c(0.7, 0.3)
  if(is.null(ind)) ind = c(0.7, 0.3)
  else(n = n)

  dat <- Feature.dna (x = x, k = k)
  dat.sel <- Feature.sel (x = x, serotypes = serotyp, n = n)$x.sel
  dat.norm <- feat.norm(dat.sel)
  serotyp <- as.factor(serotyp)

  result.mod <- vector(mode="numeric", length = m)
  #Model <- NULL
  for (i in 1:m){
    ix.train <- sample(1:nrow(dat.norm), round(ind[1]*nrow(dat.norm)))
    serotyp.train <- serotyp[ix.train]
    dta.norm.train <- dat.norm[ix.train, ]
    dta.norm.test <- dat.norm[-ix.train, ]
    serotyp.test <- serotyp[-ix.train]
    model <- randomForest(serotyp.train ~ ., data=dta.norm.train, importance = TRUE,
                          proximity = TRUE, ntree = 1000, mtry = 1)
    pred <- predict(model, dta.norm.test)
    confusionmtx.mod <- table(pred, serotyp.test)
    #confusionmtx.svm
    Acc <- sum(diag(confusionmtx.mod)/sum(confusionmtx.mod))
    result.mod[i] <- Acc
    #Model <-
  }
  res.mean.mod <- mean(result.mod)
  res.Se.mod <- sd(result.mod) / sqrt(m)
  res.err.mod <- mean(1 - result.mod)
  out.test <- list(res.mean.mod = res.mean.mod, res.Se.mod = res.Se.mod, res.err.mod = res.err.mod, model = model)
  return(out.test)
}

######## serotype prediction using trained Model
#' @title Prediction of unknown serotypes of isolates.
#'
#' @description Returns the serotypes of the novel FMD virus isolates supplied in query sequence data.
#'
#' @param x Nucleotide sequence data of virus isolates in fasta format for model training.
#' @param serotyp Serotype information of the virus isolates whose sequence data is given in x.
#' @param QuerySeq Sequence data of the virus isolates in fasta format, whose serotype will be predicted.
#'
#' @return Data frame containing probabilities for the serotypes computed for each isolate and right-most column contains predicted serotypes.
#'
#' @author Samarendra Das
#'
#' @importFrom stats predict
#' @import stats
#'
#' @export

serotypPred <- function(x, serotyp, QuerySeq){
  model1 <- SerotypPredModel (x, serotyp, ind = c(1.0, 0), m = 100, k = 3, n = 8)
  model2 <- model1$model
  feat.sel1 <- Feature.sel (x, serotyp, n = 8)$feat.sel
  dat.q <- Feature.dna(QuerySeq, k=3)
  dat.q.sel <- dat.q[, match(feat.sel1, colnames(dat.q))]
  dat.q.norm <- feat.norm(dat.q.sel)
  remove(model1, feat.sel1, dat.q, dat.q.sel); gc()
  pred.prob <- predict(model2, dat.q.norm, type = "prob")
  pred.serotypes <- predict(model2, dat.q.norm, type = "response")
  out <- data.frame(pred.prob, pred.serotypes)
  colnames(out) = c(paste(colnames(pred.prob), ".Prob", ""), "Serotypes")
  return(out)
}

# Testing of trained model
#' @title Testing of trained model.
#'
#' @description Validates the trained model on independent test data and return performance metrics.
#'
#' @param x Nucleotide sequence data of virus isolates in fasta format for model training.
#' @param serotyp vector (size equal to number of virus isolates in x) containing serotype information of virus in x.
#' @param QuerySeq Sequnce data of of virus isolates in fasta format used for testing the trained model.
#' @param serotypTest vector (size equal to QuerySeq length) containing the serotype information of the isolates of QuerySeq.
#'
#' @return List containing predicted serotypes of isolates, Confusion matrix for trained model testing, accuracy of the trained model on test data.
#' \itemize{
#'   \item1 PredictedSerotypes conatins the serotype of isolates in TestSeq data using the trained model.
#'   \item2 ConfusionMatrix confusion matrix constructed by predicted serotypes vs. actual serotypes of the isolates in TestSea data.
#'   \item3 Accuracy is the trained model accuracy measure computed from the test
#'  }
#' @author Samarendra Das
#'
#' @import stats
#'
#' @export

serotypPredTest <- function(x, serotyp, QuerySeq, serotypTest) {
  reslt <- serotypPred (x, serotyp, QuerySeq)
  pred.serotyp <- reslt$Serotypes
  remove(reslt); gc()
  serotypTest <- as.factor(serotypTest)
  confusionmtx.test <- table(serotypTest, pred.serotyp)
  names(dimnames(confusionmtx.test)) <- c("Actual", "Predicted")
  Acc <- round(sum(diag(confusionmtx.test)/sum(confusionmtx.test)) * 100, 3)
  out1 <- list(PredictedSerotypes = pred.serotyp, ConfusionMatrix = confusionmtx.test, Accuracy = Acc)
  remove(confusionmtx.test, Acc, pred.serotyp)
  return(out1)
}
