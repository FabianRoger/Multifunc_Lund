#'
#' @title Modified functions from the multifunc package
#' 
#' @description Script to hold the AIC-based turnover approach functions from
#' the multifunc package which we modified to have a stronger AIC penalty term. In
#' the standard implementation of the AIC-approach in the multifunc package, only 
#' the standard AIC is allowed but Meyer et al. (2018) extended this approach
#' to allow more stringent penalty terms. We borrow this approach from
#' Meyer et al. (2018)
#'
#'
#' @title getRedundancy
#'
#' @description
#' \code{getRedundancy} examines which species have an effect on which function
#'
#' @details getRedundancy takes a matrix of 1s,0s, and -1s, and depending on whether we're
#' interested in positive, negative, or both types of interactions looks for the
#' m-wise overlap between species and returns the overlap index for each combination. For
#' species whose effect is not different from 0 at the alpha=0.05 level, a 0 is returned.
#'
#' @author Jarrett Byrnes.
#' @param vars Vector of column names of functions
#' @param species Vector of column names of species
#' @param data data frame with species presence/absence of values of functions
#' @param negVars Vector of names of species for which a negative coefficient is actually a positive effect.
#' @param method Fitting function for statistical models.  Defaults to \code{lm}.
#' @param combine How are species combined in the model? Defaults to "+" for additive combinations.
#' @param output Will the output be sign of effect or "coefficient".  Defaults to "effect"
#' @param ... Other arguments to be supplied to fitting function.
#'
#'
#' @export
#' @return Returns a matrix of functions and the effect of species on each. 1s, -1s, and 0s for "effect" or coefficients.
#'
#'
#' @examples
#' data(all_biodepth)
#' allVars <- qw(biomassY3, root3, N.g.m2, light3, N.Soil, wood3, cotton3)
#'
#' germany <- subset(all_biodepth, all_biodepth$location == "Germany")
#'
#' vars <- whichVars(germany, allVars)
#' species <- relevantSp(germany, 26:ncol(germany))
#'
#' # re-normalize N.Soil so that everything is on the same
#' # sign-scale (e.g. the maximum level of a function is
#' # the "best" function)
#' germany$N.Soil <- -1 * germany$N.Soil + max(germany$N.Soil, na.rm = TRUE)
#'
#' res.list <- lapply(vars, function(x) sAICfun(x, species, germany))
#' names(res.list) <- vars
#'
#' getRedundancy(vars, species, germany)
#' getRedundancy(vars, species, germany, output = "coef")
#'
#'
#'
#' #########
#' # takes a vector of responses, the species that may cause them
#' # and returns a table of 1s, -1s, and 0s with regards to the kind of effect
#' # or a coefficient table, if asked for.  Arugments can take the form of the fitting function
#' # how variables are combined, and additional arguments to the fitting function
#' #########
getRedundancy2 <- function (vars, species, data, negVars = NA, method = "lm", combine = "+", k, output = "effect", ...) {
  res.list <- lapply(vars, function(x) {
    positive.desired <- T
    if (x %in% negVars) 
      positive.desired <- F
    sAICfun2(response = x, species = species, data = data, 
             positive.desired, method = method, combine = combine, k=k,
             ...)
  })
  # what if they want the coefficients
  if (output == "coef") {
    # ret<-plyr::ldply(res.list, function(x) x$coefs)
    ret <- lapply(res.list, function(x) x$coefs)
    ret <- do.call(rbind, ret) %>% as.data.frame()
  } else {
    # the default of returning the 1s, -1s, and 0s
    # ret<-ldply(res.list, function(x) x$effects)
    ret <- lapply(res.list, function(x) x$effects)
    ret <- do.call(rbind, ret) %>% as.data.frame()
  }
  
  rownames(ret) <- vars
  return(ret)
}

filterOverData2 <- function (overData, type = "positive") 
{
  neg <- which(overData < 0, arr.ind = T)
  pos <- which(overData > 0, arr.ind = T)
  if (type == "positive" & nrow(neg)>0) 
    apply(neg, 1, function(x) overData[x[1], x[2]] <<- 0)
  if (type == "negative" & nrow(pos)>0) 
    apply(pos, 1, function(x) overData[x[1], x[2]] <<- 0)
  if (type == "negative" & nrow(neg)>0) 
    apply(neg, 1, function(x) overData[x[1], x[2]] <<- 1)
  if (type == "all" & nrow(neg)>0) 
    apply(neg, 1, function(x) overData[x[1], x[2]] <<- 1)
  overData
}

#########
# sAICfun2 takes a dataset, response, and function, and then uses a stepAIC approach
# to determine the best model.  From that it extracts the species with a positive,
# negative, and neutral effect on that function
#########
sAICfun2 <- function(response, species, data, positive.desired=T, method="lm", combine="+", k=2, ...){
  # first fit the model
  obj<-sAICFit2(response, species, data, method, combine, k=k, ...)
  
  # now extract the important information about positive, negative, etc.
  
  # return that info in a list
  if(positive.desired) {
    pos.sp <- names(summary(obj)[[4]][,4][(summary(obj)[[4]][,1]>0) & names(summary(obj)[[4]][,4])!="(Intercept)"])
    neg.sp <- names(summary(obj)[[4]][,4][(summary(obj)[[4]][,1]<0) & names(summary(obj)[[4]][,4])!="(Intercept)"])
  }else{
    pos.sp <- names(summary(obj)[[4]][,4][(summary(obj)[[4]][,1]<0) & names(summary(obj)[[4]][,4])!="(Intercept)"])
    neg.sp <- names(summary(obj)[[4]][,4][(summary(obj)[[4]][,1]>0) & names(summary(obj)[[4]][,4])!="(Intercept)"])
  }
  neu.sp <- species[!(species %in% pos.sp) & !(species %in% neg.sp)]
  
  # make a vector of 1s and 0s
  effects<-rep(0, length(species))
  names(effects)<-species  
  effects[which(names(effects) %in% pos.sp)]<-1
  effects[which(names(effects) %in% neg.sp)]<- -1
  
  coefs<-c(effects,0)
  names(coefs)[length(coefs)]<-"(Intercept)"
  coefs[ match(names(coef(obj)), names(coefs)) ]<-coef(obj)
  
  return(list(pos.sp=pos.sp,neg.sp=neg.sp,neu.sp=neu.sp, functions=response, coefs=coefs, effects=effects))
}

#########
# sAICFit2 does the business of fitting a model using a stepAIC approach
#########
sAICFit2 <- function(response, species, data, method="lm", combine="+", k, ...){
  # make sure the dplyr package is installed
  if(! "dplyr" %in% installed.packages()[,1]) stop(
    "this function requires the dplyr package to be installed"
  )
  f <- as.formula(paste(response, "~" , paste(species, collapse="+")))
  fit <- eval(substitute(lm(f, data=data, ...)))
  obj <- MASS::stepAIC(fit, trace=0, k=k)
  obj
}
