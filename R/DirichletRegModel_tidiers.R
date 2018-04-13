#' Tidying methods for a Dirichlet regression model
#' 
#' These methods tidy the coefficients of a Dirichlet regression models into a summary,
#' augment the original data with information on the fitted values and
#' residuals, and construct a one-row glance of the model's statistics.
#'
#' @details If you have missing values in your model data, you may need to refit
#' the model with \code{na.action = na.exclude}.
#'
#' @return All tidying methods return a \code{data.frame} without rownames.
#' The structure depends on the method chosen.
#'
#' @seealso \code{\link{summary.lm}}
#'
#' @name DirichletReg_tidiers
#' 
#' @param x DirichletRegModel object
#' @param data Original data, defaults to the extracting it from the model
#' @param newdata If provided, performs predictions on the new data
#' @param type.predict character vector. Values to predict; passed on to
#'   \code{\link{predict.DirichletReg}}. Valid options are "mu", "alpha" and "phi"
#' @param type.residuals character vector. Type of residuals to compute; passed on to
#'   \code{\link{residuals.DirichletReg}} Valid options are "standardized", "composite" and "raw"
#'
#' @examples
#' library(DirichletReg)
#' ALake <- ArcticLake
#' ALake$Y <- DR_data(ALake[,1:3])
#' # fit a quadratic Dirichlet regression models ("common")
#' res1 <- DirichReg(Y ~ depth + I(depth^2), ALake)
NULL


#' @rdname DirichletReg_tidiers
#' 
#' @param conf.int whether to include a confidence interval
#' @param conf.level confidence level of the interval, used only if
#' \code{conf.int=TRUE}
#' @param exponentiate whether to exponentiate the coefficient estimates
#' and confidence intervals (typical for alternative parametrization)
#' 
#' @details If \code{conf.int=TRUE}, the confidence interval is computed with
#' the \code{\link{confint}} function.
#' 
#' 
#' @return \code{tidy.DirichletRegModel} returns one row for each coefficient, with six columns:
#'   \item{response}{Which response column the coefficients correspond to}
#'   \item{term}{The term in the linear model being estimated and tested}
#'   \item{estimate}{The estimated coefficient}
#'   \item{std.error}{The standard error from the linear model}
#'   \item{statistic}{t-statistic}
#'   \item{p.value}{two-sided p-value}
#' 
#' 
#' If \code{conf.int=TRUE}, it also includes columns for \code{conf.low} and
#' \code{conf.high}, computed with \code{\link{confint}}.
#' 
#' @export
tidy.DirichletRegModel <- function(x, 
         conf.int=FALSE, 
         conf.level=.95, 
         exponentiate=FALSE, ...) {
    
    # call summary and hide output
    hide <- capture.output(s <- DirichletReg::summary.DirichletRegModel(x))
    
    co <- s$coef.mat
    colnames(co) <- c("estimate", "std.error", "statistic", "p.value")
    ret <- data.frame(response=rep(s$varnames, s$n.vars), term=rownames(co), co, row.names = NULL)
    
    #* Confidence Interval
    if (conf.int){
        ci <- DirichletReg::confint.DirichletRegModel(x, level=conf.level)
        ci <- data.frame(do.call("rbind", ci$ci[[1]]), row.names = NULL)
        ci
        names(ci) <- c("conf.low", "conf.high")
        ret <- cbind(ret, ci)
    }
    
    #* Exponentiate (for Odds Ratios)
    if (exponentiate){
        if(s$parametrization=="common") warning("Warning the model is fitted with common parametrization.")
        exp.col <- if(conf.int) c("estimate", "conf.low", "conf.high") else "estimate"
        ret[, exp.col] <- lapply(ret[, exp.col, drop=FALSE], exp)
    }
    
    unrowname(ret)
}

#' @rdname DirichletReg_tidiers
#' 
#' @template augment_NAs
#' 
#' 
#' @return When \code{newdata} is not supplied \code{augment.DirichletRegModel} returns
#' one row for each observation, with seven columns added to the original
#' data:
#'   \item{.fitted}{Fitted values of model (mu for extected values, alpha values and phi for precision values )}
#'   \item{.resid}{Residuals of type "standardized", "composite", or "raw"}
#' 
#' 
#' When \code{newdata} is supplied, \code{augment.DirichletRegModel} returns one row for each
#' observation, with multiple columns added to the new data:
#'   \item{.fitted}{Fitted values of model}
#' 
#' @export
augment.DirichletRegModel <- function(x, data = x$data, newdata=NULL,
                       type.predict=c("mu", "alpha", "phi"), 
                       type.residuals=c("standardized", "composite", "raw"), ...) {  

    if(!any(type.predict %in% c("mu", "alpha", "phi")) | is.null(type.predict)) {
        stop('type.predict must be a vector of one or more types. Validl values are: "mu", "alpha", "phi".')
    }
    
    if(!any(type.residuals %in% c("standardized", "composite", "raw")) | is.null(type.residuals)) {
        stop('type.residuals must be a vector of one or more types. Validl values are: "standardized", "composite", "raw"')
    }
        
    ret <- fix_data_frame(data, newcol = ".rownames")
    
    # what to predict
    what_predict <- rep(FALSE, 3)
    names(what_predict) <- c("mu", "alpha", "phi")
    what_predict[names(what_predict) %in% type.predict] <- TRUE
    
    # what type of residuals
    what_residuals <- rep(FALSE, 3)
    names(what_residuals) <- c("standardized", "composite", "raw")
    what_residuals[names(what_residuals) %in% type.residuals] <- TRUE
    
    if(is.null(newdata)) {
        
        # fitted values
        if(any(what_predict)) {
        fitdf <- as.data.frame(DirichletReg::fitted.DirichletRegModel(x, 
                                                                      mu=what_predict["mu"], 
                                                                      alpha=what_predict["alpha"], 
                                                                      phi=what_predict["phi"]))
        names(fitdf) <- paste0(".fitted.", names(fitdf))
        } 
        
        # residuals
        residualdf <- data.frame(lapply(names(what_residuals)[what_residuals], function(type) {
            resid <- DirichletReg::residuals.DirichletRegModel(x, type=type)
            
            if(type=="composite") {
                res <- data.frame(".resid.composite"=resid)    
            } else {
                res <- data.frame(matrix(resid, ncol=ncol(resid)))    
                names(res) <- paste0(".resid.", type, ".", colnames(resid))
            }
            res
            }))
        
        if(any(what_predict) & any(what_residuals)) {
            return(cbind(ret, fitdf, residualdf))
        } else if(any(what_predict)) {
            return(cbind(ret, fitdf))
        } else if(any(what_residuals)) {
            return(cbind(ret, residualdf))
        } else {
            return(ret)
        }
        
    } else {
        ret <- fix_data_frame(newdata, newcol = ".rownames")
        fitdf <- as.data.frame(DirichletReg::predict.roxygen2(x, newdata, mu=what["mu"], alpha=what["alpha"], phi=what["phi"]))
        names(fitdf) <- paste0(".fitted.", names(fitdf))
        return(cbind(ret, fitdf))
    }
}


#' @rdname DirichletReg_tidiers
#' 
#' @param ... extra arguments (not used)
#' 
#' @return \code{glance.lm} returns a one-row data.frame with the columns
#'   \item{logLik}{the data's log-likelihood under the model}
#'   \item{AIC}{the Akaike Information Criterion}
#'   \item{BIC}{the Bayesian Information Criterion}
#'   \item{deviance}{deviance}
#'   \item{df.residual}{residual degrees of freedom}
#' 
#' @export
glance.DirichletRegModel <- function(x, ...) {
    # use summary.lm explicity, so that c("aov", "lm") objects can be
    # summarized and glanced at
    hide <- capture.output(s <- DirichletReg::summary.DirichletRegModel(x))
    ret <- data.frame(deviance=s$deviance, df.residual=s$df)
    ret <- cbind(finish_glance(ret, x), ret)
    ret
}
