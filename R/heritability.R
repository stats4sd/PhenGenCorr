#' Heritability
#'
#' Heritability as calculated in a multi environment trial (as implemented in META-R)
#' @param genotypes vector containing genotypes
#' @param environments vector containing environments
#' @param outcome vector containing outcome
#' @param rep vector containing reps
#' @param block vector containing blocks
#' @param return TRUE/FALSE return the heritability stat in the console
#' @keywords heritability
#' @importFrom lme4 lmer VarCorr
#' @export
#' @examples
#' library(agridat)
#' data(besag.met)
#' with(besag.met,heritability(genotypes=gen,environments=county,outcome=yield,rep=rep,block=block))
heritability<-
  function (genotypes, environments = NULL, outcome, rep,block=NULL,return=TRUE)
  {
    require(lme4)
    if (length(environments) > 0) {
      if(is.null(block)){
        data <- data.frame(variety = factor(genotypes), environment = factor(environments),
                           y = outcome, rep = rep)
        mod1 <- lmer(y ~ (1 | environment) + (1 | environment:rep) +
                       (1 | variety) + (1 | variety:environment), data)
      }
      else{
        data <- data.frame(variety = factor(genotypes), environment = factor(environments),
                           y = outcome, rep = rep,block=block)
        mod1 <- lmer(y ~ (1 | environment) + (1 | environment:rep) + (1 | environment:rep:block)+
                       (1 | variety) + (1 | variety:environment), data)



      }

      varcorr <- VarCorr(mod1)
      varU1 <- as.vector(varcorr$'variety')
      varErr <- attr(varcorr,'sc')^2
      varGE <- as.vector(varcorr$'variety:environment')
      varLoc <- as.vector(varcorr$'environment')
      nRep <- mean(table(data$variety,data$environment))
      nLoc <- length(unique(data$environmen))
      h2 <- varU1/(varU1 + varGE/nLoc + varErr/(nRep*nLoc))
      if(return==TRUE){
 print(h2)
      }
      output<-list(h2 = h2, comps=data.frame(varU1,varErr,varGE,varLoc,nRep,nLoc),model = mod1)
      return(invisible(output))
    }
  }
