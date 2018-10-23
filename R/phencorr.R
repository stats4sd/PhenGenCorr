#' Phenotypic Correlation
#'
#' Phenotypic Correlation (as implemented in META-R with error related to missing the reference level fixed )
#' @param genotypes vector containing genotypes
#' @param environments vector containing environments
#' @param outcomes vector of names of outcome variables to calculate correlation on
#' @param data dataframe containing outcomes
#' @param rep vector containing reps
#' @param block vector containing blocks
#' @keywords heritability
#' @importFrom lme4 lmer fixef
#' @export
#' @examples
#' library(agridat)
#' data(yacon)
#' yaconF0<-subset(yacon,dose=="F0")
#' with(yaconF0,phencorr(genotypes=entry,environments=locality,outcomes=c("height","stalks","wfr","wff","wfk","roots","FOS","glucose",
#' "fructose","brix","foliage","dry","IH"),data=yaconF0,rep=replication))
phencorr<-function (genotypes, environments = NULL, outcomes,data, rep,block=NULL)
{
  OutMat<-NULL
  if(length(outcomes)>1){
    for(i in 1:length(outcomes)){
      if(is.null(block)){
        data1 <- data.frame(variety = factor(genotypes), environment = factor(environments),
                          y = scale(data[,outcomes[i]]), rep = rep)
        mt<-lmer(y ~ -1+variety+(1 | environment) + (1 | environment:rep) +
                   (1 | variety:environment), data1)
      }
      else{
        data1 <- data.frame(variety = factor(genotypes), environment = factor(environments),
                            y = scale(data[,outcomes[i]]), rep = rep,block=block)
        mt<-lmer(y ~ -1+variety+(1 | environment) + (1 | environment:rep) +  (1 | environment:rep:block) +
                   (1 | variety:environment), data1)





      }
      OutMat<-cbind(OutMat,fixef(mt))
    }

    colnames(OutMat)<-outcomes
    return(cor(OutMat))
  }
}
