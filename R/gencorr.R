#' Genotypic Correlation
#'
#' Genotypic Correlation
#' @param genotypes vector containing genotypes
#' @param environments vector containing environments
#' @param outcomes vector of names of outcome variables to calculate correlation on
#' @param data dataframe containing outcomes
#' @param rep vector containing reps
#' @param block vector containing blocks
#' @keywords heritability
#' @importFrom lme4 lmer fixef
#' @importFrom reshape2 dcast
#' @export
#' @examples
#' library(agricolae)
#' data(yacon)
#' yaconF150<-subset(yacon,dose=="F150")
#' with(yaconF150,gencorr(genotypes=entry,environments=locality,outcomes=c("height","stalks","wfr","wff","wfk","roots","FOS","glucose",
#' "fructose","brix","foliage","dry","IH"),data=yaconF0,rep=replication))
gencorr<-function (genotypes, environments = NULL, outcomes,data, rep,block=NULL,threshold=0.05)
{
  if(length(outcomes)>1){
    if(is.null(block)){
      data1 <- data.frame(variety = factor(genotypes), environment = factor(environments),
                          data[,outcomes], rep = rep)
    }
    else{
      data1 <- data.frame(variety = genotypes, environment = environments,
                          data[,outcomes], rep = rep,block=block)
    }
    OutFrame<-expand.grid(v1=outcomes,v2=outcomes)
    OutFrame$v1<-as.character(OutFrame$v1)
    OutFrame$v2<-as.character(OutFrame$v2)
    OutFrame$gencor<-NA

    pb <- progress_bar$new(total = nrow(OutFrame))

    for(i in 1:nrow(OutFrame)){
      if(OutFrame$v1[i]==OutFrame$v2[i]){
        OutFrame$gencor[i]<-1
      }
      else{
        tmp<-data1
        tmp$v1<-scale(data1[,OutFrame$v1[i]])
        tmp$v2<-scale(data1[,OutFrame$v2[i]])
        tmp$v3<-tmp$v1+tmp$v2
        h1<-heritability(genotypes=genotypes, environments = environments, outcome=tmp$v1, rep=rep,block=block,return=FALSE)
        h2<-heritability(genotypes=genotypes, environments = environments, outcome=tmp$v2, rep=rep,block=block,return=FALSE)
        h3<-heritability(genotypes=genotypes, environments = environments, outcome=tmp$v3, rep=rep,block=block,return=FALSE)
if(h1$h2>threshold&h2$h2>threshold){
        cov1.2 <- (h3$comps$varU1-h1$comps$varU1-h2$comps$varU1)/2
        OutFrame$gencor[i] <- cov1.2/(sqrt(h1$comps$varU1)*sqrt(h2$comps$varU1))
}
        else{
          OutFrame$gencor[i] <- NA
        }
      }
      pb$tick()

    }
    OutFrame$gencor<-ifelse( OutFrame$gencor<(-1),-0.9999,ifelse( OutFrame$gencor>1,0.9999, OutFrame$gencor))
    gencor<-dcast(v1~v2,data=OutFrame,value.var="gencor")
    rownames(gencor)<-gencor$v1
    gencor<-gencor[,-1]
    gencor<-gencor[outcomes,outcomes]

    return(gencor)
  }
}
