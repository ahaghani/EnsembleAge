
# Author: Amin Haghani
# date: 06/01/2023

# Clock coefficients
epiclocks <- readRDS("manifests/MouseClocks/Clock coefficients.RDS")

#




# transformation functions
trafo= function(x,offset=0.06,adult.age=1.2) {
  y=ifelse(x<=adult.age, log(x+offset),x/(adult.age+offset) +log(adult.age+offset)- adult.age/(adult.age+offset)  );
  y }

anti.trafo= function(x,offset=0.06,adult.age=1.2) {
  ifelse(x<=log(adult.age+offset), exp(x)-offset, (adult.age +offset)*x-log(adult.age+offset)*( adult.age+offset) +adult.age) }

## Universal clock functions
F2_antitrans<-function(y,y.maxAge,y.gestation,const=1){
  x0=const*exp(-exp(-1*y))
  x1=x0*(y.maxAge+y.gestation)
  x=x1-y.gestation
  x
}

#===================================
#Funtions for clock3
#===================================
F1_logli <- function(age1, m1, m2 = m1, c1=1){
  ifelse(age1 >= m1, (age1-m1)/m2 , c1*log((age1-m1)/m2/c1 +1) )
}

F2_revtrsf <- function(y.pred, m1, m2 = m1, c1=1){
  ifelse(y.pred<0, (exp(y.pred/c1)-1)*m2*c1 + m1, y.pred*m2+m1 )
}

############# 
# The `loglifn` function shows how to calculate m1 for the transformation
# It is the `a_Logli` in the function

F3_loglifn = function(dat1,b1=1,max_tage = 4,
                      c1=5, c2 = 0.38, c0=0){
  n=nrow(dat1)
  
  age1 = (dat1$maxAgeCaesar+dat1$GestationTimeInYears)/(dat1$averagedMaturity.yrs+dat1$GestationTimeInYears)
  
  a1 = age1/(1+max_tage)
  dat1$a1_Logli = a1 #x/m1 in manuscript
  
  a2 = (dat1$GestationTimeInYears + c0)/(dat1$averagedMaturity.yrs)
  dat1$a_Logli = a_Logli = c1*a2^c2
  #m=5*(G/ASM)^0.38 from regression analysis/formula(7)
  
  
  x = dat1$Age + dat1$GestationTimeInYears
  t2 = dat1$averagedMaturity.yrs*b1 + dat1$GestationTimeInYears
  x2 = x/t2 #### log(x/t2)
  y = F1_logli(x2, a_Logli, a_Logli)
  
  dat1$LogliAge <- y
  return(dat1)
}


# 

predictAgeAndAgeAcc <- function(dat0sesame, samps){
  
  species <-  read.csv("manifests/MouseClocks/anAgeUpdatedCaesarVersion51.csv")
  species$SpeciesLatinName <- sapply(1:nrow(species), function(x){
    if(is.na(species$SpeciesLatinName[x])){
      species$SpeciesLatinName[x] <- stringr::str_replace(species$labelsCaesar[x], "_", " ")
    } else {species$SpeciesLatinName[x] <- as.character(species$SpeciesLatinName[x])}
  })
  
  species <- species %>% mutate(GestationTimeInYears = Gestation.Incubation..days./365)%>%
    mutate(logGest = log(Gestation.Incubation..days.), logLifespan = log(maxAgeCaesar), logWeight = log(weightCaesar))%>%
    dplyr::select(SpeciesLatinName, GestationTimeInYears, averagedMaturity.yrs, maxAgeCaesar, weightCaesar, logGest, logLifespan, logWeight, MammalNumberHorvath)
  
  n <- which(names(samps)%in%names(species)[-1])
  if(sum(n)>0){
    samps <- samps[,-n]
  }
  
  
  if(!"SpeciesLatinName"%in%names(samps)){
    samps <- samps%>% mutate(SpeciesLatinName="Mus musculus")
  }
  
  if(!"Age"%in%names(samps)){
    samps <- samps%>% mutate(Age=0)
  }
  
  if(!"Female"%in%names(samps)){
    samps <- samps%>% mutate(Female=NA)
  }
  
  if(!"Tissue"%in%names(samps)){
    samps <- samps%>% mutate(Tissue=NA)
  }
  
  samps <- samps %>% mutate(Age=ifelse(is.na(Age), 0, Age))%>% 
    mutate(SpeciesLatinName=ifelse(is.na(SpeciesLatinName), "Mus musculus", SpeciesLatinName))
  
  
  dat0sesame <- dat0sesame%>% tibble::column_to_rownames("CGid") %>% t(.) %>% as.data.frame(.)%>%mutate('Intercept'=1)
  samps <- samps%>% left_join(species, by="SpeciesLatinName")
  
  ageResults <- lapply(1:length(epiclocks), function(j){
    
    clocks <- epiclocks[[j]]
    
    # a loop for the main category of clocks
    ageAccel <- plyr::llply(1:length(clocks), function(i){
      #cat(paste0("Predicting age for ", names(epiclocks)[j], " clock ", names(clocks)[i], "\n"))
      clock <- clocks[[i]]
      
      dat1 <- dat0sesame %>% dplyr::select(clock$CGid)
      
      if(grepl("AgeTraf", names(epiclocks)[j])){
        samp <- samps%>%
          # predicting age
          mutate(epiAge = as.numeric(as.matrix(dat1)%*%clock$Coef))%>% mutate(AgeAccelation = as.vector(residuals(lm(epiAge~Age))))%>% dplyr::select(Basename, Age, AgeAccelation, epiAge, Female, Tissue) %>%
          # age tranformation
          mutate(epiAge = anti.trafo(epiAge))%>% mutate(AgeAccelation = as.vector(residuals(lm(epiAge~Age))))
      } else if(grepl("uniClocks", names(epiclocks)[j])&grepl("(UniClock3)|(UniBloodClock3)|(UniSkinClock3)", names(clocks)[i])){
        
        samp <- F3_loglifn(samps)#to compute m estimate for tuning point in the log-linear transformation
        samp <- samp%>%
          # predicting age
          mutate(epiAge = as.numeric(as.matrix(dat1)%*%clock$Coef))  %>% mutate(m1=a_Logli)%>%
          mutate(epiAge=F2_revtrsf(epiAge, m1) *(averagedMaturity.yrs + GestationTimeInYears) -GestationTimeInYears)%>%
          dplyr::select(-m1, -a1_Logli, -a_Logli, -LogliAge)%>% mutate(AgeAccelation = as.vector(residuals(lm(epiAge~Age))))%>% dplyr::select(Basename,epiAge, AgeAccelation, Age, Female, Tissue)
      }else if(grepl("uniClocks", names(epiclocks)[j])&grepl("(UniClock2)|(Ake.DuoHumanMouse)|(UniBloodClock2)|(UniSkinClock2)|(Ake.DNAmDuoGrimAge411)", names(clocks)[i])){
        
        samp <- samps%>%
          # predicting age
          mutate(epiAge = as.numeric(as.matrix(dat1)%*%clock$Coef)) %>% mutate(HighmaxAgeCaesar= maxAgeCaesar*1.3) %>%
          mutate(HighmaxAgeCaesar=ifelse(SpeciesLatinName%in%c("Homo sapiens", "Mus musculus"), maxAgeCaesar, HighmaxAgeCaesar))%>%
          # age transformation
          mutate(epiAge = F2_antitrans(epiAge, y.maxAge = HighmaxAgeCaesar, y.gestation = GestationTimeInYears, const = 1)) %>% mutate(AgeAccelation = as.vector(residuals(lm(epiAge~Age))))%>% dplyr::select(Basename,epiAge, AgeAccelation, Age, Female, Tissue)
      }else{
        samp <- samps %>%
          # predicting age
          mutate(epiAge = as.numeric(as.matrix(dat1)%*%clock$Coef))%>% mutate(AgeAccelation = as.vector(residuals(lm(epiAge~Age))))%>% dplyr::select(Basename, Age, AgeAccelation, epiAge, Female, Tissue)
      }
      
      return(samp)
      
    })
    names(ageAccel) <- names(clocks)
    ageAccel <- rbindlist(ageAccel, idcol = "epiClock")
  })
  
  names(ageResults) <- names(epiclocks)
  
  ageResults <- rbindlist(ageResults, idcol = "clockFamily", fill = T) %>% 
    mutate(clockFamily= ifelse(grepl("(UniClock2)|(UniBloodClock2)|(UniSkinClock2)",epiClock), "Ake.UniClock2", 
                         ifelse(grepl("(UniClock3)|(UniBloodClock3)|(UniSkinClock3)",epiClock), "Ake.UniClock3",
                          ifelse(grepl("(Ake.DuoHumanMouse)",epiClock),"Ake.DuoHumanMouse", 
                           ifelse(grepl("(Ake.DNAmDuoGrimAge411)",epiClock),"Ake.DNAmDuoGrimAge411",      
                             clockFamily)))))%>%
    mutate(epiClock = ifelse(grepl("Skin", epiClock), "Skin", 
                             ifelse(grepl("Blood", epiClock), "Blood", 
                                ifelse(grepl("(UniClock)|(DuoHumanMouse)|(Ake.DNAmDuoGrimAge)", epiClock), 
                                                                            "panTissue", epiClock)))) %>%
    mutate(epiClock=factor(epiClock, levels = c("Blood","Liver","Heart","Kidney","Muscle","Cortex","Striatum","Cerebellum", "Brain", "Fibroblast", "Tail", "panTissue", "Skin"),
                           labels = c("Blood","Liver","Heart","Kidney","Muscle","Cortex","Striatum","Cerebellum", "Brain", "Fibroblast", "Tail", "panTissue", "Skin")))
  
  # convert to wide format
  targetClocks <- c("Amin.LUC.v1", "Steve.ElasticEpi.AgeTraf", 
                    "Steve.InterventionClock.AgeTraf", "Steve.DevelopmentClock.AgeTraf", 
                    "Ake.UniClock2", "Ake.UniClock3", "Ake.DNAmDuoGrimAge411", "Ensemble.Static", "Ensemble.Static.Top")
  newNames <- c("LifespanUberClock", "DNAmAgeFinal", 
                "DNAmAgeInterventionFinal", "DNAmAgeDevelopmentFinal",
                "UniversalClock2", "UniversalClock3", "DuoGrimAge411", "Ensemble.Static", "Ensemble.Static.Top")
  
  dat1 <- ageResults %>% mutate(clockFamily=factor(clockFamily, levels=targetClocks, labels=newNames)) %>% 
    mutate(epiClock = paste(clockFamily, epiClock,"clock", "epiAge", sep=".")) %>% 
    dplyr::select(Basename, epiClock, epiAge) %>% spread(key = "epiClock", value = "epiAge") 
  
  dat2 <- ageResults %>% mutate(clockFamily=factor(clockFamily, levels=targetClocks, labels=newNames)) %>% 
    mutate(epiClock = paste(clockFamily, epiClock,"clock", "AgeAcceleration", sep=".")) %>% 
    dplyr::select(Basename, epiClock, AgeAccelation) %>% spread(key = "epiClock", value = "AgeAccelation") %>% 
    left_join(dat1, by="Basename")
  
  return(dat2)
}




