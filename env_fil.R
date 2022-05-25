env_fil=function(traits,otu,take_rows=nrow(otu), quant=5, simul=499, correction="sqrt", mc=FALSE){ #start of the function
  #Calculate original functional richness (trait hypervolume) and levels of species richness across the dataset
  Fric_origin<-dbFD(traits,otu,calc.FGR=FALSE, calc.FDiv=FALSE, calc.CWM=FALSE, messages=FALSE, corr=correction)
  spp_sel<-c(unique(specnumber(otu)))[c(unique(specnumber(otu)))>=3]
  #Start the randomization for each levels of species number. Minimum spp number=3. The function includes an if/else statement to compare the PCoA dimension reduction with the original. Not being the same would render non comparable results 
  rand<-function(var1, otu,traits, simul, corr, take_rows){
    temp_vec<-c()
    for (k in 1:simul){
      rand_col<-sample(1:nrow(traits), var1, replace=F)#Create a random subset of species 
      otu_temp<-otu
      otu_temp[1,rand_col]<-1#The random subset is the first row of the original dataset
      otu_temp[1,-c(rand_col)]<-0
      otu_temp[2,colSums(otu_temp[1:take_rows,])==0]<-1 #Avoid the columns with total zeroes that would stop the calculation.  
      Fric_rand<-dbFD(traits, otu_temp[1:take_rows,],calc.FGR=FALSE, calc.FDiv=FALSE, calc.CWM=FALSE, message=FALSE, corr=correction)#Calculate the functional richness of a part or all the samples of the dataset. 
      if(Fric_origin$qual.FRic==Fric_rand$qual.FRic) {#Compare the original reduction of dimensions to that of the obtained in line 15. If the same, keep the Functional richness of the sample with randomized species. 
        temp_vec[k]<-Fric_rand$FRic[1] 
      } else { 
        stop("increase value of take_rows")
      }
    }
    temp_vec
  }
  if(mc==TRUE) {#Implement the function on the dataset. It allows multicore processing if the argument of env_fil function mc=TRUE.
    null_vectors<-mclapply(spp_sel,mc.cores=detectCores(), rand, otu=otu, traits=traits, simul=simul, corr=correction, take_rows=take_rows)
  } else {null_vectors<-lapply(spp_sel,rand, otu=otu, traits=traits, simul=simul, corr=correction, take_rows=take_rows)
  }
  names(null_vectors)<-c(paste(spp_sel,"spp", sep="_"))
  #Summarize the distribution values of the randomized functional richness
  null_model<-lapply(null_vectors, quantile, probs=c(0.025,0.05, 0.10, 0.5, 0.90, 0.95,0.975))
  
  #Compare the observed functional richness to the specified quantile value in the "quant" argument. It assigns the putative process 
  r_fric<-cbind(specnumber(otu), Fric_origin$FRic)
  compare<-function(r_fric,null_model, quant){
    if(r_fric[1]<3){
      "NA"
    } else if (r_fric[2] < null_model[[paste(r_fric[1], "spp", sep="_")]][paste0(quant,"%")]){
      "env_fil"
      #print("env_fil")
    } else if (r_fric[2] > null_model[[paste(r_fric[1], "spp", sep="_")]][paste0(100-quant,"%")]){
      "comp_excl"
      #print("comp_excl")
    } else {
      "neutral"
      #print("neutral")
    }
  }
  comparison<-apply(r_fric,1, compare, null_model=null_model, quant=quant)#Implement function
  list("null_model"=null_model,"observed_FRic"=Fric_origin$FRic,"observed_sppnumber"=specnumber(otu), "process"=comparison)#Output randomized distribution, observed richness, functional richness and process
}#End of function
