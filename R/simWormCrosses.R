simWormCrosses = function(SimWormParams, sexChr=F) {

    #instead of attach()
    starting.sample.size=SimWormParams$starting.sample.size
    #brood size if it selfs
    brood.size.selfing=SimWormParams$brood.size.selfing
    #brood size for mating 
    brood.size.mating=SimWormParams$brood.size.mating
    #fraction of hermaphrodites that self
    selfing.rate=SimWormParams$selfing.rate
    #bottleneck per generation   
    max.per.gen=SimWormParams$max.per.gen
    #how many total generations
    max.gen=SimWormParams$max.gen
    #Founder population genotypes
    FN=SimWormParams$FN
    #simulation parameters 
    SP=SimWormParams$SP
    mating.matrix=SimWormParams$mating.matrix
    selfings=SimWormParams$selfings
    
    genMap=SimWormParams$genMap
    QTL.sims=SimWormParams$QTL.sims

    #setup of F0
    nProgenyF1.mated=starting.sample.size*(1-selfing.rate)*brood.size.mating
    nProgenyF1.selfed=starting.sample.size*(selfing.rate)*brood.size.selfing

    #third columnm, 0 indicates mating, 1 indicates selfing
    mating.matrix=cbind(mating.matrix,0)
    #generate progeny from matings
    mating.matrix=do.call(rbind, replicate(nProgenyF1.mated, mating.matrix, simplify=F))

    #generate progeny from selfing
    self.cnt=nProgenyF1.selfed
    if(self.cnt>0) {
           selfing.matrix=cbind(selfings, selfings,1)
           selfing.matrix=do.call(rbind, replicate(nProgenyF1.selfed, selfing.matrix, simplify=F))
           mating.matrix=rbind(mating.matrix, selfing.matrix)
           mating.matrix=mating.matrix[order(mating.matrix[,3]),]
    }
    # downsample, a fourth column could be added prior to downsampling to propagate through an expected genetic effect,
    # then a hard threshold or weighted sampling could occur to pick who crosses for the next generation 
    if(nrow(mating.matrix)>max.per.gen){
            mating.matrix=mating.matrix[sample(1:nrow(mating.matrix), max.per.gen),]
            mating.matrix=mating.matrix[order(mating.matrix[,3]),]
    }
    #-----------------------------------------------------------------------------------------------

    # current generation 
    current.gen=1
    while(current.gen<=max.gen) {
       
       mated.index=which(mating.matrix[,3]==0)
       selfed.index=which(mating.matrix[,3]==1)
       print(paste('progeny from mating', length(mated.index)))
       print(paste('progeny from selfing', length(selfed.index)))

       #pilot code to deal with X chromosome goes here ----------------
       if(sexChr==T) {
           matings=mating.matrix[mated.index,-3]
           #reduce ploidy for males 
           rG=reduceGenome(FN, simRecomb=F, simParam=SP)
           # hack so we can cross 
           rGr=reduceGenome(FN,simRecomb=T, simParam=SP)
           #FNmx=doubleGenome(rG)
           #FN=makeCross2(FN,FNmx, matings ,simParam=SP)
           FNr=mergeGenome(rG, rGr, matings, simParam=SP)

           if(length(selfed.index)>0) {
                selfings=mating.matrix[selfed.index,-3]
    #        #       #fake X as doubled haploid
                 print(nrow(selfings))
                 print((FN))
                selfings.cross=makeCross(FN,selfings, nProgeny=1,simParam=SP)
                FN=c(FNr,selfings.cross)
           } else {
                FN=FNr
           } 
       }
       #-------------------------------------------------------------
       if(sexChr==F) {
            # for normal autosomes 
            FN=makeCross(FN, mating.matrix[,-3], nProgeny=1,simParam=SP)
       }
       
       ######---calculate phenotype effect for each individual in FN here ----------------------
        if( !is.null(QTL.sims$sim.fitness) & QTL.sims$sim.fitness) {
            X_Q=pullMarkerGeno(FN, QTL.sims$f.add.qtl.ind, asRaw=F)
            X_Beta=QTL.sims$f.add.qtl.eff
            if(length(X_Beta)==1) {
                XB=X_Q*X_Beta
            } else {XB=X_Q%*%X_Beta  }
            simy=XB+rnorm(nrow(X_Q), mean=0, sd=QTL.sims$f.error.sd) 
            h2=var(XB)/(var(XB)+QTL.sims$f.error.sd^2)
            print(h2)
            mprob=pnorm(scale(simy))
        }
       else{ mprob=rep(1,nInd(FN))  }#; print(mprob)      }
       #----------------------------------------------------------------------------------------


       #assume half of the progeny from matings will be hermaphrodites
       herm.index=sort(sample(mated.index,size=length(mated.index)/2))


       #the rest of the mated progeny will be males 
       male.index=which(!((mated.index %in% herm.index)))
       #all of the selfed progeny will by hermaphrodites
       if(length(selfed.index)>0) {
            herm.index=c(herm.index,selfed.index)
        }
       print(paste('total herm', length(herm.index)))
       print(paste('total male', length(male.index)))
       
       
       #for all the males, sample from hermaphrodites for them to mate with
       #probability of mating 
       #mated.herm=sample(herm.index, size=min(length(herm.index), length(male.index)),prob=mprob[herm.index])
       #mated.male=sample(male.index, size=min(length(herm.index), length(male.index)),prob=mprob[male.index])

       #mating.matrix=cbind(mated.herm, mated.male, 0)
       mated.herm=sample(herm.index, size=min(length(herm.index), length(male.index)))
       mating.matrix=cbind(mated.herm,male.index, 0)

       #generate cross progeny for matings, with brood size dependent on expected brood size if mating 
       mating.matrix=do.call(rbind, replicate(brood.size.mating, mating.matrix, simplify=F))
      
       #note which hermaphrodites aren't mated 
       unmated.herm=herm.index[!(herm.index %in% mated.herm)]

       #how many do we expected to self given the total number of hermaphrodites and the selfing rate 
       self.cnt= round(length(herm.index)*selfing.rate)
       #if we expect any selfers then ... 
       if(self.cnt>0) {
           #figure out which ones self 
           #selfings=sample(unmated.herm, min(length(unmated.herm),self.cnt), prob=mprob[unmated.herm])
           selfings=sample(unmated.herm, min(length(unmated.herm),self.cnt))
           #make a selfing matrix
           selfing.matrix=cbind(selfings, selfings,1)
           #make brood size reflect brood size per selfing
           selfing.matrix=do.call(rbind, replicate(brood.size.selfing, selfing.matrix, simplify=F))
           mating.matrix=rbind(mating.matrix, selfing.matrix)
        }
       #total worms at this point z
       print(nrow(mating.matrix))

       #bottleneck
       if(nrow(mating.matrix)>max.per.gen){
            
            mprob.lookup=mprob[mating.matrix[,1]]
            #simulated brood is being downsampled here based on fitness effect of the hermaphrodite parent
            mating.matrix=mating.matrix[sample(1:nrow(mating.matrix), max.per.gen, prob=mprob.lookup),]
            mating.matrix=mating.matrix[order(mating.matrix[,3]),]
       }
       current.gen=current.gen+1
       FN@sex[male.index]='M'
    }
    return(FN)
}

