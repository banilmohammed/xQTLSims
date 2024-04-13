# forward stepwise procedure with FDR control
# G'sell 2013 procedure to detect QTL per trait

#trait and genos need to be scaled for this to work, default nperm is dangerous for very large datasets 
doTraitFDR=function(trait, genos, genos.full, FDR_thresh=.05, nperm=1e4, doLODdrop=T) {
    f.found=c()
    p.found=c()
    q.found=c()
    m.found=c()

    n=length(trait)
    L= (crossprod(trait,genos)/(n-1))^2 
    mLi=which.max(L)
    mL=max(L)
    
    yperm=replicate(nperm, sample(trait))
    nullD=(crossprod(yperm,genos)/(n-1))^2
    
    permMax=Rfast::rowMaxs(nullD,value=T)
    pNull=1-ecdf(permMax)(mL)
    if(pNull==0) {pNull=1/nperm}
    
    step=1
    
    repeat{
       p.temp=c(p.found, pNull)
       q=-mean(log(1-p.temp))
       if(q>FDR_thresh) {break;}
       p.found=c(p.found, pNull)
       q.found=c(q.found, q)
       m.found=c(m.found, colnames(genos)[mLi])
       f.found=c(f.found, mLi)
       print(paste('step=', step, 'max index=', colnames(genos)[mLi], 'max r^2=', mL, 'pnull=', pNull, 'fdr=', q))
       yr=scale(residuals(lm(trait~genos[,f.found]) ))
       L=(crossprod(yr,genos)/(n-1))^2 
       mLi=which.max(L)
       mL=max(L)
       yperm=replicate(nperm, sample(yr))
       nullD=(crossprod(yperm,genos)/(n-1))^2
       permMax=Rfast::rowMaxs(nullD, value=T) 
       pNull=1-ecdf(permMax)(mL)
       if(pNull==0) {pNull=1/nperm}
       step=step+1
   }
   results=data.frame(fscan.markers=m.found, index=f.found, p=p.found, q=q.found, stringsAsFactors=F) 
   if(doLODdrop) {
       drops=doLODdrop(trait, genos.full, results$fscan.markers)
       results=cbind(results,drops)
   }
   return(results)
}
fasterLOD=function(n.pheno, pheno.s,gdata.s, betas=FALSE, sdx=1, pheno=NULL){
   r=crossprod(pheno.s, gdata.s)/(n.pheno-1)
   LOD=(-n.pheno*log(1-r^2))/(2*log(10))
   if(betas==FALSE) {
       return(LOD)
   } else {
      # beta=r*apply(cbind(pheno),2, sd,na.rm=T)/sdx
       return(list(r=r, LOD=LOD))
   }
}

# calculate 1.5 LOD drop confidence intervals
doLODdrop=function(trait, genos.full, f.found) {
    ys=trait
    gs=genos.full
    nsegs=length(ys)
    s=-nsegs/2
   
    #print(f.found)
    doMC::registerDoMC(cores=length(f.found))
    #for likelihood calc compute shortcut see eq 11 
    # ll=function(n,RSS) {-(n/2)*log(RSS)-(n/2)*log((2*pi)/n)-(n/2) }
    #https://statproofbook.github.io/P/rsq-mll.html
    located=foreach::foreach(j=1:length(f.found), .combine='rbind') %do% { 
            print(j)
            # in 1:nrow(zf5)) { 
            #jb mod 12/13/23 speedup-------------------------
            #nm=lm(ys~gs[,f.found[-j]]-1)
            #nllik=logLik(nm)/(log(10))
            if(length(f.found)==1){
                mm = matrix(1, nrow=nsegs)
            } else {
                mm=gs[,f.found[-j]]
            }
            lNSS=log(sum((lm.fit(mm,ys))$residuals^2))
            #------------------------------------------------

            coi=strsplit(f.found[j], '_')[[1]][1]
            gcoi=gs[,grep(paste0('^', coi,'_'), colnames(gs))]
            mnames=colnames(gcoi)
            LOD=rep(0, ncol(gcoi))
            for(g in 1:ncol(gcoi)) { #ncol(gcoi)){
                #if(g%%100==0) {print(g)}
                #jb mod 12/13/23 speedup -------------------
                mm2=cbind(mm, gcoi[,g])
                lRSS=log(sum((lm.fit(mm2,ys))$residuals^2))
                LOD[g]= s*(lRSS-lNSS)/log(10)
                #system.time({
                #LOD[g]=(logLik(lm(ys~gs[,f.found[-j]]+gcoi[,g]-1))/log(10))-nllik
            }
           return(data.frame(LOD=max(LOD), pmarker=mnames[which.max(LOD)],
                             CI.l=mnames[min(which(LOD>max(LOD)-1.5))],
                             CI.r=mnames[max(which(LOD>max(LOD)-1.5))], stringsAsFactors=F))
    }
    return(located)
}

#add jitter to a genetic map 
jitterGmapVector=function(themap, amount=1e-6) {
    for (i in 1:length(themap)) {
         n <- length(themap[[i]])
         themap[[i]] <- themap[[i]] + c(0, cumsum(rep(amount, n - 1)))
    }
    return(themap)
}

#convert physical position to genetic map position 
getGmapPositions=function(vcf.cross, gmap, uchr) {
    #get physical position, split by chromosome
    p.by.chr=split(vcfR::getPOS(vcf.cross),vcfR::getCHROM(vcf.cross))
    #keep things sorted (yay yeast chr names with roman numerals)
    p.by.chr=p.by.chr[uchr]

   #where to put the variant sites, impute onto gmap
    imputed.positions=mapply( 
           function(x, y){
                approxfun(y$ppos, y$map, rule=2)(x)
            },
            x=p.by.chr, y=gmap,
            SIMPLIFY=F)

}

#input is rds object of Rockman N2 x Hawaii F10 AIL rils genetic map
#output is genetic map for F2 progeny
#note units here are centimorgans
restructureGeneticMap=function(gmap.rds, expansion.factor=1/5, expandX=T){
    #genetic map from F10 AIL 
    gmapRef=readRDS(gmap.rds) #'/data0/elegans/xQTLSims/geneticMapXQTLsnplist.rds')
    #normalize map to what's expected from F2
    #https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000419
    gmapRef$map=gmapRef$map*expansion.factor#5.3
    #reorder
    gmapRef=gmapRef[,c('map', 'chrom', 'pos')]
    #rename columns (ppos = physical pos)
    names(gmapRef)[3]='ppos'
    gmapRef[,'chrom']=as.character(gmapRef[,'chrom'])
    gmapRef=split(gmapRef, gmapRef$chrom)

    gmap=gmapRef
    #sapply(gmap, function(x) max(x$map))
    #force obligate chiasma on X
    if(expandX) {
    gmap$X$map=gmap$X$map*(51/max(gmap$X$map)) 
    }

    return(gmap)
}

#this function takes a vcf, filters out biallelic sites, removes heterozygous sites and saves qs objects of the vcf and a numeric matrix of genotypes
# assuming homozygous diploids, 0 = 0/0 = homozygous ref, 1 = 1/1 = homozygous alt 
preprocessVCF=function(elegans.isotypes.vcf,elegans.isotypes.vcf.qs,elegans.isotypes.vcf.gt.qs) {
    vcf=vcfR::read.vcfR(elegans.isotypes.vcf)
    vcf=vcf[vcfR::is.biallelic(vcf),]

    gt=vcfR::extract.gt(vcf, as.numeric=T)
    #recode hets as NA
    gt[gt=='0|1']=NA
    gt[gt=='1|0']=NA
    gt[gt=='0|0']=0
    gt[gt=='1|1']=1
    #this conversion should force anything that isn't homozygous to NA
    gt2=matrix(as.numeric(gt),ncol=ncol(gt))
    rownames(gt2)=rownames(gt)
    colnames(gt2)=colnames(gt)
    gt=gt2
    rm(gt2)
    
    qsave(vcf, file=elegans.isotypes.vcf.qs)
    qsave(gt, file=elegans.isotypes.vcf.gt.qs)
    
    return(NULL) #list(vcf=vcf, gt=gt))
}







# take the larger vcf,  genotype calls, and a subset of parents 
# extracts segregating sites 
# return an alphaSimR founder population
createFounderPop=function(vcf, gt, p.names, X.only=F, X.drop=T) { 
    gt.sub=gt[,colnames(gt) %in% p.names]

    #monomorphic=apply(gt.sub, 1, function(x) all.equal(x))
    #monomorphic sites 
    #faster to do this with math
    rSg=rowSums(gt.sub)
    #sites with hets 
    sum(is.na(rSg))
    #sites all ref
    sum(rSg==0, na.rm=T)
    #sites all alt
    sum(rSg==length(p.names), na.rm=T)
    #sites mito
    sum(grepl('MtDNA', rownames(gt.sub)))

    bad.sites= is.na(rSg) | rSg==0 | rSg==length(p.names)  | grepl('MtDNA', rownames(gt.sub))
    if(X.only) {  bad.sites = bad.sites | !(grepl('X_', rownames(gt.sub))) }
    if(X.drop) {  bad.sites = bad.sites | (grepl('X_', rownames(gt.sub))) }
    gt.sub=gt.sub[-which(bad.sites),]
    vcf.cross=vcf[match(rownames(gt.sub), rownames(gt)), samples=colnames(gt.sub)]
    #generate sample ID
    vcf.cross=vcfR::addID(vcf.cross)


    uchrU=unique(getCHROM(vcf.cross))

    imputed.positions=jitterGmapVector(getGmapPositions(vcf.cross, gmap[uchrU], uchrU)) 
    

    #genetic map positions must be in Morgans
    genMap=data.frame(markerName=paste0(getCHROM(vcf.cross),'_',getPOS(vcf.cross)), chromosome=getCHROM(vcf.cross), position=unlist(imputed.positions)/100)

    teg.GT=t(gt.sub)
    #recode
    teg.GT[teg.GT==0]=-1
    colnames(teg.GT)=paste0(getCHROM(vcf.cross),'_',getPOS(vcf.cross))
    ped=data.frame(id=rownames(teg.GT), mother=rep(0, nrow(teg.GT)), father=rep(0,nrow(teg.GT)) ) #c(0,0), father=c(0,0))
    return(importInbredGeno(geno=teg.GT, genMap=genMap, ped=ped))
}



#calc alt-based AF in chunks
calcAltAf=function(FR, markers, split.factor=2e4) {
    
    sm=split(markers, ceiling(seq_along(markers)/split.factor))
    afr=list()
    for(i in 1:length(sm)) {
        print(i)
            x=sm[[i]]
          G= pullMarkerGeno(FR, x, asRaw=T)
          afr[[as.character(i)]]=Rfast::colsums(G)/(nrow(G)*2)
    }
    afr=do.call('c', afr)
    names(afr)=markers
    return(afr)
}


#===========generate ref and alt counts given a simulated phenotype, a genotype matrix, a selection strength, 
# sequencing depth, and whether we are selectign lower tail or upper tail =============================================
simSequencingRefAlt=function(y, FR, markers, depth, sel.frac, lower.tail=F) {
    
    #sel.frac should be a number between 0 and 1 
    #fraction of the tail you're selecting (1 = the whole distribution)
    #sel.frac=.1
    if(sel.frac==1) {
        #fraction of alt alleles across pop
        #alt.af=colSums(G)/(nrow(G)*2)
        alt.af=calcAltAf(FR,markers)
        #fraction of ref alleles across pop
        ref.af=1-alt.af
    
        sel.indv.af=alt.af
        nIndiv=nInd(FR)
    } else{
        if(lower.tail==F) {
            sel.indv=which(y>quantile(y,1-sel.frac))
        }
        else{
            sel.indv=which(y< quantile(y,sel.frac))
        }
       sel.indv.af=calcAltAf(FR[sel.indv],markers)
       #sel.indv.x=G[sel.indv,]
       #sel.indv.af=colSums(sel.indv.x)/(nrow(sel.indv.x)*2)
       nIndiv=nInd(FR[sel.indv])
    }

    #freq of alt
    a=rbinom(n=length(sel.indv.af),size=depth, prob=sel.indv.af)
    #freq of ref
    r=rbinom(n=length(sel.indv.af),size=depth, prob=1-sel.indv.af)
    countdf=data.frame(ID=markers, expected=1-sel.indv.af, ref=r, alt=a)
    attr(countdf, 'nInd')=nIndiv
    return(countdf)
}
#========================================================================================================================


#=======================================================================================================================
# take the larger vcf,  genotype calls, and a subset of parents 
# extracts segregating sites 
# return an alphaSimR founder population
createFounderPop=function(vcf, gt, p.names, gmap, X.only=F, X.drop=T) { 
    gt.sub=gt[,colnames(gt) %in% p.names]

    #monomorphic=apply(gt.sub, 1, function(x) all.equal(x))
    #monomorphic sites 
    #faster to do this with math
    rSg=rowSums(gt.sub)
    #sites with hets 
    sum(is.na(rSg))
    #sites all ref
    sum(rSg==0, na.rm=T)
    #sites all alt
    sum(rSg==length(p.names), na.rm=T)
    #sites mito
    sum(grepl('MtDNA', rownames(gt.sub)))

    bad.sites= is.na(rSg) | rSg==0 | rSg==length(p.names)  | grepl('MtDNA', rownames(gt.sub))
    if(X.only) {  bad.sites = bad.sites | !(grepl('X_', rownames(gt.sub))) }
    if(X.drop) {  bad.sites = bad.sites | (grepl('X_', rownames(gt.sub))) }
    gt.sub=gt.sub[-which(bad.sites),]
    vcf.cross=vcf[match(rownames(gt.sub), rownames(gt)), samples=colnames(gt.sub)]
    #generate sample ID
    vcf.cross=vcfR::addID(vcf.cross)


    uchrU=unique(getCHROM(vcf.cross))

    imputed.positions=jitterGmapVector(getGmapPositions(vcf.cross, gmap[uchrU], uchrU)) 
    

    #genetic map positions must be in Morgans
    genMap=data.frame(markerName=paste0(getCHROM(vcf.cross),'_',getPOS(vcf.cross)), 
                      chromosome=getCHROM(vcf.cross), 
                      position=unlist(imputed.positions)/100)

    teg.GT=t(gt.sub)
    #recode
    teg.GT[teg.GT==0]=-1
    colnames(teg.GT)=paste0(getCHROM(vcf.cross),'_',getPOS(vcf.cross))
    ped=data.frame(id=rownames(teg.GT), mother=rep(0, nrow(teg.GT)), father=rep(0,nrow(teg.GT)) ) #c(0,0), father=c(0,0))
    return(importInbredGeno(geno=teg.GT, genMap=genMap, ped=ped))
}




#=============simulate phenotypes given a genetic architecture on final population after crossings====================
simPheno=function(FR, genMapMarkers, QTL.sims,ds.size=NULL,returnG=T) {
  if(is.null(ds.size) | length(ds.size)>nInd(FR) ) {
      ds=seq(1,nInd(FR)) 
  }  else{
      ds=sort(sample.int(nInd(FR),ds.size))
  }

    #G=pullSegSiteGeno(FR[ds])
    #also possible that sites that aren't segregating are assigned QTL ??? check this  
    X_Q=pullMarkerGeno(FR[ds], QTL.sims$o.add.qtl.ind, asRaw=F)

    X_Beta=QTL.sims$o.add.qtl.eff
    if(length(X_Beta)==1) {
        XB=X_Q*X_Beta
    } else {XB=X_Q%*%X_Beta    }

    #two ways to 
    if(is.null(QTL.sims$o.h2.norm) | QTL.sims$o.h2.norm==F) {
        simy=XB+rnorm(nrow(G), mean=0, sd=QTL.sims$o.error.sd) 
        h2=var(XB)/(var(XB)+QTL.sims$o.error.sd^2)
    } else {
        h2=QTL.sims$o.h2
        #g=as.vector(scale(XB))
        #simy= g + rnorm(length(g), mean=0, sd=sqrt((1-h2)/(h2*(var(g))))) #h2*XB)))))

        g=as.vector(XB)
        gv=var(XB)
        
        # to derive expected error variance 
        # tv=gv+ev
        # gv/tv=h2
        # gv=h2*tv
        # gv/h2=tv
        ev=gv/h2-gv
        simy=g+rnorm(length(g), mean=0, sd=sqrt(ev))
        g=as.vector(XB)
    }
    
    print(paste( 'simulated total h^2:' , h2))

    if(returnG) {
      G=pullMarkerGeno(FR[ds],genMapMarkers,asRaw=F)
      af=(colSums(G)/(nrow(G)*2))
      plot(af, ylab='ref/(ref+alt)', xlab='marker.index')

      #fixed alt or fixed ref sites 
      f.ref=af==0
      f.alt=af==1
      G=G[,!(f.ref|f.alt)]

        return(list(G=G,
                    ind=ds,
                    X_Q=X_Q,
                    h2=h2,
                    simy=simy))
    } else{
        return(list(ind=ds,
                    h2=h2,
                    X_Q=X_Q,
                    simy=simy))
    }
}
#==================================================================================================================


phaseBiparental=function(df, p1.name, founderPop, genMap){

        gID=genMap$id
        p1.ref=pullMarkerGeno(founderPop, genMap$id)[p1.name,]==0
        p1.ref=p1.ref[gID]
        vname=names(p1.ref)

        p1=c(df$ref[p1.ref], df$alt[!p1.ref])
        vscramb=c(vname[p1.ref], vname[!p1.ref])
        names(p1)=vscramb
        p1=p1[vname]

        p2=c(df$ref[!p1.ref], df$alt[p1.ref])
        vscramb=c(vname[!p1.ref], vname[p1.ref])
        names(p2)=vscramb
        p2=p2[vname]
        if(!is.null(df$expected)) {
            expected.phased=ifelse(p1.ref, df$expected, 1-df$expected)
            df$expected.phased=expected.phased
         }
         df$p1=p1
         df$p2=p2
    return(df)
}

