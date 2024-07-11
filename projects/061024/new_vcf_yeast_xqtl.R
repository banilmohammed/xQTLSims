library(xQTLStats)
library(qs)
library(tidyverse)
library(vcfR)
library(AlphaSimR)
library(Rfast)
library(ggpubr)
library(susieR)

#pull from github
xQTLSims.dir = '/data0/xQTLSims/'

source.dir=paste0(xQTLSims.dir, 'R/')
#data.dir=paste0(xQTLSims.dir, 'data/')
#project.dir=paste0(xQTLSims.dir, 'projects/032924/')

#function to simulate crosses, treatement of X is incomplete/broken
#and additional helper functions
source(paste0(source.dir, 'simWormCrosses.R'))
source(paste0(source.dir, 'helperFxs.R'))
source(paste0(source.dir, 'makeCountTables.R'))


#unique chromosomes 
#yeast
uchr=paste0('chr', as.character(as.roman(1:16)))
#uchr=c(as.character(as.roman(1:5)), 'X') #paste0('chr', as.roman(1:16))

#elegans
#gmap.file=paste0(data.dir, 'geneticMapXQTLsnplist.rds')
#gmap=restructureGeneticMap(gmap.file)
gmaps=readRDS(system.file('reference', 'yeast_gmaps.RDS', package='xQTLStats'))
gmap=gmaps[['A']]

#------------------------------------------------------------------------------------------------------------------------------------------------------
#read new vcfs 
vcf.dir='/data1/yeast/reference/pop_vcfs/'

    #was used to retain the set of 1011 samples to keep, saved in pop_vcfs
    #samples.to.keep=c(1,which(nchar(colnames(vcf@gt))==3))
    #filter variants in the 1011 strain collection and only retain SNPs
    setwd(vcf.dir)
    #bcftools view -S samples.txt -v snps -m2 -M2 chrI.norm.vcf.gz |  bcftools filter -e 'AC==0 || AC==AN'  | bcftools view -O z -o test.vcf.gz
    for(u in uchr) {
        print(u)
        system(paste0('bcftools view -S samples.txt -v snps -m2 -M2 ', u, ".norm.vcf.gz |  bcftools filter -e 'AC==0 || AC==AN'  | bcftools view -O z -o", ' 1011/', u, '.vcf.gz'))
    }


    #now read in each vcf and save as R objects for faster access
    #also add BY (S288C) given it is missing in the collection
    #vcfS=list()
    for( u in uchr) {
    #test.vcf.gz
    #}
        print(u)
        yeast.ref.vcf=paste0(vcf.dir, '1011/', u, '.vcf.gz') #'/data1/yeast/reference/pop_vcfs/test.vcf.gz') #', u, '.norm.vcf.gz')
        vcf=vcfR::read.vcfR(yeast.ref.vcf) #, cols=samples.to.keep)
        gt=vcfR::extract.gt(vcf, as.numeric=T)

        print('# of variants')
        print(nrow(gt))

        gcnt=rowSums(!is.na(gt))
        acnt=rowSums(gt, na.rm=T)
        af=acnt/gcnt
        monomorphic= (af==0 | af==1) # ((acnt/gcnt)==0)
        
        vcf=vcf[!monomorphic,]
        gt=gt[!monomorphic,]
        print(nrow(gt))

        gcnt=gcnt[-which(monomorphic)] #rowSums(!is.na(gt))
        acnt=acnt[-which(monomorphic)] #rowSums(gt, na.rm=T)
        af=acnt/gcnt

        #simplest to just prune multi-allelic sites
        pos=getPOS(vcf)
        multiallelic.pos=unique(pos[which(duplicated(pos))])
        biallelic=!(pos %in% multiallelic.pos)
        
        vcf=vcf[biallelic,]
        gt=gt[biallelic,]

        print(nrow(gt))

        ref.hack="0/0:0,0:0:0:.:.:0,0,0:." 
        S288C=rep(ref.hack, nrow(gt))
        vcf@gt=cbind(vcf@gt, S288C)
        S288C=rep(0, nrow(gt))
        gt=cbind(gt, S288C)

        qsave(list(gt=gt, vcf=vcf), file=paste0(vcf.dir, '1011/', u, '.gs'))
        #vcfS[[u]]=vcf
    }


    #manually filter and retain the freaking snps 
    v1=qread('/data1/yeast/reference/pop_vcfs/1011/chrI.gs')
    vcf=v1$vcf
    #manually norm the supid vcf so GATK doesn't choke 
    vcf@fix[,'REF']= substr(vcf@fix[,'REF'],1,1)
    vcf@fix[,'ALT']= substr(vcf@fix[,'ALT'],1,1)
    #snpL = nchar(vcf@fix[,'REF'])==1  & nchar(vcf@fix[,'ALT'])==1 
    gt=v1$gt
    #vcf=vcf #[snpL,]
    #gt=gt #[snpL,]
    for(u in uchr[-1]) { 
        print(u)
        vtemp=qread(paste0('/data1/yeast/reference/pop_vcfs/1011/', u, '.gs'))
        vtemp$vcf@fix[,'REF']= substr(vtemp$vcf@fix[,'REF'],1,1)
        vtemp$vcf@fix[,'ALT']= substr(vtemp$vcf@fix[,'ALT'],1,1)
        #snpL = nchar(vtemp$vcf@fix[,'REF'])==1  & nchar(vtemp$vcf@fix[,'ALT'])==1 
        vcf@fix=rbind(vcf@fix, vtemp$vcf@fix) #[snpL,])
        vcf@gt=rbind(vcf@gt, vtemp$vcf@gt) #[snpL,])
        gt=rbind(gt, vtemp$gt) #[snpL,])
    }

    # save vcf and gt for future processing and simulations 
    qsave(list(gt=gt, vcf=vcf), file=paste0(vcf.dir, '1011/mega_filtered.qs'))


    #hack a vcf with het calls at each variant site for use with gatk ASEReadCounter
    vcfout=vcf
    vcfout@gt=vcfout@gt[,c(1,2)]
    vcfout@gt[,1]='GT:DP'
    vcfout@gt[,2]='0/1:.'
    write.vcf(vcfout,'/data1/yeast/reference/pop_vcfs/1011/mega_filtered.vcf')
    #then on the command line, bgzip and gatk IndexFeatureFile


    # figure out reticulate behavior
    #library(reticulate)
    #use_condaenv('genomics')
    #maybe switch to absolute python path 
    #system(paste0('
    bcftools view -I /data1/yeast/reference/pop_vcfs/1011/mega_filtered.vcf -O z -o /data1/yeast/reference/pop_vcfs/1011/mega_filtered.vcf.gz 
    #'), intern=T)

    #bcftools norm -a -m -snps -f /media/hoffman2/lcrisman/Yeast_ref/sacCer3.fasta mega_filtered.vcf.gz > mega_filtered_snps.vcf
    #use_condaenv('genomics')
    #system(paste0("
    gatk IndexFeatureFile -I  /data1/yeast/reference/pop_vcfs/1011/mega_filtered.vcf.gz
    #"), intern=T)

    #------------------------------------------------------------------------------
    #bcftools index mega_filtered.vcf.gz 
    #tabix mega_filtered.vcf.gz
    gatk ASEReadCounter --min-mapping-quality 20 -I /media/hoffman2/lcrisman/Dmagicmarker_042024/bam_files/G1_KANa_S62.bam -V /data1/yeast/reference/pop_vcfs/1011/mega_filtered.vcf.gz -R /media/hoffman2/lcrisman/Yeast_ref/sacCer3.fasta -O /data0/xQTLSims/projects/061024/data/G1_KANa_S62.txt
    gatk ASEReadCounter --min-mapping-quality 20 -I /media/hoffman2/lcrisman/Dmagicmarker_042024/bam_files/G2_NATalpha_S63.bam -V /data1/yeast/reference/pop_vcfs/1011/mega_filtered.vcf.gz -R /media/hoffman2/lcrisman/Yeast_ref/sacCer3.fasta -O /data0/xQTLSims/projects/061024/data/G2_NATalpha_S63.txt



    #manually removed structural variants, and sorted the thing with bcftools sort 
    #yeast.ref.vcf='/media/hoffman2/jsbloom/reference/rr_parents_no_svar.vcf.gz'
    #vcf=vcfR::read.vcfR(yeast.ref.vcf)
    #vcf=vcf[vcfR::is.biallelic(vcf),]
    ##vcf=vcf[-which(duplicated(paste0(getCHROM(vcf),':', getPOS(vcf)))),]
    #vcf=vcf[!getCHROM(vcf)=='chrM',]

    #gt=vcfR::extract.gt(vcf, as.numeric=T)

#-----------------------------------------------------------------------------------------------------------------------

mega_filtered=qread(paste0(vcf.dir, '1011/mega_filtered.qs'))
vcf=mega_filtered$vcf
gt=mega_filtered$gt
rm(mega_filtered)


#sample.key=as_tibble(data.frame(s=c('G1_KANa_S62', 'G2_NATalpha_S63'), p1=rep('BYa',2), p2=rep('CBS2888a',2)))
sample.key=as_tibble(data.frame(s=c('G1_KANa_S62', 'G2_NATalpha_S63'), p1=rep('S288C',2), p2=rep('ABL',2)))
names(sample.key)=c('sample name', 'parent 1', 'parent 2')
sample.dir='/data0/xQTLSims/projects/061024/data/' #~/Desktop/' #/media/hoffman2/lcrisman/Dmagicmarker_042024/count_files/'


countdfs=makeCountTables(sample.key,sample.dir, vcf,gt, sample.suffix='.txt')



#should parallelize this .... yawn 
afds=lapply(names(countdfs), function(snn) {
       calcAFD(countdfs[[snn]], experiment.name=snn,sample.size=1e4, sel.strength=.9) #, bin.width=3000, eff.length=2000, uchr=uchr)
      })
names(afds)=names(countdfs)


plots=lapply(names(afds), function(snn) {
    plotIndividualExperiment(afds[[snn]], snn) 
      })
a=ggarrange(plots[[1]],plots[[2]],  nrow=2)

plot(a)

set.seed(1)
p.names=c('S288C', sample(colnames(gt),96))
#p.names=c('S288C', 'AAA')
mating.matrix=matrix(c(1,2), nrow=1)
# in this case, given one-way cross, potentially allow N2 herm to self 
selfings=1
founderPop = createFounderPop(vcf,gt, p.names, gmap, X.drop=F) #c('N2', 'XZ1516'))
SP=SimParam$new(founderPop)
FN=newPop(founderPop, simParam=SP)
#there's a structure here for keeping track of males that we are only partially leveraging
#FN@sex[2]='M'

genMap=getGenMap(founderPop)

#possible QTL architectures for fitness effects ---------------------
nmarker=nrow(genMap)
#for example a TA element
f.nQTL=1
#sample(genMap$id, nQTL)
f.nadditive=f.nQTL
f.add.qtl.ind  = genMap$id[sort(sample(nmarker, f.nadditive))]
f.add.qtl.eff =  sample(ifelse(runif(f.nadditive)>.5,-1,1),replace=T)
#---------------------------------------------------------------------

#possible QTL architectures for traits orthogonal to fitness 
o.nQTL=200
o.nadditive=o.nQTL
o.add.qtl.ind  = genMap$id[sort(sample(nmarker, o.nadditive))]
o.add.qtl.eff =  sample(ifelse(runif(o.nadditive)>.5,-1,1),replace=T)


#f designates QTL effects for hermaphrodites 
QTL.sims=list(#simulate fitness effects during panel construction ------------------------
              sim.fitness=F,
              #positions of fitness effect QTL
              f.add.qtl.ind=f.add.qtl.ind,
              #fitness QTL effect magnitudes
              f.add.qtl.eff=f.add.qtl.eff,
              #residual error sd of fitness effect
              f.error.sd=10,
              #---------------------------------------------------------------------------
              #simulate additional trait effects that don't manifest as fitness effects
              sim.orthogonal=T,
              # positions of trait effect QTL
              o.add.qtl.ind=o.add.qtl.ind,
              # trait effect QTL effect magnitudes
              o.add.qtl.eff=o.add.qtl.eff,
              # residual error sd of trait effect 
              o.error.sd=2,
              # toggle to ignore residual error sd of trait effect and instead normalize
              # trait variance such that residual error is 1-h^2
              o.h2.norm=F,
              o.h2=.5)
qsave(QTL.sims, file='/data0/xQTLSims/projects/yeast_sims/metaQTL.qs')
QTL.sims=qread('/data0/xQTLSims/projects/yeast_sims/metaQTL.qs')



max.per.gen=1e5
meta.results=list()
for(i in 1:96) {
    #setup sims using mega pop then mini sims per biparental
    FB=createFounderPop(vcf,gt,p.names[c(1,i+1)],gmap, X.drop=F)
    SP=SimParam$new(FB)
    genMap=getGenMap(FB)
    FB=newPop(FB, simParam=SP)
    f1=makeCross(FB, matrix(rep(c(1,2), each=1), ncol=2) , nProgeny=1, simParam=SP)
    f2=makeDH(f1, nDH=max.per.gen, simParam=SP) #matrix(rep(c(1,2), each=1), ncol=2) ,
    FR=f2
    g=pullMarkerGeno(FR[1:1000], genMap$id[genMap$chr=='chrV'])
    sg=scale(g)
    LD= crossprod(sg)/(nrow(sg)-1)

    simFR=simPheno(FR, genMapMarkers=genMap$id, QTL.sims=QTL.sims, returnG=F)
    countdf.h=simSequencingRefAlt(simFR$simy,FR, genMap$id, depth=50, sel.frac=.1, lower.tail=F)
    #if you set sel.frac=1 sample the population of existing genotypes without QTL effects 
    countdf.l=simSequencingRefAlt(y=NULL, FR,genMap$id, depth=50, sel.frac=1 , lower.tail=F)

    countdf.h=phaseBiparental(countdf.h, p.names[1], FB, genMap)
    countdf.l=phaseBiparental(countdf.l, p.names[1], FB, genMap)
  #  plot(countdf.h$p1/(countdf.h$p1+countdf.h$p2))
  #  points(countdf.h$expected, col='red') #alt/(countdf$al

    test  = calcAFD(countdf.h, experiment.name='high1',sample.size=1e4, sel.strength=.1, uchr=unique(genMap$chr) ) 
    test2 = calcAFD(countdf.l, experiment.name='unsel1',sample.size=1e4, sel.strength=1, uchr=unique(genMap$chr) )
    results=calcContrastStats(results=list(test, test2), L='_high1', R='_unsel1')
    
    plotSummary(results)

    meta.results[[p.names[i+1]]]=results
}
qsave(meta.results, file='/data0/xQTLSims/projects/yeast_sims/metaResults.qs')
meta.results=qread('/data0/xQTLSims/projects/yeast_sims/metaResults.qs')


bhat=meta.results[[1]]$contrast.beta[meta.results[[1]]$chrom=='chrV'] #[genMap$chr=='chrV']
shat=meta.results[[1]]$contrast.beta.se[meta.results[[1]]$chrom=='chrV'] # wm.se[genMap$chr=='chrV']
trunc.ind=which(is.na(bhat))
bhat=bhat[-trunc.ind]
shat=shat[-trunc.ind]
ids= meta.results[[1]]$ID[meta.results[[1]]$chrom=='chrV'][-trunc.ind]
plot(bhat/shat)
abline(v=match(QTL.sims$o.add.qtl.ind,ids), col='red')
test=susie_rss(bhat=bhat,shat=shat, n=1e5, R=LD[-trunc.ind,-trunc.ind], verbose=T)
x11()
plot(test$pip)
abline(v=match(QTL.sims$o.add.qtl.ind,ids), col='red')







#then basically full pop
founderPop = createFounderPop(vcf,gt, p.names, gmap, X.drop=F) #c('N2', 'XZ1516'))
SP=SimParam$new(founderPop)
FN=newPop(founderPop, simParam=SP)
#there's a structure here for keeping track of males that we are only partially leveraging
#FN@sex[2]='M'
genMap=getGenMap(founderPop)
mini=makeCross(FN, cbind(1, seq(1:96)), nProgeny=1, simParam=SP)
#100 progeny biparental cross to BY
mini.sample=makeDH(mini, nDH=100, simParam=SP)
#then one gen intercross 
max.per.gen=9600
m1=(sample.int(max.per.gen, round(max.per.gen/2)))
m2=which(!(seq(1,max.per.gen) %in% m1))
mmat=cbind(m1,m2)
f3=makeCross(mini.sample, mmat, nProgeny=1,simParam=SP)
#end pop about 50k
f4=makeDH(f3, nDH=10, simParam=SP) #matrix(rep(c(1,2), each=1), ncol=2) , nProgeny=1, simParam=SP)

max.per.gen=48000
m1=(sample.int(max.per.gen, round(max.per.gen/2)))
m2=which(!(seq(1,max.per.gen) %in% m1))
mmat=cbind(m1,m2)
f5=makeCross(f4, mmat, nProgeny=1,simParam=SP)
#end pop about 50k
f6=makeDH(f5, nDH=2, simParam=SP) #matrix(rep(c(1,2), each=1), ncol=2) , nP

max.per.gen=48000
m1=(sample.int(max.per.gen, round(max.per.gen/2)))
m2=which(!(seq(1,max.per.gen) %in% m1))
mmat=cbind(m1,m2)
f7=makeCross(f6, mmat, nProgeny=1,simParam=SP)
#end pop about 50k
f8=makeDH(f7, nDH=2, simParam=SP) #matrix(rep(c(1,2),

#ds=sample.int(48000,1e4)
#[ds]
g=pullMarkerGeno(f8, genMap$id[genMap$chr=='chrVI'])
sg=standardise(g)
LD= crossprod(sg)/(nrow(sg)-1)
simFR=simPheno(f8, genMapMarkers=genMap$id, QTL.sims=QTL.sims, returnG=F)

lm(simFR$simy[,1]~g[,na.omit(match(QTL.sims$o.add.qtl.ind, colnames(g)))])


test2=susie(sg, simFR$simy[,1], verbose=T)
plot(test2$pip)
abline(v=match(QTL.sims$o.add.qtl.ind, colnames(g))) #names(Zsub)))
test2$pip[na.omit(match(QTL.sims$o.add.qtl.ind, colnames(g)))]

countdf.h=simSequencingRefAlt(simFR$simy,f8, genMap$id, depth=2500, sel.frac=.1, lower.tail=F)
#if you set sel.frac=1 sample the population of existing genotypes without QTL effects 
countdf.l=simSequencingRefAlt(y=NULL, f8,genMap$id, depth=2500, sel.frac=1 , lower.tail=F)

#sp flow cell
#25000x
#s1 flow cell
#50000x

#g=pullMarkerGeno(mini.sample, genMap$id[genMap$chr=='chrV'])
#sg=scale(g)
#LD= crossprod(sg)/(nrow(sg)-1)
#LD=cora(g, large=T)
# and then fill in summary stats when variants are segregating



acnt.h=countdf.h$alt
rcnt.h=countdf.h$ref #(nrow(G.h)*2)-acnt.h

acnt.c=countdf.l$alt #colSums(G.c)
rcnt.c=countdf.l$ref #(nrow(G.c)*2)-acnt.c

#145
#library(RVAideMemoire)
#i=171974

chisq=chisq.p=rep(NA, length(acnt.h))
for(i in 1:length(acnt.h)){
    if(i%%10000==0) { print(i)} 
    # high tail alt, unselected alt
    # high tail ref, unselected ref 
    tm=cbind(c(acnt.h[i],rcnt.h[i]), 
             c(acnt.c[i],rcnt.c[i]))
    ctm=chisq.test(tm)
    chisq[i]=as.numeric(ctm$statistic)
    chisq.p[i]=as.numeric(ctm$p.value)
}
sgn=(((acnt.h/rcnt.h)/(acnt.c/rcnt.c)>1)*2)-1
plot(sgn*chisq)
abline(v=match(QTL.sims$o.add.qtl.ind, genMap$id), col=(QTL.sims$o.add.qtl.eff)+2)

Z=sgn*sqrt(chisq)
names(Z)=genMap$id
Zsub=Z[genMap$chr=='chrVI']
plot(Zsub)
abline(v=match(QTL.sims$o.add.qtl.ind, names(Zsub)))

zr=regress(Zsub~1,~LD,verbose=T)
lambda=estimate_s_rss(Zsub, LD, n=5e4)
#.721

condz_in=kriging_rss(Zsub, LD, n=5e4, s=lambda)
zcond=condz_in$conditional_dist$condmean/sqrt(condz_in$conditional_dist$condvar)

plot(Zsub, zcond)

#test=susie_rss(z=zcond, R=LD,5e3, L=10, verbose=T)
test3=susie_rss(z=zcond, R=LD,5e4, L=10, verbose=T)
test=susie_rss(z=Zsub, R=LD,5e4, L=10, verbose=T)
plot(test3$pip)
abline(v=match(QTL.sims$o.add.qtl.ind, names(Zsub)))









#then finemap the summary stats 
meta.b=meta.se=matrix(NA, nrow=nrow(genMap), ncol=length(meta.results))
rownames(meta.b)=rownames(meta.se)=genMap$id
colnames(meta.b)=colnames(meta.se)=names(meta.results)
for(n in names(meta.results)){
    #n=names(meta.results)[1]
    mr=meta.results[[n]]
    rmatch=match(mr$ID, rownames(meta.b))
    mrc=mr$contrast.beta[which(!is.na(rmatch))]
    mre=mr$contrast.beta.se[which(!is.na(rmatch))]
    meta.b[rmatch[!is.na(rmatch)],match(n,colnames(meta.b))]=mrc #mr$contrast.beta
    meta.se[rmatch[!is.na(rmatch)],match(n,colnames(meta.se))]=mre #mr$contrast.beta

}

w=1/(meta.se^2)
sw=apply(w,1,sum,na.rm=T)
wm=apply(meta.b*w, 1, sum, na.rm=T)/sw
wm.se=sqrt(1/sw)
plot(wm/wm.se)
abline(v=match(QTL.sims$o.add.qtl.ind, rownames(meta.b)), col='red')


#g=pullMarkerGeno(FN, genMap$id[genMap$chr=='chrV'])
z=wm/wm.se
#LD=cor(g)
#z=z[genMap$chr=='chrV']

bhat=wm[genMap$chr=='chrV']
bhat[is.na(bhat)]=0
shat=wm.se[genMap$chr=='chrV']
test=susie_rss(bhat=bhat[1:5000],shat=shat[1:5000], n=1e5, R=LD[1:5000,1:5000], verbose=T)

#abline(v
z=bhat/shat
plot(z)
abline(v=na.omit(match(QTL.sims$o.add.qtl.ind, names(z))))
test=susie_rss(bhat=bhat,shat=shat, n=1e6, R=LD, L=20, verbose=T, max_iter=500)

x11()
plot(test$pip>0)

#new simulation engine for yeast
max.per.gen=1e5
m1=(sample.int(max.per.gen, round(max.per.gen/2)))
m2=which(!(seq(1,max.per.gen) %in% m1))
mmat=cbind(m1,m2)

f1=makeCross(FN, matrix(rep(c(1,2), each=1), ncol=2) , nProgeny=1, simParam=SP)
f2=makeDH(f1, nDH=max.per.gen, simParam=SP) #matrix(rep(c(1,2), each=1), ncol=2) , nProgeny=1, simParam=SP)

m1=(sample.int(max.per.gen, round(max.per.gen/2)))
m2=which(!(seq(1,max.per.gen) %in% m1))
mmat=cbind(m1,m2)

f3=makeCross(f2, mmat, nProgeny=1,simParam=SP)
f4=makeDH(f3, nDH=2, simParam=SP) #matrix(rep(c(1,2), each=1), ncol=2) , nProgeny=1, simParam=SP)
#x=pullMarkerGeno(f4, genMap$id)

m1=(sample.int(max.per.gen, round(max.per.gen/2)))
m2=which(!(seq(1,max.per.gen) %in% m1))
mmat=cbind(m1,m2)

f5=makeCross(f4, mmat, nProgeny=1,simParam=SP)
f6=makeDH(f5, nDH=2, simParam=SP) #matrix(rep(c(1,2), each=1), ncol=2) , nProgeny=1, simParam=SP)
x=pullMarkerGeno(f6, genMap$id)

f1=makeCross(FN, matrix(rep(c(1,2), each=1), ncol=2) , nProgeny=1, simParam=SP)
f2=makeDH(f1, nDH=max.per.gen, simParam=SP) 

FR=f2
#see rewrite below
simFR=simPheno(FR, genMapMarkers=genMap$id, QTL.sims=QTL.sims, returnG=F)
countdf.h=simSequencingRefAlt(simFR$simy,FR, genMap$id, depth=50, sel.frac=.1, lower.tail=F)
#if you set sel.frac=1 sample the population of existing genotypes without QTL effects 
countdf.l=simSequencingRefAlt(y=NULL, FR,genMap$id, depth=50, sel.frac=1 , lower.tail=F)





#sanity checks 
plot(countdf.h$alt/(countdf.h$alt+countdf.h$ref))
#points(countdf.h$expected, col='red') #alt/(countdf$alt+countdf$ref))
countdf.h=phaseBiparental(countdf.h, p.names[1], founderPop, genMap)
countdf.l=phaseBiparental(countdf.l, p.names[1], founderPop, genMap)
plot(countdf.h$p1/(countdf.h$p1+countdf.h$p2))
points(countdf.h$expected, col='red') #alt/(countdf$al

test  = calcAFD(countdf.h, experiment.name='high1',sample.size=1e4, sel.strength=.1, uchr=unique(genMap$chr) ) 
test2 = calcAFD(countdf.l, experiment.name='unsel1',sample.size=1e4, sel.strength=1, uchr=unique(genMap$chr) )
results=calcContrastStats(results=list(test, test2), L='_high1', R='_unsel1')

# test$expected.phased_high1=1-test$expected.phased_high1

h1=plotIndividualExperiment(test, 'high1')
u1=plotIndividualExperiment(test2, 'unsel1')
c1=plotContrast(results, 'high1', 'unsel1')
ggpubr::ggarrange(h1, u1, c1, nrow=3) 


simulatedQTL=data.frame(data.table::tstrsplit(QTL.sims$o.add.qtl.ind,'_'), stringsAsFactors=F)
names(simulatedQTL)=c('chrom', 'physical.position')
simulatedQTL[,2]=as.numeric(simulatedQTL[,2])
plotSummary(results,simulatedQTL=simulatedQTL)













acnt.h=countdf.h$alt
rcnt.h=countdf.h$ref #(nrow(G.h)*2)-acnt.h

acnt.c=countdf.l$alt #colSums(G.c)
rcnt.c=countdf.l$ref #(nrow(G.c)*2)-acnt.c

#145
#library(RVAideMemoire)
#i=171974

chisq=chisq.p=rep(NA, length(acnt.h))
for(i in 1:length(acnt.h)){
    if(i%%10000==0) { print(i)} 
    # high tail alt, unselected alt
    # high tail ref, unselected ref 
    tm=cbind(c(acnt.h[i],rcnt.h[i]), 
             c(acnt.c[i],rcnt.c[i]))
    ctm=chisq.test(tm)
    chisq[i]=as.numeric(ctm$statistic)
    chisq.p[i]=as.numeric(ctm$p.value)
}
sgn=(((acnt.h/rcnt.h)/(acnt.c/rcnt.c)>1)*2)-1
plot(sgn*chisq)
abline(v=match(QTL.sims$o.add.qtl.ind, genMap$id), col=(QTL.sims$o.add.qtl.eff)+2)

n=sum(ctm$observed)
sr <-rowSums(ctm$observed) # rowSums(x)
sc <- colSums(ctm$observed)
E <- outer(sr, sc)/n
#equivalent to 
tcrossprod(sr,sc)




x1=pullMarkerGeno(test, genMap$id)

x=pullMarkerGeno(testd, genMap$id)


pedigreeCross(FN, id=c(1,2), c(1,1,1,1,1), c(2,2,2,2,2), matchID=T, simParam=SP)

SimWormParams=list(

    #how many individual crossed worms per row in mating.matrix
    starting.sample.size=1e4 ,
    #brood size if it selfs
    brood.size.selfing=20 ,
    #brood size for mating 
    brood.size.mating=20*2 ,
    #fraction of hermaphrodites that self
    selfing.rate=.90 ,
    #bottleneck per generation   
    max.per.gen=2e4,
    #how many total generations
    max.gen=12,
    #Founder population genotypes
    FN=FN ,
    #simulation parameters 
    SP=SP ,
    mating.matrix=mating.matrix ,
    selfings=selfings,
    genMap=genMap,
    QTL.sims=QTL.sims
)









##e.g.  to visualize just call
##plots[[1]]
#
##dump individual plots somewhere ---------------
##plot.dir= '/home/jbloom/Dropbox/code/xQTL/elegans/plots/'
##for(snn in names(plots)) {
##    ggsave(paste0( plot.dir, snn, '.png'), plots[[snn]], width=16)
##}
#
##I have to restructure this stuff  ---------------------------------------------
#S2_N2_XZ1516_7_RNAi=afds[[2]]
#LWM1_S56_N2_XZ1516_3_RNAi=afds[[15]]
#results=calcContrastStats(results=list(S2_N2_XZ1516_7_RNAi, LWM1_S56_N2_XZ1516_3_RNAi),
#                          L=paste0('_', names(afds)[2]), R=paste0('_', names(afds)[15]) ) #'_high1', R='_unsel1')
#
#sc=plotContrast(results, 'S2_N2_XZ1516_7_RNAi', 'LWM1_S56_N2_XZ1516_3_RNAi') #S3_N2_XZ1516_7_KO')
#
#s=plotSummary(results, effective.n.tests=2000)
#
#a=ggarrange(plots[[2]],plots[[15]], s,  nrow=3)
#a
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#    str(g2
#    gt2=gt
#    gt=cbind(gt, 0)
#    gt=vcfR::extract.gt(vcf, as.numeric=T)
#
#    df=data.frame(ind=1:nrow(gt), pos=getPOS(vcf), af=af)
#
#
#
#    monomorphic= rowSums(gt, na.rm=T)==0 |  rowSums(gt, na.rm=T)==1
#
#    vcf=vcf[!monomorphic,]
#    
#
#for( u in uchr) {
#    yeast.ref.vcf=paste0('/data1/yeast/reference/pop_vcfs/', u, '.norm.vcf.gz')
#    vcf=vcfR::read.vcfR(yeast.ref.vcf, cols=samples.to.keep)
#
#    c(1,which(nchar(colnames(vcf@gt))==3))
#    gt=vcfR::extract.gt(vcf, as.numeric=T)
#   
#   
#   
#    vcf=vcf[vcfR::is.biallelic(vcf),]
#
#    #can we hack together a humongous vcf from the consituent chromosomes
#test=vcf
#R> test@fix=rbind(test@fix, vcf2@fix)
#R> test@gt=rbind(test@gt, vcf2@gt)
#R> test
#   
#   
#   
#    yeast.ref.vcf=paste0('/data1/yeast/reference/pop_vcfs/', uchr[2], '.norm.vcf.gz')
#    vcf2=vcfR::read.vcfR(yeast.ref.vcf)
#
#}
#
#'ABL' 
#'S288C'
#
#
##yeast.ref.vcf=system.file('reference', 'parents_w_svar_sorted.vcf.gz', package='xQTLStats')
##manually removed structural variants, and sorted the thing with bcftools sort 
#yeast.ref.vcf='/media/hoffman2/jsbloom/reference/rr_parents_no_svar.vcf.gz'
#vcf=vcfR::read.vcfR(yeast.ref.vcf)
#vcf=vcf[vcfR::is.biallelic(vcf),]
##vcf=vcf[-which(duplicated(paste0(getCHROM(vcf),':', getPOS(vcf)))),]
#vcf=vcf[!getCHROM(vcf)=='chrM',]
#
#gt=vcfR::extract.gt(vcf, as.numeric=T)
#
#
#
#sample.key=as_tibble(data.frame(s=c('G1_KANa_S62', 'G2_NATalpha_S63'), p1=rep('BYa',2), p2=rep('CBS2888a',2)))
#names(sample.key)=c('sample name', 'parent 1', 'parent 2')
#sample.dir='/media/hoffman2/lcrisman/Dmagicmarker_042024/count_files/'
#
#countdfs=makeCountTables(sample.key,sample.dir, vcf,gt, sample.suffix='.txt')
#
# saveRDS(countdfs, file='/media/hoffman2/jsbloom/xQTL/yeast/050124/countdfs.RDS')
#
#
#
#                      sample.suffix='.txt'
#
##recode hets as NA
##gt[gt=='0|1']=NA
##gt[gt=='1|0']=NA
##gt[gt=='0|0']=0
##gt[gt=='1|1']=1
##this conversion should force anything that isn't homozygous to NA
##gt2=matrix(as.numeric(gt),ncol=ncol(gt))
##rownames(gt2)=rownames(gt)
##colnames(gt2)=colnames(gt)
##gt=gt2
##rm(gt2)
#
#    founderPop = createFounderPop(vcf,gt, p.names, gmap, X.drop=F)
#
#    p1='BYa'
#    p2='CBS2888a'
#    p.names=c(p1,p2)
##test=createFounderPop(vcf, gt, p.names, gmap, X.drop=F)
#
#    gt.sub=gt[,colnames(gt) %in% p.names]
#
#    #monomorphic=apply(gt.sub, 1, function(x) all.equal(x))
#    #monomorphic sites 
#    #faster to do this with math
#    rSg=rowSums(gt.sub)
#    #sites with hets 
#    sum(is.na(rSg))
#    #sites all ref
#    sum(rSg==0, na.rm=T)
#    #sites all alt
#    sum(rSg==length(p.names), na.rm=T)
#    #sites mito
#    sum(grepl('MtDNA', rownames(gt.sub)))
#
#    bad.sites= is.na(rSg) | rSg==0 | rSg==length(p.names)  | grepl('MtDNA', rownames(gt.sub))
#    gt.sub=gt.sub[-which(bad.sites),]
# 
#    vcf.cross=vcf[match(rownames(gt.sub), rownames(gt)), samples=colnames(gt.sub)]
#    #generate sample ID
#    vcf.cross=vcfR::addID(vcf.cross)
#
#    uchrU=unique(getCHROM(vcf.cross))
#
#    #get physical position, split by chromosome
#    p.by.chr=split(vcfR::getPOS(vcf.cross),vcfR::getCHROM(vcf.cross))
#    #keep things sorted (yay yeast chr names with roman numerals)
#    p.by.chr=p.by.chr[uchr]
#
#    imputed.positions=list()
#    for(chr in uchr) {
#        print(chr)
#        y=gmap[[chr]]
#        x=p.by.chr[[chr]]
#        imputed.positions[[chr]]=approxfun(y$ppos, y$map, rule=2)(x)
#    }
#     #   plot(x, approxfun(y$ppos, y$map, rule=2)(x),main=chr)
#     #   readline()
#    #}
#
#
#
#   #where to put the variant sites, impute onto gmap
#    #genetic map positions must be in Morgans
#    genMap=data.frame(markerName=paste0(getCHROM(vcf.cross),'_',getPOS(vcf.cross)), 
#                      chromosome=getCHROM(vcf.cross), 
#                      position=unlist(imputed.positions)/100)
#
#    teg.GT=t(gt.sub)
#    #recode
#    teg.GT[teg.GT==0]=-1
#    colnames(teg.GT)=paste0(getCHROM(vcf.cross),'_',getPOS(vcf.cross))
#    ped=data.frame(id=rownames(teg.GT), mother=rep(0, nrow(teg.GT)), father=rep(0,nrow(teg.GT)) ) #c(0,0), father=c(0,0))
#    return(importInbredGeno(geno=teg.GT, genMap=genMap, ped=ped))
##}
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#imputed.gmap=getGmapPositions(vcf,gmap,uchr)
#
##qsave(vcf, file=elegans.isotypes.vcf.qs)
##qsave(gt, file=elegans.isotypes.vcf.gt.qs)
#
##pretty intensive memory usage here
##include web link to vcf file 
##elegans.isotypes.vcf=paste0(data.dir,'WI.20220216.impute.isotype.vcf.gz')
#
##filtered vcf as qsave objects
## !!!! find the premade objects here folks !!!! :
## /u/project/kruglyak/jsbloom/xQTL/elegans/ref/
##and place in your data.dir
##elegans.isotypes.vcf.qs=paste0(data.dir,'WI.20220216.vcf.qs')
#
##filtered numeric gt calls from vcf  as qsave object
##elegans.isotypes.vcf.gt.qs=paste0(data.dir,'WI.20220216.vcf.GT.qs')
#
##run once, Laura skip this ============================================================
##use premade objects to save yourself the memory related headache, but this is how the objects are made
##preprocessVCF(elegans.isotypes.vcf,elegans.isotypes.vcf.qs,elegans.isotypes.vcf.gt.qs)
##======================================================================================
#
##Laura start here ======================================================================
##vcf=qread(elegans.isotypes.vcf.qs)
##gt=qread(elegans.isotypes.vcf.gt.qs)
#
#sample.key.file=paste0(project.dir, 'samplekey_032924.txt')
#sample.key=read_tsv(sample.key.file)
#sample.dir='/media/hoffman2/thatguy0/nextflow_pipes/bulkGWA/BulkGWA-20240405/undup_alignments/'
##only have tables for a subset of the samples
#sample.key=sample.key[c(1:12,21,22),]
#countdfs=makeCountTables(sample.key,sample.dir, vcf,gt)
#
#sample.key.file=paste0(project.dir, 'samplekey_earlier.txt')
#sample.key=read_tsv(sample.key.file)
#sample.key$'sample name'=paste0('LWM', 1:6, '_S', 56:61)
#sample.dir='/media/hoffman2/lwalterm/032924/count/' #media/hoffman2/thatguy0/nextflow_pipes/bulkGWA/BulkGWA-20240405/undup_alignments/'
#countdfs2=makeCountTables(sample.key,sample.dir, vcf,gt, sample.suffix='.txt')
#
##all the filtered and phased count tables, glued together
#countdfs=c(countdfs,countdfs2)
#rm(countdfs2)
#
##calculate the allele frequencies and SEs, params here need some more dialing in 
##should parallelize this .... yawn 
#afds=lapply(names(countdfs), function(snn) {
#       calcAFD(countdfs[[snn]], experiment.name=snn,sample.size=1e4, sel.strength=.95, bin.width=3000, eff.length=2000, uchr=uchr)
#      })
#names(afds)=names(countdfs)
#
#
#plots=lapply(names(afds), function(snn) {
#    plotIndividualExperiment(afds[[snn]], snn) 
#      })
#
##e.g.  to visualize just call
##plots[[1]]
#
##dump individual plots somewhere ---------------
##plot.dir= '/home/jbloom/Dropbox/code/xQTL/elegans/plots/'
##for(snn in names(plots)) {
##    ggsave(paste0( plot.dir, snn, '.png'), plots[[snn]], width=16)
##}
#
##I have to restructure this stuff  ---------------------------------------------
#S2_N2_XZ1516_7_RNAi=afds[[2]]
#LWM1_S56_N2_XZ1516_3_RNAi=afds[[15]]
#results=calcContrastStats(results=list(S2_N2_XZ1516_7_RNAi, LWM1_S56_N2_XZ1516_3_RNAi),
#                          L=paste0('_', names(afds)[2]), R=paste0('_', names(afds)[15]) ) #'_high1', R='_unsel1')
#
#sc=plotContrast(results, 'S2_N2_XZ1516_7_RNAi', 'LWM1_S56_N2_XZ1516_3_RNAi') #S3_N2_XZ1516_7_KO')
#
#s=plotSummary(results, effective.n.tests=2000)
#
#a=ggarrange(plots[[2]],plots[[15]], s,  nrow=3)
#a
#
#
#
#
#
#
#S2_N2_XZ1516_7_RNAi=afds[[2]]
#S3_N2_XZ1516_7_KO=afds[[3]]
#results=calcContrastStats(results=list(S2_N2_XZ1516_7_RNAi, S3_N2_XZ1516_7_KO),
#                          L=paste0('_', names(afds)[2]), R=paste0('_', names(afds)[3]) ) #'_high1', R='_unsel1')
#
#sc=plotContrast(results, 'S2_N2_XZ1516_7_RNAi', 'S3_N2_XZ1516_7_KO')
#
#s=plotSummary(results, effective.n.tests=2000)
#
#a=ggarrange(plots[[2]],plots[[3]], s,  nrow=3)
#a
#
#S21_N2_XZ1516_5_RNAi=afds[[13]]
#S22_N2_XZ1516_5_KO=afds[[14]] #S2_N2_XZ1516_7_RNAi=afds[[2]]
#results=calcContrastStats(results=list(S21_N2_XZ1516_5_RNAi, S22_N2_XZ1516_5_KO), #S3_N2_XZ1516_7_KO),
#                          L=paste0('_', names(afds)[13]), R=paste0('_', names(afds)[14]) ) #'_high1', R='_unsel1')
##sc=plotContrast(results, 'S2_N2_XZ1516_7_RNAi', 'S3_N2_XZ1516_7_KO')
#s=plotSummary(results, effective.n.tests=2000)
#b=ggarrange(plots[[13]],plots[[14]], s,  nrow=3)
#b
#
#
#S12_N2_XZ1516_6_RNAi=afds[[12]]
#S1_N2_XZ1516_6_KO=afds[[1]]
#results=calcContrastStats(results=list(S12_N2_XZ1516_6_RNAi, S1_N2_XZ1516_6_KO), #S3_N2_XZ1516_7_KO),
#                          L=paste0('_', names(afds)[12]), R=paste0('_', names(afds)[1]) ) 
#s=plotSummary(results, effective.n.tests=2000)
#cc=ggarrange(plots[[12]],plots[[1]], s,  nrow=3)
#cc
#
#S12_N2_XZ1516_6_RNAi=afds[[12]]
#S1_N2_XZ1516_6_KO=afds[[1]]
#results=calcContrastStats(results=list(S12_N2_XZ1516_6_RNAi, S1_N2_XZ1516_6_KO), #S3_N2_XZ1516_7_KO),
#                          L=paste0('_', names(afds)[12]), R=paste0('_', names(afds)[1]) ) 
#s=plotSummary(results, effective.n.tests=2000)
#cc=ggarrange(plots[[12]],plots[[1]], s,  nrow=3)
#cc



simPheno2=function(FR, genMapMarkers, QTL.sims,ds.size=NULL,returnG=T) {
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
        simy=XB+rnorm(nrow(X_Q), mean=0, sd=QTL.sims$o.error.sd) 
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

