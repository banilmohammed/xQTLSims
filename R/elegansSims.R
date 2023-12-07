#some more complex simulations here 
library(xQTLStats)
library(vcfR)
library(qs)
#switch to AlphaSimR for simulation backend 
library(AlphaSimR)
#library(Meiosis)

#unique chromosomes 
uchr=c(as.character(as.roman(1:5)), 'X') #paste0('chr', as.roman(1:16))


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

gmap.file='/data0/elegans/xQTLSims/geneticMapXQTLsnplist.rds'


gmap=restructureGeneticMap(gmap.file)


#sapply(gmap, function(x) max(x$map))
#data(crosses.to.parents)
#genetic maps 
#gmaps=readRDS(system.file('reference', 'yeast_gmaps.RDS', package='xQTLStats'))

#reference vcf (assume ploidy=1)
#ref.vcf=system.file('reference', 'parents_w_svar_sorted.vcf.gz', package='xQTLStats')

#pretty intensive memory usage here
#include web link to vcf file 
elegans.isotypes.vcf='/data0/elegans/xQTLSims/WI.20220216.impute.isotype.vcf.gz'

#filtered vcf as qsave object
elegans.isotypes.vcf.qs='/data0/elegans/xQTLSims/WI.20220216.vcf.qs'

#filtered numeric gt calls from vcf  as qsave object
elegans.isotypes.vcf.gt.qs='/data0/elegans/xQTLSims/WI.20220216.vcf.GT.qs'

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

#run once, Laura skip this 
preprocessVCF(elegans.isotype.vcf,elegans.isotypes.vcf.qs,elegans.isotypes.vcf.gt.qs)
                           , 
#Laura start here  
vcf=qread(elegans.isotypes.vcf.qs)
gt=qread(elegans.isotypes.vcf.gt.qs)
                      


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

#existing fog2 ko strains 
p.names=c('N2', 'ECA191', 'QG2832' ,'NIC195' ,'XZ1516' ,'QX1791' ,'QX1211','ECA369' ,'XZ1514','ECA738' ,'ECA760','ECA1255')

p.names=c('N2', 'CB4856') #XZ1516')
#founderPop = createFounderPop(vcf,gt, p.names, uchr) #c('N2', 'XZ1516'))
founderPop = createFounderPop(vcf,gt, p.names, X.drop=T) #c('N2', 'XZ1516'))
#sexChr=F
#founderPop = createFounderPop(vcf,gt, p.names, X.only=T, X.drop=F) #c('N2', 'XZ1516'))
#sexChr=T

#cChr() is how we can eventually deal with X properly 



source('/data0/elegans/xQTLSims/R/simWormCrosses.R')

SP=SimParam$new(founderPop)
#these functions to provide sufficient flexibility and are behaving unexpectedly
#SP$addTraitA(nQtlPerChr=1)
#SP$setVarE(H2=0.4)

#create founder population 
FN=newPop(founderPop, simParam=SP)
genMap=getGenMap(founderPop)

#possible QTL architectures -----------------------------------------
nmarker=nrow(genMap)
nQTL=1
#sample(genMap$id, nQTL)
nadditive=nQTL
add.qtl.ind  = genMap$id[sort(sample(nmarker, nadditive))]
add.qtl.eff =  sample(ifelse(runif(nadditive)>.5,-1,1),replace=T)
#---------------------------------------------------------------------

#f designates QTL effects for hermaphrodites 
QTL.sims=list(
              f.add.qtl.ind=add.qtl.ind,
              f.add.qtl.eff=add.qtl.eff,
              f.error.sd=10
)

# Setup for F1, deliberate crossings-------------------------------------------------------------
# what we're calling 'mating.matrix' is really a matrix of who is crossing with who and will have replicated entries given expected broods 
# female, male, mating=0/selfing=1
# herm.index, male/herm index, mating=0/herm=1
#biparental, reciprocal
#mating.matrix=matrix(rbind(c(1,2),c(2,1)),nrow=2)
#selfings=c(1,2)
#biparental, non-reciprocal
#in this case, mate herm N2 with male CB
mating.matrix=matrix(c(1,2), nrow=1)
FN@sex[2]='M'
# in this case, given one-way cross, potentially allow N2 herm to self 
selfings=1

#multiparental setup
#mating.matrix=matrix(rbind(c(1,2),c(3,4), c(5,6), c(7,8),c(9,10), c(11,12)), nrow=6)
#selfings=c(1,3,5,7,9,11)


SimWormParams=list(

    #how many individual crossed worms per row in mating.matrix
    starting.sample.size=10 ,
    #brood size if it selfs
    brood.size.selfing=20 ,
    #brood size for mating 
    brood.size.mating=20*2 ,
    #fraction of hermaphrodites that self
    selfing.rate=0.1 ,
    #bottleneck per generation   
    max.per.gen=1e4,
    #how many total generations
    max.gen=5,
    #Founder population genotypes
    FN=FN ,
    #simulation parameters 
    SP=SP ,
    mating.matrix=mating.matrix ,
    selfings=selfings,
    genMap=genMap,
    QTL.sims=QTL.sims
)

FR=simWormCrosses(SimWormParams)

#get full set of marker genotypes 
G=pullMarkerGeno(FR,genMap$id,asRaw=F)
plot(colSums(G)/(nrow(G)*2))

X_Q=pullMarkerGeno(FR, QTL.sims$add.qtl.ind, asRaw=F)

X_Beta=QTL.sims$add.qtl.eff
if(length(X_Beta)==1) {
    XB=X_Q*X_Beta
} else {XB=X_Q%*%X_Beta 
}

error.sd=1
h2norm=F
h2=.5
#two ways to 
if(h2norm==F) {
    simy=X_Beta+rnorm(nrow(G), mean=0, sd=error.sd) 
    h2=var(XB)/(var(XB)+error.sd^2)
} else {
    simy= sqrt(h2)*scale(XB) + rnorm(nrow(G), mean=0, sd=sqrt((1-h2)/(h2*var(sqrt(h2)*scale(XB)))))
}




plot(colSums(G)/(nrow(G)*2))















extractGenoMatrix(



geno=pullMarkerGeno(FN[male.index[3000:4000]],simRecomb=F),genMap$id)





#run once with X and save as FNX
#FNX=FN
#then merge #fix this
#FN@geno[[6]]=FNX@geno[[1]]
#FN@nLoci=c(FN@nLoci, FNX@nLoci)
#FN@nChr=6
#FNm=cChr(FN, FNX)

#get alt counts per site
#0/1/2 for hets

#get just hets
hi=FR@sex=='H'


t2=pullMarkerGeno(FR[which(hi)[1:1000]],getGenMap(founderPop)$id)
rcnt=colSums(t2)
plot((rcnt/(2*nrow(t2))))

qm=getQtlMap(trait=1, simParam=SP)
p=pheno(FR[1:1000])
summary(lm(p~t2[,qm$id]))

#0/1 for X
#t3=pullMarkerGeno(reduceGenome(FN[male.index[3000:4000]],simRecomb=F),genMap$id)
#t2=rbind(t2,t3)
#rcnt2=colSums(t3)
#plot(rcnt2/(nrow(t3)*2))
#plot((rcnt+rcnt2)/((nrow(t2)*2)+(nrow(t3)*2)))













#SP$addTraitA(nQtlPerChr=1)
fn2=setPheno(FN, h2=.75, simParam=SP)

t2=pullMarkerGeno(fn2[1:1000],getGenMap(founderPop)$id)

varA(fn2)
varP(fn2)
rcnt=colSums(t2)
plot(rcnt/(nrow(t2)*2))

 abline(v=cumsum(c(0,rle(data.table::tstrsplit(colnames(t2), '_')[[1]])$lengths)))

r=cor(p,t2)
plot(r[1,]^2)



   #hermaphrodites
       FN.herm.index=unique(matings[,1])
       FN.male.index=unique(matings[,2])

       FN.males=reduceGenome(FN[FN.male.index], simRecomb=F, simParam=SP)
       t2=pullMarkerGeno(FN.males, getGenMap(founderPopX)$id) #colnames(teg.GT))
       t3=pullMarkerGeno(FN.males, getGenMap(founderPopX)$id) 

#Pilot code to deal with worm XX XO thing --------------------------------------------
       
#no recombination on X for male 
#recomb for hermX
#FNhx=reduceGenome(FN,simRecomb=T, simParam=SP)
#maleX

#rG=reduceGenome(FN,nProgeny=1)
#t2=pullMarkerGeno(rG,getGenMap(founderPop)$id)


#       matings=mating.matrix[mated.index,-3]
#       rG=reduceGenome(FN, simRecomb=F, simParam=SP)
#       #capture genotypes here 
      
#       FNmx=doubleGenome(rG)
#       FN=makeCross2(FN,FNmx, matings ,simParam=SP)
#       print(current.gen) 

#       if(length(selfed.index)>0) {
#            selfings=mating.matrix[selfed.index,-3]
#        #       #fake X as doubled haploid
#            selfings.cross=makeCross(FN,selfings, nProgeny=1,simParam=SP)
#            FN=c(FN,selfings.cross)
#       }
#       
#       #t2=pullMarkerGeno(FN[1:100], getGenMap(founderPopX)$id) #colnames(teg.GT))
#       #rcnt=colSums(t2)
#       #plot(rcnt/(nrow(t2)*2))
#
#       #------------------------------------------------------------------
#mg=pullMarkerGeno(reduceGenome(FN[male.index]),getGenMap(founderPopX)$id)
#mg=mg*2
#hg=pullMarkerGeno(FN[herm.index], getGenMap(founderPopX)$id)
#t2=rbind(hg,mg)

#rcnt=colSums(t2)
#plot(rcnt/(nrow(t2)*2))
#abline(v=cumsum(c(0,rle(data.table::tstrsplit(colnames(t2), '_')[[1]])$lengths)))


