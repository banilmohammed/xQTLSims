library(xQTLStats)
library(vcfR)
library(qs)
#switch to AlphaSimR for simulation backend 
#library(Meiosis)
library(AlphaSimR)
library(Rfast)


source.dir='/data0/elegans/xQTLSims/R/'
data.dir='/data0/elegans/xQTLSims/'

#function to simulate crosses, treatement of X is incomplete/broken
source(paste0(source.dir, 'simWormCrosses.R'))
#additional helper functions 
source(paste0(source.dir, 'helperFxs.R'))

#unique chromosomes 
uchr=c(as.character(as.roman(1:5)), 'X') #paste0('chr', as.roman(1:16))

gmap.file=paste0(data.dir, 'geneticMapXQTLsnplist.rds')

gmap=restructureGeneticMap(gmap.file)

#sapply(gmap, function(x) max(x$map))
#data(crosses.to.parents)
#genetic maps 
#gmaps=readRDS(system.file('reference', 'yeast_gmaps.RDS', package='xQTLStats'))

#reference vcf (assume ploidy=1)
#ref.vcf=system.file('reference', 'parents_w_svar_sorted.vcf.gz', package='xQTLStats')

#pretty intensive memory usage here
#include web link to vcf file 
elegans.isotypes.vcf=paste0(data.dir,'WI.20220216.impute.isotype.vcf.gz')

#filtered vcf as qsave object
elegans.isotypes.vcf.qs=paste0(data.dir,'WI.20220216.vcf.qs')

#filtered numeric gt calls from vcf  as qsave object
elegans.isotypes.vcf.gt.qs=paste0(data.dir,'WI.20220216.vcf.GT.qs')


#run once, Laura skip this ============================================================
preprocessVCF(elegans.isotype.vcf,elegans.isotypes.vcf.qs,elegans.isotypes.vcf.gt.qs)
#======================================================================================

, 
#Laura start here ======================================================================
vcf=qread(elegans.isotypes.vcf.qs)
gt=qread(elegans.isotypes.vcf.gt.qs)
                      
#existing fog2 ko strains 
p.names=c('N2', 'ECA191', 'QG2832' ,'NIC195' ,'XZ1516' ,'QX1791' ,'QX1211','ECA369' ,'XZ1514','ECA738' ,'ECA760','ECA1255')

p.names=c('ECA191', 'CB4856') #XZ1516')
#founderPop = createFounderPop(vcf,gt, p.names, uchr) #c('N2', 'XZ1516'))
founderPop = createFounderPop(vcf,gt, p.names, gmap, X.drop=T) #c('N2', 'XZ1516'))
#sexChr=F
#founderPop = createFounderPop(vcf,gt, p.names, X.only=T, X.drop=F) #c('N2', 'XZ1516'))
#sexChr=T

#cChr() is how we can eventually deal with X properly 



SP=SimParam$new(founderPop)
#unfortunately these functions fail to provide sufficient flexibility wrt trait architectures
#and are behaving unexpectedly
#SP$addTraitA(nQtlPerChr=1)
#SP$setVarE(H2=0.4)

#create founder population 
FN=newPop(founderPop, simParam=SP)
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
o.nQTL=5
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
              o.error.sd=1,
              # toggle to ignore residual error sd of trait effect and instead normalize
              # trait variance such that residual error is 1-h^2
              o.h2.norm=F,
              o.h2=.5)

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

FR=simWormCrosses(SimWormParams)

#get full set of marker genotypes 



#when functionalizing, return simy 

#sanity check 
summary(lm(simy~X_Q))

rG=cor(simy, scale(G))

plot(rG[1,])
abline(v=match(QTL.sims$o.add.qtl.ind, colnames(G)))
#plot(colSums(G)/(nrow(G)*2))
#======================================================================================================================




countdf.h=simSequencingRefAlt(simy,G,depth=50, sel.frac=.1, lower.tail=F)
countdf.l=simSequencingRefAlt(simy,G,depth=50, sel.frac=1 , lower.tail=F)

plot(countdf.h$alt/(countdf.h$alt+countdf.h$ref))
points(countdf.h$expected, col='red') #alt/(countdf$alt+countdf$ref))

plot(countdf.l$alt/(countdf.l$alt+countdf.l$ref))
points(countdf.l$expected, col='red') 


countdf.h=phaseBiparental(countdf.h, p.names[1], founderPop, genMap)
countdf.l=phaseBiparental(countdf.l, p.names[1], founderPop, genMap)

#plot(df$p1/(df$p1+df$p2))
#points(df$expected.phased, col='red') #alt/(countdf$alt+countdf$ref))

#QTL.sims$o.add.qtl.ind, asRaw=F)


test  = calcAFD(countdf.h, experiment.name='high1',sample.size=1e4, sel.strength=.1, bin.width=1000, eff.length=300, uchr=as.roman(1:5))
test2 = calcAFD(countdf.l, experiment.name='unsel1',sample.size=1e4, sel.strength=.1, bin.width=1000, eff.length=300, uchr=as.roman(1:5))
results=calcContrastStats(results=list(test, test2), L='_high1', R='_unsel1')

# test$expected.phased_high1=1-test$expected.phased_high1

h1=plotIndividualExperiment(test, 'high1')
u1=plotIndividualExperiment(test2, 'unsel1')
c1=plotContrast(results, 'high1', 'unsel1')
ggpubr::ggarrange(h1, u1, c1, nrow=3) 



QS=simPheno(FR,genMap$id,QTL.sims, ds.size=2e3)


#standardize genotypes 
sG=Rfast::standardise(QS$G)
#map QTL 
qtl.detected=doTraitFDR(scale(QS$simy), sG, QS$G, nperm=1e2)




























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


