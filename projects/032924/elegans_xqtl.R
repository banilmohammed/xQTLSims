library(xQTLStats)
library(qs)
library(tidyverse)
library(vcfR)
library(AlphaSimR)
library(Rfast)
library(ggpubr)

#pull from github
xQTLSims.dir = '/data0/elegans/xQTLSims/'

source.dir=paste0(xQTLSims.dir, 'R/')
data.dir=paste0(xQTLSims.dir, 'data/')
project.dir=paste0(xQTLSims.dir, '032924/')

#function to simulate crosses, treatement of X is incomplete/broken
#and additional helper functions
source(paste0(source.dir, 'simWormCrosses.R'))
source(paste0(source.dir, 'helperFxs.R'))
source(paste0(source.dir, 'makeCountTables.R'))

#unique chromosomes 
uchr=c(as.character(as.roman(1:5)), 'X') #paste0('chr', as.roman(1:16))

gmap.file=paste0(data.dir, 'geneticMapXQTLsnplist.rds')
gmap=restructureGeneticMap(gmap.file)

#pretty intensive memory usage here
#include web link to vcf file 
#elegans.isotypes.vcf=paste0(data.dir,'WI.20220216.impute.isotype.vcf.gz')

#filtered vcf as qsave objects
# !!!! find the premade objects here folks !!!! :
# /u/project/kruglyak/jsbloom/xQTL/elegans/ref/
#and place in your data.dir
elegans.isotypes.vcf.qs=paste0(data.dir,'WI.20220216.vcf.qs')

#filtered numeric gt calls from vcf  as qsave object
elegans.isotypes.vcf.gt.qs=paste0(data.dir,'WI.20220216.vcf.GT.qs')

#run once, Laura skip this ============================================================
#use premade objects to save yourself the memory related headache, but this is how the objects are made
#preprocessVCF(elegans.isotypes.vcf,elegans.isotypes.vcf.qs,elegans.isotypes.vcf.gt.qs)
#======================================================================================

#Laura start here ======================================================================
vcf=qread(elegans.isotypes.vcf.qs)
gt=qread(elegans.isotypes.vcf.gt.qs)

sample.key.file=paste0(project.dir, 'samplekey_032924.txt')
sample.key=read_tsv(sample.key.file)
sample.dir='/media/hoffman2/thatguy0/nextflow_pipes/bulkGWA/BulkGWA-20240405/undup_alignments/'
#only have tables for a subset of the samples
sample.key=sample.key[c(1:12,21,22),]
countdfs=makeCountTables(sample.key,sample.dir, vcf,gt)

sample.key.file=paste0(project.dir, 'samplekey_earlier.txt')
sample.key=read_tsv(sample.key.file)
sample.key$'sample name'=paste0('LWM', 1:6, '_S', 56:61)
sample.dir='/media/hoffman2/lwalterm/032924/count/' #media/hoffman2/thatguy0/nextflow_pipes/bulkGWA/BulkGWA-20240405/undup_alignments/'
countdfs2=makeCountTables(sample.key,sample.dir, vcf,gt, sample.suffix='.txt')

#all the filtered and phased count tables, glued together
countdfs=c(countdfs,countdfs2)
rm(countdfs2)

#calculate the allele frequencies and SEs, params here need some more dialing in 
#should parallelize this .... yawn 
afds=lapply(names(countdfs), function(snn) {
       calcAFD(countdfs[[snn]], experiment.name=snn,sample.size=1e4, sel.strength=.95, bin.width=3000, eff.length=2000, uchr=uchr)
      })
names(afds)=names(countdfs)


plots=lapply(names(afds), function(snn) {
    plotIndividualExperiment(afds[[snn]], snn) 
      })

#e.g.  to visualize just call
#plots[[1]]

#dump individual plots somewhere ---------------
#plot.dir= '/home/jbloom/Dropbox/code/xQTL/elegans/plots/'
#for(snn in names(plots)) {
#    ggsave(paste0( plot.dir, snn, '.png'), plots[[snn]], width=16)
#}

#I have to restructure this stuff  ---------------------------------------------
S2_N2_XZ1516_7_RNAi=afds[[2]]
LWM1_S56_N2_XZ1516_3_RNAi=afds[[15]]
results=calcContrastStats(results=list(S2_N2_XZ1516_7_RNAi, LWM1_S56_N2_XZ1516_3_RNAi),
                          L=paste0('_', names(afds)[2]), R=paste0('_', names(afds)[15]) ) #'_high1', R='_unsel1')

sc=plotContrast(results, 'S2_N2_XZ1516_7_RNAi', 'LWM1_S56_N2_XZ1516_3_RNAi') #S3_N2_XZ1516_7_KO')

s=plotSummary(results, effective.n.tests=2000)

a=ggarrange(plots[[2]],plots[[15]], s,  nrow=3)
a






S2_N2_XZ1516_7_RNAi=afds[[2]]
S3_N2_XZ1516_7_KO=afds[[3]]
results=calcContrastStats(results=list(S2_N2_XZ1516_7_RNAi, S3_N2_XZ1516_7_KO),
                          L=paste0('_', names(afds)[2]), R=paste0('_', names(afds)[3]) ) #'_high1', R='_unsel1')

sc=plotContrast(results, 'S2_N2_XZ1516_7_RNAi', 'S3_N2_XZ1516_7_KO')

s=plotSummary(results, effective.n.tests=2000)

a=ggarrange(plots[[2]],plots[[3]], s,  nrow=3)
a

S21_N2_XZ1516_5_RNAi=afds[[13]]
S22_N2_XZ1516_5_KO=afds[[14]] #S2_N2_XZ1516_7_RNAi=afds[[2]]
results=calcContrastStats(results=list(S21_N2_XZ1516_5_RNAi, S22_N2_XZ1516_5_KO), #S3_N2_XZ1516_7_KO),
                          L=paste0('_', names(afds)[13]), R=paste0('_', names(afds)[14]) ) #'_high1', R='_unsel1')
#sc=plotContrast(results, 'S2_N2_XZ1516_7_RNAi', 'S3_N2_XZ1516_7_KO')
s=plotSummary(results, effective.n.tests=2000)
b=ggarrange(plots[[13]],plots[[14]], s,  nrow=3)
b


S12_N2_XZ1516_6_RNAi=afds[[12]]
S1_N2_XZ1516_6_KO=afds[[1]]
results=calcContrastStats(results=list(S12_N2_XZ1516_6_RNAi, S1_N2_XZ1516_6_KO), #S3_N2_XZ1516_7_KO),
                          L=paste0('_', names(afds)[12]), R=paste0('_', names(afds)[1]) ) 
s=plotSummary(results, effective.n.tests=2000)
cc=ggarrange(plots[[12]],plots[[1]], s,  nrow=3)
cc

S12_N2_XZ1516_6_RNAi=afds[[12]]
S1_N2_XZ1516_6_KO=afds[[1]]
results=calcContrastStats(results=list(S12_N2_XZ1516_6_RNAi, S1_N2_XZ1516_6_KO), #S3_N2_XZ1516_7_KO),
                          L=paste0('_', names(afds)[12]), R=paste0('_', names(afds)[1]) ) 
s=plotSummary(results, effective.n.tests=2000)
cc=ggarrange(plots[[12]],plots[[1]], s,  nrow=3)
cc







paste0('_', names(afds)[2]), paste0('_', names(afds)[3]) )

#subset shit
#used.markers=match(scounts$ID, genMap$id)
#FN[used.markers,]



#countdf.h=simSequencingRefAlt(simFR$simy,FR, genMap$id, depth=10000, sel.frac=.1, lower.tail=F)
#if you set sel.frac=1 sample the population of existing genotypes without QTL effects 
#countdf.l=simSequencingRefAlt(y=NULL, FR,genMap$id, depth=10000, sel.frac=1 , lower.tail=F)

#sanity checks 
plot(countdf.h$alt/(countdf.h$alt+countdf.h$ref))
points(countdf.h$expected, col='red') #alt/(countdf$alt+countdf$ref))

#sanity checks 
plot(countdf.l$alt/(countdf.l$alt+countdf.l$ref))
points(countdf.l$expected, col='red') 

###############Block of code to simulate sequencing data and analysis for a bi-parental x-qtl####################################################################
countdf.h=phaseBiparental(countdf.h, p.names[1], founderPop, genMap)
countdf.l=phaseBiparental(countdf.l, p.names[1], founderPop, genMap)

#plot(df$p1/(df$p1+df$p2))
#points(df$expected.phased, col='red') #alt/(countdf$alt+countdf$ref))

#QTL.sims$o.add.qtl.ind, asRaw=F)

test  = calcAFD(countdf.h, experiment.name='high1',sample.size=1e4, sel.strength=.1, bin.width=1000, eff.length=300, uchr=unique(genMap$chr) ) 
test2 = calcAFD(countdf.l, experiment.name='unsel1',sample.size=1e4, sel.strength=1, bin.width=1000, eff.length=300, uchr=unique(genMap$chr) )
results=calcContrastStats(results=list(test, test2), L='_high1', R='_unsel1')

# test$expected.phased_high1=1-test$expected.phased_high1

h1=plotIndividualExperiment(test, 'high1')
u1=plotIndividualExperiment(test2, 'unsel1')
c1=plotContrast(results, 'high1', 'unsel1')
ggpubr::ggarrange(h1, u1, c1, nrow=3) 

