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

#hack a vcf with het calls at each variant site for gatk ASEReadCounter
vcfout=vcf
vcfout@gt=vcfout@gt[,c(1,2)]
vcfout@gt[,1]='GT:DP'
vcfout@gt[,2]='0/1:.'
write.vcf(vcfout,'/data1/yeast/reference/pop_vcfs/1011/mega_filtered.vcf')
#then on the command line, bgzip and gatk IndexFeatureFile

# figure out reticulare behavior
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




