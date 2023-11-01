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
elegans.isotypes.vcf='/data0/elegans/xQTLSims/WI.20220216.impute.isotype.vcf.gz'

#filtered vcf as qsave object
elegans.isotypes.vcf.qs='/data0/elegans/xQTLSims/WI.20220216.vcf.qs'

#filtered numeric gt calls from vcf  as qsave object
elegans.isotypes.vcf.gt.qs='/data0/elegans/xQTLSims/WI.20220216.vcf.GT.qs'

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
#Laura load these 
vcf=qread(elegans.isotypes.vcf.qs)
gt=qread(elegans.isotypes.vcf.gt.qs)
                      

p.names=c(
'N2', 
'ECA191', 
'QG2832' ,
'NIC195' ,
'XZ1516' ,
'QX1791' ,
'QX1211',
'ECA369' ,
'XZ1514',
'ECA738' ,
'ECA760',
'ECA1255')

#take the larger vcf,  genotype calls, and a subset of parents 
#extract segregating sites 
#return an alphaSimR founder population
createFounderPop=function(vcf, gt, p.names) { 
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

    gt.sub=gt.sub[-which(bad.sites),]
    vcf.cross=vcf[match(rownames(gt.sub), rownames(gt)), samples=colnames(gt.sub)]
    #generate sample ID
    vcf.cross=vcfR::addID(vcf.cross)

    imputed.positions=jitterGmapVector(getGmapPositions(vcf.cross, gmap, uchr))

    #genetic map positions must be in Morgans
    genMap=data.frame(markerName=paste0(getCHROM(vcf.cross),'_',getPOS(vcf.cross)), chromosome=getCHROM(vcf.cross), position=unlist(imputed.positions)/100)

    teg.GT=t(gt.sub)
    #recode
    teg.GT[teg.GT==0]=-1
    colnames(teg.GT)=paste0(getCHROM(vcf.cross),'_',getPOS(vcf.cross))
    ped=data.frame(id=rownames(teg.GT), mother=rep(0, nrow(teg.GT)), father=rep(0,nrow(teg.GT)) ) #c(0,0), father=c(0,0))
    return(importInbredGeno(geno=teg.GT, genMap=genMap, ped=ped))
}

p.names=c('N2', 'XZ1516')
founderPop = createFounderPop(vcf,gt, p.names) #c('N2', 'XZ1516'))






#cChr() is how we can eventually deal with X properly 


starting.sample.size=10
#if mate or self, estimated brood size 
brood.size.selfing=20
brood.size.mating=brood.size.selfing*2

#fraction of hermaphrodites that self
selfing.rate=.2
#

max.per.gen=2e4

max.gen=12

SP=SimParam$new(founderPop)
FN=newPop(founderPop, simParam=SP)

nProgenyF1.mated=starting.sample.size*(1-selfing.rate)*brood.size.mating
nProgenyF1.selfed=starting.sample.size*(selfing.rate)*brood.size.selfing

# Setup for F1, deliberate crossings-------------------------------------------------------------
# mating matrix 
# herm.index, male/herm index, mating=0/herm=1
#biparental
mating.matrix=matrix(c(1,2),nrow=1)
selfings=c(1)

#multiparental
mating.matrix=matrix(rbind(c(1,2),c(3,4), c(5,6), c(7,8),c(9,10), c(11,12)), nrow=6)
selfings=c(1,3,5,7,9,11)


mating.matrix=cbind(mating.matrix,0)
mating.matrix=do.call(rbind, replicate(nProgenyF1.mated, mating.matrix, simplify=F))

self.cnt=nProgenyF1.selfed
if(self.cnt>0) {
       selfing.matrix=cbind(selfings, selfings,1)
       selfing.matrix=do.call(rbind, replicate(nProgenyF1.selfed, selfing.matrix, simplify=F))
       mating.matrix=rbind(mating.matrix, selfing.matrix)
       mating.matrix=mating.matrix[order(mating.matrix[,3]),]
}
if(nrow(mating.matrix)>max.per.gen){
        mating.matrix=mating.matrix[sample(1:nrow(mating.matrix), max.per.gen),]
        mating.matrix=mating.matrix[order(mating.matrix[,3]),]
}
#-----------------------------------------------------------------------------------------------


current.gen=1
while(current.gen<=max.gen) {
       print(current.gen) 
       #at least loop here 
       mated.index=which(mating.matrix[,3]==0)
       selfed.index=which(mating.matrix[,3]==1)
       print(paste('progeny from mating', length(mated.index)))
       print(paste('progeny from selfing', length(selfed.index)))

       FN=makeCross(FN, mating.matrix[,-3], nProgeny=1,simParam=SP)

       #sanity check selfing 
       #FN3=makeCross(FN2[1], matrix(c(1,1), nrow=1), nProgeny=1)
       #t2=pullMarkerGeno(FN3[1], colnames(teg.GT))
       #FN3=makeCross(FN3[1], matrix(c(1,1), nrow=1), nProgeny=1)
       #t2=pullMarkerGeno(FN3[1], colnames(teg.GT))
       #plot(t2[1,])

       #half of the progeny from matings will be hermaphrodites
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
       mated.herm=sample(herm.index, size=min(length(herm.index), length(male.index)))
       mating.matrix=cbind(male.index,mated.herm,0)

       #generate cross progeny for matings, with brood size dependent on expected brood size if mating 
       mating.matrix=do.call(rbind, replicate(brood.size.mating, mating.matrix, simplify=F))
      
       #note which hermaphrodites aren't mated 
       unmated.herm=herm.index[!(herm.index %in% mated.herm)]

       #how many do we expected to self given the total number of hermaphrodites and the selfing rate 
       self.cnt= round(length(herm.index)*selfing.rate)
       #if we expect any selfers then ... 
       if(self.cnt>0) {
           #figure out which ones self 
           selfings=sample(unmated.herm, min(length(unmated.herm),self.cnt))
           #make a selfing matrix
           selfing.matrix=cbind(selfings, selfings,1)
           #make brood size reflect brood size per selfing
           selfing.matrix=do.call(rbind, replicate(brood.size.selfing, selfing.matrix, simplify=F))
           mating.matrix=rbind(mating.matrix, selfing.matrix)
        }
       #total worms at this point 
       print(nrow(mating.matrix))

       #bottleneck
       if(nrow(mating.matrix)>max.per.gen){
            mating.matrix=mating.matrix[sample(1:nrow(mating.matrix), max.per.gen),]
            mating.matrix=mating.matrix[order(mating.matrix[,3]),]
       }
       current.gen=current.gen+1
   }


t2=pullMarkerGeno(FN[1:1000],getGenMap(founderPop)$id)

rcnt=colSums(t2)
plot(rcnt/(nrow(t2)*2))










#F1 ----------------------------------------------------------------------------------------------------------    
nProgenyF1.mated=starting.sample.size*(1-selfing.rate)*brood.size.mating
nProgenyF1.selfed=starting.sample.size*(selfing.rate)*brood.size.selfing
FN=makeCross(pop, crossPlan, nProgeny=nProgenyF1.mated,simParam=SP)
#t2=pullMarkerGeno(F1, colnames(teg.GT))
#t2=pullMarkerGeno(F2.selfed, colnames(teg.GT))
herm.index=sort(sample.int( nInd(FN), size=.5* nInd(FN)))
male.index=which(!((1: nInd(FN)) %in% herm.index))

 if(selfing.rate>0)  {
        F2.selfed=makeCross(pop, matrix(c(1,1),nrow=1,ncol=2), nProgeny=nProgenyF1.selfed, simParam=SP)
        herm.index=c(herm.index,nInd(FN)+seq(1,nInd(F2.selfed))) #length(F2.selfed)))
        FN=c(FN, F2.selfed)
    }
#--------------------------------------------------------------------------------------------------------------

    mated.herm=sample(herm.index, size=length(male.index))
    mating.matrix=cbind(male.index,mated.herm,0)
    mating.matrix=do.call(rbind, replicate(brood.size.mating, mating.matrix, simplify=F))

    unmated.herm=herm.index[!(herm.index %in% mated.herm)]
    self.cnt= round(length(herm.index)*selfing.rate)

    if(self.cnt>0) {
       selfings=sample(unmated.herm, self.cnt)
       selfing.matrix=cbind(selfings, selfings,1)
       selfing.matrix=do.call(rbind, replicate(brood.size.selfing, selfing.matrix, simplify=F))
       mating.matrix=rbind(mating.matrix, selfing.matrix)
    }
    if(nrow(mating.matrix)>max.per.gen){
        mating.matrix=mating.matrix[sample(1:nrow(mating.matrix), max.per.gen),]
        mating.matrix=mating.matrix[order(mating.matrix[,3]),]
    }
   #----------------------------------------------------------------------------------------------------------























#how many chromosomes
n_chr=length(imputed.positions)

#total gmap size
L=round(sapply(gmap, function(x) max(x$map)))
#setup for Meiosis package
xoparam = Meiosis::create_xoparam(L,obligate_chiasma=T) 
 
#nsegs=
#number of males and hermaphrodites at F0
starting.sample.size=10

#if mate or self, estimated brood size 
brood.size.selfing=20
brood.size.mating=brood.size.selfing*2

#fraction of hermaphrodites that self
selfing.rate=.2

max.per.gen=1e4

max.gen=12

#simcross
#p1=create_parent(L=55.56, allele=0)
#p2=create_parent(L=55.56, allele=1)
#p3=create_parent(L=55.56, allele=2)
#p4=create_parent(L=55.56, allele=3)
#p3=create_parent(L=55.56, allele=2)
#F1.1=cross(p1,p2,m=0, obligate_chiasma=T)
#F1.2=cross(p3,p4,m=0, obligate_chiasma=T)
#F2=cross(F1.1,F1.2)
#F2.pop=replicate(1000,cross(F1,F1))



 #parental haplotypes
    p1.ind=split(as.numeric(eg.GT[,1])+1, vcfR::getCHROM(vcf.cross))
    p1.ind=p1.ind[uchr]
    p2.ind=split(as.numeric(eg.GT[,2])+1, vcfR::getCHROM(vcf.cross))
    p2.ind=p2.ind[uchr]

    #creation of mated F1
    ind=list(p1.ind[1], p2.ind[1])
    
    #creation of non.mated F1   
    ind2=list(p1.ind, p1.ind)
    #ind3=list(p2.ind, p2.ind)

    #F2-------------------------------------------------------------------------------------------------------------
    #if it mates, 50/50 hermaphrodites vs males
    FN=replicate(starting.sample.size*(1-selfing.rate)*brood.size.mating,
                 Meiosis::cross_geno(father = ind, mother = ind, positions = imputed.positions, xoparam=xoparam), 
                       simplify=F)
    herm.index=sort(sample.int(length(FN), size=.5*length(FN)))
    male.index=which(!((1:length(FN)) %in% herm.index))
    
    #if it selfs, all hermaphrodites
    if(selfing.rate>0)  {
        F2.selfed=replicate(starting.sample.size*(selfing.rate)*brood.size.selfing,
                     Meiosis::cross_geno(father = ind2, mother = ind2, positions = imputed.positions, xoparam=xoparam), 
                           simplify=F)
        herm.index=c(herm.index,length(FN)+seq(1,length(F2.selfed)))
        FN=c(FN, F2.selfed)
    }
    
    mated.herm=sample(herm.index, size=length(male.index))
    mating.matrix=cbind(male.index,mated.herm,0)
    mating.matrix=do.call(rbind, replicate(brood.size.mating, mating.matrix, simplify=F))

    unmated.herm=herm.index[!(herm.index %in% mated.herm)]
    self.cnt= round(length(herm.index)*selfing.rate)

    if(self.cnt>0) {
       selfings=sample(unmated.herm, self.cnt)
       selfing.matrix=cbind(selfings, selfings,1)
       selfing.matrix=do.call(rbind, replicate(brood.size.selfing, selfing.matrix, simplify=F))
       mating.matrix=rbind(mating.matrix, selfing.matrix)
    }
    if(nrow(mating.matrix)>max.per.gen){
        mating.matrix=mating.matrix[sample(1:nrow(mating.matrix), max.per.gen),]
        mating.matrix=mating.matrix[order(mating.matrix[,3]),]
    }
   #----------------------------------------------------------------------------------------------------------

   current.gen=3
   while(current.gen<=max.gen) {
       print(current.gen) 
       #at least loop here 
       mated.index=which(mating.matrix[,3]==0)
       selfed.index=which(mating.matrix[,3]==1)
       print(paste('mated', length(mated.index)))
       print(paste('selfed', length(selfed.index)))
       FN=apply(mating.matrix, 1, function(x) { Meiosis::cross_geno(father=FN[[x[1]]], mother=FN[[x[2]]],  positions = imputed.positions, xoparam=xoparam) })
       #half of the progeny from matings will be hermaphrodites
       herm.index=sort(sample(mated.index,size=length(mated.index)/2))
       #the rest will be males 
       male.index=which(!((mated.index %in% herm.index)))
       if(length(selfed.index)>0) {
            herm.index=c(herm.index,selfed.index)
        }
       print(paste('total herm', length(herm.index)))
       print(paste('total male', length(male.index)))
       mated.herm=sample(herm.index, size=min(length(herm.index), length(male.index)))
       mating.matrix=cbind(male.index,mated.herm,0)
       mating.matrix=do.call(rbind, replicate(brood.size.mating, mating.matrix, simplify=F))
       unmated.herm=herm.index[!(herm.index %in% mated.herm)]
       self.cnt= round(length(herm.index)*selfing.rate)
       if(self.cnt>0) {
           selfings=sample(unmated.herm, min(length(unmated.herm),self.cnt))
           selfing.matrix=cbind(selfings, selfings,1)
           selfing.matrix=do.call(rbind, replicate(brood.size.selfing, selfing.matrix, simplify=F))
           mating.matrix=rbind(mating.matrix, selfing.matrix)
        }
       print(nrow(mating.matrix))
       if(nrow(mating.matrix)>max.per.gen){
            mating.matrix=mating.matrix[sample(1:nrow(mating.matrix), max.per.gen),]
            mating.matrix=mating.matrix[order(mating.matrix[,3]),]
       }
       current.gen=current.gen+1
   }


test=sapply(FN, function(z) { do.call('c', mapply(function(x,y) x+y, z$paternal, z$maternal)) })
rcnt=rowSums(test)
plot(rcnt/(ncol(test)*2))































just get biallelic variants for now 
#get genotypes

p1.name='N2'
p2.name='XZ1516'

#save as R object
#extract genos for parents 
gt.sub=gt[,c(p1.name,p2.name)]
#find variant sites
gt.sub=gt.sub[gt.sub[,1]!=gt.sub[,2],]
#dump sites with missing info for now
gt.sub=gt.sub[!( is.na(gt.sub[,1]) | is.na(gt.sub[,2])),]
gt.sub=gt.sub[ !(gt.sub[,1]=='0|1' | gt.sub[,1]=='1|0' | 
                gt.sub[,2]=='0|1' | gt.sub[,2]=='1|0' ),]
#remove more complex structural stuff or mitochondrial variants for now
gt.sub=gt.sub[!grepl('simple|deletion|complex|duplication|chrM|MtDNA', rownames(gt.sub)),]
#subset vcf ile
vcf.cross=vcf[match(rownames(gt.sub), rownames(gt)), samples=c(p1.name,p2.name)]
#generate sample ID
vcf.cross=vcfR::addID(vcf.cross)

#assuming input vcf with haploid parents and haploid specific variant calls 
# [1] "273614xa" "BYa"      "CBS2888a" "CLIB219x" "CLIB413a" "I14a"    
# [7] "M22"      "PW5a"     "RMx"      "Y10x"     "YJM145x"  "YJM454a" 
#[13] "YJM978x"  "YJM981x"  "YPS1009x" "YPS163a"

#specify sample size


#specify cross (see crosses.to.parents)
#cross='A'

#get the two parent names 
#p1.name=crosses.to.parents[[cross]][1]gg
#p2.name=crosses.to.parents[[cross]][2]

#extract parental genotypes at variant sites between those parents from a reference vcf
#vcf.cross=getCrossVCF(ref.vcf,p1.name, p2.name)
#vcf.cross
#gmap=gmaps[['A']]

eg.GT=vcfR::extract.gt(vcf.cross)
#recode
eg.GT[eg.GT=='0|0']=0
eg.GT[eg.GT=='1|1']=1
eg.GT=matrix(as.numeric(eg.GT), nrow(eg.GT), ncol(eg.GT))

