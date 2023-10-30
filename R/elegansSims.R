library(xQTLStats)
library(Meiosis)
library(vcfR)
library(qs)
#some more complex simulations here 

#optional for doing simulations 
#remotes::install_url("https://cran.microsoft.com/snapshot/2017-09-17/src/contrib/Meiosis_1.0.2.tar.gz")
jitterGmapVector=function(themap, amount=1e-6) {
    for (i in 1:length(themap)) {
         n <- length(themap[[i]])
         themap[[i]] <- themap[[i]] + c(0, cumsum(rep(amount, n - 1)))
    }
    return(themap)
}
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
#genetic map from F10 AIL 
gmapRef=readRDS('/data0/elegans/xQTLSims/geneticMapXQTLsnplist.rds')
#normalize map to what's expected from F2
#https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000419
gmapRef$map=gmapRef$map/5 #5.3
#reorder
gmapRef=gmapRef[,c('map', 'chrom', 'pos')]
#rename columns (ppos = physical pos)
names(gmapRef)[3]='ppos'
gmapRef[,'chrom']=as.character(gmapRef[,'chrom'])
gmapRef=split(gmapRef, gmapRef$chrom)

gmap=gmapRef
sapply(gmap, function(x) max(x$map))
gmap$X$map=gmap$X$map*(51/max(gmap$X$map))

#data(crosses.to.parents)
#genetic maps 
#gmaps=readRDS(system.file('reference', 'yeast_gmaps.RDS', package='xQTLStats'))

#reference vcf (assume ploidy=1)
#ref.vcf=system.file('reference', 'parents_w_svar_sorted.vcf.gz', package='xQTLStats')

#pretty intensive memory usage here 
big.vcf=vcfR::read.vcfR('/data0/elegans/xQTLSims/WI.20220216.impute.isotype.vcf.gz')
qsave(big.vcf, file='/data0/elegans/xQTLSims/WI.20220216.vcf.qs')

vcf=qread('/data0/elegans/xQTLSims/WI.20220216.vcf.qs')
#just get biallelic variants for now 
vcf=vcf[vcfR::is.biallelic(vcf),]
#get genotypes
gt=vcfR::extract.gt(vcf)

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
gt.sub=gt.sub[!grepl('simple|deletion|complex|duplication|chrM', rownames(gt.sub)),]
#subset vcf ile
vcf.cross=vcf[match(rownames(gt.sub), rownames(gt)), samples=c(p1.name,p2.name)]
#generate sno ID
vcf.cross=vcfR::addID(vcf.cross)




#assuming input vcf with haploid parents and haploid specific variant calls 
# [1] "273614xa" "BYa"      "CBS2888a" "CLIB219x" "CLIB413a" "I14a"    
# [7] "M22"      "PW5a"     "RMx"      "Y10x"     "YJM145x"  "YJM454a" 
#[13] "YJM978x"  "YJM981x"  "YPS1009x" "YPS163a"

#specify sample size


#specify cross (see crosses.to.parents)
#cross='A'

#get the two parent names 
#p1.name=crosses.to.parents[[cross]][1]
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

uchr=c(as.character(as.roman(1:5)), 'X') #paste0('chr', as.roman(1:16))


imputed.positions=jitterGmapVector(getGmapPositions(vcf.cross, gmap, uchr))


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
selfing.rate=.1

max.per.gen=1e4

max.gen=12


 #parental haplotypes
    p1.ind=split(as.numeric(eg.GT[,1]), vcfR::getCHROM(vcf.cross))
    p1.ind=p1.ind[uchr]
    p2.ind=split(as.numeric(eg.GT[,2]), vcfR::getCHROM(vcf.cross))
    p2.ind=p2.ind[uchr]

    #creation of mated F1
    ind=list(p1.ind, p2.ind)
    
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


