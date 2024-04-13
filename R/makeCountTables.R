#extract ref and alt counts for variants in bi-parental crosses listed in sample,key
#also phase the variants 
makeCountTables=function(sample.key, vcf, gt) {
    countdfs=list() 
    for(i in 1:nrow(sample.key) ) { 
        #nrow(sample.key)) {
        p1=sample.key$'parent 1'[i]
        p2=sample.key$'parent 2'[i]

        #in this case, mate herm N2 with male CB
        p.names=c(p1, p2) #B4856') #XZ1516')
        mating.matrix=matrix(c(1,2), nrow=1)
        # in this case, given one-way cross, potentially allow N2 herm to self 
        founderPop = createFounderPop(vcf,gt, p.names, gmap, X.drop=F) #c('N2', 'XZ1516'))
       
        genMap=getGenMap(founderPop)
        
        sn=sample.key$'sample name'[i]
        scounts=read_tsv(paste0(sample.dir, sn, sample.suffix)) #'.table'))
        scounts.sub=scounts[paste0(scounts$contig, '_', scounts$position) %in% genMap$id,]
        scounts=data.frame(id=paste0(scounts.sub$contig, '_', scounts.sub$position),ref=scounts.sub$refCount, alt=scounts.sub$altCount)
        scounts=left_join(genMap, scounts, by='id')
        #come on GATK aseReadCounter, wtf is up with different length output given different BAM input
        #for know fill out with 0s 
        scounts$ref[is.na(scounts$ref)]=0
        scounts$alt[is.na(scounts$alt)]=0
        names(scounts)[1]='ID'
      
        #phase it 
        countdf=phaseBiparental(scounts, p.names[1], founderPop, genMap)
        
        #note, we need a better structure for keeping track of which parent is which 
        attr(countdf, 'p1')=p.names[1]
        attr(countdf, 'p2')=p.names[2]

        snn=paste(sample.key[i,], collapse='_')
        countdfs[[snn]]=countdf
    }

    return(countdfs)
}

