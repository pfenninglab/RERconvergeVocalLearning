#rerun RERconverge
trait=read.table("Boreoeutheria_VLtrait_edited_converted.txt", sep="\t")

# make a named vector instead of a table
traitVals=trait$V2
names(traitVals)=trait$V1

trees=readRDS("zoonomiatrees.rds")
specstouse=trait$V1[!is.na(trait$V2)]
specstouse=specstouse[specstouse%in%trees$masterTree$tip.label]
#remove human
#specstouse=specstouse[!specstouse%in%c("REFERENCE")]
#saveRDS(specstouse, "specstouse.rds")
specstouse=readRDS("specstouse.rds")

#re-calculate RERs (without human)
#RERs=getAllResiduals(trees, useSpecies=specstouse)
#saveRDS(RERs, "VocalLearningRERs.rds")
RERs = readRDS("VocalLearningRERs.rds")

# VLpath=foreground2Paths(trait$V1[trait$V2==1], trees, clade="all", useSpecies = specstouse)
fgdspecs = names(traitVals[specstouse][traitVals[specstouse]==1])
VLpath=foreground2Paths(fgdspecs, trees, clade="all", useSpecies = specstouse)

#cors=correlateWithBinaryPhenotype(RERs, VLpath)
#saveRDS(cors, "VocalLearningCors.rds")
cors=readRDS("VocalLearningCors.rds")

#permulations
#my permulations
cors=cors
#bf=readRDS("VLBayesFactors.rds")

RERs=readRDS("VocalLearningRERs.rds")
#trees=readRDS("../NewComputerSpaceSaver/Zoonomia/zoonomiatrees.rds")
#trait=read.table("Boreoeutheria_VLtrait_edited_converted.txt", sep="\t")
# VLtrait=trait$V2
# names(VLtrait)=trait$V1
# VLtrait=VLtrait[!is.na(VLtrait)]

#t=foreground2Tree(names(VLtrait)[VLtrait==1], trees, plotTree=T, clade="all", useSpecies = names(VLtrait))
t=foreground2Tree(fgdspecs, trees, plotTree=T, clade="all", useSpecies = specstouse)
sum(t$edge.length) #45 total
sum(traitVals[specstouse]) #24 tip
sum(t$edge.length)-sum(traitVals[specstouse]) #21 internal

# trees=trees
#root="REFERENCE"
root = "vs_gorGor6" # closest species to human that has trait data
phenvec=traitVals[specstouse]
fgnum = 45
internal=21
drop = trees$masterTree$tip.label[!(trees$masterTree$tip.label %in% specstouse)]

statdf=NULL
pvaldf=NULL

count=1
while(count<=1000){
  
  #get phenotype:
  blsum=0
  if(is.null(fgnum)){
    fgnum=sum(phenvec)
  }
  tips=fgnum-internal
  t=drop.tip(trees$masterTree, drop)
  t=root.phylo(t, root, resolve.root = T)
  # while(blsum!=fgnum){
  while(blsum>(fgnum+5) | blsum<(fgnum-5)){
    rm=ratematrix(t, phenvec)
    sims=sim.char(t, rm, nsim = 1)
    nam=rownames(sims)
    s=as.data.frame(sims)
    simulatedvec=s[,1]
    names(simulatedvec)=nam
    top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
    tf=foreground2Tree(top, trees, clade="all", plotTree = F)
    blsum=sum(tf$edge.length)
  }
  
  #get path:
  p=tree2Paths(tf, trees, useSpecies = specstouse)
  #run correlation:
  c=correlateWithBinaryPhenotype(RERs, p)
  
  if(count==1){
    statdf=data.frame(c$Rho)
    rownames(statdf)=rownames(c)
    pvaldf=data.frame(c$P)
    rownames(pvaldf)=rownames(c)
  }else{
    temp=data.frame(c$Rho)
    rownames(temp)=rownames(c)
    statdf=merge(statdf, temp, by="row.names", all = T)
    rownames(statdf)=statdf$Row.names
    statdf=statdf[,-1]
    temp=data.frame(c$P)
    rownames(temp)=rownames(c)
    pvaldf=merge(pvaldf, temp, by="row.names", all = T)
    rownames(pvaldf)=pvaldf$Row.names
    pvaldf=pvaldf[,-1]
  }
  
  print(paste0("finished perm: ", count))
  count=count+1
}

saveRDS(statdf, "VLstatdf.rds")
saveRDS(pvaldf, "VLpvaldf.rds")


#get perm p-val
statdf=readRDS("VLstatdf.rds")
temp=data.frame(cors$Rho)
rownames(temp)=rownames(cors)
realandperm=merge(temp, statdf, by="row.names")
rownames(realandperm)=realandperm$Row.names
realandperm=realandperm[,-1]

permpval=data.frame(matrix(nrow=nrow(realandperm), ncol=1))
colnames(permpval)="permP"
rownames(permpval)=rownames(realandperm)

count=1
while(count<=nrow(realandperm)){
  real=realandperm[count,1]
  perms=as.numeric(realandperm[count,2:ncol(realandperm)])
  perms=perms[!is.na(perms)]
  p=sum(abs(perms)>abs(real))/length(perms)
  permpval[count, 1]=p
  print(count)
  count=count+1
}

VLcorswithpermp=cors
VLcorswithpermp=merge(VLcorswithpermp, permpval, by="row.names")
rownames(VLcorswithpermp)=VLcorswithpermp$Row.names
VLcorswithpermp=VLcorswithpermp[,-1]
VLcorswithpermp$permP.adj=p.adjust(VLcorswithpermp$permP, method="BH")
saveRDS(VLcorswithpermp, "VLcorswithpermp.rds")

#add Bayes factors
# allVLres=merge(VLcorswithpermp, bf, by="row.names", all=T)
# saveRDS(allVLres, "allVLres.rds")
# write.csv(allVLres, "allVLres.csv")





