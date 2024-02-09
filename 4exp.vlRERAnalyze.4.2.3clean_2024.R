##########################################
##### R program to analyze RERConverge Results #####
##########################################


####### Versions ####################

# 4.1.1 - Basic version of tree analysis
# 4.1.2 - Makes updates to add the Bayes factor
# 4.1.3 - Cleaning up analysis - protein Venn diagram
# 4.1.4 - Switch to add permulations p-value to results, BayesFactor
# 4.1.5 - Switch to add permulations p-value to results, BayesFactor
# 4.1.6 - Add in HyPhy Results
# 4.1.7 - Redo species predictions based on branch points
# 4.2.1 - Change R environment to TACIT
# 4.2.2_2024 - Minor Revision Figure Edits
# 4.2.3clean_2024 - Remove Legacy Code

##########################################
##### Functions #####
##########################################


###### Makes a matrix ########
makeMat <- function(x,y) {
	mat <- matrix(NA,length(x),length(y));
	rownames(mat) <- x;
	colnames(mat) <- y;

	return(mat);
}

#########################################################
##### Normalized by z-score #####
#########################################################

zNorm <- function(x) {
	#return(x - mean(x,na.rm=T)/sd(x,na.rm=T));
	return((x - mean(x,na.rm=T))/sd(x,na.rm=T));

}

#########################################################
##### Plot a tree for a gene, color based on vocal learner #####
#########################################################
#GeneName - string of gene symbol
#functionTr - Ape tree
#depends on filled speciesData2F with annotations

plotGene <- function(geneName,functionTr,functionAnnotF) {

	useTipColorsV <- c("black","red","gray");
	names(useTipColorsV) <- c("nonVocalLearner","vocalLearner","unsure");


	pdf(paste(outFd,"trees/testGene.",geneName,".1.pdf",sep=""),width=10,height=10);
		plot(functionTr,tip.color=useTipColorsV[functionAnnotF[functionTr$tip.label,"vlStatus"]]);

	dev.off();


}

#########################################################
##### Plot a tree for a gene, color based on vocal learner #####
#########################################################
#GeneName - string of gene symbol
#functionTr - Ape tree
#depends on filled speciesData2F with annotations

#plotGene2(curGene,curGeneTr,speciesData2F)

plotGene2 <- function(geneName,functionTr,functionAnnotF) {

	useTipColorsV <- c("black","red","gray");
	names(useTipColorsV) <- c("nonVocalLearner","vocalLearner","unsure");


	pdf(paste(outFd,"trees/testGene.",geneName,".noLabels.1.pdf",sep=""),width=10,height=10);
		#plot(functionTr,tip.color=useTipColorsV[functionAnnotF[functionTr$tip.label,"vlStatus"]],show.tip.label=F);
		plot(functionTr,show.tip.label=F);
		#plot(functionTr);
		#tiplabels(pch = 21, col='black')
		#tiplabels(pch = 21, col='black', fill=useTipColorsV[functionAnnotF[functionTr$tip.label,"vlStatus"]])
		tiplabels(pch = 19, col=useTipColorsV[functionAnnotF[functionTr$tip.label,"vlStatus"]])
		#tiplabels(pch=rep(21,length(functionTr$tip.label)), adj=c(0.6, 0.5),col='black')
		#tiplabels(pch=21, adj=c(0.6, 0.5),col=useTipColorsV[functionAnnotF[functionTr$tip.label,"vlStatus"]])

	#print(rep(21,length(functionTr$tip.label)));

	dev.off();

}


###################### Read in gene ontology - tables from EnrichR #######################

reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv,
              log_breaks(base = base),
              domain = c(1e-100, Inf))
}



##########################################
##### Load appropriate packages #####
##########################################

library(gplots);
library(ggplot2);

library(RColorBrewer);

library(ape);
library(VennDiagram);
library(DOSE);

library(scales)




#############################################################################################
##### Set directories #####
#############################################################################################

outFd <- "out416/"

int2classV <- c("integer","numeric","character");

fileSuffix = "4.1.7";


##########################################
##### Read in the relevant datasets #####
##########################################

########## Load RER Converge list for species left out ##########
#Original list without permulations

vocalLearningIdV <- c("allSpec","noHuman","noBat","noWhale","noSeal");
vocalLearningLooIdV <- c("noHuman","noBat","noWhale","noSeal");
vlCoefFL <- vector("list",0);

geneRowsV <- read.csv(paste("data/VL_","allSpec",".csv",sep=""),header=T)$Gene
geneRows2V <- unique(geneRowsV);

for(curVl in vocalLearningIdV) {

  curVlFn <- paste("data/VL_",curVl,".csv",sep="");
  vlCoefFL[[curVl]] <- read.csv(curVlFn,header=T);
  colnames(vlCoefFL[[curVl]])[1] <- "Gene";
  vlCoefFL[[curVl]] <- vlCoefFL[[curVl]][match(geneRows2V, vlCoefFL[[curVl]]$Gene),];
  rownames(vlCoefFL[[curVl]]) <- geneRows2V

}

head(vlCoefFL[[curVl]])


##### Load in Amanda's Bayes Factor Analysis #####

bayesFactorF <- read.csv("data/AK_allVLresFULL.csv");
rownames(bayesFactorF) <- bayesFactorF$Row.names
bayesFactorF$Gene <- as.vector(bayesFactorF$Row.names);


########## Read in Trees ##########

treeFn = "data/CombinedZoonomiaFilteredPrunedAATrees_VL.txt"

treeF <- read.delim(treeFn,row.names=1,stringsAsFactors=F);
treeV <- treeF[,1];
names(treeV) <- rownames(treeF);

treeL <- lapply(treeV,function(x) read.tree(text=x));
names(treeL) <- names(treeV);

fullTreeFn = "data/CombinedZoonomiaFilteredAATrees.txt"
fulltreeF <- read.delim(fullTreeFn,row.names=1,stringsAsFactors=F);
fulltreeV <- fulltreeF[,1];
names(fulltreeV) <- rownames(fulltreeF);

fulltreeL <- lapply(fulltreeV,function(x) read.tree(text=x));
names(fulltreeL) <- names(fulltreeV);


######## Read in data ##############

speciesDataFn <- "data/Zoonomia-Hiller_species_name_conversions.csv";
speciesDataF <- read.csv(speciesDataFn,stringsAsFactors=F);

speciesData2F <- speciesDataF[which(speciesDataF$TaxonGroup != ""),];
#sort(table(speciesData2F$Hiller_IDs..filtered.set.))
rownames(speciesData2F) <- speciesData2F$Hiller_IDs..full.set.;

table(speciesData2F$TaxonGroup);


######### Read in species sets #########

VL_AAinSpecies = c('vs_HLlepAme1', 'vs_HLoryCunCun4', 'vs_ochPri3', 'vs_HLcasCan3', 'vs_dipOrd2', 'vs_HLdipSte1', 'vs_HLperLonPac1', 'vs_HLfukDam2', 'vs_cavApe1', 'vs_cavPor3', 'vs_HLcavTsc1', 'vs_HLdolPat1', 'vs_chiLan1', 'vs_HLcteGun1', 'vs_HLcteSoc1', 'vs_HLdasPun1', 'vs_HLdinBra1', 'vs_HLhydHyd1', 'vs_HLhysCri1', 'vs_HLmyoCoy1', 'vs_octDeg1', 'vs_HLpetTyp1', 'vs_HLthrSwi1', 'vs_HLallBul1', 'vs_jacJac1', 'vs_HLzapHud1', 'vs_HLellLut1', 'vs_HLellTal1', 'vs_micOch1', 'vs_HLondZib1', 'vs_HLcriGri3', 'vs_mesAur1', 'vs_HLonyTor1', 'vs_HLperManBai2', 'vs_HLsigHis1', 'vs_HLacoCah1', 'vs_HLmerUng1', 'vs_HLpsaObe1', 'vs_HLmusPah1', 'vs_HLmusCar1', 'vs_mm39', 'vs_HLmusSpr1', 'vs_HLratNor7', 'vs_HLcriGam1', 'vs_nanGal1', 'vs_HLaplRuf1', 'vs_HLgliGli1', 'vs_HLgraMur1', 'vs_HLmusAve1', 'vs_speTri2', 'vs_HLmarMar1', 'vs_HLspeDau1', 'vs_HLxerIna1', 'vs_cerAty1', 'vs_HLcerNeg1', 'vs_chlSab2', 'vs_HLeryPat1', 'vs_HLmacFas6', 'vs_rheMac10', 'vs_macNem1', 'vs_manLeu1', 'vs_HLpapAnu5', 'vs_colAng1', 'vs_nasLar1', 'vs_HLpilTep2', 'vs_HLpygNem1', 'vs_rhiBie1', 'vs_HLrhiRox2', 'vs_HLsemEnt1', 'vs_gorGor6', 'REFERENCE', 'vs_ponAbe3', 'vs_HLnomLeu4', 'vs_aotNan1', 'vs_HLaloPal1', 'vs_HLateGeo1', 'vs_HLsagImp1', 'vs_HLcebAlb1', 'vs_cebCap1', 'vs_HLsaiBol1', 'vs_HLpitPit1', 'vs_HLdauMad1', 'vs_HLcheMed1', 'vs_micMur3', 'vs_HLmirCoq1', 'vs_proCoq1', 'vs_HLeulFla1', 'vs_HLeulFul1', 'vs_HLlemCat1', 'vs_otoGar3', 'vs_HLnycCou1', 'vs_tupChi1', 'vs_HLantAme1', 'vs_HLbeaHun1', 'vs_panHod1', 'vs_HLsaiTat1', 'vs_bisBis1', 'vs_HLbosInd2', 'vs_HLbosMut2', 'vs_bosTau9', 'vs_HLbubBub2', 'vs_HLammLer1', 'vs_HLcapAeg1', 'vs_HLhemHyl1', 'vs_HLoviAri5', 'vs_HLoviCan2', 'vs_HLelaDav1', 'vs_HLodoVir3', 'vs_HLranTar1', 'vs_HLgirTip1', 'vs_HLokaJoh2', 'vs_HLmosMos1', 'vs_HLtraJav1', 'vs_susScr11', 'vs_HLcatWag1', 'vs_HLcamFer3', 'vs_HLeubJap1', 'vs_orcOrc1', 'vs_HLturTru4', 'vs_HLlniGeo1', 'vs_HLdelLeu2', 'vs_HLailFul2', 'vs_HLlycPic3', 'vs_HLvulLag1', 'vs_HLspiGra1', 'vs_HLpteBra2', 'vs_HLmelCap1', 'vs_HLodoRos1', 'vs_HLzalCal1', 'vs_lepWed1', 'vs_HLmirAng2', 'vs_HLailMel2', 'vs_ursMar1', 'vs_HLcryFer2', 'vs_HLaciJub2', 'vs_HLfelNig1', 'vs_HLpumCon1', 'vs_HLpanOnc2', 'vs_HLpanPar1', 'vs_panTig1', 'vs_HLhelPar1', 'vs_HLmunMug1', 'vs_HLsurSur1', 'vs_HLhyaHya1', 'vs_HLparHer1', 'vs_pteAle1', 'vs_HLpteVam2', 'vs_HLrouAeg4', 'vs_HLhipArm1', 'vs_HLhipGal1', 'vs_HLmegLyr2', 'vs_HLtadBra1', 'vs_ptePar1', 'vs_HLcarPer3', 'vs_HLdesRot2', 'vs_HLrhiSin1', 'vs_HLmurAurFea1', 'vs_myoBra1', 'vs_HLmyoMyo6', 'vs_HLpipPip2', 'vs_eriEur2', 'vs_HLsolPar1', 'vs_sorAra2', 'vs_conCri1', 'vs_HLscaAqu1', 'vs_HLuroGra1', 'vs_equPrz1', 'vs_HLcerSimCot1', 'vs_cerSim1', 'vs_HLdicSum1', 'vs_HLdicBic1', 'vs_HLtapInd2', 'vs_HLtapTer1', 'vs_HLmanJav2', 'vs_HLmanPen2', 'vs_HLzipCav1')

VLnobats_AAinSpecies = c('vs_HLlepAme1', 'vs_HLoryCunCun4', 'vs_ochPri3', 'vs_HLcasCan3', 'vs_dipOrd2', 'vs_HLdipSte1', 'vs_HLperLonPac1', 'vs_HLfukDam2', 'vs_cavApe1', 'vs_cavPor3', 'vs_HLcavTsc1', 'vs_HLdolPat1', 'vs_chiLan1', 'vs_HLcteGun1', 'vs_HLcteSoc1', 'vs_HLdasPun1', 'vs_HLdinBra1', 'vs_HLhydHyd1', 'vs_HLhysCri1', 'vs_HLmyoCoy1', 'vs_octDeg1', 'vs_HLpetTyp1', 'vs_HLthrSwi1', 'vs_HLallBul1', 'vs_jacJac1', 'vs_HLzapHud1', 'vs_HLellLut1', 'vs_HLellTal1', 'vs_micOch1', 'vs_HLondZib1', 'vs_HLcriGri3', 'vs_mesAur1', 'vs_HLonyTor1', 'vs_HLperManBai2', 'vs_HLsigHis1', 'vs_HLacoCah1', 'vs_HLmerUng1', 'vs_HLpsaObe1', 'vs_HLmusPah1', 'vs_HLmusCar1', 'vs_mm39', 'vs_HLmusSpr1', 'vs_HLratNor7', 'vs_HLcriGam1', 'vs_nanGal1', 'vs_HLaplRuf1', 'vs_HLgliGli1', 'vs_HLgraMur1', 'vs_HLmusAve1', 'vs_speTri2', 'vs_HLmarMar1', 'vs_HLspeDau1', 'vs_HLxerIna1', 'vs_cerAty1', 'vs_HLcerNeg1', 'vs_chlSab2', 'vs_HLeryPat1', 'vs_HLmacFas6', 'vs_rheMac10', 'vs_macNem1', 'vs_manLeu1', 'vs_HLpapAnu5', 'vs_colAng1', 'vs_nasLar1', 'vs_HLpilTep2', 'vs_HLpygNem1', 'vs_rhiBie1', 'vs_HLrhiRox2', 'vs_HLsemEnt1', 'vs_gorGor6', 'REFERENCE', 'vs_ponAbe3', 'vs_HLnomLeu4', 'vs_aotNan1', 'vs_HLaloPal1', 'vs_HLateGeo1', 'vs_HLsagImp1', 'vs_HLcebAlb1', 'vs_cebCap1', 'vs_HLsaiBol1', 'vs_HLpitPit1', 'vs_HLdauMad1', 'vs_HLcheMed1', 'vs_micMur3', 'vs_HLmirCoq1', 'vs_proCoq1', 'vs_HLeulFla1', 'vs_HLeulFul1', 'vs_HLlemCat1', 'vs_otoGar3', 'vs_HLnycCou1', 'vs_tupChi1', 'vs_HLantAme1', 'vs_HLbeaHun1', 'vs_panHod1', 'vs_HLsaiTat1', 'vs_bisBis1', 'vs_HLbosInd2', 'vs_HLbosMut2', 'vs_bosTau9', 'vs_HLbubBub2', 'vs_HLammLer1', 'vs_HLcapAeg1', 'vs_HLhemHyl1', 'vs_HLoviAri5', 'vs_HLoviCan2', 'vs_HLelaDav1', 'vs_HLodoVir3', 'vs_HLranTar1', 'vs_HLgirTip1', 'vs_HLokaJoh2', 'vs_HLmosMos1', 'vs_HLtraJav1', 'vs_susScr11', 'vs_HLcatWag1', 'vs_HLcamFer3', 'vs_HLeubJap1', 'vs_orcOrc1', 'vs_HLturTru4', 'vs_HLlniGeo1', 'vs_HLdelLeu2', 'vs_HLailFul2', 'vs_HLlycPic3', 'vs_HLvulLag1', 'vs_HLspiGra1', 'vs_HLpteBra2', 'vs_HLmelCap1', 'vs_HLodoRos1', 'vs_HLzalCal1', 'vs_lepWed1', 'vs_HLmirAng2', 'vs_HLailMel2', 'vs_ursMar1', 'vs_HLcryFer2', 'vs_HLaciJub2', 'vs_HLfelNig1', 'vs_HLpumCon1', 'vs_HLpanOnc2', 'vs_HLpanPar1', 'vs_panTig1', 'vs_HLhelPar1', 'vs_HLmunMug1', 'vs_HLsurSur1', 'vs_HLhyaHya1', 'vs_HLparHer1', 'vs_eriEur2', 'vs_HLsolPar1', 'vs_sorAra2', 'vs_conCri1', 'vs_HLscaAqu1', 'vs_HLuroGra1', 'vs_equPrz1', 'vs_HLcerSimCot1', 'vs_cerSim1', 'vs_HLdicSum1', 'vs_HLdicBic1', 'vs_HLtapInd2', 'vs_HLtapTer1', 'vs_HLmanJav2', 'vs_HLmanPen2', 'vs_HLzipCav1')

VLnohuman_AAinSpecies = c('vs_HLlepAme1', 'vs_HLoryCunCun4', 'vs_ochPri3', 'vs_HLcasCan3', 'vs_dipOrd2', 'vs_HLdipSte1', 'vs_HLperLonPac1', 'vs_HLfukDam2', 'vs_cavApe1', 'vs_cavPor3', 'vs_HLcavTsc1', 'vs_HLdolPat1', 'vs_chiLan1', 'vs_HLcteGun1', 'vs_HLcteSoc1', 'vs_HLdasPun1', 'vs_HLdinBra1', 'vs_HLhydHyd1', 'vs_HLhysCri1', 'vs_HLmyoCoy1', 'vs_octDeg1', 'vs_HLpetTyp1', 'vs_HLthrSwi1', 'vs_HLallBul1', 'vs_jacJac1', 'vs_HLzapHud1', 'vs_HLellLut1', 'vs_HLellTal1', 'vs_micOch1', 'vs_HLondZib1', 'vs_HLcriGri3', 'vs_mesAur1', 'vs_HLonyTor1', 'vs_HLperManBai2', 'vs_HLsigHis1', 'vs_HLacoCah1', 'vs_HLmerUng1', 'vs_HLpsaObe1', 'vs_HLmusPah1', 'vs_HLmusCar1', 'vs_mm39', 'vs_HLmusSpr1', 'vs_HLratNor7', 'vs_HLcriGam1', 'vs_nanGal1', 'vs_HLaplRuf1', 'vs_HLgliGli1', 'vs_HLgraMur1', 'vs_HLmusAve1', 'vs_speTri2', 'vs_HLmarMar1', 'vs_HLspeDau1', 'vs_HLxerIna1', 'vs_cerAty1', 'vs_HLcerNeg1', 'vs_chlSab2', 'vs_HLeryPat1', 'vs_HLmacFas6', 'vs_rheMac10', 'vs_macNem1', 'vs_manLeu1', 'vs_HLpapAnu5', 'vs_colAng1', 'vs_nasLar1', 'vs_HLpilTep2', 'vs_HLpygNem1', 'vs_rhiBie1', 'vs_HLrhiRox2', 'vs_HLsemEnt1', 'vs_gorGor6', 'vs_ponAbe3', 'vs_HLnomLeu4', 'vs_aotNan1', 'vs_HLaloPal1', 'vs_HLateGeo1', 'vs_HLsagImp1', 'vs_HLcebAlb1', 'vs_cebCap1', 'vs_HLsaiBol1', 'vs_HLpitPit1', 'vs_HLdauMad1', 'vs_HLcheMed1', 'vs_micMur3', 'vs_HLmirCoq1', 'vs_proCoq1', 'vs_HLeulFla1', 'vs_HLeulFul1', 'vs_HLlemCat1', 'vs_otoGar3', 'vs_HLnycCou1', 'vs_tupChi1', 'vs_HLantAme1', 'vs_HLbeaHun1', 'vs_panHod1', 'vs_HLsaiTat1', 'vs_bisBis1', 'vs_HLbosInd2', 'vs_HLbosMut2', 'vs_bosTau9', 'vs_HLbubBub2', 'vs_HLammLer1', 'vs_HLcapAeg1', 'vs_HLhemHyl1', 'vs_HLoviAri5', 'vs_HLoviCan2', 'vs_HLelaDav1', 'vs_HLodoVir3', 'vs_HLranTar1', 'vs_HLgirTip1', 'vs_HLokaJoh2', 'vs_HLmosMos1', 'vs_HLtraJav1', 'vs_susScr11', 'vs_HLcatWag1', 'vs_HLcamFer3', 'vs_HLeubJap1', 'vs_orcOrc1', 'vs_HLturTru4', 'vs_HLlniGeo1', 'vs_HLdelLeu2', 'vs_HLailFul2', 'vs_HLlycPic3', 'vs_HLvulLag1', 'vs_HLspiGra1', 'vs_HLpteBra2', 'vs_HLmelCap1', 'vs_HLodoRos1', 'vs_HLzalCal1', 'vs_lepWed1', 'vs_HLmirAng2', 'vs_HLailMel2', 'vs_ursMar1', 'vs_HLcryFer2', 'vs_HLaciJub2', 'vs_HLfelNig1', 'vs_HLpumCon1', 'vs_HLpanOnc2', 'vs_HLpanPar1', 'vs_panTig1', 'vs_HLhelPar1', 'vs_HLmunMug1', 'vs_HLsurSur1', 'vs_HLhyaHya1', 'vs_HLparHer1', 'vs_pteAle1', 'vs_HLpteVam2', 'vs_HLrouAeg4', 'vs_HLhipArm1', 'vs_HLhipGal1', 'vs_HLmegLyr2', 'vs_HLtadBra1', 'vs_ptePar1', 'vs_HLcarPer3', 'vs_HLdesRot2', 'vs_HLrhiSin1', 'vs_HLmurAurFea1', 'vs_myoBra1', 'vs_HLmyoMyo6', 'vs_HLpipPip2', 'vs_eriEur2', 'vs_HLsolPar1', 'vs_sorAra2', 'vs_conCri1', 'vs_HLscaAqu1', 'vs_HLuroGra1', 'vs_equPrz1', 'vs_HLcerSimCot1', 'vs_cerSim1', 'vs_HLdicSum1', 'vs_HLdicBic1', 'vs_HLtapInd2', 'vs_HLtapTer1', 'vs_HLmanJav2', 'vs_HLmanPen2', 'vs_HLzipCav1')

VLnoseals_AAinSpecies = c('vs_HLlepAme1', 'vs_HLoryCunCun4', 'vs_ochPri3', 'vs_HLcasCan3', 'vs_dipOrd2', 'vs_HLdipSte1', 'vs_HLperLonPac1', 'vs_HLfukDam2', 'vs_cavApe1', 'vs_cavPor3', 'vs_HLcavTsc1', 'vs_HLdolPat1', 'vs_chiLan1', 'vs_HLcteGun1', 'vs_HLcteSoc1', 'vs_HLdasPun1', 'vs_HLdinBra1', 'vs_HLhydHyd1', 'vs_HLhysCri1', 'vs_HLmyoCoy1', 'vs_octDeg1', 'vs_HLpetTyp1', 'vs_HLthrSwi1', 'vs_HLallBul1', 'vs_jacJac1', 'vs_HLzapHud1', 'vs_HLellLut1', 'vs_HLellTal1', 'vs_micOch1', 'vs_HLondZib1', 'vs_HLcriGri3', 'vs_mesAur1', 'vs_HLonyTor1', 'vs_HLperManBai2', 'vs_HLsigHis1', 'vs_HLacoCah1', 'vs_HLmerUng1', 'vs_HLpsaObe1', 'vs_HLmusPah1', 'vs_HLmusCar1', 'vs_mm39', 'vs_HLmusSpr1', 'vs_HLratNor7', 'vs_HLcriGam1', 'vs_nanGal1', 'vs_HLaplRuf1', 'vs_HLgliGli1', 'vs_HLgraMur1', 'vs_HLmusAve1', 'vs_speTri2', 'vs_HLmarMar1', 'vs_HLspeDau1', 'vs_HLxerIna1', 'vs_cerAty1', 'vs_HLcerNeg1', 'vs_chlSab2', 'vs_HLeryPat1', 'vs_HLmacFas6', 'vs_rheMac10', 'vs_macNem1', 'vs_manLeu1', 'vs_HLpapAnu5', 'vs_colAng1', 'vs_nasLar1', 'vs_HLpilTep2', 'vs_HLpygNem1', 'vs_rhiBie1', 'vs_HLrhiRox2', 'vs_HLsemEnt1', 'vs_gorGor6', 'REFERENCE', 'vs_ponAbe3', 'vs_HLnomLeu4', 'vs_aotNan1', 'vs_HLaloPal1', 'vs_HLateGeo1', 'vs_HLsagImp1', 'vs_HLcebAlb1', 'vs_cebCap1', 'vs_HLsaiBol1', 'vs_HLpitPit1', 'vs_HLdauMad1', 'vs_HLcheMed1', 'vs_micMur3', 'vs_HLmirCoq1', 'vs_proCoq1', 'vs_HLeulFla1', 'vs_HLeulFul1', 'vs_HLlemCat1', 'vs_otoGar3', 'vs_HLnycCou1', 'vs_tupChi1', 'vs_HLantAme1', 'vs_HLbeaHun1', 'vs_panHod1', 'vs_HLsaiTat1', 'vs_bisBis1', 'vs_HLbosInd2', 'vs_HLbosMut2', 'vs_bosTau9', 'vs_HLbubBub2', 'vs_HLammLer1', 'vs_HLcapAeg1', 'vs_HLhemHyl1', 'vs_HLoviAri5', 'vs_HLoviCan2', 'vs_HLelaDav1', 'vs_HLodoVir3', 'vs_HLranTar1', 'vs_HLgirTip1', 'vs_HLokaJoh2', 'vs_HLmosMos1', 'vs_HLtraJav1', 'vs_susScr11', 'vs_HLcatWag1', 'vs_HLcamFer3', 'vs_HLeubJap1', 'vs_orcOrc1', 'vs_HLturTru4', 'vs_HLlniGeo1', 'vs_HLdelLeu2', 'vs_HLailFul2', 'vs_HLlycPic3', 'vs_HLvulLag1', 'vs_HLspiGra1', 'vs_HLpteBra2', 'vs_HLmelCap1', 'vs_HLailMel2', 'vs_ursMar1', 'vs_HLcryFer2', 'vs_HLaciJub2', 'vs_HLfelNig1', 'vs_HLpumCon1', 'vs_HLpanOnc2', 'vs_HLpanPar1', 'vs_panTig1', 'vs_HLhelPar1', 'vs_HLmunMug1', 'vs_HLsurSur1', 'vs_HLhyaHya1', 'vs_HLparHer1', 'vs_pteAle1', 'vs_HLpteVam2', 'vs_HLrouAeg4', 'vs_HLhipArm1', 'vs_HLhipGal1', 'vs_HLmegLyr2', 'vs_HLtadBra1', 'vs_ptePar1', 'vs_HLcarPer3', 'vs_HLdesRot2', 'vs_HLrhiSin1', 'vs_HLmurAurFea1', 'vs_myoBra1', 'vs_HLmyoMyo6', 'vs_HLpipPip2', 'vs_eriEur2', 'vs_HLsolPar1', 'vs_sorAra2', 'vs_conCri1', 'vs_HLscaAqu1', 'vs_HLuroGra1', 'vs_equPrz1', 'vs_HLcerSimCot1', 'vs_cerSim1', 'vs_HLdicSum1', 'vs_HLdicBic1', 'vs_HLtapInd2', 'vs_HLtapTer1', 'vs_HLmanJav2', 'vs_HLmanPen2', 'vs_HLzipCav1')

VLnowhales_AAinSpecies = c('vs_HLlepAme1', 'vs_HLoryCunCun4', 'vs_ochPri3', 'vs_HLcasCan3', 'vs_dipOrd2', 'vs_HLdipSte1', 'vs_HLperLonPac1', 'vs_HLfukDam2', 'vs_cavApe1', 'vs_cavPor3', 'vs_HLcavTsc1', 'vs_HLdolPat1', 'vs_chiLan1', 'vs_HLcteGun1', 'vs_HLcteSoc1', 'vs_HLdasPun1', 'vs_HLdinBra1', 'vs_HLhydHyd1', 'vs_HLhysCri1', 'vs_HLmyoCoy1', 'vs_octDeg1', 'vs_HLpetTyp1', 'vs_HLthrSwi1',
'vs_HLallBul1', 'vs_jacJac1', 'vs_HLzapHud1', 'vs_HLellLut1', 'vs_HLellTal1', 'vs_micOch1', 'vs_HLondZib1', 'vs_HLcriGri3', 'vs_mesAur1', 'vs_HLonyTor1', 'vs_HLperManBai2', 'vs_HLsigHis1', 'vs_HLacoCah1', 'vs_HLmerUng1', 'vs_HLpsaObe1', 'vs_HLmusPah1', 'vs_HLmusCar1', 'vs_mm39', 'vs_HLmusSpr1', 'vs_HLratNor7', 'vs_HLcriGam1', 'vs_nanGal1', 'vs_HLaplRuf1', 'vs_HLgliGli1', 'vs_HLgraMur1',
'vs_HLmusAve1', 'vs_speTri2', 'vs_HLmarMar1', 'vs_HLspeDau1', 'vs_HLxerIna1', 'vs_cerAty1', 'vs_HLcerNeg1', 'vs_chlSab2', 'vs_HLeryPat1', 'vs_HLmacFas6', 'vs_rheMac10', 'vs_macNem1', 'vs_manLeu1', 'vs_HLpapAnu5', 'vs_colAng1', 'vs_nasLar1', 'vs_HLpilTep2', 'vs_HLpygNem1', 'vs_rhiBie1', 'vs_HLrhiRox2', 'vs_HLsemEnt1', 'vs_gorGor6', 'REFERENCE', 'vs_ponAbe3', 'vs_HLnomLeu4', 'vs_aotNan1',
'vs_HLaloPal1', 'vs_HLateGeo1', 'vs_HLsagImp1', 'vs_HLcebAlb1', 'vs_cebCap1', 'vs_HLsaiBol1', 'vs_HLpitPit1', 'vs_HLdauMad1', 'vs_HLcheMed1', 'vs_micMur3', 'vs_HLmirCoq1', 'vs_proCoq1', 'vs_HLeulFla1', 'vs_HLeulFul1', 'vs_HLlemCat1', 'vs_otoGar3', 'vs_HLnycCou1', 'vs_tupChi1', 'vs_HLantAme1', 'vs_HLbeaHun1', 'vs_panHod1', 'vs_HLsaiTat1', 'vs_bisBis1', 'vs_HLbosInd2', 'vs_HLbosMut2',
'vs_bosTau9', 'vs_HLbubBub2', 'vs_HLammLer1', 'vs_HLcapAeg1', 'vs_HLhemHyl1', 'vs_HLoviAri5', 'vs_HLoviCan2', 'vs_HLelaDav1', 'vs_HLodoVir3', 'vs_HLranTar1', 'vs_HLgirTip1', 'vs_HLokaJoh2', 'vs_HLmosMos1', 'vs_HLtraJav1', 'vs_susScr11', 'vs_HLcatWag1', 'vs_HLcamFer3', 'vs_HLailFul2', 'vs_HLlycPic3', 'vs_HLvulLag1', 'vs_HLspiGra1', 'vs_HLpteBra2', 'vs_HLmelCap1', 'vs_HLodoRos1', 'vs_HLzalCal1',
'vs_lepWed1', 'vs_HLmirAng2', 'vs_HLailMel2', 'vs_ursMar1', 'vs_HLcryFer2', 'vs_HLaciJub2', 'vs_HLfelNig1', 'vs_HLpumCon1', 'vs_HLpanOnc2', 'vs_HLpanPar1', 'vs_panTig1', 'vs_HLhelPar1', 'vs_HLmunMug1', 'vs_HLsurSur1', 'vs_HLhyaHya1', 'vs_HLparHer1', 'vs_pteAle1', 'vs_HLpteVam2', 'vs_HLrouAeg4', 'vs_HLhipArm1', 'vs_HLhipGal1', 'vs_HLmegLyr2', 'vs_HLtadBra1', 'vs_ptePar1', 'vs_HLcarPer3',
'vs_HLdesRot2', 'vs_HLrhiSin1', 'vs_HLmurAurFea1', 'vs_myoBra1', 'vs_HLmyoMyo6', 'vs_HLpipPip2', 'vs_eriEur2', 'vs_HLsolPar1', 'vs_sorAra2', 'vs_conCri1', 'vs_HLscaAqu1', 'vs_HLuroGra1', 'vs_equPrz1', 'vs_HLcerSimCot1', 'vs_cerSim1', 'vs_HLdicSum1', 'vs_HLdicBic1', 'vs_HLtapInd2', 'vs_HLtapTer1', 'vs_HLmanJav2', 'vs_HLmanPen2');

batsV <- setdiff(VL_AAinSpecies,VLnobats_AAinSpecies);
humanV <- setdiff(VL_AAinSpecies,VLnohuman_AAinSpecies);
whalesV <- setdiff(VL_AAinSpecies,VLnowhales_AAinSpecies);
sealsV <- setdiff(VL_AAinSpecies,VLnoseals_AAinSpecies);

allVl <- c(batsV,humanV,whalesV,sealsV)

speciesData2F$vlStatus <- "nonVocalLearner"
speciesData2F[allVl,"vlStatus"] <- "vocalLearner";
allVlCladesV <- rownames(speciesData2F)[which(!is.na(match(speciesData2F$TaxonGroup,c("Bat","Cetaceans","Pinniped"))))]
speciesData2F[setdiff(allVlCladesV,allVl),"vlStatus"] <- "unsure";
table(speciesData2F$vlStatus);


#Print out specices Data
#write.csv(speciesData2F,quote=F,row.names=F,col.names=T,file=paste(outFd,"speciesAnnotations.",fileSuffix,".1.csv",sep=""))


########## Quick test of a gene - Figure S2A, Figure S2B ##########

if(F) {


	curGene <- "GRM8"; #Best Accelerated
	curGene <- "CENPC" #Best conserved?


	curGeneTr <- fulltreeL[[curGene]]; #The tree for the current gene
	curVlTr <- drop.tip(curGeneTr,setdiff(curGeneTr$tip.label,allVl),trim.internal=T);  #tree with only species of interest
	curBatTr <- drop.tip(curGeneTr,setdiff(curGeneTr$tip.label,batsV),trim.internal=T);  #tree with only species of interest
	curNonvlTr <- drop.tip(curGeneTr,allVl,trim.internal=T); #tree without species of interest or vocal learners


	plotGene(curGene,curGeneTr,speciesData2F);

	plotGene2(curGene,curGeneTr,speciesData2F)


	if(F) {

		useTipColorsV <- c("black","red","gray");
		names(useTipColorsV) <- c("nonVocalLearner","vocalLearner","unsure");


		pdf(paste(outFd,"trees/testGene.",curGene,".1.pdf",sep=""),width=10,height=10);
			plot(curGeneTr,tip.color=useTipColorsV[speciesData2F[curGeneTr$tip.label,"vlStatus"]]);

		dev.off();

	}

}


########## Histogram of p-values comparing methods ##########

curPplotF <- bayesFactorF;
curPplotF$pvalueMethod <- "Raw";
curPplotPermF <- bayesFactorF;
curPplotPermF$P <- curPplotPermF$permP;
curPplotPermF$pvalueMethod <- "Permulations";
curPplotF <- rbind(curPplotF,curPplotPermF);

if(F) {

	pdf(paste(outFd,"pvlaueCompare",fileSuffix,".2.pdf",sep=""),width=5,height=4);

		curP <- ggplot(curPplotF, aes(x=P, fill = pvalueMethod)) +
			geom_histogram(color="black",position='identity',alpha=0.5) +
		#	geom_density(color=pvalueMethod, alpha = 0) +
			#geom_density(adjust=0.5,fill="black") +
			#geom_bar(stat="identity") +
			#scale_x_continuous(trans='log2') +
			theme_bw();
			#theme(axis.text.x = element_text(angle = 90, hjust = 1));
			print(curP);
	dev.off();

}


####### Print out Bayes Factor Results ##########

bayesFactorPrintF <- bayesFactorF;
bayesFactorPrintF$X <- NULL;
bayesFactorPrintF$Row.names <- NULL;
#write.csv(bayesFactorPrintF,quote=F,row.names=F,col.names=T,file=paste(outFd,"vocalLearningRERConvergeMaster.",fileSuffix,".2.csv",sep=""))



if(F) {
	loadedDataWsFn <- paste("loadedData.4.2.1.1.ws",sep=""); #The basic normalized data, fixed names bug
	#save.image(loadedDataWsFn);
	load(loadedDataWsFn);


}

##########################################
##### Create a new list of vocal learning genes with Bayes Factor and Permulations #####
##########################################

#Predict Which Species are Vocal Learners based on scores


########### Functions to traverse tree #############

### Get the parent of a node ###
#Based on its number
#getParent(curGeneTr,149)
getParent <- function(funPhylo,funChild){
	funRow <- which(funPhylo$edge[,2] == funChild);
	return(funPhylo$edge[funRow,1]);
}

### Get all tips for an internal node ###
#Internal node and input is based on number
#Tips are based on name
getTips <- function(funPhylo,funParent){
	if(funParent <= length(funPhylo$tip.label)) {
		return(funPhylo$tip.label[funParent]);
	} else {
		furChildrenRowV <- funPhylo$edge[which(funPhylo$edge[,1] == funParent),2];
		funTipsV <- c();
		for(funChild in furChildrenRowV) {
			funTipsV <- c(funTipsV,getTips(funPhylo,funChild));
		}
		return(funTipsV);
	}
}

if(F) {
	tmpParent <- 298;
	tmpChildrenRowV <- which(curGeneTr$edge[,1] == tmpParent);
}


########### Aggregrate the evolutionary rates for the tips #############


curGene <- "NCAN"
#curGene <- "GRM8"
#curGene <- "TSHZ3";

adjThresh <- 0.01;
adjPermThresh <- 0.01;

#Set which species result to use as the basis for the set of genes to use
curRerResF <- bayesFactorF;
curRerResTopF <- curRerResF[which(curRerResF$p.adj <= adjThresh & curRerResF$permP.adj <= adjPermThresh),];

#Matrix where the rows are the top genes and the columns represent the edge to the tip of each species
tipMatrixZM <- makeMat(curRerResTopF$Gene,speciesData2F$Hiller_IDs..full.set.);
tipMatrixZvsnonM <- tipMatrixZM; #Rather than normalize again all other species, normalize relative to non-vocal learners
tipMatrixRawM <- tipMatrixZM; #Raw tip values, not just RERs

#Loop through all gene trees for the top genes based on all vocal learners
for(curGene in curRerResTopF$Gene) {
	#curGeneTr <- treeL[[curGene]]; #The tree for the current gene
	curGeneTr <- fulltreeL[[curGene]]; #The tree for the current gene

	#Loop through all tips and get tip edges
	tipEdgesV <- c();
	for(curTipNum in 1:length(curGeneTr$tip.label)) {
		curTip <- curGeneTr$tip.label[curTipNum]
		curParent <- getParent(curGeneTr,curTipNum)
		curRelativesV <- setdiff(getTips(curGeneTr,curParent),curTip);
		curEdgeName <- paste(curTip,"-",paste(curRelativesV,collapse=","),sep="");
		tipEdgesV[curEdgeName] <- curGeneTr$edge.length[which(curGeneTr$edge[,2] == curTipNum)];
	}

	names(tipEdgesV) <- curGeneTr$tip.label;
	tipEdges2V <- tipEdgesV[which(!is.na(match(names(tipEdgesV),speciesData2F$Hiller_IDs..full.set.)))]
	curSpeciesV <- names(tipEdges2V);
	tipVlZ2V <- zNorm(tipEdges2V);

	curNonVlSpeciesV <- intersect(speciesData2F[which(speciesData2F$vlStatus == "nonVocalLearner"),"Hiller_IDs..filtered.set."],curGeneTr$tip.label);
	tipVlZvsnon2V <- tipVlZ2V - mean(tipVlZ2V[curNonVlSpeciesV],na.rm=T)

	tipMatrixRawM[curGene,curSpeciesV] <- tipEdges2V; #Matrix represent the  evolutionary rate for each tip branch
	tipMatrixZM[curGene,curSpeciesV] <- tipVlZ2V; #Matrix represent the  evolutionary rate for each branch, z normalized
	tipMatrixZvsnonM[curGene,curSpeciesV] <- tipVlZvsnon2V; #Matrix represent the evolution rate of that tip relative to vocal non-learners lengths


	}

}

########### Calculate sets of vocal learning species across each clade #######################

speciesSetV <- c("bats","whales","seals","human");
speciesSetsL <- vector("list",0);
speciesSetsL[["bats"]] <- batsV;
speciesSetsL[["whales"]] <- whalesV;
speciesSetsL[["seals"]] <- sealsV;
speciesSetsL[["human"]] <- humanV;


######## Explore all human RERs ##############

humanFullRerF <- data.frame(gene=names(fulltreeL),edgeName=rep(NA,length(fulltreeL)),edgeLen=rep(NA,length(fulltreeL)),
														parentName=rep(NA,length(fulltreeL)),parentRate=rep(NA,length(fulltreeL)),stringsAsFactors=F);

rownames(humanFullRerF) <- names(fulltreeL)

#Loop through all gene trees for the top genes based on all vocal learners
for(curGene in names(fulltreeL)) {
	curGeneTr <- fulltreeL[[curGene]]; #The tree for the current gene

	if("REFERENCE" %in% curGeneTr$tip.label) {
		curTipNum <- which(curGeneTr$tip.label=="REFERENCE");
		curParent <- getParent(curGeneTr,curTipNum);
		curRelativesV <- setdiff(getTips(curGeneTr,curParent),curTip);
		curEdgeName <- paste(curTip,"-",paste(curRelativesV,collapse=","),sep="");

		#Get the length of the edge from great apes to human
		humanFullRerF[curGene,"edgeName"] <- curEdgeName;
		humanFullRerF[curGene,"edgeLen"] <- curGeneTr$edge.length[which(curGeneTr$edge[,2] == curTipNum)];

		#Get the length of the edge from great apes to human
		curParentParent <- getParent(curGeneTr,curParent);
		curParentRate <- curGeneTr$edge.length[which(curGeneTr$edge[,2] == curParent)];
		curParentParentNames <- paste(getTips(curGeneTr,curParentParent),collapse=",",sep="")
		humanFullRerF[curGene,"parentName"] <- curParentParentNames;
		humanFullRerF[curGene,"parentRate"] <- curParentRate
	}

}

#These are the ones where there are no recent ancestors missing
reliableRERV <- humanFullRerF[which(humanFullRerF$edgeName == "vs_HLmanJav2-REFERENCE,vs_panPan3,vs_panTro6" & humanFullRerF$parentName == "vs_gorGor6,REFERENCE,vs_panPan3,vs_panTro6"),"gene"]

#Human rate relative to parent
humanFullRerF$RER <- humanFullRerF$edgeLen - humanFullRerF$parentRate

positiveHgRerThresh <- 0.01;
negativeHgRerThresh <- -0.005;

humanAccelerateV <- intersect(humanFullRerF[which(humanFullRerF$RER >= positiveHgRerThresh),"gene"],reliableRERV);
humanConservedV <- intersect(humanFullRerF[which(humanFullRerF$RER <= negativeHgRerThresh),"gene"],reliableRERV);

### Incorporate Human with Bayes Factor ####

bayesFactor2F <- bayesFactorF[rownames(curRerResTopF),]; #Rownames are same as curRerResTopF
bayesFactorThresh <- 5;

bayesFactor2F["SETBP1",];

bayesFactorAccel2F <- bayesFactor2F[which(bayesFactor2F$Rho > 0),];
proteinVennAccelL <- vector("list",0);
proteinVennAccelL[["Bats"]] <- bayesFactorAccel2F[which(bayesFactorAccel2F$bat > bayesFactorThresh),"Gene"];
proteinVennAccelL[["Cetaceans"]] <- bayesFactorAccel2F[which(bayesFactorAccel2F$whales > bayesFactorThresh),"Gene"];
proteinVennAccelL[["Pinnipeds"]] <- bayesFactorAccel2F[which(bayesFactorAccel2F$seals > bayesFactorThresh),"Gene"];
proteinVennAccelL[["Human"]] <- intersect(humanAccelerateV,bayesFactorAccel2F$Gene);
lapply(proteinVennAccelL,length);

proteinAccelCladesV <- sort(table(unlist(proteinVennAccelL)));
proteinAccelCladesV
sort(table(unlist(proteinVennAccelL)));


bayesFactorCons2F <- bayesFactor2F[which(bayesFactor2F$Rho < 0),];
proteinVennConsL <- vector("list",0);
proteinVennConsL[["Bats"]] <- bayesFactorCons2F[which(bayesFactorCons2F$bat > bayesFactorThresh),"Gene"];
proteinVennConsL[["Cetaceans"]] <- bayesFactorCons2F[which(bayesFactorCons2F$whales > bayesFactorThresh),"Gene"];
proteinVennConsL[["Pinnipeds"]] <- bayesFactorCons2F[which(bayesFactorCons2F$seals > bayesFactorThresh),"Gene"];
proteinVennConsL[["Human"]] <- intersect(humanConservedV,bayesFactorCons2F$Gene);
lapply(proteinVennConsL,length);


sort(table(unlist(proteinVennConsL)));
proteinConsCladesV <- sort(table(unlist(proteinVennConsL)))
proteinConsCladesV

nonBayesfactSetV <- setdiff(rownames(bayesFactor2F),c(rownames(proteinAccelCladesV),rownames(proteinConsCladesV)));

fullAccelSetF <- bayesFactor2F[which(bayesFactor2F$Rho > 0),];
fullConsSetF <- bayesFactor2F[which(bayesFactor2F$Rho < 0),];


####### Plot Venn Diagrams - Figure 1B, 1C #########
if(F) {

	currVennL <- proteinVennConsL
	currVennFn <- paste(outFd,"venn/vennDiagram.cons.",fileSuffix,".2rename.png",sep="");
	venn.diagram(x=currVennL,
			height = 1500, width = 2000,
			category.names = names(currVennL),
			fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
			filename = currVennFn,
			imagetype ="png"
	);

	currVennL <- proteinVennAccelL
	currVennFn <- paste(outFd,"venn/vennDiagram.accel.",fileSuffix,".2rename.png",sep="");
	venn.diagram(x=currVennL,
			category.names = names(currVennL),
			height = 1500, width = 2000,
			fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
			filename = currVennFn,
			imagetype ="png"
	);

}



##########################################
##### Read in HyPhy Results #####
##########################################

########## Read in Results of HyPhy ##########

#1 ID: gene
#2 Sequences: number of sequences in the alignment
#3 Sites: length of the alignment
#4 FG q-value: q-value (Benjamini-Hochberg) of the BUSTED test on foreground branches
#5 BG q-value: q-value (Benjamini-Hochberg) of the BUSTED test on background branches
#6 DIFF q-value: q-value (Benjamini-Hochberg) for the difference between FG and BG
#7 S.Sites (FG): Sites that are showing support for positive selection in the Foreground branches with evidence ratio of 100 or higher
#8 S.Sites (BG): Sites that are showing support for positive selection in the Background branches with evidence ratio of 100 or higher
#9 L (FG): tree length  foreground
#10 L (BG): tree length background
#11 RELAX q-value:
#12 RELAX K:  the foreground is modeled as (omega)^K so it matters if K<1 or K>1
#13 -18 Foreground parameters for BUSTED-PH
#19 PH: Selection associated with phenotype
#20 TYPE: Binary summary of 4,5,6


hyphyCatV <- c("cons","accel","rand");
useHyphy <- "srv";
#useHyphy <- "nosrv";

hyphyResL <- vector("list",0);

plotHyphyF <- c();

for(curCat in hyphyCatV) {
	curHyphyFn <- paste("data/HyPhy/hyPhy.vocal.",curCat,".",useHyphy,".csv",sep="")
	hyphyResL[[curCat]] <- read.csv(curHyphyFn,header=T,stringsAsFactors=F)
	hyphyResL[[curCat]]$geneSet <- curCat;
	plotHyphyF <- rbind(plotHyphyF,hyphyResL[[curCat]]);
}
hyphyResL[[curCat]][1:5,];

unlist(lapply(hyphyResL,function(x) length(which(x$TYPE == "C000"))/dim(x)[1]))
unlist(lapply(hyphyResL,function(x) length(which(x$TYPE == "C101"))/dim(x)[1]))
unlist(lapply(hyphyResL,function(x) length(which(x$TYPE == "C011"))/dim(x)[1]))

plotHyphyF$geneSetT <- factor(plotHyphyF$geneSet,levels=c("rand","cons","accel"));

quantile(log2(plotHyphyF$RELAX.K),c(0:10)*.1)

if(F) { #Comopare HyPhy Results to RERConverge Results

	pdf(paste(outFd,"hyPhy/genesetCompare.k",fileSuffix,".1.pdf",sep=""),width=5,height=3);
		curP <- ggplot(plotHyphyF, aes(x = log2(RELAX.K), fill=geneSetT)) +
		geom_density(adjust=1,alpha=0.5) +
		#geom_smooth(method = "lm", se=FALSE, color="blue") +
		theme_bw() +
		scale_x_continuous(limits = c(-1, 1)) +
		scale_fill_manual(values=c("gray","purple","green"));
		#theme(axis.text.x = element_text(angle = 90,hjust=1)) +
		#scale_fill_manual(values=c("gray","#997570","red"));
		#scale_fill_gradient2(low="red",	mid = "grey50",high = "blue",	midpoint = 0.5,	na.value = "grey50", guide = "colourbar");
		#geom_line(data=curKnlData2F, aes(x=age, y=pred), colour="red", size=1.5);
		print(curP);
	dev.off();

}

wilcox.test(x=log2(hyphyResL[["cons"]]$RELAX.K),y=log2(hyphyResL[["rand"]]$RELAX.K),paired=F,alternative="greater");
wilcox.test(x=log2(hyphyResL[["accel"]]$RELAX.K),y=log2(hyphyResL[["rand"]]$RELAX.K),paired=F,alternative="less")


newConservedV <- hyphyResL[["cons"]][which(hyphyResL[["cons"]]$RELAX.K > 1 & hyphyResL[["cons"]]$RELAX.q.value < 0.05),"ID"]
newConservedSymV <- unique(unlist(lapply(strsplit(newConservedV,"\\."),function(x) x[2])));

#write.table(newConservedSymV,file=paste(outFd,"hyPhy/newConservedSym.RELAX",fileSuffix,".1.txt",sep=""), sep="\n",quote=F,row.names=F,col.names=F)

newAccelV <- hyphyResL[["accel"]][which(hyphyResL[["accel"]]$RELAX.K < 1 & hyphyResL[["accel"]]$RELAX.q.value < 0.05),"ID"]
newAccelSymV <- unique(unlist(lapply(strsplit(newAccelV,"\\."),function(x) x[2])));

#write.table(newAccelSymV,file=paste(outFd,"hyPhy/newAccelSym.RELAX",fileSuffix,".1.txt",sep=""), sep="\n",quote=F,row.names=F,col.names=F)

######## HyPhy Augmented Gene Ontology from EnrichR ##############

image.scale <- function(z, zlim, col = heat.colors(12),
breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
 if(!missing(breaks)){
  if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
 }
 if(missing(breaks) & !missing(zlim)){
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 if(missing(breaks) & missing(zlim)){
  zlim <- range(z, na.rm=TRUE)
  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 poly <- vector(mode="list", length(col))
 for(i in seq(poly)){
  poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
 }
 xaxt <- ifelse(horiz, "s", "n")
 yaxt <- ifelse(horiz, "n", "s")
 if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
 if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
 if(missing(xlim)) xlim=XLIM
 if(missing(ylim)) ylim=YLIM
 plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)
 for(i in seq(poly)){
  if(horiz){
   polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
  }
  if(!horiz){
   polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
  }
 }
}


goResFn <- c("GO_Biological_Process_2023_HyPhy.Conserved.txt","GO_Biological_Process_2023_HyPhyAccel.txt","Human_Phenotype_Ontology_HyPhy.Conserved.txt","Human_Phenotype_Ontology_HyPhyAccel.txt");
names(goResFn) <- c("Conserved_GoBP","Accelerated_GoBP","Conserved_HumanPhenotype","Accelerated_HumanPhenotype");
goResL <- vector("list",0);

curGo <- "Conserved_GoBP";

for(curGO in names(goResFn)) {


	curGeneGoF <- read.delim(paste("data/",goResFn[curGO],sep=""),sep="\t",header=T,stringsAsFactors=F);
	curGeneGenesL <- strsplit(curGeneGoF[,"Genes"],";");
	curGeneGoF$termT <- factor(curGeneGoF$Term,levels=rev(curGeneGoF$Term))
	curGeneGoF$logAdjustedP <- -log10(curGeneGoF$Adjusted.P.value);
	curGeneGoF$numGenes <- unlist(lapply(curGeneGenesL,length));

	curGeneGo2F <- curGeneGoF[which(curGeneGoF$Adjusted.P.value < 0.05 & curGeneGoF$numGenes > 4),]

	if(length(curGeneGo2F$numGenes) >= 7) {
		curGeneGo2F <- curGeneGo2F[1:6,];
	}

	goResL[[curGO]] <- curGeneGo2F;

}


allGoResF <- do.call(rbind,goResL);
allGoResF$adjOddsRatio <- allGoResF$Odds.Ratio/15; #manually adjust scale to 10-5
tmpColorM <- colorRamp(c("blue","red"),space= "rgb")(allGoResF$adjOddsRatio)
allGoResF$plotColor <- apply(tmpColorM,1,function(x) rgb(x[1],x[2],x[3],maxColorValue=255));

for(curGO in names(goResL)) {
	goResL[[curGO]]$plotColor <- allGoResF[match(goResL[[curGO]]$Term,allGoResF$Term),"plotColor"];

}

if(F) { #Plot color legend

	pdf(paste(outFd,"goPlot.legend.2024.4oddratio.pdf",sep=""),width=12,height=3);
		image.scale(0:50, col=colorRampPalette(c("blue", "red"))(50))
		box();
	dev.off();

}


curGO <- "Conserved_GoBP";
curGO <- "Conserved_HumanPhenotype";


### Plot gene ontology results - Figure 1D, E ####

for(curGO in names(goResL)) {

	#print(curGO)

	curGeneGo2F <- goResL[[curGO]];
	#print(curGeneGo2F);

	#pdf(paste(outFd,"goPlot.",curGO,".",fileSuffix,".2024.4hyphy.pdf",sep=""),width=12,height=2);
	pdf(paste(outFd,"goPlot.",curGO,".",fileSuffix,".2024.4hyphy.pdf",sep=""),width=12,height=1.5);

		ggplot(data = curGeneGo2F, aes(x = Adjusted.P.value, y = termT, color=termT)) +
		  geom_point(size=3) +
			#scale_x_log10() +
			#scale_x_reverse() +
			scale_x_continuous(trans=reverselog_trans(base=10)) +
		  scale_colour_manual(values=rev(curGeneGo2F$plotColor),guide = "none") +
			#scale_colour_manual(values=rev(curGeneGo2F$plotColor)) +
		  theme_bw() +
		  ylab("") +
		  xlab("") +
			ggtitle(paste("EnrichR - ",curGO,sep=""));

	dev.off();

}















#
