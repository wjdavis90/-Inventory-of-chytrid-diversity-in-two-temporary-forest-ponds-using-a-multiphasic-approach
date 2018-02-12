## OTU alpha diversity
The first set of scripts analyzes the alpha diversity indices calculated from the OTU data.

First, the data was read into R.

`>TVP<-read.delim("Alpha_indices_151026.txt", header= TRUE)`

The ANOVAs were run using the built-in aov in R. The models were created using the example below as a template. 

`>TVP.aov1=aov([alpha_diversity_index]~Location*Month, data=TVP)`

The results were then viewed and converted into a text file.

`>summary(TVP.aov1)`
`> capture.output(summary(TVP.aov1), file="all_otu_anova_[alpha_diversity_index].txt")`

The same was done with the alpha diversity indices calculated with only chytrid OTUs.

`>chytrid<-read.csv("chytrid_alpha_indices_151026.csv")`

`>chytrid.[alpha_diversity_index.aov=aov([alpha_diversity_index]~Location*Month,data=chytrid)`

`>summary(chytrid.[alpha_diversity_index].aov)`

`>capture.output(summary(chytrid.[alpha_diversity_index].aov), file = “chytrid_otu_anova_[alpha_diversity_index].txt” )`

# Presence/Absence data
## Alpha diversity & Individual taxa
The following scripts require the packages lme4 and picante.

Read the data file into R

`> TVP <- read.delim("TVP_data.txt", header=TRUE)`

This first set of scripts is for the mixed effect ANOVA analyses of the individual species using #lme4. In the scripts that follow, [species_name] is a place holder. First, a null model of each species was generated using only the random effect block.

`>TVP.[species_name].null <- lmer([species_name] ~ (1|Pond), data=TVP,  REML=FALSE)`

Next, the fixed effect variables are added to the model one at a time.

`>TVP. [species_name].month <- lmer([species_name] ~ Month + (1|Pond), data=TVP, REML=FALSE)`

`>TVP. [species_name].location <- lmer([species_name] ~ Location + (1|Pond), data=TVP, REML=FALSE)`

`>TVP. [species_name].interaction <- lmer([species_name] ~ Month*Location + (1|Pond), data=TVP, REML=FALSE)`

Each model with a fixed effect was then tested against the null model to determine the #significance of the fixed effect.

`>anova(TVP. [species_name].null, TVP. [species_name].month)`

`>anova(TVP. [species_name].null, TVP. [species_name].location)`

`>anova(TVP. [species_name].null, TVP. [species_name].interaction)`

The same scripts were run to build and test models of species richness (Sr).

`>TVP.Sr.null <- lmer(Species_richness ~ (1|Pond), data=TVP, REML=FALSE)`

`>TVP.Sr.month <- lmer(Species_richness ~ Month + (1|Pond), data=TVP, REML=FALSE)`

`>TVP.Sr.location <- lmer(Species_richness ~ Location + (1|Pond), data=TVP, REML=FALSE)`

`>TVP.Sr.interaction <- lmer(Species_richness ~ Month*Location + (1|Pond), data=TVP, REML=FALSE)`

Phylogenetic diversity (PD) was calculated with the package picante. The raw presence/absence data was loaded into R.

`>TVP.1<-read.delim("TVP_Data_151204.txt", header=TRUE, row.names=1)`

The tree was loaded into R.

`>Tree<-read.nexus(file="picante_3.tre")`

A reduced dataset with only the taxa was generated.

`>myvars.1<-names(TVP.1) %in% c("Pond","Month","Location","Species_richness")`

`>TVP.2<-TVP.1[!myvars.1]`

The taxa names in the tree need to match the taxa names in the dataset. Some of the taxa in the dataset needed to be renamed.

`>library(plyr)`

`TVP.3<-rename(TVP.2, c("Powellomycetacase.sp_"="Powellomycetaceae_sp", "Lobulomyces_pocularus" ="Lobulomyces_poculatus", "Rhizoclosmatum_globosum"="Rhizoclosmatium_globosum", "Rhizoclosmatum_aurantiacum" ="Rhizoclosmatium_aurantiacum", "Odontochytrium_millerii" ="Odontochytrium_milleri"))`

The dataset was reduced to only those taxa included in the tree.

`>myvars.2<-c("Boothiomyces_macroporosum","Cladochytrium_replicatum","Chytriomyces_hyalinus", "Fayochytriomyces_spinosus", "Entophlyctis_sp", "Lobulomyces_poculatus", "Lobulomyces_angularis", "Nowakowskiella_sp", "Rhizoclosmatium_globosum", "Rhizoclosmatium_aurantiacum", "Odontochytrium_milleri", "Powellomycetaceae_sp")`

`>TVP.4<-TVP.3[myvars.2]

To be sure there was an exact match of taxa in the tree to taxa in the dataset, the tree was pruned.

`>pruneTree<-prune.sample(TVP.4, Tree)` 

PD was calculated and stored as an object.
`>PD <- pd(TVP.4, pruneTree, include.root=TRUE)`

An object with the variables Pond, Month, and Location was created.
`>TVP.5<-TVP.1[myvars.1]`

Models were built and tested as above using lme4.

`>TVP.PD.null <- lmer(PD$PD ~ (1|TVP.5$Pond), data=TVP.5, REML=FALSE)`

`>TVP.PD.month <- lmer(PD$PD ~ TVP.5$Month + (1|TVP.5$Pond), data=TVP.5, REML=FALSE)`

`>TVP.PD.Location <- lmer(PD$PD ~ TVP.5$Location + (1|TVP.5$Pond), data=TVP.5, REML=FALSE)`

`>TVP.PD.interaction <- lmer(PD$PD ~ TVP.5$Month*TVP.5$Location + (1|TVP.5$Pond), data=TVP.5, REML=FALSE)`

`>anova(TVP.PD.null, TVP.PD.month)`

`>anova(TVP.PD.null, TVP.PD.Location)`

`>anova(TVP.PD.null, TVP.PD.interaction)`

Results were stored as text files using the capture.output function.

## Beta diversity

