source('~/Documents/GitHub/R-General/.Rprofile')
source('~/Documents/GitHub/R-LampPrimerAnalyzer/LAMPUtils.R')

library(zoo)
library(primer3)
library(seqinr)
# library(reticulate)
# source /Users/jaywarrick/.pyenv/versions/3.10.2/envs/PyRStudio/bin/activate
# pip3 install seqfold
# pip3 install viennarna
# use_python("")

# Calculate Different Score Metrics
scoreAlgorithm <- function(primerSet, DimerTarget=0, HpTarget=0, HpTmTarget=55, HpDBTarget=0, power=2)
{
	pSet <- copy(primerSet)
	pSet[is.na(HpDeltaG), HpDeltaG:=0]
	primerScore <- with(pSet, {
		DimerScore <- sum(sapply(DimerDeltaG[Primer %in% c('FIP','BIP')]-DimerTarget, min, 0)^power)
		HpPrimerScore <- sum(sapply(HpDeltaG[Primer %in% c('FIP','BIP')]-HpTarget, min, 0)^power)
		HpDBScore <- sum(sapply(HpDeltaG[Primer %in% c('DBF','DBB')]-HpDBTarget, min, 0)^power)
		TmScore1 <- sum(sapply(Tm[Primer %in% c('F2','B2')]-defaultTemps$F2, max, 0)^power)
		TmScore2 <- sum(sapply(Tm[Primer %in% c('F1c','B1c')]-defaultTemps$F1, max, 0)^power)
		TmScore3 <- sum(sapply((HpTm-HpTmTarget)/2.5, max, 0)^power)
		primerScore <- DimerScore + HpPrimerScore + HpDBScore + TmScore1 + TmScore2 + TmScore3
		return(primerScore)
	})
	return(primerScore)
}

plotEnergies <- function(primerSets)
{
	primerNames <- c('F3','B3','F2','B2','F1c','B1c','FIP','BIP','LF','LB','PNAF','PNAB','DBF','DBB')
	ret <- data.table(Primer=primerNames,
							Seq=as.character(sapply(lapply(primerNames, function(x){vals$NTs[[x]]}), paste, collapse='')),
							Sense=c('Sense','Antisense','Sense','Antisense','Sense','Antisense','NA','NA','Antisense','Sense','Sense','Antisense','Sense','Sense'))
	ret[, Stability3p:=end_stability(Seq, n=vals$stabilityN, threePrimeEnd=T, mv=MvConc, dv=DvConc, dntp=DntpConc, temp_c=RxnTemp), by='Primer']
	ret[, Stability5p:=end_stability(Seq, n=vals$stabilityN, threePrimeEnd=F, mv=MvConc, dv=DvConc, dntp=DntpConc, temp_c=RxnTemp), by='Primer']
	ret[, Start5p:='NA']
	ret[Primer == 'F3', Start5p:=vals$Start$F3]
	ret[Primer == 'B3', Start5p:=getEndSetting('B3c', vals)]
	ret[Primer == 'F2', Start5p:=vals$Start$F2]
	ret[Primer == 'F1c', Start5p:=getEndSetting('F1', vals)]
	ret[Primer == 'B2', Start5p:=getEndSetting('B2c', vals)]
	ret[Primer == 'B1c', Start5p:=vals$Start$B1c]
	ret[Primer == 'LF', Start5p:=getEndSetting('LFc', vals)]
	ret[Primer == 'LB', Start5p:=vals$Start$LB]
	ret[Primer == 'PNAF', Start5p:=vals$Start$PNAF]
	ret[Primer == 'PNAB', Start5p:=getEndSetting('PNABc', vals)]
	ret[, Len:=length(s2c(Seq)), by='Primer']
	ret[, c('Tm','HpTm','HpDeltaG','HpStruct'):=getPrimerStats(Seq), by='Primer']
	checked <- gsub('c','',names(vals$Check)[as.logical(vals$Check)])
	checked <- checked[checked %!in% c('F1','F2','B1','B2')]
	primerBaseNames <- gsub('c','',ret$Primer)
	primersToCalc <- primerBaseNames %in% c(checked, 'FIP', 'BIP')
	# browser()
	ret[primersToCalc, c('P2','DimerTm','DimerDeltaG','DimerStruct'):=getDimerStats(Primer, Seq, func=getDimer, mv=MvConc, dv=DvConc, dntp=DntpConc, temp_c=RxnTemp)]
	# ret[Primer %in% c('FIP','BIP'), Tm:='NA']
	ret[, KeyEndStability:=ifelse(Primer %in% c('F1c','B1c'), Stability5p, Stability3p)]
	
	# Data for energy plot
	ret2 <- melt.data.table(data=ret, id.vars='Primer', measure.vars=c('KeyEndStability','HpDeltaG','DimerDeltaG'))
	ret2[, value:=suppressWarnings(as.numeric(value))]
	ret2[value > 0, value:=0]
	ret2[variable=='KeyEndStability', value:=-1*value]
	ret2[, Primer:=factor(ret2$Primer, levels=primerNames)]
	
	ret <- ggplot(data=results[!is.na(value) & variable != 'KeyEndStability'], aes( x=Primer, y=value, fill=variable)) +
		geom_col() +
		geom_col(data=results[!is.na(value) & variable == 'KeyEndStability']) +
		geom_hline(yintercept=c(4,-4)) +
		labs(x='Sequence', y='Energy [kJ/mol]') +
		scale_y_continuous(limits=c(-25,10),oob = rescale_none) +
		theme(axis.text=element_text(size=rel(2.0)),
				axis.title=element_text(size=rel(2.0)),
				legend.text=element_text(size=rel(2.0)),
				legend.title=element_text(size=rel(2.0)))
	return(ret)
}

InitList <- data.table.expand.grid(Len=11:28)
InitList <- InitList[, getSeqsWithLen(rpoB_500_c, len=.BY[[1]], align='left'), by='Len']
InitList[, Tm:=getTm(Seq), by=c('Start','Len')]
InitList[, HmdDeltaG:=getHmd(Seq), by=c('Start','Len')]
InitList[, c('HpTm','HpDeltaG','HpStruct'):=getHp(Seq), by=c('Start','Len')]
InitList[is.na(HpDeltaG), HpDeltaG:=999]
InitList[, End:=Start+Len-1]
InitList[, Id:=as.character(1:.N)]

# Generate Candidate Primer Sets
# primerSets <- getPossibleWindowsOverlappingIndex(primer='B2c', TmDelta=1, TmOffset=0, Index=260, HmdLimit=-10, HpLimit=-2)
# primerSets[, c('HpTm','HpDeltaG','HpStruct'):=getHp(Seq)[c('temp','dg','structure')], by='primerId']

primerSets <- getPossibleWindowsFromEnd(primer='B2c', primerList=InitList, TmDelta=1, TmOffset=0, endMin=260, endMax=264, HmdLimit=-2, HpLimit=-1)
if(nrow(primerSets)==0){stop('No candidate primers produced')}
primerSets <- tryPrimersWithPrimerSets(primerSets, getPossibleWindowsFromBounds(primer='B1c', primerList=InitList, TmDelta=1, TmOffset=0, bMin=175, bMax=210, HmdLimit=0, HpLimit=-1))
if(nrow(primerSets)==0){stop('No candidate primers produced')}
primerSets <- tryPrimersAdjacentToPrimerSets(primerSets, primerList=InitList, right=F, relToPrimer='B1c', newPrimer='F1', TmDelta=1, minSpace=0, maxSpace=25, TmOffset=0, HmdLimit=-6, HpLimit=-1)
if(nrow(primerSets)==0){stop('No candidate primers produced')}
primerSets <- tryPrimersWithPrimerSets(primerSets, getPossibleWindowsFromEnd(primer='F2', primerList=InitList, TmDelta=1, TmOffset=0, endMin=93, endMax=123, HmdLimit=0))
if(nrow(primerSets)==0){stop('No candidate primers produced')}
# lapply(c('F2','F1','B2c','B1c'), FUN=function(x){
# 	getUniqueCombos(primerSets, by=c('Primer','Id','HmdDeltaG'))
# })
# primerSets <- tryPrimersWithPrimerSets(getPossibleWindowsFromStart(primer='F1', TmDelta=1, TmOffset=1, startMin=105, startMax=115), primerSets)
# primerSets <- tryPrimersAdjacentToPrimerSets(primerSets, right=F, relToPrimer='F1', newPrimer='F2', TmDelta=1, minSpace=20, maxSpace=30, TmOffset=0, HmdLimit=-7)


# Make FIP, BIP, DBF, and DBB
primerSets2 <- primerSets[, assemblePrimers(.SD, polyNT = 'a'), by='primerId']
hist(primerSets2[Primer %in% c('BIP')]$HmdDeltaG, main='BIP')
primerSets3 <- primerSets2[, if(min(HmdDeltaG[Primer %in% c('BIP')]) >= -4) .SD, by='primerId']
hist(primerSets3[Primer %in% c('FIP')]$HmdDeltaG, main='FIP')
primerSets3 <- primerSets3[, if(min(HmdDeltaG[Primer %in% c('FIP')]) >= 0) .SD, by='primerId']

if(nrow(primerSets3)==0){stop('No candidate primers produced')}

# Make Copy and Calculate Hairpin Energy
primerSets3[Primer %!in% c('DBB','DBF'), c('HpTm','HpDeltaG','HpStruct'):=getHp(Seq), by=c('primerId','Primer')]
primerSets3[, c('HpTm','HpDeltaG','HpStruct'):=getDBHps(.SD, width=25), by='primerId']
# write(primerSets3, 'PrimerSets3.csv')
primerSets3[, PrimerFactor:=factor(Primer, levels=finalPrimerNames)]
primerSets3[, PrimerNum:=as.numeric(PrimerFactor)]
# data.table.plot.all(primerSets3, xcol='TypeNum', ycol='HpDeltaG', by='primerId', legend.plot = F)
levels(primerSets3$PrimerFactor)
primerSets4 <- copy(primerSets3) #primerSets3[, .SD[all(HpDeltaG >= -3)], by='primerId']
if(nrow(primerSets4)==0){stop('No candidate primers produced')}

# Calculate Dimer and KeyEndStability Energy
primerSets4[Primer %in% c('FIP','BIP'), c('P2','DimerTm','DimerDeltaG'):=getDimerComboStats(primers=Primer, seqs=Seq), by='primerId']
primerSets4[Primer %!in% c('DBF', 'DBB'), KeyEndStability:=ifelse(Primer %in% c('F1c','B1c'), getStability(Seq, threePrimeEnd = F), getStability(Seq, threePrimeEnd = T)), by=c('primerId','Primer')]
data.table.plot.all(primerSets4, xcol='PrimerNum', ycol='DimerDeltaG', by='primerId', legend.plot = F)

# Score, Rank, and Filter
primerSets4[, Score1:=scoreAlgorithm(.SD, DimerTarget=max(primerSets4$DimerDeltaG, na.rm=T), HpTarget=max(primerSets4$HpDeltaG[primerSets4$Primer %!in% c('DBF','DBB')], na.rm=T), HpTmTarget=55, HpDBTarget=max(primerSets4$HpDeltaG[primerSets4$Primer %in% c('DBF','DBB')], na.rm=T)), by='primerId']
setorder(primerSets4, Score1, primerId)
primerSets4[, FIPBIPRank:=.GRP, by='primerId']
primerSets5 <- primerSets4[FIPBIPRank <= 10]

# Add F3 and B3 to best sets
primerSets5 <- tryPrimersAdjacentToPrimerSets(primerSets5, primerList=InitList, right=F, relToPrimer='F2', newPrimer='F3', TmDelta=1, minSpace=0, maxSpace=25, TmOffset=0, HmdLimit=-4, HpLimit=-1)
primerSets5 <- tryPrimersAdjacentToPrimerSets(primerSets5, primerList=InitList, right=T, relToPrimer='B2', newPrimer='B3c', TmDelta=1, minSpace=0, maxSpace=50, TmOffset=0, HmdLimit=-4, HpLimit=-1)
if(nrow(primerSets)==0){stop('No candidate primers produced')}

# Calculate added Dimer and EndStabilityEnergies
primerSets5[Primer %in% c('FIP','BIP'), c('P2','DimerTm','DimerDeltaG'):=getDimerComboStats(primers=Primer, seqs=Seq), by='primerId']
primerSets5[Primer %!in% c('DBF', 'DBB'), KeyEndStability:=ifelse(Primer %in% c('F1c','B1c'), getStability(Seq, threePrimeEnd = F), getStability(Seq, threePrimeEnd = T)), by=c('primerId','Primer')]
data.table.plot.all(primerSets5, xcol='PrimerNum', ycol='DimerDeltaG', by='primerId', legend.plot = F)

# Rescore, Rerank, and Refilter
primerSets5[, Score1:=scoreAlgorithm(.SD, DimerTarget=max(primerSets5$DimerDeltaG, na.rm=T), HpTarget=max(primerSets5$HpDeltaG[primerSets5$Primer %!in% c('DBF','DBB')], na.rm=T), HpTmTarget=55, HpDBTarget=max(primerSets5$HpDeltaG[primerSets4$Primer %in% c('DBF','DBB')], na.rm=T)), by='primerId']
setorder(primerSets5, Score1, primerId)
primerSets5[, F3B3Rank:=.GRP, by='primerId']
primerSets6 <- primerSets5[F3B3Rank <= 10]

# Add Loop Primers

# Add PNA Primers

# Calculate Hairpin Energy

# Calculate Dimer Energy

# Calculate Different Score Metrics

# Sort and Rank

# Output Files/Plots/Tables for Top Ranked Primer Sets

