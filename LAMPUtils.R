#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

# install.packages("devtools")
# library(devtools)
# devtools::install_github("SalusDiscovery/primer3")
# devtools::install_github("JohnCoene/marker")

# list.of.packages <- c('shiny','zoo','XNAString','Biostrings','seqinr','data.table','primer3')
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("XNAString")

library(data.table)
library(zoo)
library(grDevices)
library(XNAString)
library(Biostrings)
library(seqinr)
library(primer3)
library(pwalign)
library(ggplot2)
library(scales)

##### Constants #####
sensePrimerNames <- c('F3','F2','F1','LFc','B1c','B2c','B3c','LB','PNAF','PNABc')
dependentPrimerNames <- c('F1c', 'B2', 'B3', 'LF', 'LB', 'PNAB', 'FIP', 'BIP', 'DBF', 'DBB')
finalPrimerNames <- c('F3','B3','F2','B2','F1c','B1c','FIP','BIP','LF','LB','PNAF','PNAB','DBF','DBB')
allPrimerNames <- c('F3','F3c','F2','F2c','F1','F1c','B1','B1c','B2','B2c','B3','B3c','LF','LFc','LB','LBc','FIP','BIP','DBF','DBB','PNAF','PNAFc','PNAB','PNABc')

# Per PrimerExplorerV5 defaults
# 5' Stability: < -3
# 3' Stability: < -4
# Dimer Stability: > -2.5
# F1c B1c Tm: 66 ± 2
# Other Tm: 61 ± 2
# F1c-F2 Space: 40-60
# F3-F2 Space: 0-20
# F1c-B1c Space: 0-100
defaultTemps <- list(F3=61, F3c=61, F2=61, F2c=61, B3=61, B3c=61, B2=61, B2c=61, LF=66, LFc=66, LB=66, LBc=66, PNAF=61, PNAFc=61, PNAB=61, PNABc=61, F1=66, F1c=66, B1=66, B1c=66)
defaultWeights <- list(F1=3, F1c=3, B1=3, B1c=3, FIP=5, BIP=5, F3=1, F3c=1, F2=2, F2c=2, B3=1, B3c=1, B2=2, B2c=2, LF=1.5, LFc=1.5, LB=1.5, LBc=1.5, PNAF=1.5, PNAFc=1.5, PNAB=1.5, PNABc=1.5, DBB=2, DBF=2)

# Settings to match PrimerExplorerV5 (Eiken) melt temperature and end stability calculations and mFold deltaG at Na=50, temp=65
paramsDNA <- list(
	Tm = list(
		salt_method = 'Schildkraut',
		tm_method = 'SantaLucia',
		mv = 50,       # [mM]
		dv = 4.45,     # [mM]
		dntp = 0, # [mM]
		dna = 100 # [mM]
	),
	Dimer = list(
		temp_c = 65,   # [°C]
		mv = 50,       # [mM]
		dv = 7,        # [mM]
		dntp = 0, # [mM]
		dna = 100 # [mM]
	),
	Stability = list(
		n=6,           #
		temp_c = 32,   # [°C]
		mv = 50,       # [mM]
		dv = 7,        # [mM]
		dntp = 0, # [mM]
		dna = 100 # [mM]
	),
	Hp = list(
		temp_c = 66.5,  # [°C]
		mv = 50,   # [mM]
		dv = 7, # [mM]
		dntp = 0,  # [mM]
		dna = 100  # [mM]
	)
)

##### TB Constants #####
rpoB_500 <- 'atcgaccacttcggcaaccgccgcctgcgtacggtcggcgagctgatccaaaaccagatccgggtcggcatgtcgcggatggagcgggtggtccgggagcggatgaccacccaggacgtggaggcgatcacaccgcagacgttgatcaacatccggccggtggtcgccgcgatcaaggagttcttcggcaccagccagctgagccaattcatgGACcagaacaacccgctgtcggggttgaccCACaagcgccgactgTCGgcgctggggcccggcggtctgtcacgtgagcgtgccgggctggaggtccgcgacgtgcacccgtcgcactacggccggatgtgcccgatcgaaacccctgaggggcccaacatcggtctgatcggctcgctgtcggtgtacgcgcgggtcaacccgttcgggttcatcgaaacgccgtaccgcaaggtggtcgacggcgtggttagcgacgagatcgtgtacctgaccgccgacgagga'
rpoB_500_c <- s2c(rpoB_500)

##### Helper Functions #####

getDefault <- function(x, default, test=is.null)
{
	ret <- copy(x)
	result <- test(x)
	if(length(result) > 1)
	{
		ret[result] <- default
	}
	else
	{
		if(result)
		{
			ret <- default
		}
	}
	# ret[test(x)] <- default
	return(ret)
}

appendToList <- function(items, ...)
{
	# Default values in list1
	# Overriding values in list2
	args <- list(...)
	if(any(names(args) %in% names(items)))
	{
		stop('At least one of the items already exists in the list provided. Use updateList instead.')
	}
	for(name in names(args))
	{
		items[[name]] <- args[[name]]
	}
	return(items)
}

updateList <- function(items, ..., allowNew=F)
{
	# Default values in list1
	# Overriding values in list2
	args <- list(...)
	if(!allowNew && !all(names(args) %in% names(items)))
	{
		stop('Item to update does not exist in list.')
	}
	for(name in names(args))
	{
		items[[name]] <- args[[name]]
	}
	return(items)
}

loadRData <- function(fileName)
{
	#loads an RData file, and returns it
	load(fileName)
	get(ls()[ls() != "fileName"])
}

is.even <- function(x)
{
	return(x %% 2 == 0)
}

is.odd <- function(x)
{
	x %% 2 != 0
}

revC <- function(x, keepCase=F)
{
	if(any(is.na(x)) || any(grepl( '[^AGTCagtcNn]', x)))
	{
		print(x)
		stop("Found a non coding character in revC.")
	}
	if(length(x) == 1)
	{
		ret <- toString(reverseComplement(DNAString(toUpper(x))))
		if(keepCase)
		{
			ret <- matchCase(template=paste(rev(s2c(x)), collapse=''), target=ret)
		}
	}
	else
	{
		ret <- s2c(toString(reverseComplement(DNAString(toUpper(paste(x, collapse=''))))))
		if(keepCase)
		{
			ret <- matchCase(template=rev(x), target=ret)
		}
	}
	return(ret)
}

toLower <- function(x)
{
	# print(x)
	return(tolower(enc2utf8(x)))
}

toUpper <- function(x)
{
	# print(x)
	return(toupper(enc2utf8(x)))
}

'%!in%' <- function(x,y)!('%in%'(x,y))

calcLen <- function(oligo)
{
	return(ifelse(length(oligo)==1, nchar(oligo), length(oligo)))
}

calcEnd <- function(start, len)
{
	return(start+len-1)
}

getCombos <- function(x, colnames=c('Var1','Var2'))
{
	Var1 <- c()
	Var2 <- c()
	for(i in seq_along(x))
	{
		for(j in i:length(x))
		{
			Var1 <- c(Var1, x[i])
			Var2 <- c(Var2, x[j])
		}
	}
	ret <- data.table(Var1=Var1, Var2=Var2)
	setnames(ret, old=c('Var1','Var2'), new=colnames)
	return(ret)
}

setColor <- function(my.colors, alpha)
{
	x <- col2rgb(my.colors)/255
	x <- data.table(t(x), alpha=alpha)
	ret <- rgb(x$red,x$green,x$blue,x$alpha)
	return(ret)
}

sig.digits <- function(x, nSig=2, trim.spaces=T, trim.zeros=F)
{
	ret <- getPrettyNum(x, sigFigs = nSig, dropTrailingZeros = trim.zeros)
	if(trim.spaces)
	{
		ret <- trimws(ret)
	}
	return(ret)
}

getPrettyNum <- function(x, sigFigs=3, dropTrailingZeros=F)
{
	ret <- formatC(signif(x,digits=sigFigs), digits=sigFigs,format="fg", flag="#", drop0trailing = dropTrailingZeros)
	if(any(endsWith(ret, '.')))
	{
		for(i in 1:length(ret))
		{
			if(endsWith(ret[i], '.'))
			{
				ret[i] <- substr(ret[i], 1, nchar(ret[i])-1)
			}
		}
	}
	return(ret)
}


#' Read table from the clipboard
#'
#' This is a cool way to import data using the clipboard. The clipboard table
#' is typically copied as a tab delimited text 'file' connection
#'
#' @param os - c('mac','win'), string value indicating the platform to use
#' @param header - TRUE or FALSE, whether the header is included in the copied table
#' @param sep - text, defining the separator character used between values (needs to be in double quotes)
#' @param use.data.table - TRUE or FALSE, whether to return a data.table (default)
#' @param ... - additional variables supplied are passed onto the underlying read.table function (e.g., stringsAsFactors, comment.char, col.names)
#'
#' @export
read.clipboard <- function(os=c('mac','win'), header=T, sep="\t", use.data.table=T, ...)
{
	if(os[1] %in% c('Darwin','mac'))
	{
		ret <- read.table(pipe('pbpaste'), header=header, sep=sep, ...) # Mac
	}
	else
	{
		ret <- read.table('clipboard', header=header, sep=sep, ...) # Windows
	}
	if(use.data.table)
	{
		return(data.table(ret))
	}
	else
	{
		return(ret)
	}
}

copyPrimerSet <- function(primerSets, rank, scorecol='Score')
{
	nameSubs <- list(F1c='F1', F2c='F2', F3c='F3', B1='B1c', B2='B2c', B3='B3c', PNAFc='PNAF', PNAB='PNABc', LF='LFc', LBc='LB')
	temp <- copy(primerSets)
	setorder(temp, primerId)
	setorderv(temp, scorecol)
	temp[, myRank:=.GRP, by=c(scorecol, 'primerId')]
	temp[Primer %in% names(nameSubs), c('Primer', 'Seq'):=list(nameSubs[[.BY[[2]]]], revC(Seq, keepCase=T)), by=c('primerId','Primer')]
	ret <- temp[myRank==rank & Primer %in% nameSubs, c('Primer','Start','End')]
	clip <- pipe("pbcopy", "w")                       
	write.table(ret, file=clip, col.names = F, row.names = F, quote = F, sep='\t')                               
	close(clip)
	print(ret)
}

getSettingsFile <- function(username='jay@salusdiscovery.com', project='mTB', batch=1, set=1)
{
	path <- paste('~/Library/CloudStorage/GoogleDrive-', 'jay@salusdiscovery.com', '/Shared drives/Salus Discovery/Projects/3STEP/3STEP Salus/Primer ReferenceMaterials/', project, ' Primers/Batch ', batch, '/Completed Sets/Primer Set ', set, '/Primer Set ', set, '/Settings.RData', sep='')
	ret <- loadRData(path)$results
}

getEditablePrimerSet <- function(settings, primerList)
{
	# Example: b2s1 <- getEditablePrimerSet(getSettingsFile(batch='2', set='1'), InitList)
	settings[, primerId:=1]
	settings[, Id:=as.character(-1*(1:.N))]
	settings[Primer %!in% c('FIP','BIP','DBF','DBB'), Id:=ifelse(Seq %in% InitList$Seq, InitList$Id[match(Seq, InitList$Seq)], InitList$Id[match(revC(Seq, keepCase = T), InitList$Seq)]), by='Primer']
	settings[Primer %in% c('FIP','DBF'), Id:=paste(settings$Id[settings$Primer=='F1c'], settings$Id[settings$Primer=='F2'], sep='.')]
	settings[Primer %in% c('BIP','DBB'), Id:=paste(settings$Id[settings$Primer=='B1c'], settings$Id[settings$Primer=='B2'], sep='.')]
	settings <- settings[primerList, c('HmdTm','HmdDeltaG','HmdStruct'):=list(HmdTm, HmdDeltaG, HmdStruct), on='Id']
	settings[, primerId:=1]
	return(settings[, c('Primer', names(primerList), 'primerId'), with=F])
}

data.table.expand.grid <- function(..., BaseTable=NULL, KEEP.OUT.ATTRS=T, stringsAsFactors = T)
{
	if(!is.null(BaseTable))
	{
		return(BaseTable[, data.table.expand.grid(..., BaseTable=NULL, KEEP.OUT.ATTRS=KEEP.OUT.ATTRS, stringsAsFactors=stringsAsFactors), by=names(BaseTable)])
	}
	else
	{
		return(data.table(expand.grid(..., KEEP.OUT.ATTRS=KEEP.OUT.ATTRS, stringsAsFactors=stringsAsFactors)))
	}
}

#' rbind.results
#'
#' @param dt We have to have a separate parameter for the table to be worked on
#' @param dt.expression unquoted data.table expression that uses 'dt' as the table to operate on
#' @param ... parameters to expand.grid on and run the data.table expression
#'
#' @return the aggregated results of the dt expression for each possible
#' value combination of the arguments. rbindlist is used to aggregate the
#' results. Argument values are added as new columns to the results table.
#'
#' @examples
#' duh <- data.table(a=1:4, b=c(1,2,1,2), c=c(1,1,2,2))
#' zeta <- 25
#' rbind.results(duh, duh[, c('d','e','f'):=.(a+alpha, a+beta, a+zeta), by=c('a')], alpha=c(1,2,3), beta=c(1))
#'
rbind.results <- function(dt.expression, grpColName='paramSet', ...)
{
	args <- list(...)
	args.table <- data.table.expand.grid(args)
	rets <- list()
	# copy the table to the current environment to perform the calculation
	table.name <- strsplit(deparse(match.call()$dt.expression, width.cutoff = 500), '[', fixed=T)[[1]][1]
	assign(table.name, get(table.name, envir=parent.frame()), envir=environment())
	for(i in 1:nrow(args.table))
	{
		args.list <- as.list(args.table[i])
		make.vars(args.list, env=environment())
		rets[[i]] <- copy(eval(match.call()$dt.expression, envir=environment()))
		rets[[i]][, names(args.table):=args.list]
		rets[[i]][[grpColName]] <- i
	}
	return(rbindlist(rets))
}

##### DNA Calc/Manipulation Functions #####

getSeqFromStartAndLen <- function(start, len, seq)
{
	return(seq[start:(start+len-1)])
}

getStartFromSeq <- function(primer, target)
{
	temp <- pwalign::pairwiseAlignment(toUpper(primer), toUpper(target), type='local', gapOpening=1000000, gapExtension=1000000)
	return(start(subject(temp)))
}

getTm <- function(x, params=paramsDNA$Tm)
{
	return(calculate_tm(toUpper(paste(x, collapse='')), salt_conc=params$mv, divalent_conc=params$dv, dntp_conc=params$dntp, dna_conc=params$dna, tm_method=params$tm_method, salt_correction=params$salt_method))
}

getHp <- function(x, useXNAString=F, params=paramsDNA$Hp)
{
	if(useXNAString)
	{
		library(XNAString)
		ret1 <- predictMfeStructure(XNAString(base = toUpper(paste(x, collapse=''))))
		ret2 <- list(structure_found = T,
						 temp = 0,
						 dg = ret1$mfe,
						 structure = ret1$structure)
	}
	else
	{
		ret1 <- toUpper(paste(x, collapse=''))
		ret2 <- calculate_hairpin(ret1, mv=params$mv, dv=params$dv, dntp=params$dntp, temp_c=params$temp_c, dna=params$dna)
		ret2$dg <- ret2$dg/1000
	}
	if(is.na(ret2$dg)){ret2$structure<-''}
	ret3 <- list(HpTm=ret2$temp,
					 HpDeltaG=ret2$dg,
					 HpStruct=ifelse(!is.na(ret2$dg), getHpCamelCase(paste(x, collapse=''), ret2$structure), ''))
	return(ret3)
}

getHmd <- function(x, params=paramsDNA$Dimer)
{
	seq <- toUpper(paste(x, collapse=''))
	ret <- getDimer(seq, seq, params=params)
	names(ret) <- c('HmdTm','HmdDeltaG','HmdStruct')
	return(ret)
}

#' getDimer
#' 
#' Function that is used with getDimerComboStats in a data.table. Returns the dimer
#' melt temperature, deltaG energy in kJ/mol, and struct string where capitals
#' indicate bonded NTs and lower case as unbonded NTs with ampersand between the
#' two sequences
#' 
#' @param x single string or character vector of individual nucleotides (lower and/or uppcase)
#' @param y single string or character vector of individual nucleotides (lower and/or uppcase)
#'
#' @return list(DimerTm, DimerDeltaG, DimerStruct)
#' @export
#'
#' @examples
getDimer <- function(x, y, params=paramsDNA$Dimer)
{
	hetd <- calculate_dimer(toUpper(paste(x, collapse='')), toUpper(paste(y, collapse='')), mv=params$mv, dv=params$dv, dntp=params$dntp, temp_c=params$temp_c, dna=params$dna)
	hetd$dg <- hetd$dg/1000
	ret <- list(DimerTm=hetd$temp,
					DimerDeltaG=hetd$dg,
					DimerStruct=ifelse(!is.na(hetd$dg), hetd$structure, ''))
	return(ret)
}

getStability <- function(x, threePrimeEnd=TRUE, params=paramsDNA$Stability)
{
	if(length(x) == 1)
	{
		x <- s2c(x)
	}
	if(threePrimeEnd)
	{
		x <- c(rep('N',20-params$n), x[(length(x)-params$n):length(x)])
	}
	else
	{
		x <- c(x[1:params$n], rep('N',20-params$n))
	}
	return(getDimer(x, revC(x), params=params)$DimerDeltaG)
}

getDimerComboStats <- function(primers, seqs)
{
	duh2 <- data.table(getCombos(primers, colnames=c('P1','P2')), getCombos(seqs, c('Seq1','Seq2')))
	duh2 <- duh2[duh2[, getDimer(Seq1, Seq2), by=c('P1','P2')], on=c('P1','P2')]
	valNames <- names(duh2)
	valNames <- valNames[valNames %!in% c('P1','P2','Seq1','Seq2')]
	blah <- rbindlist(list(duh2, data.table(Seq1=duh2$Seq2, Seq2=duh2$Seq1, P1=duh2$P2, P2=duh2$P1, duh2[, mget(valNames)])), use.names=T)
	blah <- unique(blah)
	blah[, max.i:=which.min(DimerDeltaG)[1], by='P1']
	blah <- blah[, .SD[max.i[1]], by='P1'][,c('P1','P2',valNames), with=F]
	return(blah[match(primers, P1), mget(c('P2', valNames))])
}

matchCase <- function(template, target)
{
	origLength <- length(target)
	# e.g., matchCase('HHHhhh','agtCCt') -> "AGTcct"
	target <- toLower(target)
	if(origLength == 1)
	{
		target <- s2c(target)
	}
	if(length(template) == 1)
	{
		template <- s2c(template)
	}
	if(length(target) != length(template))
	{
		# browser()
		stop("Target and template must be same length.")
	}
	uppers <- grepl("^[[:upper:]]+$", template)
	target[uppers] <- toUpper(target)[uppers]
	if(origLength == 1)
	{
		target <- paste(target, collapse='')
	}
	return(target)
}

getHpTable <- function(seq_c, width=25, align='left', trim=F)
{
	ret <- data.table(bp=seq_along(seq_c))
	ret[, Seq:=paste(getSeqAtIndex(seq_c=seq_c, pos=.BY[[1]], width=width, align=align), collapse=''), by='bp']
	ret[, c('HpTm','HpDeltaG','HpStruct'):=getHp(Seq), by='bp']
	if(trim)
	{
		ret <- ret[bp > width%/%2 & bp <= length(seq_c)-width%/%2]
	}
	return(ret[])
}

getHpCamelCase <- function(seq, brackets)
{
	ret <- toLower(s2c(seq))
	bonds <- s2c(brackets) %in% c('(',')')
	ret[bonds] <- toUpper(s2c(seq))[bonds]
	return(paste(ret, collapse=''))
}

getDBHps <- function(primerSet, width=25)
{
	if(is.even(width))
	{
		warning('Width must be odd. Adding 1 to width provided.')
		width <- width + 1
	}
	F1cr <- primerSet[Primer=='F1c']
	F2r <- primerSet[Primer=='F2']
	B1cr <- primerSet[Primer=='B1c']
	B2r <- primerSet[Primer=='B2']
	FIPr <- primerSet[Primer=='FIP']
	BIPr <- primerSet[Primer=='BIP']
	DBFr <- primerSet[Primer=='DBF']
	DBBr <- primerSet[Primer=='DBB']
	DBFSeqLin <- c(s2c(F1cr$Seq)[(F1cr$Len+1-width%/%2):F1cr$Len], s2c(DBFr$Seq), revC(s2c(F1cr$Seq), keepCase=T)[1:(width%/%2)]) # Linear version DB padded by F1/F1c
	DBFSeqCirc <- c(s2c(DBFr$Seq)[(DBFr$Len+1-width%/%2):DBFr$Len], s2c(DBFr$Seq), s2c(DBFr$Seq)[1:(width%/%2)]) # Circular version DB padded by DB loop
	DBFSeqLinHp <- getHpTable(DBFSeqLin, width=width, align='center', trim=T)
	DBFSeqLinHp[, bp:=bp-min(bp)+1]
	DBFSeqCircHp <- getHpTable(DBFSeqCirc, width=width, align='center', trim=T)
	DBFSeqCircHp[, bp:=bp-min(bp)+1]
	DBFSeqLinHp[, HpType:='Lin']
	DBFSeqCircHp[, HpType:='Circ']
	DBFTable <- rbindlist(list(DBFSeqCircHp, DBFSeqLinHp))
	DBFI <- which.min(DBFTable$HpDeltaG)
	DBFTable[HpStruct != '', HpStruct:=paste(HpType, ':', bp-width%/%2, '-', bp+width%/%2, ':', HpStruct, '\n', getHpCamelCase(Seq, HpStruct), sep=''), by=c('bp','HpType')]
	ret <- copy(primerSet)
	if(length(DBFI) >= 1)
	{
		ret[Primer=='DBF', c('HpTm','HpDeltaG','HpStruct'):=DBFTable[DBFI, c('HpTm','HpDeltaG','HpStruct')]]
	}
	else
	{
		ret[Primer=='DBF', c('HpTm','HpDeltaG','HpStruct'):=data.table(HpTm=as.numeric(0.0), HpDeltaG=as.numeric(NA), HpStruct='')]
	}
	
	DBBSeqLin <- c(s2c(B1cr$Seq)[(B1cr$Len+1-width%/%2):B1cr$Len], s2c(DBBr$Seq), revC(s2c(B1cr$Seq), keepCase=T)[1:(width%/%2)]) # Linear version DB padded by F1/F1c
	DBBSeqCirc <- c(s2c(DBBr$Seq)[(DBBr$Len+1-width%/%2):DBBr$Len], s2c(DBBr$Seq), s2c(DBBr$Seq)[1:(width%/%2)]) # Circular version DB padded by DB loop
	DBBSeqLinHp <- getHpTable(DBBSeqLin, width=width, align='center', trim=T)
	DBBSeqLinHp[, bp:=bp-min(bp)+1]
	DBBSeqCircHp <- getHpTable(DBBSeqCirc, width=width, align='center', trim=T)
	DBBSeqCircHp[, bp:=bp-min(bp)+1]
	DBBSeqLinHp[, HpType:='Lin']
	DBBSeqCircHp[, HpType:='Circ']
	DBBTable <- rbindlist(list(DBBSeqCircHp, DBBSeqLinHp))
	DBBTable[HpStruct != '', HpStruct:=paste(HpType, ':', bp-width%/%2, '-', bp+width%/%2, ':', HpStruct, sep='')]
	DBBI <- which.min(DBBTable$HpDeltaG)
	if(length(DBBI) >= 1)
	{
		ret[Primer=='DBB', c('HpTm','HpDeltaG','HpStruct'):=DBBTable[DBBI, c('HpTm','HpDeltaG','HpStruct')]]
	}
	else
	{
		ret[Primer=='DBB', c('HpTm','HpDeltaG','HpStruct'):=data.table(HpTm=as.numeric(0.0), HpDeltaG=as.numeric(NA), HpStruct='')]
	}
	return(ret[, c('HpTm','HpDeltaG','HpStruct')])
}

getSeqsWithLen <- function(seq_c, len, align=c('center','left','right'))
{
	ret <- rollapply(seq_along(rpoB_500_c),
						  FUN=function(x){paste(rpoB_500_c[x], collapse='')},
						  width=len,
						  align=align[1],
						  partial=FALSE)
	return(list(Start=1:length(ret), Seq=ret))
}

getPrimerStats <- function(x)
{
	hp <- getHp(x)
	Tm <- getTm(x)
	ret <- list(Tm=sig.digits(Tm, 3),
					HpTm=sig.digits(hp$temp, 3),
					HpDeltaG=sig.digits(hp$dg, 3),
					HpStruct=ifelse(hp$structure_found, hp$structure, ''))
	return(ret)
}

allLegalStartLenSettings <- function(settings)
{
	getEndSetting('F3', settings) <= length(settings$seq) &&
		getEndSetting('F2', settings) <= length(settings$seq) &&
		getEndSetting('F1', settings) <= length(settings$seq) &&
		getEndSetting('LFc', settings) <= length(settings$seq) &&
		getEndSetting('B3c', settings) <= length(settings$seq) &&
		getEndSetting('B2c', settings) <= length(settings$seq) &&
		getEndSetting('B1c', settings) <= length(settings$seq) &&
		getEndSetting('LB', settings) <= length(settings$seq) &&
		getEndSetting('PNAF', settings) <= length(settings$seq) &&
		getEndSetting('PNABc', settings) <= length(settings$seq)
}

getSeqAtIndex <- function(seq_c, pos, width, align=c('center','left','right'))
{
	if(is.even(width))
	{
		warning("Window width is even. Adding one to make it odd.")
		width <- width + 1
	}
	if(align[1]=='center')
	{
		window <- c(max(1,pos-(width%/%2)), min(length(seq_c), pos+(width%/%2)))
	}
	else if(align[1]=='left')
	{
		window <- c(max(1,pos), min(length(seq_c), pos+width-1))
	}
	else if(align[1]=='right')
	{
		window <- c(max(1,pos-width+1), min(length(seq_c), pos))
	}
	else
	{
		stop('Unrecognized alignment argument.')
	}
	return(seq_c[window[1]:window[2]])
}

##### Plotting/Displate Functions #####
getPrimerColor <- function(...)
{
	args <- list(...)
	if(length(args) == 1 && length(args[[1]])>1)
	{
		return(getAllPrimerColors()[match(as.character(args[[1]]), Primer)]$plotColor)
	}
	else
	{
		return(getAllPrimerColors()[match(args, Primer)]$plotColor)
	}
	
}

plotHairpin <- function(seq, fullHairpinEnergy, dbSeq, dbStart, dbEnd, starts, ends, sensePrimerNames)
{
	bp <- dbStart:dbEnd
	seq.new <- s2c(dbSeq)
	# bpRange <- c(max(1,start.new), min(length(seq), start.new + length(seq.new) - 1))
	# bp <- bp[bp > 0 & bp < calcLen(seq)]
	hairpin.new <- rollapply(seq.new, width=25, FUN=function(x){return(min(0, getDefault(getHp(x)$HpDeltaG, 0, test=is.na)))}, partial=T, align='center')
	hairpin.new[is.na(hairpin.new) | hairpin.new > 0] <- 0
	# validBp <- bp[bp > 0 & bp < length(seq)]
	hairpin <- fullHairpinEnergy # rollapply(seq, width=25, FUN=function(x){return(getHp(x)$HpDeltaG/1000)}, partial=T, align='center')
	# hairpin <- rollapply(seq[validBp], width=25, FUN=function(x){return(getHp(x)$HpDeltaG/1000)}, partial=T, align='center')
	hairpin[is.na(hairpin) | hairpin > 0] <- 0
	if(any(seq %in% c('A','G','T','C')))
	{
		codons <- which(seq %in% c('A','G','T','C'))
	}
	else
	{
		codons <- c()
	}
	ret <- ggplot(NULL) +
		geom_rect(data=data.table(Primer=factor(sensePrimerNames, levels=sensePrimerNames), xmin=starts, xmax=ends), aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill=Primer), color='black')	+
		scale_fill_manual(values = setColor(getPrimerColor(sensePrimerNames), 0.3)) +
		geom_path(data=data.table(x=1:calcLen(seq), y=hairpin), aes(x=x, y=y)) +
		geom_path(data=data.table(x=bp, y=hairpin.new), aes(x=x, y=y), linewidth=3) +
		geom_hline(yintercept=-2, col=setColor('red', 0.2)) +
		geom_vline(xintercept=codons, col=setColor('red', 0.2)) +
		scale_x_continuous(limits=c(1,calcLen(seq))) +
		scale_y_continuous(limits=c(-6,0.1)) + 
		xlab('Nucleotide Position') +
		ylab('Energy [kJ/mol]') +
		theme_classic() + 
		theme(axis.text=element_text(size=rel(2.0)),
				axis.title=element_text(size=rel(2.0)),
				legend.text=element_text(size=rel(2.0)),
				legend.title=element_text(size=rel(2.0)),
				legend.position='none',
				panel.border = element_rect(colour = "black", fill=NA, size=1))
	return(ret)
}

plotTm <- function(results)
{
	# browser()
	# setorder(results, PrimerFactor)
	results2 <- copy(results)
	results2[, PrimerFactor:=factor(Primer, levels=finalPrimerNames)]
	setorder(results2, PrimerFactor)
	ret <- ggplot(data=results2, aes(x=PrimerFactor, y=Tm, fill=PrimerFactor)) + 
		geom_col() +
		scale_fill_manual(values = getPrimerColor(as.character(results2$PrimerFactor))) + 
		geom_hline(yintercept=c(defaultTemps$F1, defaultTemps$F2)) +
		scale_y_continuous(limits=c(50,82),oob = rescale_none) +
		theme(axis.text=element_text(size=rel(2.0)),
				axis.title=element_text(size=rel(2.0)),
				legend.text=element_text(size=rel(2.0)),
				legend.title=element_text(size=rel(2.0)))
	return(ret)
}

plotGC <- function(seq, fullGC, dbSeq, dbStart, dbEnd, starts, ends, sensePrimerNames)
{
	gc.frac <- fullGC # rollapply(gc, width=20, FUN=mean, partial=T, align='center')
	start.new <- dbStart
	seq.new <- s2c(dbSeq)
	gc.new <- ifelse(seq.new %in% c("g","G","c","C"), 1, 0)
	gc.frac.new <- rollapply(gc.new, width=20, FUN=mean, partial=T, align='center')
	bp <- start.new:(start.new + length(seq.new) - 1)
	bp <- bp[bp > 0 & bp < calcLen(seq)]
	if(any(seq[bp] %in% c('A','G','T','C')))
	{
		codons <- bp[which(seq[bp] %in% c('A','G','T','C'))]
	}
	else
	{
		codons <- c()
	}
	par(mar = c(4, 5, 1, 1), mgp=c(2.7, 1, 0))
	
	ret <- ggplot(NULL) +
		geom_rect(data=data.table(Primer=factor(sensePrimerNames, levels=sensePrimerNames), xmin=starts, xmax=ends), aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill=Primer), color='black')	+
		scale_fill_manual(values = setColor(getPrimerColor(sensePrimerNames), 0.3)) +
		geom_path(data=data.table(x=seq_along(gc.frac), y=gc.frac), aes(x=x, y=y)) +
		geom_path(data=data.table(x=start.new:(start.new+(length(gc.frac.new)-1)), y=gc.frac.new), aes(x=x, y=y), linewidth=3) +
		geom_hline(yintercept=0.6, col=setColor('red', 0.2)) +
		geom_vline(xintercept=codons, col=setColor('red', 0.2)) +
		scale_x_continuous(limits=c(1,calcLen(seq))) +
		scale_y_continuous(limits=c(0,1)) + 
		xlab('GC Fraction') +
		ylab('Energy [kJ/mol]') +
		theme_classic() + 
		theme(axis.text=element_text(size=rel(2.0)),
				axis.title=element_text(size=rel(2.0)),
				legend.text=element_text(size=rel(2.0)),
				legend.title=element_text(size=rel(2.0)),
				legend.position='none',
				panel.border = element_rect(colour = "black", fill=NA, size=1))
	return(ret)
}

plotPrimerEnergies <- function(primerSet)
{
	# Data for energy plot
	ret2 <- melt.data.table(data=primerSet, id.vars='Primer', measure.vars=c('KeyEndStability','HpDeltaG','DimerDeltaG'))
	ret2[, value:=suppressWarnings(as.numeric(value))]
	ret2[value > 0, value:=0]
	ret2[variable=='KeyEndStability', value:=-1*value]
	ret2[, Primer:=factor(Primer, levels=finalPrimerNames)]
	plotEnergies(ret2)
}

plotEnergies <- function(results)
{
	ret <- ggplot(data=results[!is.na(value) & variable == 'KeyEndStability'], aes( x=Primer, y=value, fill=variable)) +
		geom_col() +
		geom_col(data=results[!is.na(value) & variable != 'KeyEndStability']) +
		scale_fill_manual(values = hue_pal()(3)[c(3,1,2)]) + 
		geom_hline(yintercept=c(4,-2)) +
		labs(x='Sequence', y='Energy [kJ/mol]') +
		scale_y_continuous(limits=c(-15,10),oob = rescale_none) +
		theme(axis.text=element_text(size=rel(2.0)),
				axis.title=element_text(size=rel(2.0)),
				legend.text=element_text(size=rel(2.0)),
				legend.title=element_text(size=rel(2.0)))
	return(ret)
}

getAllPrimerColors <- function()
{
	return(data.table(Primer=allPrimerNames,
							plotColor=c('#99ff99', # F3
											'#99ff99', # F3
											'#66ffff', # F2
											'#66ffff', # F2
											'#ffff66', # F1
											'#ffff66', # F1
											'#ffff66', # B1c
											'#ffff66', # B1c
											'#66ffff', # B2c
											'#66ffff', # B2c
											'#99ff99', # B3c
											'#99ff99', # B3c
											'#cc99ff', # LFc
											'#cc99ff', # LFc
											'#cc99ff', # LB
											'#cc99ff', # LB
											'#003399', # FIP
											'#003399', # BIP
											'#ff0066', # DBF
											'#ff0066', # DBB
											'#ff3300', # PNAF
											'#ff3300', # PNAF
											'#ff3300', # PNABc
											'#ff3300') # PNABc
	))
}

sensePrimerColors <- function()
{
	getPrimerColor(sensePrimerNames)
}

getPossibleWindowsFromBounds <- function(primer, primerList, TmDelta=1, bMin, bMax, TmOffset=0, TmTarget=defaultTemps[[primer]], HmdTmLimit=35, HmdLimit=-7, HpTmLimit=35, HpLimit=-1, idName='primerId')
{
	ret <- data.table(Primer=primer, primerList[abs(Tm-TmTarget) <= TmDelta & Start>=bMin & End<=bMax & HmdDeltaG>=HmdLimit & HpDeltaG >= HpLimit])
	if(nrow(ret)==0)
	{
		stop("No possible windows found with provided filters.")
	}
	ret[, c(idName):=1:.N]
	if(primer %in% c('F1c','F2c','F3c','B1','B2','B3','LF','LBc','PNAFc','PNAB'))
	{
		ret[, Seq:=revC(Seq, keepCase=T), by=c(idName)]
	}
	return(ret[!is.na(Start)])
}

getPossibleWindowsOverlappingIndex <- function(primer, primerList, TmDelta=1, Index, TmOffset=0, TmTarget=defaultTemps[[primer]], HmdTmLimit=35, HmdLimit=-7, HpTmLimit=35, HpLimit=-1, idName='primerId')
{
	ret <- data.table(Primer=primer, primerList[abs(Tm-TmTarget) <= TmDelta & Start<=Index & End>=Index & HmdDeltaG>=HmdLimit & HpDeltaG >= HpLimit])
	if(nrow(ret)==0)
	{
		stop("No possible windows found with provided filters.")
	}
	ret[, c(idName):=1:.N]
	if(primer %in% c('F1c','F2c','F3c','B1','B2','B3','LF','LBc','PNAFc','PNAB'))
	{
		ret[, Seq:=revC(Seq, keepCase=T), by=c(idName)]
	}
	return(ret[!is.na(Start)])
}

getPossibleWindowsFromStart <- function(primer, primerList, TmDelta=1, startMin, startMax, TmOffset=0, TmTarget=defaultTemps[[primer]], HmdTmLimit=35, HmdLimit=-7, HpTmLimit=35, HpLimit=-1, idName='primerId')
{
	ret <- primerList[abs(Tm-TmTarget) <= TmDelta & Start>=startMin & Start<=startMax & HmdDeltaG>=HmdLimit & HpDeltaG >= HpLimit]
	if(nrow(ret)==0)
	{
		stop("No possible windows found with provided filters.")
	}
	ret <- data.table(Primer=primer, ret)
	ret[, c(idName):=1:.N]
	if(primer %in% c('F1c','F2c','F3c','B1','B2','B3','LF','LBc','PNAFc','PNAB'))
	{
		ret[, Seq:=revC(Seq, keepCase=T), by=c(idName)]
	}
	return(ret[!is.na(Start)])
}

getPossibleWindowsFromEnd <- function(primer, primerList, TmDelta=1, endMin, endMax, TmOffset=0, TmTarget=defaultTemps[[primer]], HmdTmLimit=35, HmdLimit=-7, HpTmLimit=35, HpLimit=-1, idName='primerId')
{
	ret <- data.table(Primer=primer, primerList[abs(Tm-TmTarget) <= TmDelta & End>=endMin & End<=endMax & HmdDeltaG>=HmdLimit & HpDeltaG >= HpLimit])
	if(nrow(ret)==0)
	{
		stop("No possible windows found with provided filters.")
	}
	ret[, c(idName):=1:.N]
	if(primer %in% c('F1c','F2c','F3c','B1','B2','B3','LF','LBc','PNAFc','PNAB'))
	{
		ret[, Seq:=revC(Seq, keepCase=T), by=c(idName)]
	}
	return(ret[!is.na(Start)])
}

# primerSetIsLegal <- function(primerSet)
# {
# 	
# 	if(any(duplicated(primerSet$Id))){return(FALSE)}
# 	if(any(duplicated(primerSet$Primer))){return(FALSE)}
# 	if(nrow(primerSet[Primer=='F3'])==1 && any(primerSet[Primer=='F3', list(End=Start+Len)] > ){return(FALSE)}
# }

tryPrimersWithPrimerSet <- function(primers, primerSet)
{
	primers <- primers[, data.table(rbindlist(list(primerSet[, names(primerSet) != 'primerId', with=F], .SD), use.names=T, fill=T)), by='primerId']
	return(primers[])
}

tryPrimersWithPrimerSets <- function(primers, primerSets)
{
	primerSets <- primerSets[Primer %!in% unique(primers$Primer)]
	if('primerId' %in% names(primerSets))
	{
		setnames(primerSets, old='primerId', new='primerId_temp')
	}
	ret <- primerSets[, tryPrimersWithPrimerSet(primers, .SD), by='primerId_temp']
	ret[, primerId:=.GRP, by=c('primerId','primerId_temp')]
	ret[, primerId_temp:=NULL]
	return(ret[])
}

# makeIdString <- function(primers)
# {
# 	paste.cols()
# }

tryPrimersAdjacentToPrimerSets <- function(primerSets, primerList, right=T, relToPrimer, newPrimer, TmDelta=1, minSpace, maxSpace, TmOffset=0, TmTarget=defaultTemps[[newPrimer]], HmdTmLimit=35, HmdLimit=-7, HpTmLimit=35, HpLimit=-4)
{
	if('primerId' %in% names(primerSets))
	{
		setnames(primerSets, old='primerId', new='primerId_temp')
	}
	if(length(relToPrimer) > 1)
	{
		ret <- primerSets[, tryPrimersWithPrimerSet(getPossibleWindowsFromBounds(primer=newPrimer,
																										primerList=primerList, 
																										TmDelta=TmDelta,
																										bMin=(min(End[Primer %in% relToPrimer]) + 1) + minSpace[1],
																										bMax=(max(Start[Primer %in% relToPrimer]) + 1) - ifelse(length(minSpace)>1, minSpace[2], minSpace[1]),
																										TmOffset=TmOffset,
																										TmTarget=TmTarget,
																										HmdTmLimit=HmdTmLimit,
																										HmdLimit=HmdLimit,
																										HpTmLimit=HpTmLimit,
																										HpLimit = HpLimit),
																  .SD), by='primerId_temp']
	}
	else if(right)
	{
		ret <- primerSets[, tryPrimersWithPrimerSet(getPossibleWindowsFromStart(primer=newPrimer,
																										primerList=primerList, 
																										TmDelta=TmDelta,
																										startMin=(End[Primer==relToPrimer] + 1) + minSpace,
																										startMax=(End[Primer==relToPrimer] + 1) + maxSpace,
																										TmOffset=TmOffset,
																										TmTarget=TmTarget,
																										HmdTmLimit=HmdTmLimit,
																										HmdLimit=HmdLimit,
																										HpTmLimit=HpTmLimit,
																										HpLimit = HpLimit),
																  .SD), by='primerId_temp']
	}
	else
	{
		ret <- primerSets[, tryPrimersWithPrimerSet(getPossibleWindowsFromEnd(primer=newPrimer,
																									 primerList=primerList, 
																									 TmDelta=TmDelta,
																									 endMin=(Start[Primer==relToPrimer] - 1) - maxSpace,
																									 endMax=(Start[Primer==relToPrimer] - 1) - minSpace,
																									 TmOffset=TmOffset,
																									 TmTarget=TmTarget,
																									 HmdTmLimit=HmdTmLimit,
																									 HmdLimit=HmdLimit,
																									 HpTmLimit=HpTmLimit,
																									 HpLimit=HpLimit),
																  .SD), by='primerId_temp']
	}
	ret[, primerId:=.GRP, by=c('primerId','primerId_temp')]
	ret[, primerId_temp:=NULL]
	return(ret[])
}

assemblePrimers <- function(primerSet, polyNT='a', n=3)
{
	F1r <- primerSet[Primer=='F1']
	F1cr <- primerSet[Primer=='F1c']
	F2r <- primerSet[Primer=='F2']
	B1cr <- primerSet[Primer=='B1c']
	B2cr <- primerSet[Primer=='B2c']
	B2r <- primerSet[Primer=='B2']
	if(nrow(B2cr)==0)
	{
		B2cr <- copy(B2r)
		B2cr[, Seq:=revC(Seq, keepCase = T)]
	}
	else
	{
		B2r <- copy(B2cr)
		B2r[, Seq:=revC(Seq, keepCase = T)]
	}
	if(nrow(F1r)==0)
	{
		F1r <- F1cr
		F1r[, Seq:=revC(Seq, keepCase = T)]
	}
	else
	{
		F1cr <- F1r
		F1cr[, Seq:=revC(Seq, keepCase = T)]
	}
	F1c <- s2c(F1cr$Seq)
	F2 <- s2c(F2r$Seq)
	B1c <- s2c(B1cr$Seq) # rpoB_500_c[primerSet[Primer=='B1c']$Start:primerSet[Primer=='B1c']$End]
	B2 <- s2c(B2r$Seq)
	B2c <- s2c(B2cr$Seq) # rpoB_500_c[primerSet[Primer=='B2c']$Start:primerSet[Primer=='B2c']$End]
	
	FIP <- paste(c(F1c, rep(polyNT, n), F2), collapse='') 
	BIP <- paste(c(B1c, rep(polyNT, n), B2), collapse='')
	DBF <- paste(c(rep(polyNT, n), F2, rpoB_500_c[(F2r$End+1):(F1r$Start-1)]), collapse='')
	DBB <- paste(c(rpoB_500_c[(B1cr$End+1):(B2cr$Start-1)], B2c, revC(c(rep(polyNT, n)), keepCase=T)), collapse='')
	ret <- data.table(Primer=c('F1c','F2','B1c','B2','FIP','BIP','DBF','DBB'), 
							Len=c(F1r$Len, F2r$Len, B1cr$Len, B2cr$Len, calcLen(FIP), calcLen(BIP), calcLen(DBF), calcLen(DBB)), 
							Tm=c(F1r$Tm, F2r$Tm, B1cr$Tm, B2cr$Tm, getTm(FIP), getTm(BIP), NA, NA), 
							HmdDeltaG=c(F1r$HmdDeltaG, F2r$HmdDeltaG, B1cr$HmdDeltaG, B2cr$HmdDeltaG, getHmd(FIP)$HmdDeltaG, getHmd(BIP)$HmdDeltaG, NA, NA), 
							Start=c(F1r$Start, F2r$Start, B1cr$Start, B2cr$Start, F2r$End-calcLen(FIP)+1, B2cr$Start, F2r$Start-n, B1cr$End+1), 
							End=c(F1r$End, F2r$End, B1cr$End, B2cr$End, F2r$End, B2cr$Start+calcLen(BIP)-1, F1r$Start-1, B2cr$End+n), 
							Seq=c(paste(F1c, collapse=''), paste(F2, collapse=''), paste(B1c, collapse=''), paste(B2, collapse=''), FIP, BIP, DBF, DBB),
							Id=c(F1r$Id,
								  F2r$Id,
								  B1cr$Id,
								  B2cr$Id,
								  paste(F1r$Id, '.', F2r$Id, sep=''), 
								  paste(B1cr$Id, '.', B2cr$Id, sep=''),
								  paste(F1r$Id, '.', F2r$Id, sep=''), 
								  paste(B1cr$Id, '.', B2cr$Id, sep=''))
	)
	ret <- rbindlist(list(ret, primerSet[Primer %!in% unique(ret$Primer), names(ret), with=F]))
	return(ret)
}

##### Plotting #####

getGCFractionPlot <- function(settings)
{
	plotGC(seq = settings$seq,
			 fullGC = calcGCFraction(),
			 dbSeq = settings$DBAll,
			 dbStart = settings$DBStart,
			 dbEnd = settings$DBEnd,
			 starts = starts(), 
			 ends   = stops(),
			 sensePrimerNames = sensePrimerNames)
}

getHairpinPlot <- function(settings)
{
	plotHairpin(settings$seq,
					fullHairpinEnergy=calcHairpinEnergy(),
					dbSeq = settings$DBAll,
					dbStart = settings$DBStart,
					dbEnd = settings$DBEnd,
					starts = starts(), 
					ends   = stops(),
					sensePrimerNames = sensePrimerNames)
}

getEnergiesPlot <- function(settings)
{
	plotEnergies(settings$results2)
}

getTmPlot <- function(settings)
{
	plotTm(settings$results3)
}

##### Stats/Table Generating Functions #####

#' calcHairpinEnergy
#'
#' @param seq characgter vector of individual nucleotides
#'
#' @return
#' @export
#'
#' @examples
calcHairpinEnergy <- function(seq, params=paramsDNA$Hp)
{
	rollapply(seq, width=25, FUN=function(x){return(getHp(x, params=params)$HpDeltaG)}, partial=T, align='center')
}
	
#' calcGCFraction
#'
#' @param seq character vector of individual nucleotides
#'
#' @return
#' @export
calcGCFraction <- function(seq)
{
	gc <- seq %in% c("g","G","c","C")
	rollapply(gc, width=20, FUN=mean, partial=T, align='center')
}

getPrimerStats <- function(x)
{
	hp <- getHp(x)
	Tm <- getTm(x)
	ret <- list(Tm=sig.digits(Tm, 3),
					HpTm=sig.digits(hp$temp, 3),
					HpDeltaG=sig.digits(hp$dg, 3),
					HpStruct=ifelse(hp$structure_found, hp$structure, ''))
	return(ret)
}

getDimerComboStats <- function(primers, seqs)
{
	duh2 <- data.table(getCombos(primers, colnames=c('P1','P2')), getCombos(seqs, c('Seq1','Seq2')))
	duh2 <- duh2[duh2[, getDimer(Seq1, Seq2), by=c('P1','P2')], on=c('P1','P2')]
	valNames <- names(duh2)
	valNames <- valNames[valNames %!in% c('P1','P2','Seq1','Seq2')]
	blah <- rbindlist(list(duh2, data.table(Seq1=duh2$Seq2, Seq2=duh2$Seq1, P1=duh2$P2, P2=duh2$P1, duh2[, mget(valNames)])), use.names=T)
	blah <- unique(blah)
	blah[, max.i:=which.min(DimerDeltaG)[1], by='P1']
	blah <- blah[, .SD[max.i[1]], by='P1'][,c('P1','P2',valNames), with=F]
	return(blah[match(primers, P1), mget(c('P2', valNames))])
}

getUniqueDimerComboStats <- function(primers, seqs, func=daFunc, mv=50.0, dv=4.45, dntp=0, temp_c=65, dna=100, tm_method='SantaLucia', salt_method='Schildkraut')
{
	duh2 <- data.table(getCombos(primers, colnames=c('P1','P2')), getCombos(seqs, c('Seq1','Seq2')))
	duh2 <- duh2[duh2[, func(Seq1, Seq2, mv=mv, dv=dv, dntp=dntp, temp_c=temp_c, dna=dna, tm_method=tm_method, salt_method=salt_method), by=c('P1','P2')], on=c('P1','P2')]
	valNames <- names(duh2)
	valNames <- valNames[valNames %!in% c('P1','P2','Seq1','Seq2')]
	# blah <- rbindlist(list(duh2, data.table(Seq1=duh2$Seq2, Seq2=duh2$Seq1, P1=duh2$P2, P2=duh2$P1, duh2[, mget(valNames)])), use.names=T)
	# blah <- unique(blah)
	blah[, max.i:=which.min(DimerDeltaG)[1], by='P1']
	blah <- blah[, c(list(P2=P2[max.i]), lapply(mget(valNames), function(x){x[max.i[1]]})), by='P1']
	return(blah[match(primers, P1), mget(c('P2', valNames))])
}

calcResultsTables <- function(settings)
{
	ret <- data.table(Primer=finalPrimerNames,
							Seq=as.character(sapply(lapply(finalPrimerNames, function(x){settings$NTs[[x]]}), paste, collapse='')),
							Sense=c('Sense','Antisense','Sense','Antisense','Sense','Antisense','NA','NA','Antisense','Sense','Sense','Antisense','Sense','Sense'))
	ret[, Stability3p:=getStability(Seq, threePrimeEnd=T), by='Primer']
	ret[, Stability5p:=getStability(Seq, threePrimeEnd=F), by='Primer']
	ret[, Start:=as.numeric(NA)]
	ret[, End:=as.numeric(NA)]
	ret[Primer == 'F3', Start:=settings$Start$F3]
	ret[Primer == 'F3', End:=getEndSetting('F3', settings)]
	ret[Primer == 'B3', Start:=settings$Start$B3c]
	ret[Primer == 'B3', End:=getEndSetting('B3c', settings)]
	ret[Primer == 'F2', Start:=settings$Start$F2]
	ret[Primer == 'F2', End:=getEndSetting('F2', settings)]
	ret[Primer == 'F1c', Start:=settings$Start$F1]
	ret[Primer == 'F1c', End:=getEndSetting('F1', settings)]
	ret[Primer == 'B2', Start:=settings$Start$B2c]
	ret[Primer == 'B2', End:=getEndSetting('B2c', settings)]
	ret[Primer == 'B1c', Start:=settings$Start$B1c]
	ret[Primer == 'B1c', End:=getEndSetting('B1c', settings)]
	ret[Primer == 'LF', Start:=settings$Start$LFc]
	ret[Primer == 'LF', End:=getEndSetting('LFc', settings)]
	ret[Primer == 'LB', Start:=settings$Start$LB]
	ret[Primer == 'LB', End:=getEndSetting('LB', settings)]
	ret[Primer == 'PNAF', Start:=settings$Start$PNAF]
	ret[Primer == 'PNAF', End:=getEndSetting('PNAF', settings)]
	ret[Primer == 'PNAB', Start:=settings$Start$PNABc]
	ret[Primer == 'PNAB', End:=getEndSetting('PNABc', settings)]
	ret[Primer == 'FIP', Start:=getEndSetting('F2', settings) - calcLen(settings$Seq$FIP) + 1] # F2r$End-calcLen(FIP)+1
	ret[Primer == 'FIP', End:=getEndSetting('F2', settings)] # F2r$End
	ret[Primer == 'BIP', Start:=settings$Start$B2c] # B2cr$Start
	ret[Primer == 'BIP', End:=settings$Start$B2c + calcLen(settings$Seq$BIP) - 1] # B2cr$Start+calcLen(BIP)-1
	ret[Primer == 'DBF', Start:=settings$Start$F2 - calcLen(settings$FIPlinker)] # F2r$Start-n
	ret[Primer == 'DBF', End:=Start+calcLen(Seq) - 1] # F1r$Start-1
	ret[Primer == 'DBB', Start:=getEndSetting('B1c', settings) + 1] # B1cr$End+1
	ret[Primer == 'DBB', End:=Start + calcLen(Seq) - 1] # B2cr$End+n
	ret[, Start:=as.integer(Start)]
	ret[, End:=as.integer(End)]
	
	ret[, Len:=length(s2c(Seq)), by='Primer']
	
	ret[, Tm:=getTm(Seq), by='Primer']
	ret[Primer %!in% c('DBB','DBF'), c('HpTm','HpDeltaG','HpStruct'):=getHp(Seq), by=c('Primer')]
	ret[, c('HpTm','HpDeltaG','HpStruct'):=getDBHps(.SD, width=25)]
	
	checked <- gsub('c','',names(settings$Check)[as.logical(settings$Check)])
	checked <- checked[checked %!in% c('F1','F2','B1','B2')]
	primerBaseNames <- gsub('c','',ret$Primer)
	primersToCalc <- primerBaseNames %in% c(checked, 'FIP', 'BIP')
	ret[primersToCalc, c('P2','DimerTm','DimerDeltaG','DimerStruct'):=getDimerComboStats(Primer, Seq)]
	ret[, KeyEndStability:=ifelse(Primer %in% c('F1c','B1c'), Stability5p, Stability3p)]
	
	# Data for energy plot
	ret2 <- melt.data.table(data=ret, id.vars='Primer', measure.vars=c('KeyEndStability','HpDeltaG','DimerDeltaG'))
	ret2[, value:=suppressWarnings(as.numeric(value))]
	ret2[value > 0, value:=0]
	ret2[variable=='KeyEndStability', value:=-1*value]
	ret2[, Primer:=factor(Primer, levels=finalPrimerNames)]
	
	# Data for Tm plot
	ret3 <- ret[, c('Primer','Tm','HpTm')]
	ret3[Primer %in% c('DBB','DBF'), Tm:=HpTm]
	ret3[, HpTm:=NULL]
	# ret3[, Tm:=suppressWarnings(as.numeric(Tm))]
	primerColors <- getAllPrimerColors()
	setkey(ret3, Primer)
	setkey(primerColors, Primer)
	ret3 <- primerColors[ret3]
	ret3[, Primer:=factor(Primer, levels=finalPrimerNames)]
	return(list(results=ret, results2=ret2, results3=ret3))
}

##### Settings Functions #####

updateSettingsGroupItem <- function(settings, group, id, val)
{
	if(length(val) > 0)
	{
		if(length(settings[[group]]) == 0)
		{
			settings[[group]] <- list()
		}
		if(length(settings[[group]][[id]])==0 || any(settings[[group]][[id]] != val))
		{
			settings[[group]][[id]] <- val
			return(TRUE)
		}
	}
	return(FALSE)
}

updateSettingsItem <- function(settings, id, val)
{
	if(length(settings)>0 && length(val)>0)
	{
		if(length(settings[[id]])==0 || length(settings[[id]]) != length(val) || any(settings[[id]] != val))
		{
			settings[[id]] <- val
			return(TRUE)
		}
	}
	return(FALSE)
}

getStartSetting <- function(id, settings)
{
	return(settings[['Start']][[id]])
}

getLenSetting <- function(id, settings)
{
	return(settings[['Len']][[id]])
}

getEndSetting <- function(id, settings)
{
	return(calcEnd(settings[['Start']][[id]], settings[['Len']][[id]]))
}

initLocs <- function(settings)
{
	return(round(seq(1, length(vals$seq)-20, length.out=8)))
}

starts <- function(settings)
{
	as.numeric(sapply(sensePrimerNames, function(x){settings$Start[[x]]}))
}

stops <- function(settings)
{
	as.numeric(sapply(sensePrimerNames, function(x){settings$Start[[x]]}))+as.numeric(lapply(sensePrimerNames, function(x){settings$Len[[x]]}))-1
}

InitList <- data.table.expand.grid(Len=11:28)
InitList <- InitList[, getSeqsWithLen(rpoB_500_c, len=.BY[[1]], align='left'), by='Len']
InitList[, End:=Start+Len-1]
InitList[, Tm:=getTm(Seq), by=c('Start','Len')]
InitList[, c('HmdTm','HmdDeltaG','HmdStruct'):=getHmd(Seq), by=c('Start','Len')]
InitList[, c('HpTm','HpDeltaG','HpStruct'):=getHp(Seq), by=c('Start','Len')]
InitList[, HpDeltaG:=ifelse(is.na(HpDeltaG), max(HpDeltaG, na.rm=T), HpDeltaG)]
InitList[, Id:=as.character(1:.N)]