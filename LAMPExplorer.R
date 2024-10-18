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


library(shiny)
library(shinythemes)
library(data.table)
library(zoo)
library(grDevices)
library(XNAString)
library(Biostrings)
library(seqinr)
library(primer3)
library(pwalign)
library(marker)
library(ggplot2)
library(scales)
library(purrr)

##### Constants #####
rpoB_500 <- 'atcgaccacttcggcaaccgccgcctgcgtacggtcggcgagctgatccaaaaccagatccgggtcggcatgtcgcggatggagcgggtggtccgggagcggatgaccacccaggacgtggaggcgatcacaccgcagacgttgatcaacatccggccggtggtcgccgcgatcaaggagttcttcggcaccagccagctgagccaattcatgGACcagaacaacccgctgtcggggttgaccCACaagcgccgactgTCGgcgctggggcccggcggtctgtcacgtgagcgtgccgggctggaggtccgcgacgtgcacccgtcgcactacggccggatgtgcccgatcgaaacccctgaggggcccaacatcggtctgatcggctcgctgtcggtgtacgcgcgggtcaacccgttcgggttcatcgaaacgccgtaccgcaaggtggtcgacggcgtggttagcgacgagatcgtgtacctgaccgccgacgagga'
rpoB_500_c <- s2c(rpoB_500)

# Per PrimerExplorerV5 defaults
# 5' Stability: < -3
# 3' Stability: < -4
# Dimer Stability: > -2.5
# F1c B1c Tm: 66 ± 2
# Other Tm: 61 ± 2
# F1c-F2 Space: 40-60
# F3-F2 Space: 0-20
# F1c-B1c Space: 0-100
defaultTemps <- list(F3=61, F2=61, B3c=61, B2c=61, LFc=61, LB=61, PNAF=61, PNABc=61, F1=66, B1c=66)

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
##### Helper Functions #####

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

updateList <- function(items, ...)
{
	# Default values in list1
	# Overriding values in list2
	args <- list(...)
	if(!all(names(args) %in% names(items)))
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
		browser()
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
	temp <- pwalign::pairwiseAlignment(primer, target, type='local', gapOpening=1000000, gapExtension=1000000)
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
					 HpStruct=ifelse(!is.na(ret2$dg), ret2$structure, ''))
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
	DBFSeqLinHp <- getHpTable(DBFSeqLin, width=width, align='center', trim=F)
	DBFSeqLinHp[, bp:=bp-min(bp)+1]
	DBFSeqCircHp <- getHpTable(DBFSeqCirc, width=width, align='center', trim=T)
	DBFSeqCircHp[, bp:=bp-min(bp)+1]
	DBFSeqLinHp[, HpType:='Lin']
	DBFSeqCircHp[, HpType:='Circ']
	DBFTable <- rbindlist(list(DBFSeqCircHp, DBFSeqLinHp))
	DBFI <- which.min(DBFTable$HpDeltaG)
	DBFTable[HpStruct != '', HpStruct:=paste(HpType, ':', bp-width%/%2, '-', bp+width%/%2, ':', HpStruct, sep='')]
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

getAllPrimerColors <- function()
{
	return(data.table(Primer=c('F3',
										'F3c',
										'F2',
										'F2c',
										'F1',
										'F1c',
										'B1',
										'B1c',
										'B2',
										'B2c',
										'B3',
										'B3c',
										'LF',
										'LFc',
										'LB',
										'LBc',
										'FIP',
										'BIP',
										'DBF',
										'DBB',
										'PNAF',
										'PNAFc',
										'PNAB',
										'PNABc'),
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

plotEnergies <- function(results)
{
	ret <- ggplot(data=results[variable != 'KeyEndStability'], aes( x=Primer, y=value, fill=variable)) +
		geom_col() +
		geom_col(data=results[variable == 'KeyEndStability']) +
		geom_hline(yintercept=c(4,-2)) +
		labs(x='Sequence', y='Energy [kJ/mol]') +
		scale_y_continuous(limits=c(-15,10),oob = rescale_none) +
		theme(axis.text=element_text(size=rel(2.0)),
				axis.title=element_text(size=rel(2.0)),
				legend.text=element_text(size=rel(2.0)),
				legend.title=element_text(size=rel(2.0)))
	return(ret)
}

plotTm <- function(results)
{
	ret <- ggplot(data=results, aes(x=Primer, y=Tm, fill=Primer)) + 
		geom_col() +
		scale_fill_manual(values = getPrimerColor(as.character(levels(results$Primer)))) + 
		geom_hline(yintercept=c(defaultTemps$F1, defaultTemps$F2)) +
		scale_y_continuous(limits=c(40,80),oob = rescale_none) +
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

##### Calculated Constants #####

InitList <- data.table.expand.grid(Len=11:28)
InitList <- InitList[, getSeqsWithLen(rpoB_500_c, len=.BY[[1]], align='left'), by='Len']
InitList[, Tm:=getTm(Seq), by=c('Start','Len')]
InitList[, HmdDeltaG:=getHmd(Seq), by=c('Start','Len')]
InitList[, c('HpTm','HpDeltaG','HpStruct'):=getHp(Seq), by=c('Start','Len')]
InitList[is.na(HpDeltaG), HpDeltaG:=999]
InitList[, End:=Start+Len-1]
InitList[, Id:=1:.N]

# Define UI for application that draws a histogram
ui <- fluidPage(
	theme = shinytheme("darkly"),
	tags$head(
		# Note the wrapping of the string in HTML()
		tags$style(HTML(paste("
      markF3	{background-color: ", getPrimerColor('F3'), "; 	color: black;}
      markF2	{background-color: ", getPrimerColor('F2'), ";	color: black;}
      markLFc	{background-color: ", getPrimerColor('LFc'), ";	color: black;}
      markF1	{background-color: ", getPrimerColor('F1'), ";	color: black;}
      markB1c	{background-color: ", getPrimerColor('B1c'), ";	color: black;}
      markLB	{background-color: ", getPrimerColor('LB'), ";	color: black;}
      markB2c	{background-color: ", getPrimerColor('B2c'), ";	color: black;}
      markB3c	{background-color: ", getPrimerColor('B3c'), ";	color: black;}
      markPNAF {background-color: ", getPrimerColor('PNAF'), ";	color: black;}
      markPNABc{background-color: ", getPrimerColor('PNABc'), ";	color: black;}
      ", collapse='')))
	),
	tags$script("
            Shiny.addCustomMessageHandler('txt', function (txt) {
                navigator.clipboard.writeText(txt);
            });
        "), # this is new
	
	useMarker(),
	
	# Application title
	titlePanel("LAMP Primer Optimizer"),
	
	# Sidebar with a slider input for number of bins 
	sidebarLayout(
		sidebarPanel(width=5,
						 fluidPage(fluidRow(column(5, h4("Target Sequence:")),
						 						 column(7, fileInput("importSettings", '', placeholder='Settings.RData',
						 						 						  multiple = FALSE,
						 						 						  accept = c(".RData"))))
						 ),
						 tags$textarea(id="Seq", rows=10, cols=60, 
						 				  'atcgaccacttcggcaaccgccgcctgcgtacggtcggcgagctgatccaaaaccagatccgggtcggcatgtcgcggatggagcgggtggtccgggagcggatgaccacccaggacgtggaggcgatcacaccgcagacgttgatcaacatccggccggtggtcgccgcgatcaaggagttcttcggcaccagccagctgagccaattcatgGACcagaacaacccgctgtcggggttgaccCACaagcgccgactgTCGgcgctggggcccggcggtctgtcacgtgagcgtgccgggctggaggtccgcgacgtgcacccgtcgcactacggccggatgtgcccgatcgaaacccctgaggggcccaacatcggtctgatcggctcgctgtcggtgtacgcgcgggtcaacccgttcgggttcatcgaaacgccgtaccgcaaggtggtcgacggcgtggttagcgacgagatcgtgtacctgaccgccgacgagga'
						 ),
						 hr(),
						 fluidPage(fluidRow(column(4, selectInput("linker", "Linker NT", choices=c('a','g','t','c'), selected='a')),
						 						 column(4, numericInput("polyT", "Linker Length", 3)),
						 						 column(4, numericInput("stabilityN", "End Stability Bp's:", 5))),
						 			 textInput('revCTool', 'Rev. Compl. Tool', placeholder='Auto-copy RevC to Clipboard'),
						 			 textOutput('revCOutput'),
						 			 actionButton("importButton", "Import Clipboard")),
						 hr(),
						 # textOutput("Warnings"),
						 # actionButton("Update", "Update Plots"),
						 uiOutput("F3"),
						 uiOutput("F2"),
						 uiOutput("F1"),
						 uiOutput("B1c"),
						 uiOutput("B2c"),
						 uiOutput("B3c"),
						 uiOutput("LFc"),
						 uiOutput("LB"),
						 uiOutput("PNAF"),
						 uiOutput("PNABc")
		),
		
		# Show a plot of the generated distribution
		mainPanel(width=7,
					 uiOutput("coloredSeq"),
					 wellPanel(style="background-color: white;", 
					 			 plotOutput("GCPlot"),
					 			 plotOutput("HairpinPlot"),
					 			 plotOutput("EnergyPlot"),
					 			 plotOutput("TmPlot"),
					 ),
					 hr(),
					 fluidPage(fluidRow(column(2, downloadButton('download',"Download Table")),
					 						 column(3, textInput('downloadName',NULL, value='Primer Set 1')),
					 						 column(1, h4(".csv")),
					 						 column(6, ))),
					 tableOutput("ResultsTable"),
		)
	)
)

renderPrimerControl <- function(seq, controlId, initLocs=c(1), Loc=1, initLen=20)
{
	initStart <- ifelse(Loc > 0, initLocs[Loc], 1)
	fluidPage({
		fluidRow(style = "height:20px;",
					column(1, HTML(paste(tags$h4(tags$span(style=paste("color: ", getPrimerColor(controlId), sep=''), controlId), sep='')))),
					column(1, checkboxInput(paste0(controlId, 'Check'), '', value=ifelse(grepl("PNA", controlId), F, T))), #value=grepl('[12]', controlId))),
					column(2, numericInput(paste0(controlId, 'Start'), 'Start', value=initStart)),
					column(2, numericInput(paste0(controlId, 'Len'), 'Length', value=initLen)),
					column(6, textInput(paste0(controlId, 'NTs'), 'Seq', value=paste(seq[initStart:(initStart+initLen-1)], collapse='')))
		)
	})
}

getInputs <- function(prefixes, suffix, input)
{
	lapply(prefixes, getInput, suffix=suffix, input=input)
	# input[paste(prefixes, suffix, sep='')]
}

getInput <- function(prefix, suffix, input)
{
	input[[paste(prefix, suffix, sep='')]]
}

updateSettingsGroupItem <- function(id, group, input, settings)
{
	# print(paste(c(id, group, input[[paste(id, group, sep='')]], settings[[group]])))
	val <- input[[paste(id, group, sep='')]]
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

updateSettingsItem <- function(id, val, settings, group=NULL)
{
	if(is.null(group))
	{
		if(length(settings)>0 && length(val)>0)
		{
			if(length(settings[[id]])==0 || length(settings[[id]]) != length(val) || any(settings[[id]] != val))
			{
				settings[[id]] <- val
				return(TRUE)
			}
		}
	}
	else
	{
		if(length(settings)>0 && length(val)>0)
		{
			if(length(settings[[group]])==0)
			{
				settings[[group]] <- list()
			}
			if(length(settings[[group]][[id]])==0 || any(settings[[group]][[id]] != val))
			{
				settings[[group]][[id]] <- val
				return(TRUE)
			}
		}
	}
	return(FALSE)
}

getStartSetting <- function(id, settings)
{
	return(settings[[Start]][[id]])
}

getEndSetting <- function(id, settings)
{
	return(settings[['Start']][[id]] + settings[['Len']][[id]] - 1)
}

getStart <- function(input, controlId)
{
	return(input[[paste0(controlId, 'Start')]])
}

getEnd <- function(input, controlId)
{
	return(input[[paste0(controlId, 'Start')]] + (input[[paste0(controlId, 'Len')]] - 1))
}

silenceNTUpdate <- 0
silenceStartUpdate <- 0
silenceLenUpdate <- 0

# Define server logic required to draw a histogram
server <- function(input, output, session) {
	
	settings <- reactiveValues(render=1)
	
	mySeqHTML <- reactive({
		req(settings$seq)
		markup=p(id='text-to-mark', style = "word-wrap: break-word;", paste(settings$seq, collapse=''))
	})
	
	fullHairpinEnergy <- reactive({
		req(settings$seq, allLegal())
		print("Updating hairpin energy plot.")
		rollapply(settings$seq, width=25, FUN=function(x){return(min(0, getDefault(getHp(x)$HpDeltaG, 0, test=is.na)))}, partial=T, align='center')
	})
	
	fullGC <- reactive({
		req(settings$seq)
		print("Updating fullGC.")
		gc <- settings$seq %in% c("g","G","c","C")
		rollapply(gc, width=20, FUN=mean, partial=T, align='center')
	})
	
	initLocs <- reactive({
		req(settings$seq)
		print("Updating initLocs.")
		return(round(seq(1, length(settings$seq)-20, length.out=8)))
	})
	
	output$coloredSeq <- renderUI({
		req(mySeqHTML())
		print("Rendering HTML Seq.")
		mySeqHTML()
	})
	
	observeEvent( list(settings$seq, settings$render > 0),{
		req(settings$seq)
		if(calcLen(settings$seq)==0){settings$render <- settings$render + 1}
		print("Rendering primer controls.")
		output$F3 <- renderUI(isolate(renderPrimerControl(settings$seq, 'F3', initLocs=initLocs(), Loc=2)))
		output$F2 <- renderUI(isolate(renderPrimerControl(settings$seq, 'F2', initLocs=initLocs(), Loc=3)))
		output$LFc <- renderUI(isolate(renderPrimerControl(settings$seq, 'LFc', initLocs=initLocs(), Loc=0)))
		output$F1 <- renderUI(isolate(renderPrimerControl(settings$seq, 'F1', initLocs=initLocs(), Loc=4)))
		output$B1c <- renderUI(isolate(renderPrimerControl(settings$seq, 'B1c', initLocs=initLocs(), Loc=5)))
		output$LB <- renderUI(isolate(renderPrimerControl(settings$seq, 'LB', initLocs=initLocs(), Loc=0)))
		output$B2c <- renderUI(isolate(renderPrimerControl(settings$seq, 'B2c', initLocs=initLocs(), Loc=6)))
		output$B3c <- renderUI(isolate(renderPrimerControl(settings$seq, 'B3c', initLocs=initLocs(), Loc=7)))
		output$PNAF <- renderUI(isolate(renderPrimerControl(settings$seq, 'PNAF', initLocs=initLocs(), Loc=0)))
		output$PNABc <- renderUI(isolate(renderPrimerControl(settings$seq, 'PNABc', initLocs=initLocs(), Loc=0)))
		# print(settings$render)
		if(all(sapply(sensePrimerNames, function(x){calcLen(input[[paste(x, 'NTs', sep='')]]) > 0}))){
			settings$render <- 0
			silenceLenUpdate <<- 0
			silenceStartUpdate <<- 0
			silenceNTUpdate <<- 0
		}
	})
	
	observeEvent(list(mySeqHTML(), settings$NTs, settings$Check), {
		# Highlight output
		# print(paste("Marking ", controlId, ": ", input[[paste0(controlId, "NTs")]]))
		req(mySeqHTML(), all(sapply(sensePrimerNames, function(x){!is.null(settings$NTs[[x]])})), all(sapply(sensePrimerNames, function(x){!is.null(settings$Check[[x]])})))
		print("Updating markers")
		my_marker <- marker$new("#text-to-mark")
		my_marker$unmark()
		lapply(names(settings$Check), function(x){if(settings$Check[[x]]){my_marker$mark(settings$NTs[[x]], element=paste0("mark", x))}})
	}, ignoreInit = F)
	
	sensePrimerNames <- c('F3','F2','F1','LFc','B1c','B2c','B3c','LB','PNAF','PNABc')
	
	dependentPrimerNames <- c('F1c', 'B2', 'B3', 'LF', 'LB', 'PNAB', 'FIP', 'BIP', 'DBF', 'DBB')
	
	starts <- reactive({
		as.numeric(sapply(sensePrimerNames, function(x){settings$Start[[x]]}))
	})
	
	stops <- reactive({
		as.numeric(sapply(sensePrimerNames, function(x){settings$Start[[x]]}))+as.numeric(lapply(sensePrimerNames, function(x){settings$Len[[x]]}))-1
	})
	
	primerColors <- reactive({
		getPrimerColor(sensePrimerNames)
	})
	
	observeEvent(list(input$Seq, settings$render > 0), {
		updateSettingsItem('seq', {
			temp <- s2c(input$Seq)
			temp[temp %in% c('a','g','t','c','A','G','T','C')]
			if(input$Seq != paste(temp, collapse=''))
			{
				warning('Non A, G, T, C, a, g, t, c, characters found. Repairing sequence')
				if(settings$seq != temp)
				{
					settings$seq <- temp
				}
				if(input$Seq == '')
				{
					updateTextAreaInput(session, inputId='Seq', value=paste(temp, collapse=''))
				}
			}
			# print(temp)
			temp
		}, settings)
	})
	
	observeEvent(getInputs(sensePrimerNames, 'Check', input), {
		do.call('req', lapply(getInputs(sensePrimerNames, 'Check', input), function(x){length(x)>0}))
		lapply(sensePrimerNames, updateSettingsGroupItem, group='Check', input=input, settings)
	})
	
	observeEvent(getInputs(sensePrimerNames, 'Start', input), {
		do.call('req', getInputs(sensePrimerNames, 'Start', input))
		lapply(sensePrimerNames, updateSettingsGroupItem, group='Start', input=input, settings)
	})
	
	observeEvent(getInputs(sensePrimerNames, 'Len', input), {
		do.call('req', getInputs(sensePrimerNames, 'Len', input))
		lapply(sensePrimerNames, updateSettingsGroupItem, group='Len', input=input, settings)
	})
	
	observeEvent(input$linker, {
		req(input$linker)
		updateSettingsItem('linker', input$linker, settings)
	})
	
	# Watch all the inputs and create a ground truth stored version of everything
	observeEvent(list(settings$seq, settings$polyT, settings$linker, getInputs(sensePrimerNames, 'NTs', input)), {
		tempPNABc <- input$PNABc
		do.call('req', getInputs(sensePrimerNames, 'NTs', input))
		req(settings$seq, settings$polyT)
		lapply(sensePrimerNames, updateSettingsGroupItem, group='NTs', input=input, settings=settings)
		updateSettingsItem('F1c', revC(settings$NTs$F1, keepCase=T), settings, group='NTs')
		updateSettingsItem('B2', revC(settings$NTs$B2c, keepCase=T), settings, group='NTs')
		updateSettingsItem('B3', revC(settings$NTs$B3c, keepCase=T), settings, group='NTs')
		updateSettingsItem('LF', revC(settings$NTs$LFc, keepCase=T), settings, group='NTs')
		updateSettingsItem('LB', settings$NTs$LB, settings, group='NTs')
		updateSettingsItem('FIP', {
			paste(c(settings$NTs$F1c, rep(settings$linker, settings$polyT), settings$NTs$F2), collapse='')	
		}, settings, group='NTs')
		updateSettingsItem('BIP', {
			paste(c(settings$NTs$B1c, rep(settings$linker, settings$polyT), settings$NTs$B2), collapse='')	
		}, settings, group='NTs')
		updateSettingsItem('DBF', paste(c(rep(settings$linker, settings$polyT), s2c(settings$NTs$F2), settings$seq[(getEndSetting('F2', settings)+1):(settings$Start$F1-1)]), collapse=''), settings, group='NTs')
		updateSettingsItem('DBB', paste(c(settings$seq[(getEndSetting('B1c', settings)+1):(settings$Start$B2c-1)], s2c(settings$NTs$B2c), revC(rep(settings$linker, settings$polyT), keepCase=T)), collapse=''), settings, group='NTs')
		updateSettingsItem('PNAF', {
			# If the PNA straddles the start of F2
			ifelse(settings$Start$PNAF < settings$Start$F2 && getEndSetting('PNAF', settings) >= settings$Start$F2, 
					 paste(
					 	{
					 		# Paste together the upstream portion of FIP with the downstream portion of the PNAF
					 		FIPN <- length(settings$NTs$FIP)
					 		F2N <- settings$Len$F2
					 		Offset <- settings$Start$F2-settings$Start$PNAF
					 		if(getEndSetting('PNAF', settings) > getEndSetting('F2', settings))
					 		{
					 			warning('Having a PNAF longer than F2 is not supported yet.')
					 			updateTextInput(session, 'PNAFNTs', settings$NTs$PNAF)
					 			s2c(settings$NTs$PNAF) # Don't change anything
					 		}
					 		else
					 		{
					 			c(settings$NTs$FIP[(FIPN-F2N-Offset):(FIPN-F2N)], # Chunk from start of PNAF up to start of F2 in FIP
					 			  settings$NTsPNAF[(Offset+1):settings$Len$PNAF]) # Remaining part of PNAF
					 		}
					 	}, 
					 	collapse=''),
					 settings$NTs$PNAF
			)
		}, settings, group='NTs')
		updateSettingsItem('PNAB', {
			# If the PNA straddles the end of B2c
			ifelse(settings$Start$PNABc < getEndSetting('B2c', settings) && getEndSetting('PNABc', settings) >= getEndSetting('B2c', settings), 
					 paste(
					 	{
					 		# Paste together the upstream portion of FIP with the downstream portion of the PNAF
					 		BIPN <- length(settings$NTs$BIP)
					 		B2cN <- settings$Len$B2c
					 		Offset <- getEndSetting('PNABc', settings)-getEndSetting('B2c', settings)
					 		if(settings$Start$PNABc < settings$Start$B2c)
					 		{
					 			stop('Having a PNAB longer than B2 is not supported yet.')
					 			# Reset the source PNABc value to what it was before this observed event
					 			settings$NTs$PNABc <- tempPNABc
					 			# Reset the input to the same value
					 			updateTextInput(session, 'PNABcNTs', tempPNABc)
					 			s2c(settings$NTs$PNAB) # Don't change anything by returning the same value
					 		}
					 		else
					 		{
					 			c(settings$NTs$BIP[(BIPN-B2cN-Offset):(BIPN-B2cN)], # Chunk from end of PNABc up to end of B2c in BIP
					 			  revC(settings$NTs$PNABc, keepCase=T)[1:(settings$Len$PNABc-Offset)], keepCase=T) # Remaining (i.e., beginning) part of revC(PNABc)
					 		}
					 	}, 
					 	collapse=''),
					 revC(settings$NTs$PNABc, keepCase=T))
		}, settings, group='NTs')
		updateSettingsItem('DBAll', {
			paste(c(settings$NTs$F1c, rep(settings$linker, input$polyT), settings$seq[settings$Start$F2:getEndSetting('B2c', settings)], rep(settings$linker, input$polyT), revC(settings$NTs$B1c, keepCase=T)), collapse='')
		}, settings)
		updateSettingsItem('DBStart', settings$Start$F2-settings$polyT-settings$Len$F1, settings)
		updateSettingsItem('DBEnd', (settings$Start$B2c + settings$Len$B2c - 1)+settings$polyT+settings$Len$B1c, settings)
	})
	
	observeEvent(getInputs(sensePrimerNames, 'Check', input), {
		do.call('req', getInputs(sensePrimerNames, 'Check', input))
		lapply(sensePrimerNames, updateSettingsGroupItem, group='Check', input=input, settings)
	})
	
	observeEvent(input$polyT, {
		updateSettingsItem('polyT', input$polyT, settings)
	})
	
	observeEvent(list(settings$Check, settings$NTs), {
		# Highlight output
		# print(paste("Marking ", controlId, ": ", input[[paste0(controlId, "NTs")]]))
		
		req(!is.null(settings$seq), !is.null(settings$Check), length(settings$Check) > 0)
		my_marker <- marker$new("#text-to-mark")
		my_marker$unmark()
		if(input$F3Check){ my_marker$mark(input[[paste0("F3", "NTs")]], element=paste0("markF3")) }
		if(input$F2Check){ my_marker$mark(input[[paste0("F2", "NTs")]], element=paste0("markF2")) }
		if(input$LFcCheck){ my_marker$mark(input[[paste0("LFc", "NTs")]], element=paste0("markLFc")) }
		if(input$F1Check){ my_marker$mark(input[[paste0("F1", "NTs")]], element=paste0("markF1")) }
		if(input$B1cCheck){ my_marker$mark(input[[paste0("B1c", "NTs")]], element=paste0("markB1c")) }
		if(input$LBCheck){ my_marker$mark(input[[paste0("LB", "NTs")]], element=paste0("markLB")) }
		if(input$B2cCheck){ my_marker$mark(input[[paste0("B2c", "NTs")]], element=paste0("markB2c")) }
		if(input$B3cCheck){ my_marker$mark(input[[paste0("B3c", "NTs")]], element=paste0("markB3c")) }
		if(input$PNAFCheck){ my_marker$mark(input[[paste0("PNAF", "NTs")]], element=paste0("markPNAF")) }
		if(input$PNABcCheck){ my_marker$mark(input[[paste0("PNABc", "NTs")]], element=paste0("markPNABc")) }
	}, ignoreInit = F)
	
	observeEvent(input$stabilityN, {
		
		updateSettingsItem('stabilityN', input$stabilityN, settings)	
	})
	
	# # Dumbbell regions
	# DBF <- reactive(paste(c(rep('t', settings$polyT), mySeq()[getStart(input, 'F2'):(getStart(input, 'F1')-1)]), collapse=''))
	# DBB <- reactive(paste(c(rep('t', input$polyT), mySeq()[getEnd(input, 'B2c'):(getEnd(input, 'B1c')-1)]), collapse=''))
	
	# # Combo primers
	# FIP <- reactive(paste(c(F1c(), rep('t', input$polyT), F2()), collapse=''))
	# BIP <- reactive(paste(c(B1c(), rep('t', input$polyT), B2()), collapse=''))
	
	# # PNAs
	# PNAF <- reactive({
	# 	# If the PNA straddles the start of F2
	# 	ifelse(getStart(input, 'F2') %in% c(getStart(input, 'PNAF'):getEnd(input, 'PNAF')),
	# 			 paste(
	# 			 	{
	# 			 		# Paste together the upstream portion of FIP with the downstream portion of the PNAF
	# 			 		FIPN <- length(s2c(FIP()))
	# 			 		F2N <- input$F2Len
	# 			 		Offset <- getStart(input, 'F2')-getStart(input, 'PNAF')
	# 			 		if(getEnd(input, 'PNAF') > getEnd(input, 'F2'))
	# 			 		{
	# 			 			stop('Having a PNAF longer than F2 is not supported yet.')
	# 			 		}
	# 			 		c(s2c(FIP())[(FIPN-F2N-Offset):(FIPN-F2N)], # Chunk from start of PNAF up to start of F2 in FIP
	# 			 		  s2c(input$PNAFNTs)[(Offset+1):input$PNAFLen]) # Remaining part of PNAF
	# 			 	}, 
	# 			 	collapse=''),
	# 			 input$PNAFNTs
	# 	)
	# })
	# PNAB <- reactive({
	# 	# If the PNA straddles the end of B2c
	# 	ifelse(getEnd(input, 'B2c') %in% c(getStart(input, 'PNABc'):getEnd(input, 'PNABc')),
	# 			 paste(
	# 			 	{
	# 			 		# Paste together the upstream portion of FIP with the downstream portion of the PNAF
	# 			 		BIPN <- length(s2c(BIP()))
	# 			 		B2cN <- input$B2cLen
	# 			 		Offset <- getEnd(input, 'PNABc')-getEnd(input, 'B2c')
	# 			 		if(getStart(input, 'PNABc') < getStart(input, 'B2c'))
	# 			 		{
	# 			 			stop('Having a PNAB longer than B2 is not supported yet.')
	# 			 		}
	# 			 		c(s2c(BIP())[(BIPN-B2cN-Offset):(BIPN-B2cN)], # Chunk from end of PNABc up to end of B2c in BIP
	# 			 		  revC(s2c(input$PNABcNTs), keepCase=T)[1:(input$PNAFLen-Offset)], keepCase=T) # Remaining (i.e., beginning) part of revC(PNABc)
	# 			 	}, 
	# 			 	collapse=''),
	# 			 revC(input$PNABcNTs, keepCase=T)
	# 	)
	# })
	
	# Now feed the ground truth out to inputs as necessary (if value in values in settings already matches, no additional triggers will setoff)
	# Feed NT values if Start input values change
	observeEvent(settings$Start, {
		req(all(sapply(sensePrimerNames, function(x){!is.null(settings$Start[[x]])})))
		if(silenceStartUpdate == 0)
		{
			silenceNTUpdate <<- silenceNTUpdate + 1
			ntUpdated <- FALSE
			for(controlId in sensePrimerNames)
			{
				bpStart <- settings$Start[[controlId]]
				bpEnd <- (settings$Start[[controlId]]+settings$Len[[controlId]]-1)
				if(length(bpStart) != 1)
				{
					browser()
				}
				bp <- seq(from=bpStart, to=bpEnd)
				temp <- paste(settings$seq[bp], collapse='')
				if(settings$NTs[[controlId]] != temp) # && silenceNTUpdate == 0)
				{
					# print(settings$render)
					ntUpdated <- TRUE
					print(paste("Updating NTs from Start:", silenceNTUpdate))
					updateTextInput(session, paste0(controlId, 'NTs'), value=temp)
				}
			}
			if(!ntUpdated){silenceNTUpdate <<- silenceNTUpdate - 1}
		}
		else if(silenceStartUpdate > 0)
		{
			silenceStartUpdate <<- silenceStartUpdate - 1
			print(paste("Silencing Start Update:", silenceStartUpdate))
		}
	})
	
	# Feed NT values if Len input values change
	observeEvent(settings$Len, {
		req(all(sapply(sensePrimerNames, function(x){!is.null(settings$Len[[x]])})))
		if(silenceLenUpdate == 0)
		{
			silenceNTUpdate <<- silenceNTUpdate + 1
			ntUpdated <- FALSE
			for(controlId in sensePrimerNames)
			{
				bpStart <- settings$Start[[controlId]]
				bpEnd <- (settings$Start[[controlId]]+settings$Len[[controlId]]-1)
				if(length(bpStart) != 1)
				{
					browser()
				}
				bp <- seq(from=bpStart, to=bpEnd)
				temp <- paste(settings$seq[bp], collapse='')
				if(settings$NTs[[controlId]] != temp) # && silenceNTUpdate == 0)
				{
					ntUpdated <- TRUE
					print(paste("Updating NTs from Len:", silenceNTUpdate))
					updateTextInput(session, paste0(controlId, 'NTs'), value=temp)
				}
			}
			if(!ntUpdated){silenceNTUpdate <<- silenceNTUpdate - 1}
		}
		else if(silenceLenUpdate > 0)
		{
			silenceLenUpdate <<- silenceLenUpdate - 1
			print(paste("Silencing Len Update:", silenceLenUpdate))
		}
	})
	
	# Feed Start and Len input values if NT values change
	observeEvent(settings$NTs, {
		req(all(sapply(sensePrimerNames, function(x){!is.null(settings$NTs[[x]])})))
		if(silenceNTUpdate > 0)
		{
			silenceNTUpdate <<- silenceNTUpdate - 1
			print(paste("Silencing NT Update:", silenceNTUpdate))
		}
		else
		{
			silenceStartUpdate <<- silenceStartUpdate + 1
			silenceLenUpdate <<- silenceLenUpdate + 1
			startUpdated <- FALSE
			lenUpdated <- FALSE
			for(controlId in sensePrimerNames)
			{
				temp <- settings$NTs[[controlId]]
				start <- getStartFromSeq(temp, paste(settings$seq, collapse=''))
				len <- calcLen(temp)
				if(settings$Start[[controlId]] != start) # && silenceStartUpdate == 0)
				{
					startUpdated <- TRUE
					print(paste("Update Start UI from NTs change:", silenceStartUpdate))
					updateNumericInput(session, paste0(controlId, 'Start'), value=start)
				}
				if(settings$Len[[controlId]] != len) # && silenceLenUpdate == 0)
				{
					lenUpdated <- TRUE
					print(paste("Update Len UI from NTs change", silenceLenUpdate))
					updateNumericInput(session, paste0(controlId, 'Len'), value=len)
				}
			}
			if(!startUpdated){
				# browser()
				silenceStartUpdate <<- silenceStartUpdate - 1
			}
			if(!lenUpdated){
				# browser()
				silenceLenUpdate <<- silenceLenUpdate - 1
			}
		}
	})
	
	allLegal <- reactive({
		req(all(sapply(sensePrimerNames, function(x){!is.null(settings$Start[[x]])})), all(sapply(sensePrimerNames, function(x){!is.null(settings$Len[[x]])})))
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
	})
	
	observeEvent(list(settings$NTs, settings$stabilityN, settings$Check), {
		# toInclude <- c(input$F3Check, input$B3Check, T, T, T, T, T, T, input$LFcCheck, input$LBCheck, input$PNACheck, input$PNAcCheck, T, T)
		primerNames <- c('F3','B3','F2','B2','F1c','B1c','FIP','BIP','LF','LB','PNAF','PNAB','DBF','DBB')
		req(settings$seq, allLegal(), all(as.logical(settings$NTs != '')))
		ret <- data.table(Primer=primerNames,
								Seq=as.character(sapply(lapply(primerNames, function(x){settings$NTs[[x]]}), paste, collapse='')),
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
		ret[Primer == 'DBF', Start:=settings$Start$F2 - settings$polyT] # F2r$Start-n
		ret[Primer == 'DBF', End:=settings$Start$F1 - 1] # F1r$Start-1
		ret[Primer == 'DBB', Start:=getEndSetting('B1c', settings) + 1] # B1cr$End+1
		ret[Primer == 'DBB', End:=getEndSetting('B2c', settings) + settings$polyT] # B2cr$End+n
		
		ret[, Len:=length(s2c(Seq)), by='Primer']
		
		# Old
		# ret[, c('Tm','HpTm','HpDeltaG','HpStruct'):=getPrimerStats(Seq), by='Primer']
		# New
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
		ret2[, Primer:=factor(Primer, levels=primerNames)]
		
		# Data for Tm plot
		ret3 <- ret[, c('Primer','Tm')]
		ret3[, Tm:=suppressWarnings(as.numeric(Tm))]
		primerColors <- getAllPrimerColors()
		setkey(ret3, Primer)
		setkey(primerColors, Primer)
		ret3 <- primerColors[ret3]
		ret3[, Primer:=factor(Primer, levels=primerNames)]
		settings$results2 <- ret2
		settings$results3 <- ret3
		settings$results <- ret
	})
	
	GCPlot <- reactive({
		req(settings$seq, settings$DBAll, settings$results2, calcLen(settings$DBAll) == (settings$DBEnd-settings$DBStart+1))
		plotGC(seq = settings$seq,
				 fullGC = isolate(fullGC()),
				 dbSeq = isolate(settings$DBAll),
				 dbStart = isolate(settings$DBStart),
				 dbEnd = isolate(settings$DBEnd),
				 starts = starts(), 
				 ends   = stops(),
				 sensePrimerNames = sensePrimerNames)
	})
	
	HairpinPlot <- reactive({
		req(settings$seq, settings$DBAll, calcLen(settings$DBAll) == (settings$DBEnd-settings$DBStart+1))
		plotHairpin(settings$seq,
						fullHairpinEnergy=isolate(fullHairpinEnergy()),
						dbSeq = isolate(settings$DBAll),
						dbStart = isolate(settings$DBStart),
						dbEnd = isolate(settings$DBEnd),
						starts = starts(), 
						ends   = stops(),
						sensePrimerNames = sensePrimerNames)
	})
	
	EnergiesPlot <- reactive({
		req(settings$results2) # 
		plotEnergies(settings$results2)
	})
	
	TmPlot <- reactive({
		req(settings$results3)
		plotTm(settings$results3)
	})
	
	observeEvent(list(settings$seq, settings$NTs), {
		
		output$GCPlot <- renderPlot({
			GCPlot()
		})
		
		output$HairpinPlot <- renderPlot({
			HairpinPlot()
		})
		
		output$EnergyPlot <- renderPlot({
			EnergiesPlot()
		})
		
		output$TmPlot <- renderPlot({
			TmPlot()
		})
	})
	
	observeEvent(input$revCTool,{
		req(input$revCTool != '')
		output$revCOutput <- renderText({
			revc <- revC(input$revCTool, keepCase=T)
			session$sendCustomMessage("txt", revc)
			paste('COPIED:', revc)
		})
	})
	
	output$ResultsTable <- renderTable({
		req(settings$results)
		settings$results
	})
	
	observeEvent(input$importButton, {
		x <- read.clipboard(os=Sys.info()['sysname'], header=F)
		x <- x[x[[1]] %in% c('F1','F1c','F2','F2c','F3','F3c','B1','B1c','B2','B2c','B3','B3c','PNAF','PNAFc','PNAB','PNABc','LF','LFc','LB','LBc')]
		if(nrow(x) == 0)
		{
			txt <- "First column should contain the name of the primer type (e.g., F3 etc). Re-copy primer information."
			warning(txt)
			output$revCOutput <- renderText(txt)
		}
		else if(ncol(x) != 3)
		{
			txt <- "Please copy just the primer name column and the start and end position columns (3 columns total). Re-copy primer information."
			warning(txt)
			output$revCOutput <- renderText(txt)
		}
		else
		{
			primerColName <- names(x)[1]
			startColName <- names(x)[2]
			endColName <- names(x)[3]
			x[, isRC:=grepl('[c]', get(primerColName))]
			x[, isSense:=get(primerColName) %in% c('F1','F1c','F2','F2c','F3','F3c','PNAF','PNAFc','LB','LBc')]
			# x[, seq:=paste(vars$seq[get(startColName):get(endColName)], collapse=''), by=c(primerColName)]
			# x[isRC & isSense, seq:=revC(seq, keepCase=T), by=c(primerColName)]
			x[isRC & isSense, c(primerColName):=gsub('[c]', '', get(primerColName))]
			# x[!isRC & !isSense, seq:=revC(seq, keepCase=T), by=c(primerColName)]
			x[!isRC & !isSense, c(primerColName):=paste(get(primerColName), 'c', sep='')]
			# browser()
			x[get(primerColName) %in% c('F1','F2','F3','B1c','B2c','B3c','PNAF','PNABc','LFc','LB')]
			for(i in 1:nrow(x))
			{
				startSetting <- x[i][[startColName]]
				lenSetting <- x[i][[endColName]]-x[i][[startColName]]+1
				NTSetting <- paste(settings$seq[seq(startSetting, startSetting + lenSetting -1)], collapse='')
				if(settings$Start[[x[i][[primerColName]]]] != startSetting) #updateSettingsItem(x[i][[primerColName]], startSetting, settings, group='Start'))
				{
					# Only mark down 1 update for an import event instead one for each change made.
					if(silenceStartUpdate == 0){silenceStartUpdate <<- silenceStartUpdate + 1}
					updateNumericInput(session, paste(x[i][[primerColName]], 'Start', sep=''), value=startSetting)
				}
				if(settings$Len[[x[i][[primerColName]]]] != lenSetting) # updateSettingsItem(x[i][[primerColName]], lenSetting, settings, group='Len'))
				{
					# Only mark down 1 update for an import event instead one for each change made.
					if(silenceLenUpdate == 0){silenceLenUpdate <<- silenceLenUpdate + 1}
					updateNumericInput(session, paste(x[i][[primerColName]], 'Len', sep=''), value=lenSetting)
				}
				if(settings$NTs[[x[i][[primerColName]]]] != NTSetting) # updateSettingsItem(x[i][[primerColName]], NTSetting, settings, group='NTs'))
				{
					# Only mark down 1 update for an import event instead one for each change made.
					if(silenceNTUpdate == 0){silenceNTUpdate <<- silenceNTUpdate + 1}
					updateTextInput(session, paste(x[i][[primerColName]], 'NTs', sep=''), value=NTSetting)
				}
			}
			output$revCOutput <- renderText("Clipboard imported.")
		}
	})
	
	observeEvent(input$importSettings, {
		# browser()
		settingsToSave <- loadRData(input$importSettings$datapath) # previously saved as 'settingsToSave' during download
		for(item in names(settingsToSave))
		{
			print(paste("Importing:", item))
			settings[[item]] <- settingsToSave[[item]]
		}
		updateTextAreaInput(session, 'Seq', value=paste(settings$seq, collapse=''))
		for(name in sensePrimerNames)
		{
			updateNumericInput(session, inputId=paste(name, 'Start', sep=''), value=settings$Start[[name]])
			updateNumericInput(session, inputId=paste(name, 'Len', sep=''), value=settings$Len[[name]])
		}
	})
	
	output$download <- downloadHandler(
		# filename = function(){paste(input$downloadName, '.csv')},
		# content = function(fname){
		# 	fwrite(settings$results, fname)
		# }

		filename = function(){paste(input$downloadName, '.zip')},
		content = function(fname){

			# Set temporary working directory
			owd <- setwd( tempdir())
			print(getwd())
			on.exit( setwd( owd))
			
			temp <- list(GCPlot=GCPlot(), HairpinPlot=HairpinPlot(), EnergiesPlot=EnergiesPlot(), TmPlot=TmPlot())
			temp <- temp[sapply(temp, isTruthy)]
			lapply(names(temp), function(x){
				ggsave(paste(x, '.png', sep=''), plot = temp[[x]], device = "png", width = 18, height = 6, dpi = 150, units = "in")
			})

			# print(settings$results)
			fwrite(settings$results, 'SummaryTable.csv')
			
			# Save the inputs so they can be reloaded later if desired.
			settingsToSave <- list()
			for(item in names(settings))
			{
				settingsToSave[[item]] <- settings[[item]]
			}
			save(settingsToSave, file='Settings.RData')
			
			# Zip them up
			# print(fname)
			return(zip( zipfile=fname, files=c('Settings.RData', 'SummaryTable.csv', paste(names(temp), '.png', sep=''))))
		},
		contentType = "application/zip"
	)
}

# Run the application 
shinyApp(ui = ui, server = server)

