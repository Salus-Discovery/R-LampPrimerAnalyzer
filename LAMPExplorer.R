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

sink(stdout(), type="message")

getMultiVals = function(field, vals){
	isolate(purrr::pluck(vals, !!!field))
}

setMultiVals = function(rv,field,value){
	isolate(purrr::pluck(rv, !!!field) <- value )
}

fold <- function(x)
{
	if(calcLen(x) > 47)
	{
		if(any(grepl('[Nn]',x)))
		{
			browser()
		}
		duh <- predictMfeStructure(XNAString(base = toUpper(paste(x, collapse=''))))
		ret <- list(structure_found = T, 
						temp = 0, 
						dg = duh$mfe,
						structure = duh$structure)
	}
	else
	{
		duh <- toUpper(paste(x, collapse=''))
		ret <- calculate_hairpin(duh)
		ret$dg <- ret$dg/1000
	}
	# obj1 <- XNAString(base = duh)
	# ret <- predictMfeStructure(obj1)
	return(ret)
}

homodimer <- function(x)
{
	duh <- toUpper(paste(x, collapse=''))
	if(is.null(duh))
	{
		browser()
	}
	ret <- calculate_homodimer(duh)
	ret$dg <- ret$dg/1000
	return(ret)
}

heterodimer <- function(x, y)
{
	ret <- calculate_dimer(toUpper(paste(x, collapse='')), toUpper(paste(y, collapse='')))
	ret$dg <- ret$dg/1000
	return(ret)
}

end_stability <- function(x, n=5, threePrimeEnd=TRUE)
{
	if(length(x) == 1)
	{
		x <- s2c(x)
	}
	if(threePrimeEnd)
	{
		x <- c(rep('N',20), x[(length(x)-n):length(x)])
		# x[1:(length(x)-n)] <- 'N'
	}
	else
	{
		x <- c(x[1:n], rep('N',20))
		# x[(n+1):length(x)] <- 'N'
	}
	# print(x)
	# print(revC(x))
	return(sig.digits(heterodimer(x, revC(x))$dg, nSig=2))
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

plotHairpin <- function(seq, fullHairpinEnergy, dbSeq, dbStart, dbEnd, starts, ends, sensePrimerNames)
{
	bp <- dbStart:dbEnd
	seq.new <- s2c(dbSeq)
	# bpRange <- c(max(1,start.new), min(length(seq), start.new + length(seq.new) - 1))
	# bp <- bp[bp > 0 & bp < calcLen(seq)]
	hairpin.new <- rollapply(seq.new, width=25, FUN=function(x){return(fold(x)$dg)}, partial=T, align='center')
	hairpin.new[is.na(hairpin.new) | hairpin.new > 0] <- 0
	# validBp <- bp[bp > 0 & bp < length(seq)]
	hairpin <- fullHairpinEnergy # rollapply(seq, width=25, FUN=function(x){return(fold(x)$dg/1000)}, partial=T, align='center')
	# hairpin <- rollapply(seq[validBp], width=25, FUN=function(x){return(fold(x)$dg/1000)}, partial=T, align='center')
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
		scale_fill_manual(values = setColor(getColor(sensePrimerNames), 0.3)) +
		geom_path(data=data.table(x=1:calcLen(seq), y=hairpin), aes(x=x, y=y)) +
		geom_path(data=data.table(x=bp, y=hairpin.new), aes(x=x, y=y), linewidth=3) +
		geom_hline(yintercept=-4, col=setColor('red', 0.2)) +
		geom_vline(xintercept=codons, col=setColor('red', 0.2)) +
		scale_x_continuous(limits=c(1,calcLen(seq))) +
		scale_y_continuous(limits=c(-15,0.1)) + 
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

plotTm <- function(results)
{
	ret <- ggplot(data=results, aes(x=Primer, y=Tm, fill=Primer)) + 
		geom_col() +
		scale_fill_manual(values = getColor(as.character(levels(results$Primer)))) + 
		geom_hline(yintercept=c(50,60,70)) +
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
		scale_fill_manual(values = setColor(getColor(sensePrimerNames), 0.3)) +
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

shade <- function(start, end, color, alpha=1)
{
	ymin <- par('usr')[3]
	ymax <- par('usr')[4]
	polygon(c(start, end, end, start), c(ymin, ymin, ymax, ymax), col=setColor(color, alpha))
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

getPrimerStats <- function(x)
{
	hp <- fold(x)
	# homd <- homodimer(x)
	Tm <- calculate_tm(toUpper(x))
	ret <- list(Tm=sig.digits(Tm, 3),
					HpTm=sig.digits(hp$temp, 3),
					HpDeltaG=sig.digits(hp$dg, 3),
					HpStruct=ifelse(hp$structure_found, hp$structure, ''))
	# HmdTm=sig.digits(homd$temp, 3),
	# HmdDeltaG=sig.digits(homd$dg, 3),
	# HmdStruct=ifelse(homd$structure_found, homd$structure, ''))
	return(ret)
}

getColor <- function(...)
{
	args <- list(...)
	if(length(args) == 1 && length(args[[1]])>1)
	{
		return(getPrimerColors()[match(as.character(args[[1]]), Primer)]$plotColor)
	}
	else
	{
		return(getPrimerColors()[match(args, Primer)]$plotColor)
	}
	
}

getPrimerColors <- function()
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

getHtdStats <- function(primers, seqs, func=daFunc)
{
	duh2 <- data.table(getCombos(primers, colnames=c('P1','P2')), getCombos(seqs, c('Seq1','Seq2')))
	duh2 <- duh2[duh2[, func(Seq1, Seq2), by=c('P1','P2')], on=c('P1','P2')]
	valNames <- names(duh2)
	valNames <- valNames[valNames %!in% c('P1','P2','Seq1','Seq2')]
	blah <- rbindlist(list(duh2, data.table(Seq1=duh2$Seq2, Seq2=duh2$Seq1, P1=duh2$P2, P2=duh2$P1, duh2[, mget(valNames)])), use.names=T)
	blah <- unique(blah)
	blah[, max.i:=which.min(HtdDeltaG)[1], by='P1']
	blah <- blah[, c(list(P2=P2[max.i]), lapply(mget(valNames), function(x){x[max.i[1]]})), by='P1']
	return(blah[match(primers, P1), mget(c('P2', valNames))])
}


getDimerStats <- function(x, y)
{
	hetd <- heterodimer(x, y)
	ret <- list(HtdTm=sig.digits(hetd$temp, 3),
					HtdDeltaG=sig.digits(hetd$dg, 3),
					HtdStruct=ifelse(hetd$structure_found, hetd$structure, ''))
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
	if(os[1]=='mac')
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

# Define UI for application that draws a histogram
ui <- fluidPage(
	theme = shinytheme("darkly"),
	tags$head(
		# Note the wrapping of the string in HTML()
		tags$style(HTML(paste("
      markF3	{background-color: ", getColor('F3'), "; 	color: black;}
      markF2	{background-color: ", getColor('F2'), ";	color: black;}
      markLFc	{background-color: ", getColor('LFc'), ";	color: black;}
      markF1	{background-color: ", getColor('F1'), ";	color: black;}
      markB1c	{background-color: ", getColor('B1c'), ";	color: black;}
      markLB	{background-color: ", getColor('LB'), ";	color: black;}
      markB2c	{background-color: ", getColor('B2c'), ";	color: black;}
      markB3c	{background-color: ", getColor('B3c'), ";	color: black;}
      markPNAF {background-color: ", getColor('PNAF'), ";	color: black;}
      markPNABc{background-color: ", getColor('PNABc'), ";	color: black;}
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
						 fluidPage(fluidRow(column(6, numericInput("polyT", "FIP/BIP PolyT-Spacer Length:", 3)),
						 						 column(6, numericInput("stabilityN", "End Stability Bp's:", 5))),
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

getSeqFromStartAndLen <- function(start, len, mySeq)
{
	return(mySeq[start:(start+len-1)])
}

getStartFromSeq <- function(primer, target)
{
	temp <- pwalign::pairwiseAlignment(primer, target, type='local', gapOpening=1000000, gapExtension=1000000)
	return(start(subject(temp)))
}

renderPrimerControl <- function(seq, controlId, initLocs=c(1), Loc=1, initLen=20)
{
	initStart <- ifelse(Loc > 0, initLocs[Loc], 1)
	fluidPage({
		fluidRow(style = "height:20px;",
					column(1, HTML(paste(tags$h4(tags$span(style=paste("color: ", getColor(controlId), sep=''), controlId), sep='')))),
					column(1, checkboxInput(paste0(controlId, 'Check'), '', value=T)), #value=grepl('[12]', controlId))),
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

updateValsGroupItem <- function(id, group, input, vals)
{
	# print(paste(c(id, group, input[[paste(id, group, sep='')]], vals[[group]])))
	val <- input[[paste(id, group, sep='')]]
	if(length(val) > 0)
	{
		if(length(vals[[group]]) == 0)
		{
			vals[[group]] <- list()
		}
		if(length(vals[[group]][[id]])==0 || any(vals[[group]][[id]] != val))
		{
			vals[[group]][[id]] <- val
		}
	}
}

updateValsItem <- function(id, val, vals, group=NULL)
{
	if(is.null(group))
	{
		if(length(vals)>0 && length(val)>0)
		{
			if(length(vals[[id]])==0 || length(vals[[id]]) != length(val) || any(vals[[id]] != val))
			{
				vals[[id]] <- val
				return(TRUE)
			}
		}
	}
	else
	{
		if(length(vals)>0 && length(val)>0)
		{
			if(length(vals[[group]])==0)
			{
				vals[[group]] <- list()
			}
			if(length(vals[[group]][[id]])==0 || any(vals[[group]][[id]] != val))
			{
				vals[[group]][[id]] <- val
				return(TRUE)
			}
		}
	}
	return(FALSE)
}

getValStart <- function(id, vals)
{
	return(vals[[Start]][[id]])
}

getValEnd <- function(id, vals)
{
	return(vals[['Start']][[id]] + vals[['Len']][[id]] - 1)
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
	
	vals <- reactiveValues(render=1)
	
	mySeqHTML <- reactive({
		req(vals$seq)
		markup=p(id='text-to-mark', style = "word-wrap: break-word;", paste(vals$seq, collapse=''))
	})
	
	fullHairpinEnergy <- reactive({
		req(vals$seq, allLegal())
		print("Updating hairpin energy plot.")
		rollapply(vals$seq, width=25, FUN=function(x){return(fold(x)$dg)}, partial=T, align='center')
	})
	
	fullGC <- reactive({
		req(vals$seq)
		print("Updating fullGC.")
		gc <- vals$seq %in% c("g","G","c","C")
		rollapply(gc, width=20, FUN=mean, partial=T, align='center')
	})
	
	initLocs <- reactive({
		req(vals$seq)
		print("Updating initLocs.")
		return(round(seq(1, length(vals$seq)-20, length.out=8)))
	})
	
	output$coloredSeq <- renderUI({
		req(mySeqHTML())
		print("Rendering HTML Seq.")
		mySeqHTML()
	})
	
	observeEvent( list(vals$seq, vals$render > 0),{
		req(vals$seq)
		if(calcLen(vals$seq)==0){vals$render <- vals$render + 1}
		print("Rendering primer controls.")
		output$F3 <- renderUI(isolate(renderPrimerControl(vals$seq, 'F3', initLocs=initLocs(), Loc=2)))
		output$F2 <- renderUI(isolate(renderPrimerControl(vals$seq, 'F2', initLocs=initLocs(), Loc=3)))
		output$LFc <- renderUI(isolate(renderPrimerControl(vals$seq, 'LFc', initLocs=initLocs(), Loc=0)))
		output$F1 <- renderUI(isolate(renderPrimerControl(vals$seq, 'F1', initLocs=initLocs(), Loc=4)))
		output$B1c <- renderUI(isolate(renderPrimerControl(vals$seq, 'B1c', initLocs=initLocs(), Loc=5)))
		output$LB <- renderUI(isolate(renderPrimerControl(vals$seq, 'LB', initLocs=initLocs(), Loc=0)))
		output$B2c <- renderUI(isolate(renderPrimerControl(vals$seq, 'B2c', initLocs=initLocs(), Loc=6)))
		output$B3c <- renderUI(isolate(renderPrimerControl(vals$seq, 'B3c', initLocs=initLocs(), Loc=7)))
		output$PNAF <- renderUI(isolate(renderPrimerControl(vals$seq, 'PNAF', initLocs=initLocs(), Loc=0)))
		output$PNABc <- renderUI(isolate(renderPrimerControl(vals$seq, 'PNABc', initLocs=initLocs(), Loc=0)))
		# print(vals$render)
		if(all(sapply(sensePrimerNames, function(x){calcLen(input[[paste(x, 'NTs', sep='')]]) > 0}))){
			vals$render <- 0
			silenceLenUpdate <<- 0
			silenceStartUpdate <<- 0
			silenceNTUpdate <<- 0
		}
	})
	
	observeEvent(list(mySeqHTML(), getInputs(sensePrimerNames, 'NTs', input), getInputs(sensePrimerNames, 'Check', input)), {
		# Highlight output
		# print(paste("Marking ", controlId, ": ", input[[paste0(controlId, "NTs")]]))
		req(mySeqHTML(), all(sapply(sensePrimerNames, function(x){!is.null(vals$NTs[[x]])})), all(sapply(sensePrimerNames, function(x){!is.null(vals$Check[[x]])})))
		print("Updating markers")
		my_marker <- marker$new("#text-to-mark")
		my_marker$unmark()
		# browser()
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
	
	# setUpControlLinks('F3', mySeq, initLocs)
	# setUpControlLinks('F2', mySeq, initLocs)
	# setUpControlLinks('LFc', mySeq, initLocs)
	# setUpControlLinks('F1', mySeq, initLocs)
	# setUpControlLinks('B1c', mySeq, initLocs)
	# setUpControlLinks('LB', mySeq, initLocs)
	# setUpControlLinks('B2c', mySeq, initLocs)
	# setUpControlLinks('B3c', mySeq, initLocs)
	# setUpControlLinks('PNAF', mySeq, initLocs)
	# setUpControlLinks('PNABc', mySeq, initLocs)
	
	sensePrimerNames <- c('F3','F2','F1','LFc','B1c','B2c','B3c','LB','PNAF','PNABc')
	
	dependentPrimerNames <- c('F1c', 'B2', 'B3', 'LF', 'LB', 'PNAB', 'FIP', 'BIP', 'DBF', 'DBB')
	
	starts <- reactive({
		as.numeric(sapply(sensePrimerNames, function(x){vals$Start[[x]]}))
	})
	
	stops <- reactive({
		as.numeric(sapply(sensePrimerNames, function(x){vals$Start[[x]]}))+as.numeric(lapply(sensePrimerNames, function(x){vals$Len[[x]]}))-1
	})
	
	primerColors <- reactive({
		getColor(sensePrimerNames)
	})
	
	observeEvent(list(input$Seq, vals$render > 0), {
		updateValsItem('seq', {
			temp <- s2c(input$Seq)
			temp[temp %in% c('a','g','t','c','A','G','T','C')]
			if(input$Seq != paste(temp, collapse=''))
			{
				warning('Non A, G, T, C, a, g, t, c, characters found. Repairing sequence')
				if(vals$seq != temp)
				{
					vals$seq <- temp
				}
				if(input$Seq == '')
				{
					updateTextAreaInput(session, inputId='Seq', value=paste(temp, collapse=''))
				}
			}
			# print(temp)
			temp
		}, vals)
	})
	
	observeEvent(getInputs(sensePrimerNames, 'Start', input), {
		do.call('req', getInputs(sensePrimerNames, 'Start', input))
		lapply(sensePrimerNames, updateValsGroupItem, group='Start', input=input, vals)
	})
	
	observeEvent(getInputs(sensePrimerNames, 'Len', input), {
		do.call('req', getInputs(sensePrimerNames, 'Len', input))
		lapply(sensePrimerNames, updateValsGroupItem, group='Len', input=input, vals)
	})
	
	# Watch all the inputs and create a ground truth stored version of everything
	observeEvent(list(vals$seq, vals$polyT, getInputs(sensePrimerNames, 'NTs', input)), {
		tempPNABc <- input$PNABc
		do.call('req', getInputs(sensePrimerNames, 'NTs', input))
		req(vals$seq, vals$polyT)
		lapply(sensePrimerNames, updateValsGroupItem, group='NTs', input=input, vals=vals)
		updateValsItem('F1c', revC(vals$NTs$F1, keepCase=T), vals, group='NTs')
		updateValsItem('B2', revC(vals$NTs$B2c, keepCase=T), vals, group='NTs')
		updateValsItem('B3', revC(vals$NTs$B3c, keepCase=T), vals, group='NTs')
		updateValsItem('LF', revC(vals$NTs$LFc, keepCase=T), vals, group='NTs')
		updateValsItem('LB', vals$NTs$LB, vals, group='NTs')
		updateValsItem('FIP', {
			paste(c(vals$NTs$F1c, rep('t', vals$polyT), vals$NTs$F2), collapse='')	
		}, vals, group='NTs')
		updateValsItem('BIP', {
			paste(c(vals$NTs$B1c, rep('t', vals$polyT), vals$NTs$B2), collapse='')	
		}, vals, group='NTs')
		updateValsItem('DBF', paste(c(rep('t', vals$polyT), vals$seq[vals$Start$F2:(vals$Start$F1-1)]), collapse=''), vals, group='NTs')
		updateValsItem('DBB', paste(c(rep('t', vals$polyT), vals$seq[getValEnd('B2c', vals):(getValEnd('B1c', vals)+1)]), collapse=''), vals, group='NTs')
		updateValsItem('PNAF', {
			# If the PNA straddles the start of F2
			ifelse(vals$Start$PNAF < vals$Start$F2 && getValEnd('PNAF', vals) >= vals$Start$F2, 
					 paste(
					 	{
					 		# Paste together the upstream portion of FIP with the downstream portion of the PNAF
					 		FIPN <- length(vals$NTs$FIP)
					 		F2N <- vals$Len$F2
					 		Offset <- vals$Start$F2-vals$Start$PNAF
					 		if(getValEnd('PNAF', vals) > getValEnd('F2', vals))
					 		{
					 			warning('Having a PNAF longer than F2 is not supported yet.')
					 			updateTextInput(session, 'PNAFNTs', vals$NTs$PNAF)
					 			s2c(vals$NTs$PNAF) # Don't change anything
					 		}
					 		else
					 		{
					 			c(vals$NTs$FIP[(FIPN-F2N-Offset):(FIPN-F2N)], # Chunk from start of PNAF up to start of F2 in FIP
					 			  vals$NTsPNAF[(Offset+1):vals$Len$PNAF]) # Remaining part of PNAF
					 		}
					 	}, 
					 	collapse=''),
					 vals$NTs$PNAF
			)
		}, vals, group='NTs')
		updateValsItem('PNAB', {
			# If the PNA straddles the end of B2c
			ifelse(vals$Start$PNABc < getValEnd('B2c', vals) && getValEnd('PNABc', vals) >= getValEnd('B2c', vals), 
					 paste(
					 	{
					 		# Paste together the upstream portion of FIP with the downstream portion of the PNAF
					 		BIPN <- length(vals$NTs$BIP)
					 		B2cN <- vals$Len$B2c
					 		Offset <- getValEnd('PNABc', vals)-getValEnd('B2c', vals)
					 		if(vals$Start$PNABc < vals$Start$B2c)
					 		{
					 			stop('Having a PNAB longer than B2 is not supported yet.')
					 			# Reset the source PNABc value to what it was before this observed event
					 			vals$NTs$PNABc <- tempPNABc
					 			# Reset the input to the same value
					 			updateTextInput(session, 'PNABcNTs', tempPNABc)
					 			s2c(vals$NTs$PNAB) # Don't change anything by returning the same value
					 		}
					 		else
					 		{
					 			c(vals$NTs$BIP[(BIPN-B2cN-Offset):(BIPN-B2cN)], # Chunk from end of PNABc up to end of B2c in BIP
					 			  revC(vals$NTs$PNABc, keepCase=T)[1:(vals$Len$PNABc-Offset)], keepCase=T) # Remaining (i.e., beginning) part of revC(PNABc)
					 		}
					 	}, 
					 	collapse=''),
					 revC(vals$NTs$PNABc, keepCase=T))
		}, vals, group='NTs')
		updateValsItem('DBAll', {
			paste(c(vals$NTs$F1c, rep('t', input$polyT), vals$seq[vals$Start$F2:getValEnd('B2c', vals)], rep('t', input$polyT), revC(vals$NTs$B1c, keepCase=T)), collapse='')
		}, vals)
		updateValsItem('DBStart', vals$Start$F2-vals$polyT-vals$Len$F1, vals)
		updateValsItem('DBEnd', (vals$Start$B2c + vals$Len$B2c - 1)+vals$polyT+vals$Len$B1c, vals)
	})
	
	observeEvent(getInputs(sensePrimerNames, 'Check', input), {
		do.call('req', getInputs(sensePrimerNames, 'Check', input))
		lapply(sensePrimerNames, updateValsGroupItem, group='Check', input=input, vals)
	})
	
	observeEvent(input$polyT, {
		updateValsItem('polyT', input$polyT, vals)
	})
	
	observeEvent(list(vals$Check, vals$NTs), {
		# Highlight output
		# print(paste("Marking ", controlId, ": ", input[[paste0(controlId, "NTs")]]))
		
		req(!is.null(vals$seq), !is.null(vals$Check), length(vals$Check) > 0)
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
		
		updateValsItem('stabilityN', input$stabilityN, vals)	
	})
	
	# # Dumbbell regions
	# DBF <- reactive(paste(c(rep('t', vals$polyT), mySeq()[getStart(input, 'F2'):(getStart(input, 'F1')-1)]), collapse=''))
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
	
	# Now feed the ground truth out to inputs as necessary (if value in values in vals already matches, no additional triggers will setoff)
	# Feed NT values if Start input values change
	observeEvent(vals$Start, {
		req(all(sapply(sensePrimerNames, function(x){!is.null(vals$Start[[x]])})))
		if(silenceStartUpdate == 0)
		{
			silenceNTUpdate <<- silenceNTUpdate + 1
			ntUpdated <- FALSE
			for(controlId in sensePrimerNames)
			{
				bpStart <- vals$Start[[controlId]]
				bpEnd <- (vals$Start[[controlId]]+vals$Len[[controlId]]-1)
				if(length(bpStart) != 1)
				{
					browser()
				}
				bp <- seq(from=bpStart, to=bpEnd)
				temp <- paste(vals$seq[bp], collapse='')
				if(vals$NTs[[controlId]] != temp) # && silenceNTUpdate == 0)
				{
					# print(vals$render)
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
	observeEvent(vals$Len, {
		req(all(sapply(sensePrimerNames, function(x){!is.null(vals$Len[[x]])})))
		if(silenceLenUpdate == 0)
		{
			silenceNTUpdate <<- silenceNTUpdate + 1
			ntUpdated <- FALSE
			for(controlId in sensePrimerNames)
			{
				bpStart <- vals$Start[[controlId]]
				bpEnd <- (vals$Start[[controlId]]+vals$Len[[controlId]]-1)
				if(length(bpStart) != 1)
				{
					browser()
				}
				bp <- seq(from=bpStart, to=bpEnd)
				temp <- paste(vals$seq[bp], collapse='')
				if(vals$NTs[[controlId]] != temp) # && silenceNTUpdate == 0)
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
	observeEvent(vals$NTs, {
		req(all(sapply(sensePrimerNames, function(x){!is.null(vals$NTs[[x]])})))
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
				temp <- vals$NTs[[controlId]]
				start <- getStartFromSeq(temp, paste(vals$seq, collapse=''))
				len <- calcLen(temp)
				if(vals$Start[[controlId]] != start) # && silenceStartUpdate == 0)
				{
					startUpdated <- TRUE
					print(paste("Update Start UI from NTs change:", silenceStartUpdate))
					updateNumericInput(session, paste0(controlId, 'Start'), value=start)
				}
				if(vals$Len[[controlId]] != len) # && silenceLenUpdate == 0)
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
		req(all(sapply(sensePrimerNames, function(x){!is.null(vals$Start[[x]])})), all(sapply(sensePrimerNames, function(x){!is.null(vals$Len[[x]])})))
		getValEnd('F3', vals) <= length(vals$seq) &&
			getValEnd('F2', vals) <= length(vals$seq) &&
			getValEnd('F1', vals) <= length(vals$seq) &&
			getValEnd('LFc', vals) <= length(vals$seq) &&
			getValEnd('B3c', vals) <= length(vals$seq) &&
			getValEnd('B2c', vals) <= length(vals$seq) &&
			getValEnd('B1c', vals) <= length(vals$seq) &&
			getValEnd('LB', vals) <= length(vals$seq) &&
			getValEnd('PNAF', vals) <= length(vals$seq) &&
			getValEnd('PNABc', vals) <= length(vals$seq)
	})
	
	observeEvent(list(vals$NTs, vals$stabilityN, vals$Check), {
		# toInclude <- c(input$F3Check, input$B3Check, T, T, T, T, T, T, input$LFcCheck, input$LBCheck, input$PNACheck, input$PNAcCheck, T, T)
		primerNames <- c('F3','B3','F2','F1c','B2','B1c','FIP','BIP','LF','LB','PNAF','PNAB','DBF','DBB')
		req(vals$seq, allLegal(), all(as.logical(vals$NTs != '')))
		ret <- data.table(Primer=primerNames,
								Seq=as.character(sapply(lapply(primerNames, function(x){vals$NTs[[x]]}), paste, collapse='')),
								Sense=c('Sense','Antisense','Sense','Antisense','Sense','Antisense','NA','NA','Antisense','Sense','Sense','Antisense','Sense','Sense'))
		ret[, Stability3p:=end_stability(Seq, n=vals$stabilityN, threePrimeEnd=T), by='Primer']
		ret[, Stability5p:=end_stability(Seq, n=vals$stabilityN, threePrimeEnd=F), by='Primer']
		ret[, Start5p:='NA']
		ret[Primer == 'F3', Start5p:=vals$Start$F3]
		ret[Primer == 'B3', Start5p:=getValEnd('B3c', vals)]
		ret[Primer == 'F2', Start5p:=vals$Start$F2]
		ret[Primer == 'F1c', Start5p:=getValEnd('F1', vals)]
		ret[Primer == 'B2', Start5p:=getValEnd('B2c', vals)]
		ret[Primer == 'B1c', Start5p:=vals$Start$B1c]
		ret[Primer == 'LF', Start5p:=getValEnd('LFc', vals)]
		ret[Primer == 'LB', Start5p:=vals$Start$LB]
		ret[Primer == 'PNAF', Start5p:=vals$Start$PNAF]
		ret[Primer == 'PNAB', Start5p:=getValEnd('PNABc', vals)]
		ret[, Len:=length(s2c(Seq)), by='Primer']
		ret[, c('Tm','HpTm','HpDeltaG','HpStruct'):=getPrimerStats(Seq), by='Primer']
		ret[Primer %!in% c('DBF','DBB'), c('P2','DimerTm','DimerDeltaG','DimerStruct'):=getHtdStats(Primer, Seq, func=getDimerStats)]
		# ret[Primer %in% c('FIP','BIP'), Tm:='NA']
		ret[, KeyEndStability:=ifelse(Primer %in% c('F1c','B1c'), Stability5p, Stability3p)]
		
		# Data for energy plot
		ret2 <- melt.data.table(data=ret, id.vars='Primer', measure.vars=c('KeyEndStability','HpDeltaG','DimerDeltaG'))
		ret2[, value:=suppressWarnings(as.numeric(value))]
		ret2[value > 0, value:=0]
		ret2[variable=='KeyEndStability', value:=-1*value]
		ret2[, Primer:=factor(ret2$Primer, levels=c('F3','B3','F2','B2','LF','LB','F1c','B1c','FIP','BIP','PNAF','PNAB','DBF','DBB'))]
		
		# Data for Tm plot
		ret3 <- ret[, c('Primer','Tm')]
		ret3[, Tm:=suppressWarnings(as.numeric(Tm))]
		primerColors <- getPrimerColors()
		setkey(ret3, Primer)
		setkey(primerColors, Primer)
		ret3 <- primerColors[ret3]
		ret3[, Primer:=factor(ret3$Primer, levels=c('F3','B3','F2','B2','LF','LB','F1c','B1c','PNAF','PNAB','FIP','BIP','DBF','DBB'))]
		vals$results2 <- ret2
		vals$results3 <- ret3
		vals$results <- ret[Primer %in% primerNames[
			c(vals$Check$F3,
			  vals$Check$B3c,
			  vals$Check$F2,
			  vals$Check$F1,
			  vals$Check$B2c,
			  vals$Check$B1c,
			  vals$Check$LFc,
			  vals$Check$LB,
			  vals$Check$PNAF,
			  vals$Check$PNABc,
			  vals$Check$F1 && vals$Check$F2,
			  vals$Check$B1c && vals$Check$B2c,
			  T,
			  T)]]
	})
	
	GCPlot <- reactive({
		req(vals$seq, vals$DBAll, vals$results2, calcLen(vals$DBAll) == (vals$DBEnd-vals$DBStart+1))
		plotGC(seq = vals$seq,
				 fullGC = isolate(fullGC()),
				 dbSeq = isolate(vals$DBAll),
				 dbStart = isolate(vals$DBStart),
				 dbEnd = isolate(vals$DBEnd),
				 starts = starts(), 
				 ends   = stops(),
				 sensePrimerNames = sensePrimerNames)
	})
	
	HairpinPlot <- reactive({
		req(vals$seq, vals$DBAll, calcLen(vals$DBAll) == (vals$DBEnd-vals$DBStart+1))
		plotHairpin(vals$seq,
						fullHairpinEnergy=isolate(fullHairpinEnergy()),
						dbSeq = isolate(vals$DBAll),
						dbStart = isolate(vals$DBStart),
						dbEnd = isolate(vals$DBEnd),
						starts = starts(), 
						ends   = stops(),
						sensePrimerNames = sensePrimerNames)
	})
	
	EnergiesPlot <- reactive({
		req(vals$results2) # 
		plotEnergies(vals$results2)
	})
	
	TmPlot <- reactive({
		req(vals$results3)
		plotTm(vals$results3)
	})
	
	observeEvent(list(vals$seq, vals$NTs), {
		
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
		req(vals$results)
		vals$results
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
				startVal <- x[i][[startColName]]
				lenVal <- x[i][[endColName]]-x[i][[startColName]]+1
				NTVal <- paste(vals$seq[seq(startVal, startVal + lenVal -1)], collapse='')
				if(vals$Start[[x[i][[primerColName]]]] != startVal) #updateValsItem(x[i][[primerColName]], startVal, vals, group='Start'))
				{
					# Only mark down 1 update for an import event instead one for each change made.
					if(silenceStartUpdate == 0){silenceStartUpdate <<- silenceStartUpdate + 1}
					updateNumericInput(session, paste(x[i][[primerColName]], 'Start', sep=''), value=startVal)
				}
				if(vals$Len[[x[i][[primerColName]]]] != lenVal) # updateValsItem(x[i][[primerColName]], lenVal, vals, group='Len'))
				{
					# Only mark down 1 update for an import event instead one for each change made.
					if(silenceLenUpdate == 0){silenceLenUpdate <<- silenceLenUpdate + 1}
					updateNumericInput(session, paste(x[i][[primerColName]], 'Len', sep=''), value=lenVal)
				}
				if(vals$NTs[[x[i][[primerColName]]]] != NTVal) # updateValsItem(x[i][[primerColName]], NTVal, vals, group='NTs'))
				{
					# Only mark down 1 update for an import event instead one for each change made.
					if(silenceNTUpdate == 0){silenceNTUpdate <<- silenceNTUpdate + 1}
					updateTextInput(session, paste(x[i][[primerColName]], 'NTs', sep=''), value=NTVal)
				}
			}
			output$revCOutput <- renderText("Clipboard imported.")
		}
	})
	
	observeEvent(input$importSettings, {
		# browser()
		load(input$importSettings$datapath)
		for(item in names(settings))
		{
			print(paste("Importing:", item))
			vals[[item]] <- settings[[item]]
		}
		updateTextAreaInput(session, 'Seq', value=paste(vals$seq, collapse=''))
		for(name in sensePrimerNames)
		{
			updateNumericInput(session, inputId=paste(name, 'Start', sep=''), value=vals$Start[[name]])
			updateNumericInput(session, inputId=paste(name, 'Len', sep=''), value=vals$Len[[name]])
		}
	})
	
	output$download <- downloadHandler(
		# filename = function(){paste(input$downloadName, '.csv')},
		# content = function(fname){
		# 	fwrite(vals$results, fname)
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

			# print(vals$results)
			fwrite(vals$results, 'SummaryTable.csv')
			
			# Save the inputs so they can be reloaded later if desired.
			settings <- list()
			for(item in names(vals))
			{
				settings[[item]] <- vals[[item]]
			}
			save(settings, file='Settings.RData')
			
			# Zip them up
			# print(fname)
			return(zip( zipfile=fname, files=c('Settings.RData', 'SummaryTable.csv', paste(names(temp), '.png', sep=''))))
		},
		contentType = "application/zip"
	)
}

# Run the application 
shinyApp(ui = ui, server = server)

