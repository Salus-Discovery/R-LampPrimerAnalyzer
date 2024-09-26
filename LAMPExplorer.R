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
	return(ret)
}

heterodimer <- function(x, y)
{
	ret <- calculate_dimer(toUpper(paste(x, collapse='')), toUpper(paste(y, collapse='')))
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
	return(sig.digits(heterodimer(x, revC(x))$dg/1000, nSig=2))
}

revC <- function(x, keepCase=F)
{
	if(any(is.na(x)))
	{
		browser()
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

plotHairpin <- function(seq, dbSeq, dbStart, dbEnd, starts, ends, colors)
{
	start.new <- dbStart
	seq.new <- s2c(dbSeq)
	bp <- start.new:(start.new + length(seq.new) - 1)
	bp <- bp[bp > 0 & bp < calcLen(seq)]
	hairpin.new <- rollapply(seq.new[bp], width=25, FUN=function(x){return(fold(x)$dg/1000)}, partial=T, align='center')
	hairpin <- rollapply(seq[bp], width=25, FUN=function(x){return(fold(x)$dg/1000)}, partial=T, align='center')
	if(any(seq[bp] %in% c('A','G','T','C')))
	{
	  codons <- bp[which(seq[bp] %in% c('A','G','T','C'))]
	}
	else
	{
	  codons <- c()
	}
	
	if(length(codons) > 0)
	{
		plot(bp, hairpin,
			  type='l',
			  # xlim=c(max(0, min(c(codons, starts, start.new))-25), min(length(seq), max(c(codons, ends, start.new+(length(hairpin.new)-1)))+25)),
			  main='',
			  ylab='Hairpin Energy [kJ/mol]',
			  xlab='Position [bp]',
			  cex.lab=1.5)
		lines(x=bp, y=hairpin.new, col='gray', lwd=4)
		abline(v=codons, h=0.6, col=setColor('red', 0.2))
	}
	else
	{
		plot(bp, hairpin,
			  type='l',
			  # xlim=c(max(0, min(c(starts, start.new))-25), min(length(seq), max(c(ends, start.new+(length(hairpin.new)-1)))+25)),
			  main='',
			  ylab='Hairpin Energy [kJ/mol]',
			  xlab='Position [bp]',
			  cex.lab=1.5)
		lines(x=bp, y=hairpin.new, col='gray', lwd=4)
	}
	for(i in seq_along(starts))
	{
		shade(starts[i], ends[i], colors[i], alpha=0.2)
	}
}

plotGC <- function(seq, dbSeq, dbStart, dbEnd, starts, ends, colors)
{
	gc <- ifelse(seq %in% c("g","G","c","C"), 1, 0)
	gc.frac <- rollapply(gc, width=20, FUN=mean, partial=T, align='center')
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
	
	if(length(codons) > 0)
	{
		plot(seq_along(gc.frac), gc.frac,
			  type='l',
			  xlim=c(max(0, min(c(codons, starts, start.new))-25), min(length(seq), max(c(codons, ends, start.new+(length(gc.frac.new)-1)))+25)),
			  main='',
			  ylim=c(0.35,1),
			  ylab='GC Fraction',
			  xlab='Position [bp]',
			  cex.lab=1.5)
		lines(x=start.new:(start.new+(length(gc.frac.new)-1)), y=gc.frac.new, col='gray', lwd=4)
		abline(v=codons, h=0.6, col=setColor('red', 0.2))
	}
	else
	{
		plot(seq_along(gc.frac), gc.frac,
			  type='l',
			  xlim=c(max(0, min(c(starts, start.new))-25), min(length(seq), max(c(ends, start.new+(length(gc.new)-1)))+25)),
			  main='',
			  ylim=c(0.35,1),
			  ylab='GC Fraction',
			  xlab='Position [bp]',
			  cex.lab=1.5)
		lines(x=start.new:(start.new+(length(gc.frac.new)-1)), y=gc.frac.new, col='gray', lwd=4)
	}
	for(i in seq_along(starts))
	{
		shade(starts[i], ends[i], colors[i], alpha=0.2)
	}
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

getPrimerStats <- function(x, type)
{
	hp <- fold(x)
	homd <- homodimer(x)
	ret <- list(Tm=sig.digits(calculate_tm(toUpper(x)), 3),
					HpTm=sig.digits(hp$temp, 3),
					HpDeltaG=sig.digits(hp$dg/1000, 3),
					HpStruct=ifelse(hp$structure_found, hp$structure, ''),
					HmdTm=sig.digits(homd$temp, 3),
					HmdDeltaG=sig.digits(homd$dg/1000, 3),
					HmdStruct=ifelse(homd$structure_found, homd$structure, ''))
	return(ret)
}

getDimerStats <- function(x, y)
{
	hetd <- heterodimer(x, y)
	ret <- list(HtdTm=sig.digits(hetd$temp, 3),
					HtdDeltaG=sig.digits(hetd$dg/1000, 3),
					HtdStruct=ifelse(hetd$structure_found, hetd$structure, ''))
	return(ret)
}

# Define UI for application that draws a histogram
ui <- fluidPage(
	tags$head(
		# Note the wrapping of the string in HTML()
		tags$style(HTML("
      markF3	{background-color: #99ff99; 	color: black;}
      markF2	{background-color: #66ffff;	color: black;}
      markLFc	{background-color: #ffff66;	color: black;}
      markF1	{background-color: #cc99ff;	color: black;}
      markB1c	{background-color: #cc99ff;	color: black;}
      markLB	{background-color: #ffff66;	color: black;}
      markB2c	{background-color: #66ffff;	color: black;}
      markB3c	{background-color: #99ff99;	color: black;}
      markPNAF	{background-color: #ff3300;	color: black;}
      markPNABc{background-color: #ff3300;	color: black;}
		"))
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
		sidebarPanel(
			h4("Target Sequence:"),
			tags$textarea(id="Seq", rows=10, cols=60, 
							  'atcgaccacttcggcaaccgccgcctgcgtacggtcggcgagctgatccaaaaccagatccgggtcggcatgtcgcggatggagcgggtggtccgggagcggatgaccacccaggacgtggaggcgatcacaccgcagacgttgatcaacatccggccggtggtcgccgcgatcaaggagttcttcggcaccagccagctgagccaattcatgGACcagaacaacccgctgtcggggttgaccCACaagcgccgactgTCGgcgctggggcccggcggtctgtcacgtgagcgtgccgggctggaggtccgcgacgtgcacccgtcgcactacggccggatgtgcccgatcgaaacccctgaggggcccaacatcggtctgatcggctcgctgtcggtgtacgcgcgggtcaacccgttcgggttcatcgaaacgccgtaccgcaaggtggtcgacggcgtggttagcgacgagatcgtgtacctgaccgccgacgagga'
			),
			hr(),
			fluidPage(fluidRow(column(6, numericInput("polyT", "FIP/BIP PolyT-Spacer Length:", 3)),
									 column(6, numericInput("stabilityN", "End Stability Bp's:", 5))),
						 textInput('revCTool', 'Rev. Compl. Tool', placeholder='Auto-copy RevC to Clipboard'),
						 textOutput('revCOutput')),
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
		mainPanel(
			plotOutput("GCPlot"),
			plotOutput("HairpinPlot"),
			plotOutput("EnergyPlot"),
			plotOutput("TmPlot"),
			uiOutput("coloredSeq"),
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

getStartAndLenFromSeq <- function(primer, mySeq)
{
	temp <- pwalign::pairwiseAlignment(primer, mySeq, type='local', gapOpening=1000000, gapExtension=1000000)
	return(start(subject(temp)))
}

renderPrimerControl <- function(input, output, session, vals, controlId, mySeq, initLocs=c(1), Loc=1, initLen=20, initNTs='')
{
	fluidPage({
		fluidRow(
			column(1, h4(controlId)),
			column(1, checkboxInput(paste0(controlId, 'Check'), '', value=T)), #value=grepl('[12]', controlId))),
			column(2, numericInput(paste0(controlId, 'Start'), 'Start', value=ifelse(Loc > 0, initLocs[Loc], 1))),
			column(2, numericInput(paste0(controlId, 'Len'), 'Length', value=initLen)),
			column(6, textInput(paste0(controlId, 'NTs'), 'Seq', value=initNTs))
		)
	})
}

setUpControlLinks <- function(input, output, session, vals, controlId, mySeq, initLocs)
{
	observeEvent(input[[paste0(controlId, 'Start')]], {
		req(mySeq(), input[[paste0(controlId, 'Start')]])
		if(is.null(vals[[paste0(controlId, 'Start')]]) || vals[[paste0(controlId, 'Start')]] != input[[paste0(controlId, 'Start')]]) # Then go ahead and update
		{
			# Reset if the new start is less than 1
			if(input[[paste0(controlId, 'Start')]] <= 0)
			{
				# print(paste("Updating ", controlId, "Start1", sep=''))
				updateNumericInput(session, inputId=paste0(controlId, 'Start'), value=1)
			}
			else
			{
				# Reset the start if the new Start and Len are too big.
				if(input[[paste0(controlId, 'Start')]] > (length(mySeq()) - (input[[paste0(controlId, 'Len')]] - 1)))
				{
					# print(paste("Updating ", controlId, "Start2", sep=''))
					updateNumericInput(session, inputId=paste0(controlId, 'Start'), value=(length(mySeq()) - (input[[paste0(controlId, 'Len')]] - 1)))
				}
				else
				{
					# update reactive value for Start
					vals[[paste0(controlId, 'Start')]] <- input[[paste0(controlId, 'Start')]]
					
					# update linked controls
					newSeq <- getSeqFromStartAndLen(start=input[[paste0(controlId, 'Start')]], len=input[[paste0(controlId, 'Len')]], mySeq())
					if(input[[paste0(controlId, 'NTs')]] != paste(newSeq, collapse=''))
					{
						# Update reactive value and UI for NTs
						# print(paste("Updating ", controlId, "NTs1", sep=''))
						vals$updatingNTCount <- vals$updatingNTCount + 1
						vals[[paste0(controlId, 'NTs')]] <- paste(newSeq, collapse='')
						updateTextInput(session, inputId=paste0(controlId, 'NTs'), value=paste(newSeq, collapse=''))
					}
				}
			}
		}
	}, ignoreInit = F)
	
	observeEvent(input[[paste0(controlId, 'Len')]], {
		req(mySeq(), input[[paste0(controlId, 'Len')]])
		if(is.null(vals[[paste0(controlId, 'Len')]]) || vals[[paste0(controlId, 'Len')]] != input[[paste0(controlId, 'Len')]]) # Then go ahead and update
		{
			# Reset if the new start is less than 1
			if(input[[paste0(controlId, 'Len')]] <= 0)
			{
				# print(paste("Updating ", controlId, "Len1", sep=''))
				updateNumericInput(session, inputId=paste0(controlId, 'Len'), value=1)
			}
			else
			{
				# Reset the start if the new Start and Len are too big.
				if(input[[paste0(controlId, 'Start')]] > (length(mySeq()) - (input[[paste0(controlId, 'Len')]] - 1)))
				{
					# print(paste("Updating ", controlId, "Len2", sep=''))
					updateNumericInput(session, inputId=paste0(controlId, 'Len'), value=(length(mySeq()) - (input[[paste0(controlId, 'Start')]] + 1)))
				}
				else
				{
					# update reactive value for Len
					vals[[paste0(controlId, 'Len')]] <- input[[paste0(controlId, 'Len')]]
					
					# update linked controls
					newSeq <- getSeqFromStartAndLen(start=input[[paste0(controlId, 'Start')]], len=input[[paste0(controlId, 'Len')]], mySeq())
					if(input[[paste0(controlId, 'NTs')]] != paste(newSeq, collapse=''))
					{
						# Update reactive value and UI for NTs
						# print(paste("Updating ", controlId, "NTs2", sep=''))
						vals$updatingNTCount <- vals$updatingNTCount + 1
						vals[[paste0(controlId, 'NTs')]] <- paste(newSeq, collapse='')
						updateTextInput(session, inputId=paste0(controlId, 'NTs'), value=paste(newSeq, collapse=''))
					}
				}
			}
		}
	}, ignoreInit = F)
	
	observeEvent(input[[paste0(controlId, 'NTs')]], {
		if(!is.null(vals$updatingNTCount) && !vals$updatingNTCount && (is.null(vals[[paste0(controlId, 'NTs')]]) || vals[[paste0(controlId, 'NTs')]] != input[[paste0(controlId, 'NTs')]]))
		{
			req(mySeq(), input[[paste0(controlId, 'NTs')]])
			newPrimer <- s2c(input[[paste0(controlId, 'NTs')]])
			if(length(newPrimer) > 0 && all(newPrimer %in% c('a','g','t','c','A','G','T','C')))
			{
				primerStart <- getStartAndLenFromSeq(toUpper(paste(newPrimer, collapse='')), toUpper(paste(mySeq(), collapse='')))
				primerLen <- length(newPrimer)
				if(is.finite(primerStart) && primerStart > 0 && primerStart < (length(mySeq()) - (primerLen - 1)))
				{
					# update reactive value
					vals[[paste0(controlId, 'Start')]] <- primerStart
					vals[[paste0(controlId, 'Len')]] <- primerLen
					
					# update linked controls
					# print(paste("Updating ", controlId, "Start3", sep=''))
					updateNumericInput(session, inputId=paste0(controlId, 'Start'), value=primerStart)
					# print(paste("Updating ", controlId, "Len3", sep=''))
					updateNumericInput(session, inputId=paste0(controlId, 'Len'), value=primerLen)
				}
			}
		}
		vals$updatingNTCount <- vals$updatingNTCount - 1
	}, ignoreInit = F)
	
	observeEvent(mySeq(), {
	  req(mySeq(), input[[paste0(controlId, 'NTs')]], input[[paste0(controlId, 'Start')]], input[[paste0(controlId, 'Len')]])
	  primerStart <- input[[paste0(controlId, 'Start')]]
	  primerLen <- input[[paste0(controlId, 'Len')]]
	  if(is.finite(primerStart) && is.finite(primerLen) && primerStart > 0 && primerStart < (length(mySeq()) - (primerLen - 1)))
	  {
	    # Then start and len are still legal and just update the sequence to match the new sequence
	    newPrimer <- mySeq()[primerStart:(primerStart+primerLen-1)]
	    vals$updatingNTCount <- vals$updatingNTCount + 1
	    vals[[paste0(controlId, 'NTs')]] <- paste(newPrimer, collapse='')
	    updateTextInput(session, inputId = paste0(controlId, 'NTs'), value=paste(newPrimer, collapse=''))
	  }
	  else 
	  {
	    # Then the start or len are illegal and we need to reinitialize all the controls.
	    # We use the start and len inputs to trigger the update of the NTs box
	    indexes <- c(2,3,4,0,7,6,5,0,0,0)
	    names <- c('F3','F2','F1','LFc','B3c','B2c','B1c','LB','PNAF','PNABc')
	    thisIndex <- indexes[which(controlId == names)]
	    theStart <- ifelse(thisIndex > 0, initLocs()[thisIndex], 1)
	    theEnd <- ifelse(length(mySeq()) >= (theStart + primerLen - 1), theStart + primerLen -1, length(mySeq()))
	    theLen <- theEnd-theStart+1
	    theSeq <- mySeq()[theStart:theEnd]
	    updateNumericInput(session, inputId=paste0(controlId, 'Start'), value=theStart)
	    updateNumericInput(session, inputId=paste0(controlId, 'Len'), value=theLen)
	  }
	}, ignoreInit = F)
	
	observeEvent(list(input[[paste0(controlId, 'Check')]], vals$DBAll), {
		# Highlight output
		# print(paste("Marking ", controlId, ": ", input[[paste0(controlId, "NTs")]]))
		req(!is.null(input$Seq), !is.null(input$F3Check), !is.null(input$F2Check), !is.null(input$F1Check), !is.null(input$LFcCheck), !is.null(input$LBCheck), !is.null(input$B3cCheck), !is.null(input$B2cCheck), !is.null(input$B1cCheck), !is.null(input$PNAFCheck), !is.null(input$PNABcCheck))
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
}

getStart <- function(input, controlId)
{
	return(input[[paste0(controlId, 'Start')]])
}

getEnd <- function(input, controlId)
{
	return(input[[paste0(controlId, 'Start')]] + (input[[paste0(controlId, 'Len')]] - 1))
}

# Define server logic required to draw a histogram
server <- function(input, output, session) {
	
	vals <- reactiveValues(results=NULL,
								  DBAll=NULL,
								  updatingNTCount=0,
								  updatingStartOrLenCount=0)
	
	mySeq <- reactive({
		temp <- s2c(input$Seq)
	  temp[temp %in% c('a','g','t','c','A','G','T','C')]
	})
	
	mySeqHTML <- reactive({
		markup=p(id='text-to-mark', input$Seq)
	})
	
	initLocs <- reactive({
		req(mySeq())
		return(round(seq(1, length(mySeq())-20, length.out=8)))
	})
	
	output$coloredSeq <- renderUI(mySeqHTML())
	
	output$F3 <- renderUI(isolate(renderPrimerControl(input, output, session, vals, 'F3', mySeq, initLocs=initLocs(), Loc=2)))
	output$F2 <- renderUI(isolate(renderPrimerControl(input, output, session, vals, 'F2', mySeq, initLocs=initLocs(), Loc=3)))
	output$LFc <- renderUI(isolate(renderPrimerControl(input, output, session, vals, 'LFc', mySeq, initLocs=initLocs(), Loc=0)))
	output$F1 <- renderUI(isolate(renderPrimerControl(input, output, session, vals, 'F1', mySeq, initLocs=initLocs(), Loc=4)))
	output$B1c <- renderUI(isolate(renderPrimerControl(input, output, session, vals, 'B1c', mySeq, initLocs=initLocs(), Loc=5)))
	output$LB <- renderUI(isolate(renderPrimerControl(input, output, session, vals, 'LB', mySeq, initLocs=initLocs(), Loc=0)))
	output$B2c <- renderUI(isolate(renderPrimerControl(input, output, session, vals, 'B2c', mySeq, initLocs=initLocs(), Loc=6)))
	output$B3c <- renderUI(isolate(renderPrimerControl(input, output, session, vals, 'B3c', mySeq, initLocs=initLocs(), Loc=7)))
	output$PNAF <- renderUI(isolate(renderPrimerControl(input, output, session, vals, 'PNAF', mySeq, initLocs=initLocs(), Loc=0)))
	output$PNABc <- renderUI(isolate(renderPrimerControl(input, output, session, vals, 'PNABc', mySeq, initLocs=initLocs(), Loc=0)))
	
	setUpControlLinks(input, output, session, vals, 'F3', mySeq, initLocs)
	setUpControlLinks(input, output, session, vals, 'F2', mySeq, initLocs)
	setUpControlLinks(input, output, session, vals, 'LFc', mySeq, initLocs)
	setUpControlLinks(input, output, session, vals, 'F1', mySeq, initLocs)
	setUpControlLinks(input, output, session, vals, 'B1c', mySeq, initLocs)
	setUpControlLinks(input, output, session, vals, 'LB', mySeq, initLocs)
	setUpControlLinks(input, output, session, vals, 'B2c', mySeq, initLocs)
	setUpControlLinks(input, output, session, vals, 'B3c', mySeq, initLocs)
	setUpControlLinks(input, output, session, vals, 'PNAF', mySeq, initLocs)
	setUpControlLinks(input, output, session, vals, 'PNABc', mySeq, initLocs)
	
	starts <- reactive({
		c(getStart(input, 'F3'), 
		  getStart(input, 'F2'), 
		  getStart(input, 'F1'), 
		  getStart(input, 'LFc'), 
		  getStart(input, 'B1c'), 
		  getStart(input, 'B2c'), 
		  getStart(input, 'B3c'), 
		  getStart(input, 'LB'), 
		  getStart(input, 'PNAF'),
		  getStart(input, 'PNABc'))
	})
	
	stops <- reactive({
		c(getEnd(input, 'F3'), 
		  getEnd(input, 'F2'), 
		  getEnd(input, 'F1'), 
		  getEnd(input, 'LFc'), 
		  getEnd(input, 'B1c'), 
		  getEnd(input, 'B2c'), 
		  getEnd(input, 'B3c'), 
		  getEnd(input, 'LB'), 
		  getEnd(input, 'PNAF'),
		  getEnd(input, 'PNABc'))
	})
	
	primerColors <- reactive({
		c('#99ff99', # F3
		  '#66ffff', # F2
		  '#ffff66', # F1
		  '#cc99ff', # LFc
		  '#ffff66', # B1c
		  '#66ffff', # B2c
		  '#99ff99', # B3c
		  '#cc99ff', # LB
		  '#ff3300', # PNAF
		  '#ff3300') # PNABc
	})
	
	# Sense primers
	F3 <- reactive(input$F3NTs)
	F2 <- reactive(input$F2NTs)
	B1c <- reactive(input$B1cNTs)
	LB <- reactive(input$LBNTs)
	
	# Antisense primers
	B3 <- reactive(revC(input$B3cNTs, keepCase=T))
	B2 <- reactive(revC(input$B2cNTs, keepCase=T))
	F1c <- reactive(revC(input$F1NTs, keepCase=T))
	LF <- reactive(revC(input$LFcNTs, keepCase=T))
	
	# Combo primers
	FIP <- reactive(paste(c(F1c(), rep('t', input$polyT), F2()), collapse=''))
	BIP <- reactive(paste(c(B1c(), rep('t', input$polyT), B2()), collapse=''))
	
	# Amplicon
	observeEvent(list(input$F3NTs, input$F2NTs, input$F1NTs, input$LFcNTs, input$LBNTs, input$B3cNTs, input$B2cNTs, input$B1cNTs, input$PNAFNTs, input$PNABcNTs), 
					 {
					 	req(mySeq(), 
					 		 input$F3NTs != '', 
					 		 input$F2NTs != '', 
					 		 input$F1NTs != '', 
					 		 input$LFcNTs != '', 
					 		 input$LBNTs != '', 
					 		 input$B1cNTs != '', 
					 		 input$B2cNTs != '', 
					 		 input$B3cNTs != '', 
					 		 input$PNAFNTs != '',
					 		 input$PNABcNTs != '')
					 	vals$DBAll <- paste(c(F1c(), rep('t', input$polyT), mySeq()[getStart(input, 'F2'):(getEnd(input, 'B2c')-1)], rep('t', input$polyT), revC(B1c(), keepCase=T)), collapse='')
					 }
	)
	DBStart <- reactive(getStart(input, 'F2')-input$polyT-(getEnd(input, 'F1')-getStart(input, 'F1')+1))
	DBEnd <- reactive(getEnd(input, 'B2c')+(getEnd(input, 'B1c')-getStart(input, 'B1c')+1))
	
	# PNAs
	PNAF <- reactive({
		# If the PNA straddles the start of F2
		ifelse(getStart(input, 'F2') %in% c(getStart(input, 'PNAF'):getEnd(input, 'PNAF')),
				 paste(
				 	{
				 		# Paste together the upstream portion of FIP with the downstream portion of the PNAF
				 		FIPN <- length(s2c(FIP()))
				 		F2N <- input$F2Len
				 		Offset <- getStart(input, 'F2')-getStart(input, 'PNAF')
				 		if(getEnd(input, 'PNAF') > getEnd(input, 'F2'))
				 		{
				 			stop('Having a PNAF longer than F2 is not supported yet.')
				 		}
				 		c(s2c(FIP())[(FIPN-F2N-Offset):(FIPN-F2N)], # Chunk from start of PNAF up to start of F2 in FIP
				 		  s2c(input$PNAFNTs)[(Offset+1):input$PNAFLen]) # Remaining part of PNAF
				 	}, 
				 	collapse=''),
				 input$PNAFNTs
		)
	})
	PNAB <- reactive({
		# If the PNA straddles the end of B2c
		ifelse(getEnd(input, 'B2c') %in% c(getStart(input, 'PNABc'):getEnd(input, 'PNABc')),
				 paste(
				 	{
				 		# Paste together the upstream portion of FIP with the downstream portion of the PNAF
				 		BIPN <- length(s2c(BIP()))
				 		B2cN <- input$B2cLen
				 		Offset <- getEnd(input, 'PNABc')-getEnd(input, 'B2c')
				 		if(getStart(input, 'PNABc') < getStart(input, 'B2c'))
				 		{
				 			stop('Having a PNAB longer than B2 is not supported yet.')
				 		}
				 		c(s2c(BIP())[(BIPN-B2cN-Offset):(BIPN-B2cN)], # Chunk from end of PNABc up to end of B2c in BIP
				 		  revC(s2c(input$PNABcNTs), keepCase=T)[1:(input$PNAFLen-Offset)], keepCase=T) # Remaining (i.e., beginning) part of revC(PNABc)
				 	}, 
				 	collapse=''),
				 revC(input$PNABcNTs, keepCase=T)
		)
	})
	
	# Dumbbell regions
	DBF <- reactive(paste(c(rep('t', input$polyT), mySeq()[getStart(input, 'F2'):(getStart(input, 'F1')-1)]), collapse=''))
	DBB <- reactive(paste(c(rep('t', input$polyT), mySeq()[getEnd(input, 'B2c'):(getEnd(input, 'B1c')-1)]), collapse=''))
	
	allLegal <- reactive({
	  getEnd(input, 'F3') <= length(mySeq()) &&
	    getEnd(input, 'F2') <= length(mySeq()) &&
	    getEnd(input, 'F1') <= length(mySeq()) &&
	    getEnd(input, 'LFc') <= length(mySeq()) &&
	    getEnd(input, 'B3c') <= length(mySeq()) &&
	    getEnd(input, 'B2c') <= length(mySeq()) &&
	    getEnd(input, 'B1c') <= length(mySeq()) &&
	    getEnd(input, 'LB') <= length(mySeq()) &&
	    getEnd(input, 'PNAF') <= length(mySeq()) &&
	    getEnd(input, 'PNABc') <= length(mySeq())
	})
	
	observeEvent(list(input$Seq, input$F3NTs, input$F2NTs, input$F1NTs, input$LFcNTs, input$LBNTs, input$B3cNTs, input$B2cNTs, input$B1cNTs, input$PNAFNTs, input$PNABcNTs, input$polyT, input$stabilityN,
							input$F3Check, input$F2Check, input$F1Check, input$LFcCheck, input$LBCheck, input$B3cCheck, input$B2cCheck, input$B1cCheck, input$PNAFCheck, input$PNABcCheck), {
								
								# toInclude <- c(input$F3Check, input$B3Check, T, T, T, T, T, T, input$LFcCheck, input$LBCheck, input$PNACheck, input$PNAcCheck, T, T)
								req(mySeq(), allLegal(), input$F3NTs != '', input$F2NTs != '', input$LFcNTs != '', input$F1NTs != '', input$B1cNTs != '', input$LBNTs != '', input$B2cNTs != '', input$B3cNTs != '', input$PNAFNTs != '', input$PNABcNTs != '')
								ret <- data.table(Primer=c('F3','B3','F2','F1c','B2','B1c','FIP','BIP','LF','LB','PNAF','PNAB','DBF','DBB'),
														Seq=c(F3(), B3(), F2(), F1c(), B2(), B1c(), FIP(), BIP(), LF(), LB(), PNAF(), PNAB(), DBF(), DBB()),
														Sense=c('Sense','Antisense','Sense','Antisense','Sense','Antisense','NA','NA','Antisense','Sense','Sense','Antisense','Sense','Sense'))
								ret[, Stability3p:=end_stability(Seq, n=input$stabilityN, threePrimeEnd=T), by='Primer']
								ret[, Stability5p:=end_stability(Seq, n=input$stabilityN, threePrimeEnd=F), by='Primer']
								ret[, Start5p:='NA']
								ret[Primer == 'F3', Start5p:=getStart(input, 'F3')]
								ret[Primer == 'B3', Start5p:=getEnd(input, 'B3c')]
								ret[Primer == 'F2', Start5p:=getStart(input, 'F2')]
								ret[Primer == 'F1c', Start5p:=getEnd(input, 'F1')]
								ret[Primer == 'B2', Start5p:=getEnd(input, 'B2c')]
								ret[Primer == 'B1c', Start5p:=getStart(input, 'B1c')]
								ret[Primer == 'LF', Start5p:=getEnd(input, 'LFc')]
								ret[Primer == 'LB', Start5p:=getStart(input, 'LB')]
								ret[Primer == 'PNAF', Start5p:=getStart(input, 'PNAF')]
								ret[Primer == 'PNAB', Start5p:=getEnd(input, 'PNABc')]
								ret[, Len:=length(s2c(Seq)), by='Primer']
								ret[, c('Tm','HpTm','HpDeltaG','HpStruct','HmdTm','HmdDeltaG','HmdStruct'):=getPrimerStats(Seq, .BY[[1]]), by='Primer']
								ret[Primer %in% c('FIP','BIP'), c('HtdTm','HtdDeltaG','HtdStruct'):=getDimerStats(Seq, ifelse(.BY[[1]]=='FIP', BIP(), FIP())), by='Primer']
								ret[Primer %in% c('FIP','BIP'), Tm:='NA']
								ret[, KeyEndStability:=ifelse(Primer %in% c('F1c','B1c'), Stability5p, Stability3p)]
								
								# Data for energy plot
								ret2 <- melt.data.table(data=ret, id.vars='Primer', measure.vars=c('KeyEndStability','HpDeltaG','HmdDeltaG','HtdDeltaG'))
								ret2[, value:=suppressWarnings(as.numeric(value))]
								ret2[value > 0, value:=0]
								ret2[variable=='KeyEndStability', value:=-1*value]
								ret2[, Primer:=factor(ret2$Primer, levels=c('F3','B3','F2','B2','LF','LB','F1c','B1c','FIP','BIP','PNAF','PNAB','DBF','DBB'))]
								
								# Data for Tm plot
								ret3 <- ret[Primer %!in% c('DBF','DBB','FIP','BIP'), c('Primer','Tm')]
								ret3[, Tm:=suppressWarnings(as.numeric(Tm))]
								primerColors <- data.table(Primer=c('F3',
																				'F2',
																				'F1c',
																				'LF',
																				'B1c',
																				'B2',
																				'B3',
																				'LB',
																				'PNAF',
																				'PNAB'),
																	plotColor=c('#99ff99', # F3
																					'#66ffff', # F2
																					'#ffff66', # F1
																					'#cc99ff', # LFc
																					'#ffff66', # B1c
																					'#66ffff', # B2c
																					'#99ff99', # B3c
																					'#cc99ff', # LB
																					'#ff3300', # PNAF
																					'#ff3300') # PNABc
								)
								setkey(ret3, Primer)
								setkey(primerColors, Primer)
								ret3 <- primerColors[ret3]
								ret3[, Primer:=factor(ret3$Primer, levels=c('F3','B3','F2','B2','LF','LB','F1c','B1c','PNAF','PNAB','FIP','BIP','DBF','DBB'))]
								vals$results2 <- ret2
								vals$results3 <- ret3
								vals$results <- ret[Primer %in% c('F3',
																			 'B3',
																			 'F2',
																			 'F1c',
																			 'B2',
																			 'B1c',
																			 'LF',
																			 'LB',
																			 'PNAF',
																			 'PNAB',
																			 'FIP',
																			 'BIP',
																			 'DBF',
																			 'DBB')[
																			 	c(input$F3Check,
																			 	  input$B3cCheck,
																			 	  input$F2Check,
																			 	  input$F1Check,
																			 	  input$B2cCheck,
																			 	  input$B1cCheck,
																			 	  input$LFcCheck,
																			 	  input$LBCheck,
																			 	  input$PNAFCheck,
																			 	  input$PNABcCheck,
																			 	  input$F1Check && input$F2Check,
																			 	  input$B1cCheck && input$B2cCheck,
																			 	  T,
																			 	  T)]]
							})
	
	observeEvent(list(mySeq(), starts(), stops()), {
		output$HairpinPlot <- renderPlot({
			req(mySeq(), vals$DBAll)
			plotHairpin(mySeq(),
							dbSeq = isolate(vals$DBAll),
							dbStart = isolate(DBStart()),
							dbEnd = isolate(DBEnd()),
							starts = isolate(starts()), 
							ends   = isolate(stops()),
							colors = isolate(primerColors()))
		})
		
		output$GCPlot <- renderPlot({
			req(mySeq(), vals$DBAll)
			plotGC(seq = mySeq(),
					 dbSeq = isolate(vals$DBAll),
					 dbStart = isolate(DBStart()),
					 dbEnd = isolate(DBEnd()),
					 starts = isolate(starts()), 
					 ends   = isolate(stops()),
					 colors = isolate(primerColors()))
		})
		
		output$EnergyPlot <- renderPlot({
			req(mySeq(), vals$results2)
			ggplot(data=vals$results2[!is.na(value) & variable != 'KeyEndStability'], aes( x=Primer, y=value, fill=variable)) +
				geom_col() +
				geom_col(data=vals$results2[!is.na(value) & variable == 'KeyEndStability']) +
				geom_hline(yintercept=c(4,-4)) +
				labs(x='Sequence', y='Energy [kJ/mol]') +
				scale_y_continuous(limits=c(-20,10),oob = rescale_none)
		})
		
		output$TmPlot <- renderPlot({
			req(mySeq(), vals$results3)
			ggplot(data=vals$results3, aes(x=Primer, y=Tm, fill=plotColor)) + 
				geom_col() +
				geom_hline(yintercept=c(50,60,70)) +
				scale_y_continuous(limits=c(40,80),oob = rescale_none)
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
		req(mySeq(), vals$results)
		return(vals$results)
	})
	
	output$download <- downloadHandler(
		filename = function(){paste(input$downloadName, '.csv')}, 
		content = function(fname){
			fwrite(vals$results, fname)
		}
	)
	
	observe({
		updateCheckboxInput(session, inputId='F3Check', value=F)
		updateCheckboxInput(session, inputId='F3Check', value=T)
	})
}

# Run the application 
shinyApp(ui = ui, server = server)

