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
library(marker)
library(curl)

sourceGitHubFile <- function(user, repo, branch='main', file)
{
	destfile <- tempfile(fileext='.R')
	fileToGet <- paste0("https://raw.githubusercontent.com/", user, "/", repo, "/", branch, "/", file)
	curl_download(url=fileToGet, destfile)
	print(destfile)
	source(destfile)
}

sourceGitHubFile('Salus-Discovery','R-LampPrimerAnalyzer','main','LAMPUtils.R')
# source('~/Documents/GitHub/R-LampPrimerAnalyzer/LAMPUtils.R')

##### Calculated Constants #####

InitList <- data.table.expand.grid(Len=11:28)
InitList <- InitList[, getSeqsWithLen(rpoB_500_c, len=.BY[[1]], align='left'), by='Len']
InitList[, Tm:=getTm(Seq), by=c('Start','Len')]
InitList[, HmdDeltaG:=getHmd(Seq), by=c('Start','Len')]
InitList[, c('HpTm','HpDeltaG','HpStruct'):=getHp(Seq), by=c('Start','Len')]
InitList[is.na(HpDeltaG), HpDeltaG:=999]
InitList[, End:=Start+Len-1]
InitList[, Id:=1:.N]

##### Define UI for application #####
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

##### Server Helper Functions #####
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

##### Define server logic #####
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
	observeEvent(list(settings$seq, settings$polyT, settings$linker, getInputs(sensePrimerNames, 'NTs', input)), {# getInputs(sensePrimerNames, 'NTs', input)), {
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
			ntUpdated <- FALSE
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
				if(temp != input[[paste0(controlId, 'NTs')]]) # && silenceLenUpdate == 0)
				{
					ntUpdated <- TRUE
					print(paste("Update NTs UI from NTs change (e.g., seq change)", silenceLenUpdate))
					updateTextInput(session, paste0(controlId, 'NTs'), value=temp)
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
			if(!ntUpdated){
				# browser()
				silenceNTUpdate <<- silenceNTUpdate - 1
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
		req(settings$seq, allLegal(), all(as.logical(settings$NTs != '')))
		tables <- calcResultsTables(settings)
		settings$results <- tables$results
		settings$results2 <- tables$results2
		settings$results3 <- tables$results3
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

		filename = function(){paste(input$downloadName, '.zip', sep='')},
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

