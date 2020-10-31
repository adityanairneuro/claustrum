library(shiny)
library(pracma)
library(DT)
library(shinythemes)
require(neuralnet)
require(nnet)
library(httr)
library(plotly)
library(radarchart)
library(rhandsontable)
library(lime)
library(ggplot2)
library(caret)

source("Loading_Functions.R")
source("Analysis_Functions.R")

footer = h6("Nair, Graf and Augustine., Lee Kong Chian School of Medicine, Nanyang Technological University, Singapore ")




# Define UI ----
ui <- fluidPage(theme = shinytheme("yeti"),
                
                tags$head(
                  tags$style(HTML("
                                  @import url('//fonts.googleapis.com/css?family=Roboto|Cabin:400,700');
                                  "),
                             HTML(".shiny-notification {
              height: 100px;
                                  width: 800px;
                                  position:fixed;
                                  top: calc(80% - 50px);;
                                  left: calc(60% - 400px);;
                                  }
                                  "
                             )),
                  tags$link(rel = "icon", type = "image/png", href = "favicon.png")
                  ),
  
  # App title ----
    headerPanel(fluidRow(
    column(4,img(src="CLA_Logo_2.png", width = "100%")),
    column(6, img(src="CLA_12.png", width = "120%"))
    
    ), windowTitle = "Claustrum Classifier"
    
    
    ),
  
  
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select a file ----
      fileInput("file1", "Choose ABF File",
                multiple = FALSE,
                accept = c(".abf",
                           ".abfs")),
      
      # Horizontal line ----
      tags$hr(),
      
      numericInput("obs",
                   label = "Select trace to view:",
                   value = 10),
      # Horizontal line ----
      tags$hr(),
      
      selectInput("tracesel", "Criteria for trace selection",
                  choices = c("Current Threshold x 2" = "ct2")),
      
      # Horizontal line ----
      tags$hr(),
      
      
      fileInput("filebatch", "Batch Process Several ABF File",
                multiple = TRUE,
                accept = c(".abf",
                           ".abfs")),
      
      tags$hr(),
      
      uiOutput("uibatch")
      
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      tabsetPanel(
        tabPanel("Introduction",
                 br(),
                 p("The claustrum classifier is a program that allows users to upload ABF (Axon Binary File) format files to automatically extract electrophsiological information and classify cells according to the scheme written Graf & Augustine, 2017.",a("Abstract", href = "http://www.abstractsonline.com/pp8/#!/4376/presentation/33214",target = "_blank")),
                 p("It was written at the laboratory of ",a("George Augustine", href = "http://www.lkcmedicine.ntu.edu.sg/aboutus/Faculty-and-Staff/Pages/George-Augustine.aspx",target = "_blank") ," by ",a("Aditya Nair", href = "https://adityanairneuro.github.io/",target = "_blank"), " at the Lee Kong Chian School of Medicine, NTU, Singapore for the automated classification of claustral neurons."),
                 p(" The program uses several packages internally such as abf2, ggplot2, peakdet etc. A full list of packages and attributions can be found in the upcoming paper."),
                 br(),
                 img(src = "CLA_Methods_3.png", width = "100%"),
                 br(),
                 br(),
                 footer
        ),
        tabPanel("View your traces",
                 br(),
                 p("This tab allows the user to view and interact with traces after uploaded files in the ABF format."),
                 
      h4("Trace used for analysis"),
      plotlyOutput("selectedtraceplot"),
      verbatimTextOutput("trace10"),
      
      h4("Trace Output"),
      verbatimTextOutput("currentinjec"),
      plotlyOutput("traceplot")
      
      
      ),
      tabPanel("Extracted Properties",
               br(),
               p("This tab displays 14 extracted properties from the selected trace in the previous tab. The properties are described as per Graf et al., 2020"),
      h4("Result"),
      fluidRow(
      column (width = 6, DT::dataTableOutput("extractedproperties")),
      column (width = 6, h5("AP Number - Current Injection"),plotlyOutput("plotfi"),h5("Instantaneous Frequency") ,plotlyOutput("plotif") )
      )
      ),
               
      
      tabPanel("Classification",
               br(),
               p("The extracted properties are used to classify the given cell into the different cell-types described in Graf et al., 2020 using feedforward neural networks with a single hidden layer"),
                h4("Result"),
               fluidRow(
                 column( width = 10, textOutput("classification"),h5("Confidence of Classification"), DT::dataTableOutput("classficiationaccuracy"), plotOutput("plotLIME"))
               )
               
               
      )
      ,
      tabPanel("Manual Classifier",
               br(),
               p("This tab can be used to classify claustral neurons using properties calculated from other software and for instances where the ABF format is not used"),
               p("This example uses a PV interneuron"),
               h4("Enter your values"),
               fluidRow(
                 column (width = 4, rHandsontableOutput("hot")),
                 column (width = 8,textOutput("classificationmanual"),plotOutput("plotLIMEmanual"))
                 )
               
               
      )
      

      
    )
    
  )
))

# Define server logic ----
server <- function(input, output) {
  options(shiny.maxRequestSize=30*1024^2)
  
  
  response <- GET(url = "https://www.dropbox.com/s/l12ttqo86zc0xcy/Classifier_3_Level_May_2020_V1.RData?dl=1") 
  load(rawConnection(response$content))
  
  
# Required variables ------------------------------------------------------
  padNA <- function (mydata, rowsneeded, first = TRUE) 
  {
    temp1 = colnames(mydata)
    rowsneeded = rowsneeded - nrow(mydata)
    temp2 = setNames(
      data.frame(matrix(rep(NA, length(temp1) * rowsneeded), 
                        ncol = length(temp1))), temp1)
    if (isTRUE(first)) rbind(mydata, temp2)
    else rbind(temp2, mydata)
  }
  
  dotnames <- function(...) {
    vnames <- as.list(substitute(list(...)))[-1L]
    vnames <- unlist(lapply(vnames,deparse), FALSE, FALSE)
    vnames
  }
  Cbind <- function(..., first = TRUE) {
    Names <- dotnames(...)
    datalist <- setNames(list(...), Names)
    nrows <- max(sapply(datalist, function(x) 
      ifelse(is.null(dim(x)), length(x), nrow(x))))
    datalist <- lapply(seq_along(datalist), function(x) {
      z <- datalist[[x]]
      if (is.null(dim(z))) {
        z <- setNames(data.frame(z), Names[x])
      } else {
        if (is.null(colnames(z))) {
          colnames(z) <- paste(Names[x], sequence(ncol(z)), sep = "_")
        } else {
          colnames(z) <- paste(Names[x], colnames(z), sep = "_")
        }
      }
      padNA(z, rowsneeded = nrows, first = first)
    })
    do.call(cbind, datalist)
  }
  
  header2Def <- data.frame(
    field=c("headerSize", "episodes", "startDate", "startTime",
            "stopwatchTime", "fileType", "dataFormat", "nScan",
            "CRCEnable", "fileCRC", "fileGUID", "creatorVersion",
            "creatorNameIndex", "modifierVersion", "modifierNameIndex", "protocolPathIndex"),
    offset=c(8, 12, 16, 20,
             24, 28, 30, 32,
             34, 36, 40, 56,
             60, 64, 68, 72),
    type=c("uint32", "uint32", "uint32", "uint32",
           "uint32", "integer", "integer", "integer",
           "integer", "uint32", "uint32", "uint32",
           "uint32", "uint32", "uint32", "uint32"),
    bytes=c(4, 4, 4, 4,
            4, 2, 2, 2,
            2, 4, 4, 4,
            4, 4, 4, 4),
    stringsAsFactors=FALSE
  )
  
  # sections
  abfSections <- c(
    "ProtocolSection",
    "ADCSection",
    "DACSection",
    "EpochSection",
    "ADCPerDACSection",
    "EpochPerDACSection",
    "UserListSection",
    "StatsRegionSection",
    "MathSection",
    "StringSection",
    "DataSection",
    "TagSection",
    "ScopeSection",
    "DeltaSection",
    "VoiceTagSection",
    "SynchArraySection",
    "AnnotationSection",
    "StatsSection"
  )
  
  sectionInfo <- data.frame(
    field=c("blockIndex", "bytes", "numEntries"),
    type=c("uint32", "uint32", "integer"),
    bytes=c(4,4,8),
    stringsAsFactors=FALSE
  )
  
  abfProtocolInfo <- data.frame(
    field=c("operationMode", "ADCSequenceInterval", "fileCompressionEnabled", "unused",
            "fileCompressionRatio", "synchTimeUnit", "secondsPerRun", "samplesPerEpisode",
            "pretriggerSamples", "episodesPerRun", "runsPerTrial", "nTrials",
            "averagingMode", "undoRunCount", "firstEpInRun", "triggerThreshold",
            "triggerSource", "triggerAction", "triggerPolarity", "scopeOutputInterval",
            "episodeStartToStart", "runStartToStart", "averageCount", "trialStartToStart",
            "autoTriggerStrategy", "firstRunDelay", "channelStatsStrategy", "samplesPerTrace",
            "startDisplayNum", "finishDisplayNum", "showPNRawData", "statsPeriod",
            "statsMeasurements", "statsSaveStrategy", "ADCRange", "DACRange",
            "ADCResolution", "DACResolution", "experimentType", "manualInfoStrategy",
            "commentsEnable", "fileCommentIndex", "autoAnalyseEnable", "signalType",
            "digitalEnable", "activeDACChannel", "digitalHolding", "digitalInterEpisode",
            "digitalDACChannel", "digitalTrainActiveLogic", "statsEnable", "statsClearStrategy",
            "levelHysteresis", "timeHysteresis", "allowExternalTags", "averageAlgorithm",
            "averageWeighting", "undoPromptStrategy", "trialTriggerSource", "statsDisplayStrategy",
            "externalTagType", "scopeTriggerOut", "LTPType", "alternateDACOutputState",
            "alternateDigitalOutputState", "cellID1", "cellID2", "cellID3",
            "digitizerADCs", "digitizerDACs", "digitizerTotalDigOuts", "digitizerSynchDigOuts",
            "digitizerType"),
    type=c("integer", "numeric", "logical", "skip",
           "uint32", "numeric", "numeric", "integer",
           "integer", "integer", "integer", "integer",
           "integer", "integer", "integer", "numeric",
           "integer", "integer", "integer", "numeric",
           "numeric", "numeric", "integer", "numeric",
           "integer", "numeric", "integer", "integer",
           "integer", "integer", "integer", "numeric",
           "integer", "integer", "numeric", "numeric",
           "integer", "integer", "integer", "integer",
           "integer", "integer", "integer", "integer",
           "integer", "integer", "integer", "integer",
           "integer", "integer", "integer", "integer",
           "integer", "integer", "integer", "integer",
           "numeric", "integer", "integer", "integer",
           "integer", "integer", "integer", "integer",
           "integer", "numeric", "numeric", "numeric",
           "integer", "integer", "integer", "integer",
           "integer"),
    bytes=c(2, 4, 1, 3,
            4, 4, 4, 4,
            4, 4, 4, 4,
            2, 2, 2, 4,
            2, 2, 2, 4,
            4, 4, 4, 4,
            2, 4, 2, 4,
            4, 4, 2, 4,
            4, 2, 4, 4,
            4, 4, 2, 2,
            2, 4, 2, 2,
            2, 2, 2, 2,
            2, 2, 2, 2,
            2, 4, 2, 2,
            4, 2, 2, 2,
            2, 2, 2, 2,
            2, 4, 4, 4,
            2, 2, 2, 2,
            2),
    stringsAsFactors=FALSE
  )
  
  
  abfEpochDacInfo <- data.frame(
    field=c("nEpochNum", "nDACNum", "nEpochType", "fEpochInitLevel",
            "fEpochLevelInc", "lEpochInitDuration", "lEpochDurationInc", "lEpochPulsePeriod",
            "lEpochPulseWidth"),
    type=c("integer", "integer", "integer", "numeric",
           "numeric", "uint32", "uint32", "uint32",
           "uint32"),
    bytes=c(2, 2, 2, 8,
            8, 4, 4, 4,
            4),
    stringsAsFactors=FALSE
  )
  
  
  
  
  abfMathInfo <- data.frame(
    field=c("mathEnable", "mathExpression", "mathOperatorIndex", "mathUnitsIndex",
            "mathUpperLimit", "mathLowerLimit", "mathADCNum1", "mathADCNum2",
            "unused", "mathK1", "mathK2", "mathK3",
            "mathK4", "mathK5", "mathK6"),
    type=c("integer", "integer", "uint32", "uint32",
           "numeric", "numeric", "integer", "integer",
           "skip", "numeric", "numeric", "numeric",
           "numeric", "numeric", "numeric"),
    bytes=c(2, 2, 4, 4,
            4, 4, 2, 2,
            16, 4, 4, 4,
            4, 4, 4),
    stringsAsFactors=FALSE
  )
  
  abfADCInfo <- data.frame(
    field=c("ADCNum", "teleEnable", "teleInstrument", "teleAddGain",
            "teleFilter", "teleMembraneCap", "teleMode", "teleAccResist",
            "ADCPtoLChannelMap", "ADCSamplingSeq", "ADCProgGain", "ADCDispAmp",
            "ADCDispOffset", "instScaleFactor", "instOffset", "signalGain",
            "signalOffset", "signalLowpass", "signalHighpass", "lowpassType",
            "highpassType", "postprocLowpass", "postprocLowpassType", "enabledDuringPN",
            "statsChannelPolarity", "ADCChannelNameIndex", "ADCUnitsIndex"),
    type=c("integer", "integer", "integer", "numeric",
           "numeric", "numeric", "integer", "numeric",
           "integer", "integer", "numeric", "numeric",
           "numeric", "numeric", "numeric", "numeric",
           "numeric", "numeric", "numeric", "integer",
           "integer", "numeric", "integer", "logical",
           "integer", "integer", "integer"),
    bytes=c(2, 2, 2, 4,
            4, 4, 2, 4,
            2, 2, 4, 4,
            4, 4, 4, 4,
            4, 4, 4, 1,
            1, 4, 1, 1,
            2, 4, 4),
    stringsAsFactors=FALSE
  )
  
  abfTagInfo <- data.frame(
    field=c("tagTime", "comment", "tagType", "annotIndex"),
    type=c("integer", "string", "integer", "integer"),
    bytes=c(4, 56, 2, 2),
    stringsAsFactors=FALSE
  )
  abfEpochDacInfo <- data.frame(
    field=c("nEpochNum", "nDACNum", "nEpochType", "fEpochInitLevel",
            "fEpochLevelInc", "lEpochInitDuration", "lEpochDurationInc", "lEpochPulsePeriod",
            "lEpochPulseWidth"),
    type=c("integer", "integer", "integer", "numeric",
           "numeric", "uint32", "uint32", "uint32",
           "uint32"),
    bytes=c(2, 2, 2, 4,
            4, 8, 8, 8,
            8),
    stringsAsFactors=FALSE
  )
  
  
  
# Plots and Tables --------------------------------------------------------
  
  
  output$uibatch <- renderUI( {
    req(input$filebatch)
    #req(output$batch)
    tagList(
      tags$h6("Download Quantified Data"),
      downloadButton("downloadData", "Download Data")
    )
  })
  
# Obtain traces from ABF File  
  abftrace <- reactive({
    
    
    
    req(input$file1)
    req(input$obs)
    abftest <- abfloadv2(input$file1$datapath)
    numberoftraces <- abftest$episodes
    lengthofepisode <- length(abftest$traces[1,])/numberoftraces
    n <-  input$obs;
    if (n > numberoftraces){
      
      notrace <- NULL
      
    } else {
    n1L <- {n-1}*lengthofepisode+1
    n2L <- lengthofepisode*n
    tn <- abftest$traces[1,n1L:n2L];
    timeT<- abftest$s
    timeTn <- timeT[n1L:n2L]
    df.abf <- data.frame(timeTn,tn)
      
    }
    
      
      
    })
  
  # Check whether the traces has enough traces (Current Threshold x 2)
  ctx2check <- reactive({
    
    req(input$file1)
    abftest <- abfloadv2(input$file1$datapath)
    
    numberoftraces <- abftest$episodes
    lengthofepisode <- length(abftest$traces[1,])/numberoftraces
    sweep.max <- seq(1,numberoftraces,1)
    
    abftest$epoch$fEpochInitLevel
    abftest$epoch$fEpochLevelInc
    c2 = vector();
    for (i in sweep.max){
      
      n <-  i;
      n1L <- {n-1}*lengthofepisode+1
      n2L <- lengthofepisode*n
      tn <- abftest$traces[1,n1L:n2L];
      peak1list <- peakdet(tn, 20)
      example1  <- peak1list[[1]]
      if (length(example1[,1]) == 1){
        c2 <- i
        
      } else if (length(example1[,1]) > 1){
        c2<- i
        break
      }
      
      
    }
    
    creq <- c2-1;
    Ct <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(c2-1)
    
    # Ctx2
    Ctx2 <- Ct*2;
    Ctx2
    c2 <- vector();
    
    ct_times_two_c = (Ctx2 - abftest$epoch$fEpochInitLeve)/abftest$epoch$fEpochLevelInc
    ct_times_two = ct_times_two_c +1
    
    if (ct_times_two > numberoftraces ){
      
      notrace <- NULL 
      
    } else{
      
      ct_times_two
      
    }
    
    
    
  })
  
  # Get trace number for selected condition, example: Ctx2
  tracenumber <- reactive({
    
    
    
    req(input$file1)
    if (input$tracesel == "ap11"){
      c1 = vector()
      abftest <- abfloadv2(input$file1$datapath)
      numberoftraces <- abftest$episodes
      lengthofepisode <- length(abftest$traces[1,])/numberoftraces
      sweep.max <- seq(1,numberoftraces,1)
      
      withProgress(message = 'Finding the right trace to use', value = 0,{
        n <- length(sweep.max)
        for (i in sweep.max){
          n <-  i;
          n1L <- {n-1}*lengthofepisode+1
          n2L <- lengthofepisode*n
          tn <- abftest$traces[1,n1L:n2L];
          peak1list <- peakdet(tn, 20)
          example1  <- peak1list[[1]]
          if (length(example1[,1]) == 11){
            c1 <- i
            break
          } else if (length(example1[,1]) > 11){
            c1<- i
            break
          }
          incProgress(1/n, detail = paste("Doing part", i))
          
        }
        
      })
      c1
      
    } else if (input$tracesel == "ct2"){
      
      abftest <- abfloadv2(input$file1$datapath)
      numberoftraces <- abftest$episodes
      lengthofepisode <- length(abftest$traces[1,])/numberoftraces
      sweep.max <- seq(1,numberoftraces,1)
      
      abftest$epoch$fEpochInitLevel
      abftest$epoch$fEpochLevelInc
      c2 = vector();
      for (i in sweep.max){
        
        n <-  i;
        n1L <- {n-1}*lengthofepisode+1
        n2L <- lengthofepisode*n
        tn <- abftest$traces[1,n1L:n2L];
        peak1list <- peakdet(tn, 20)
        example1  <- peak1list[[1]]
        if (length(example1[,1]) == 1){
          c2 <- i
          
        } else if (length(example1[,1]) > 1){
          c2<- i
          break
        }
        
        
      }
      creq <- c2-1;
      
      
      Ct <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(c2-1)
      
      # Ctx2
      Ctx2 <- Ct*2;
      Ctx2
      c2 <- vector();
      
      ct_times_two_c = (Ctx2 - abftest$epoch$fEpochInitLeve)/abftest$epoch$fEpochLevelInc
      ct_times_two = ct_times_two_c +1
      
      
    } else if (input$tracesel == "ct1"){
      abftest <- abfloadv2(input$file1$datapath)
      numberoftraces <- abftest$episodes
      lengthofepisode <- length(abftest$traces[1,])/numberoftraces
      sweep.max <- seq(1,numberoftraces,1)
      
      abftest$epoch$fEpochInitLevel
      abftest$epoch$fEpochLevelInc
      c2 = vector();
      for (i in sweep.max){
        
        n <-  i;
        n1L <- {n-1}*lengthofepisode+1
        n2L <- lengthofepisode*n
        tn <- abftest$traces[1,n1L:n2L];
        peak1list <- peakdet(tn, 20)
        example1  <- peak1list[[1]]
        if (length(example1[,1]) == 1){
          c2 <- i
          
        } else if (length(example1[,1]) > 1){
          c2<- i
          break
        }
        
        
      }
      creq <- c2-1;
      ct_plus_one <- c2+1
      
      Ct <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(c2-1)
      
      # Ctx2
      Ctx2 <- Ct*2;
      Ct_1 <- Ct + abftest$epoch$fEpochLevelInc;
      Ct_1
      c2 <- vector();
      
      ct_plus_one
      
      
      
      
      
    }
    
    
    
  })
  


  
 # Plot traces selected  
  output$traceplot <- renderPlotly({
    
    req(input$file1)
    req(input$obs)
    if( is.null(abftrace()) ){
      
      text = paste("\n   There is no trace to display.\n",
                   "       Please try a step protocol with more current injections\n")
      ggplot() + 
        annotate("text", x = 4, y = 25, size=4, label = text) + 
        theme_bw() +
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.position="none",
              panel.background=element_blank(),
              panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank()) + labs(x=NULL, y=NULL)
      
    } else{
      traceplot <- ggplot(abftrace(),aes(x=timeTn, y=tn)) + geom_line(colour="red") + labs(x = "Time (sec)",y = "Membrane Potential (mV)")
      traceplot + labs(x = "Time (sec)",y = "Membrane Potential (mV)")
      traceplot <- plotly_build(traceplot, registerFrames = TRUE)
      traceplot$elementId <- NULL
    ggplotly(traceplot)
    }
  }
    
  )
  
  # Obtain full trace for selected condition, example: Ctx2
  
  selectedtrace <- reactive({
    
    
    
    req(input$file1)
    c1 = vector()
    abftest <- abfloadv2(input$file1$datapath)
    numberoftraces <- abftest$episodes
    lengthofepisode <- length(abftest$traces[1,])/numberoftraces
    sweep.max <- seq(1,numberoftraces,1)
    
    n <-  tracenumber();
    if (length(n) == 0){
      
      notrace <- NULL
      
    } else {
    n1L <- {n-1}*lengthofepisode+1
    n2L <- lengthofepisode*n
    tn <- abftest$traces[1,n1L:n2L];
    timeT<- abftest$s
    timeTn <- timeT[n1L:n2L]
    df.abf2 <- data.frame(timeTn,tn)
    
    }
  })
  
  
  
# Plot trace required for analysis  
  
  output$selectedtraceplot <- renderPlotly({
    
    if( is.null(selectedtrace()) ){
      
     
      text = paste("\n   There is no trace to display as the file does not have enough traces for analysis.\n",
                   "       Please choose another file with enough traces and current injection steps\n")
      ggplot() + 
        annotate("text", x = 4, y = 25, size=4, label = text) + 
        theme_bw() +
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.position="none",
              panel.background=element_blank(),
              panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank()) + labs(x=NULL, y=NULL)
      
    } else{
      selectedtraceplot <- ggplot(selectedtrace(),aes(x=timeTn, y=tn)) + geom_line(colour="blue") + labs(x = "Time (sec)",y = "Membrane Potential (mV)")
      selectedtraceplot + labs(x = "Time (sec)",y = "Membrane Potential (mV)")
      selectedtraceplot <- plotly_build(selectedtraceplot, registerFrames = TRUE)
      selectedtraceplot$elementId <- NULL
    ggplotly(selectedtraceplot)
    
    
    
    }}
  
  )
  
  output$trace10 <- renderText({ 
    
    paste("Selected trace is", as.character(tracenumber()))
    
  })
   
  output$currentinjec <- renderText({ 
    
    req(input$file1)
    req(input$obs)
    
    abftest <- abfloadv2(input$file1$datapath)
    
    CurrentInjection <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(input$obs-1)
    
    paste("Current Injection is ", as.character(CurrentInjection)," pA")
    
  })
  
## Calculation of parameters from ABF files
  
    
#Action Potential Properties
  
# Get AP properties from selected traces 
  
  apcalc <- reactive({
      
    req(input$file1)
    #req(ctx2check())
    if (is.null(ctx2check())){
      
      withProgress(message = 'Trace Analysis Begins', value = 0,{
        
        c1 = vector()
        abftest <- abfloadv2(input$file1$datapath)
        numberoftraces <- abftest$episodes
        lengthofepisode <- length(abftest$traces[1,])/numberoftraces
        sweep.max <- seq(1,numberoftraces,1)
        
        
        n <-  tracenumber();
        n1L <- {n-1}*lengthofepisode+1
        n2L <- lengthofepisode*n
        tn <- abftest$traces[1,n1L:n2L];
        timeT<- abftest$s
        timeTn <- timeT[n1L:n2L]
        #df.abf <- data.frame(timeTn,tn)
        #df <- data.frame(abftest1$dataX,abftest1$dataY)
        #ggplot(df,aes(x=abftest$dataX, y=abftest$dataY)) + geom_line()
        peak1 <- tn;
        time1 <-  timeTn*1000;
        dy <- diff(peak1)/diff(time1);
        #df1 <- data.frame(time1[2:100000],dy);
        #ggplot(df1,aes(x=time1[1:99999], y=dy)) + geom_line()
        peakMag2 <- peakdet(dy,10)
        
        peakMag2x <- peakMag2[[1]]
        
        locfx <- which.max(peakMag2x[,2])
        locfx2 <- peakMag2x[locfx,1]
        
        ds <- dy[1:locfx2]
        dsf <- rev(ds)
        
        tn2<- vector();
        
        for (n in 1:locfx2){
          if (dsf[n] < 10){
            tn2 =n
            break
          }
          
        }
        
        
        
        ds2 <- peak1[1:locfx2];
        ds3 <- rev(ds2);
        apthres <- ds3[tn2]
        
        n <-  tracenumber();
        
        n1L <- {n-1}*lengthofepisode+1
        n2L <- lengthofepisode*n
        tn <- abftest$traces[1,n1L:n2L];
        peak1list <- peakdet(tn, 10)
        example1  <- peak1list[[1]]
        
        df.abf <- data.frame(timeTn,tn)
        #df <- data.frame(abftest1$dataX,abftest1$dataY)
        ggplot(df.abf,aes(x=timeTn, y=tn)) + geom_line()
        
        
        
        APAmp <- example1[1,2] - apthres
        APAmp2 <- example1[2,2] - apthres
        APAmp3 <- example1[3,2] - apthres
        APAmpEnd <- example1[end,2] - apthres
       
        
        AP_1_2_R <- APAmp/APAmp2;
        AP_2_3_R <- APAmp2/APAmp3;
        AP_End_3_R <- APAmp/APAmpEnd;
        
        APThres_APAmp <- apthres/APAmp;
        
        #APRise and Fall
        #Rise
        tnMax <- tn[1:example1[1,1]];
        timeTnMax <- timeTn[1:example1[1,1]]
        peaktnMax <- tnMax;
        time1tnMax <-  timeTnMax*1000;
        dytnMax <- diff(peaktnMax)/diff(time1tnMax);
        peakMag2Max <- peakdet(dytnMax,10)
        #df <- data.frame(timeTnMax,tnMax)
        #ggplot(df,aes(x=timeTnMax, y=tnMax)) + geom_line()
        
        peakMag2x <- peakMag2Max[[1]]
        peakMag3x <- peakMag2x[[2]]
        Max_Rise <- (max(peakMag3x[]))
        
        
        #Fall
        tnFall <- tn[example1[1,1]:example1[2,1]];
        timeTnFall <- timeTn[example1[1,1]:example1[2,1]]
        peaktnFall <- tnFall;
        time1tnFall <-  timeTnFall*1000;
        dytnFall <- diff(peaktnFall)/diff(time1tnFall);
        peakMag2Fall <- peakdet(dytnFall,10)
        #df <- data.frame(timeTnFall,tnFall)
        #ggplot(df,aes(x=timeTnFall, y=tnFall)) + geom_line()
        
        peakMag2xF <- peakMag2Fall[[2]]
        peakMag3xF <- peakMag2xF[[2]]
        Max_Fall <- (min(peakMag3xF[]))
        
        
        peakLoc <- example1[[1]]
        
        
        
        #AHP
        peakAHP <-peakdet(tnFall,10)
        peakAHP2 <- peakAHP[[2]]
        peakAHP3 <-  peakAHP2[[2]]
        yAHP <- min(peakAHP3[])
        
        AHPAmp <- yAHP - apthres;
        
        AP_Rise_Decay_R <- -Max_Rise/Max_Fall;
        
        #HalfWidth
        
        half1 <- apthres + APAmp/2;
        yp1x <- peakLoc[1]-200
        yp1 <- peak1[yp1x:peakLoc[1]]
        tt1 <- which.min(abs(yp1 - half1)) 
        time <- timeTn
        
        timett1 <- time[yp1x:peakLoc[1]]
        
        tt1x <- timett1[tt1]
        yp2x <- peakLoc[1]+200
        yp2 <- peak1[peakLoc[1]:yp2x]
        
        tt2 <- which.min(abs(yp2 - half1)) 
        
        timett2 <- time[peakLoc[1]:yp2x]
        tt2x <- timett2[tt2]
        
        halfwidth <- (tt2x - tt1x)*1000
        
        AP_Rise_Half <- Max_Rise/halfwidth
        
        
        ##ISIRatio
        end <- length(peakLoc)
        sumisi <- (peakLoc[end]-peakLoc[end-1])+ (peakLoc[end-1]-peakLoc[end-2])+ (peakLoc[end-2]-peakLoc[end-3]);
        ISIRatio <- (peakLoc[2]-peakLoc[1])/( sumisi/3);
        
        n3<- 3
        incProgress(1/n3, detail = paste("Starting Analysis"))
        
        ##Max Adaptation
        tnl<- length(tn)/2
        
        sumadapt <- 1/((peakLoc[end]-peakLoc[end-1])/tnl)+ 1/((peakLoc[end-1]-peakLoc[end-2])/tnl)+ 1/((peakLoc[end-2]-peakLoc[end-3])/tnl);
        
        
        MaxAdaptation <- (1/((peakLoc[2]-peakLoc[1])/tnl)- (sumadapt/3) );
        
        ##Cv2 and CV1
        
        nCV <- length(peakLoc);
        nCV1 <- nCV-1
        CV= vector();
        
        for (i in 1:nCV1){
          
          CV<- append(CV, (peakLoc[i+1]-peakLoc[i])/50)
          
        }
        
        CVn <- vector();
        n2CV <- length(CV);
        n2CV1 <- n2CV-1;
        
        for (i in 1:n2CV1){
          
          CVn <- append(CVn, 2*( abs((CV[i]-CV[i+1]))/((CV[i]+CV[i+1]))))
          
        }
        
        n3CV <- length(CVn);
        CVall <- mean(CVn[1:n3CV]);
        CV1 <- mean(CVn[2:n3CV]);
        CV2 <- mean(CVn[3:n3CV]);
        
        #Rm, RMP, Ct
        abftest$epoch$fEpochInitLevel
        abftest$epoch$fEpochLevelInc
        c2 = vector();
        for (i in sweep.max){
          
          n <-  i;
          n1L <- {n-1}*lengthofepisode+1
          n2L <- lengthofepisode*n
          tn <- abftest$traces[1,n1L:n2L];
          peak1list <- peakdet(tn, 20)
          example1  <- peak1list[[1]]
          if (length(example1[,1]) == 1){
            c2 <- i
            
          } else if (length(example1[,1]) > 1){
            c2<- i
            break
          }
          
          
        }
        creq <- c2-1;
        iv.v <-  vector();
        iv.i <-  vector();
        withProgress(message = 'Analyzing traces', value = 0,{
          n2 <- creq
          for (i in 1:creq){
            n <-  i;
            n1L <- {n-1}*lengthofepisode+1
            n2L <- lengthofepisode*n
            tn <- abftest$traces[1,n1L:n2L];
            lt <- length(tn);
            peakrmp <- tn[lt/2];
            peaki <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(i-1)
            iv.v<- append(iv.v,peakrmp)
            iv.i<- append(iv.i,peaki)
            
            incProgress(1/n2, detail = paste("Doing part", i))
            
          }})
        
        P <- polyfit(iv.i,iv.v,1)
        Rm <- P[1]*1000
        RMP <- P[2]
        
        Ct <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(c2-1)
        
        
        
        n3<- 5
        incProgress(1/n3, detail = paste("Halfway there!"))
        
        
        ###Max Ap
        
        apnumberlist <- vector();
        iv.i.entire <- vector();
        for (i in sweep.max){
          
          n <-  i;
          n1L <- {n-1}*lengthofepisode+1
          n2L <- lengthofepisode*n
          tn <- abftest$traces[1,n1L:n2L];
          peak1 <- tn
          peak1list <- peakdet(peak1, 20)
          example1  <- peak1list[[1]]
          apnumber2 <- length(example1[,1]) 
          apnumberlist<- append(apnumberlist,apnumber2)
          peaki2 <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(i-1)
          iv.i.entire<- append(iv.i.entire,peaki2)
        }  
        
        n3<- 7
        incProgress(1/n3, detail = paste("Almost!"))
        
        apnumberlist  
        iv.i.entire  
        
        #dfap <- data.frame(iv.i.entire,apnumberlist)
        #ggplot(dfap,aes(x=iv.i.entire, y=apnumberlist)) + geom_line()
        
        dyap <- diff(apnumberlist)/diff(iv.i.entire)
        
        maxaploc <- which.max(apnumberlist)
        maxapsweep <- maxaploc
        
        n <-  maxapsweep;
        n1L <- {n-1}*lengthofepisode+1
        n2L <- lengthofepisode*n
        tn <- abftest$traces[1,n1L:n2L];
        peak1 <- tn
        peak1list <- peakdet(peak1, 20)
        example1  <- peak1list[[1]]
        apnumbermax <- length(example1[,1]) 
        peakLoc <- example1[,1]
        #ggplot(dfap,aes(x=iv.i.entire, y=apnumberlist)) + geom_line()
        
        ##ISIRatio
        end <- length(peakLoc)
        sumisimax <- (peakLoc[end]-peakLoc[end-1])+ (peakLoc[end-1]-peakLoc[end-2])+ (peakLoc[end-2]-peakLoc[end-3]);
        ISIRatiomax <- (peakLoc[2]-peakLoc[1])/( sumisimax/3);
        
        
        
        ##Max Adaptation
        
        tnl<- length(tn)/2
        sumadaptmax <- 1/((peakLoc[end]-peakLoc[end-1])/tnl)+ 1/((peakLoc[end-1]-peakLoc[end-2])/tnl)+ 1/((peakLoc[end-2]-peakLoc[end-3])/tnl);
        
        
        MaxAdaptationmax <- (1/((peakLoc[2]-peakLoc[1])/tnl)- (sumadaptmax/3) );
        
        
        
        
        Value <- c(RMP, Rm, Ct, Max_Fall, AHPAmp,ISIRatio,MaxAdaptation,CVall,CV2,AP_Rise_Decay_R,CV1,APThres_APAmp,AP_Rise_Half,AP_1_2_R,AP_2_3_R, apnumbermax,MaxAdaptationmax,ISIRatiomax,halfwidth)
        
        
        #Property.Name <- c("RMP","Rm","Ct","Max Fall Rate","AHP Amp","ISI Ratio","Max Adapt","CV2","CV2- 1st/2nd AP","AP Rise Decay R","CV2 1st AP","AP Thres/Amp","AP Rise/HalfWidth","1st/2nd AP","2nd/3rd AP","AP Number at Ctx2","Max Adaptation at Ctx2","ISIRatio at Ctx2","Max AP Number","Max Adaptation at Max AP No","ISIRatio at Max AP No", "Halfwidth")
        
        #data1 <- data.frame(Property.Name, Value);
        
        n3<- 10
        incProgress(1/n3, detail = paste("Done!"))
      })
      Value 
      
    }else {
    
    withProgress(message = 'Trace Analysis Begins', value = 0,{
      
      c1 = vector()
      abftest <- abfloadv2(input$file1$datapath)
      numberoftraces <- abftest$episodes
      lengthofepisode <- length(abftest$traces[1,])/numberoftraces
      sweep.max <- seq(1,numberoftraces,1)
      
      
      n <-  tracenumber();
      n1L <- {n-1}*lengthofepisode+1
      n2L <- lengthofepisode*n
      tn <- abftest$traces[1,n1L:n2L];
      timeT<- abftest$s
      timeTn <- timeT[n1L:n2L]
      #df.abf <- data.frame(timeTn,tn)
      #df <- data.frame(abftest1$dataX,abftest1$dataY)
      #ggplot(df,aes(x=abftest$dataX, y=abftest$dataY)) + geom_line()
      peak1 <- tn;
      time1 <-  timeTn*1000;
      dy <- diff(peak1)/diff(time1);
      #df1 <- data.frame(time1[2:100000],dy);
      #ggplot(df1,aes(x=time1[1:99999], y=dy)) + geom_line()
      peakMag2 <- peakdet(dy,10)
      
      peakMag2x <- peakMag2[[1]]
      
      locfx <- which.max(peakMag2x[,2])
      locfx2 <- peakMag2x[locfx,1]
      
      ds <- dy[1:locfx2]
      dsf <- rev(ds)
      
      tn2<- vector();
      
      for (n in 1:locfx2){
        if (dsf[n] < 10){
          tn2 =n
          break
        }
        
      }
      
      
      
      ds2 <- peak1[1:locfx2];
      ds3 <- rev(ds2);
      apthres <- ds3[tn2]
      
      n <-  tracenumber();
      
      n1L <- {n-1}*lengthofepisode+1
      n2L <- lengthofepisode*n
      tn <- abftest$traces[1,n1L:n2L];
      peak1list <- peakdet(tn, 20)
      example1  <- peak1list[[1]]
      
      df.abf <- data.frame(timeTn,tn)
      #df <- data.frame(abftest1$dataX,abftest1$dataY)
      #ggplot(df.abf,aes(x=timeTn, y=tn)) + geom_line()
      
    
      APAmp <- example1[1,2] - apthres
      APAmp2 <- example1[2,2] - apthres
      APAmp3 <- example1[3,2] - apthres
      APLast <- tail(example1, n = 1);
      APAmpEnd <- APLast[1,2] - apthres
      
      AP_1_2_R <- APAmp/APAmp2;
      AP_2_3_R <- APAmp2/APAmp3;
      AP_end_3_R_OG <- APAmp/APAmpEnd;
      
      APThres_APAmp <- apthres/APAmp;
      
      #APRise and Fall
      #Rise
      tnMax <- tn[1:example1[1,1]];
      timeTnMax <- timeTn[1:example1[1,1]]
      peaktnMax <- tnMax;
      time1tnMax <-  timeTnMax*1000;
      dytnMax <- diff(peaktnMax)/diff(time1tnMax);
      peakMag2Max <- peakdet(dytnMax,10)
      #df <- data.frame(timeTnMax,tnMax)
      #ggplot(df,aes(x=timeTnMax, y=tnMax)) + geom_line()
      
      peakMag2x <- peakMag2Max[[1]]
      peakMag3x <- peakMag2x[[2]]
      Max_Rise <- (max(peakMag3x[]))
      
      
      #Fall
      tnFall <- tn[example1[1,1]:example1[2,1]];
      timeTnFall <- timeTn[example1[1,1]:example1[2,1]]
      peaktnFall <- tnFall;
      time1tnFall <-  timeTnFall*1000;
      dytnFall <- diff(peaktnFall)/diff(time1tnFall);
      peakMag2Fall <- peakdet(dytnFall,10)
      #df <- data.frame(timeTnFall,tnFall)
      #ggplot(df,aes(x=timeTnFall, y=tnFall)) + geom_line()
      
      peakMag2xF <- peakMag2Fall[[2]]
      peakMag3xF <- peakMag2xF[[2]]
      Max_Fall <- (min(peakMag3xF[]))
      
      
      peakLoc <- example1[[1]]
      
      
      
      #AHP
      peakAHP <-peakdet(tnFall,10)
      peakAHP2 <- peakAHP[[2]]
      peakAHP3 <-  peakAHP2[[2]]
      yAHP <- min(peakAHP3[])
      
      AHPAmp <- yAHP - apthres;
      
      AP_Rise_Decay_R <- -Max_Rise/Max_Fall;
      
      #HalfWidth
      
      half1 <- apthres + APAmp/2;
      yp1x <- peakLoc[1]-200
      yp1 <- peak1[yp1x:peakLoc[1]]
      tt1 <- which.min(abs(yp1 - half1)) 
      time <- timeTn
      
      timett1 <- time[yp1x:peakLoc[1]]
      
      tt1x <- timett1[tt1]
      yp2x <- peakLoc[1]+200
      yp2 <- peak1[peakLoc[1]:yp2x]
      
      tt2 <- which.min(abs(yp2 - half1)) 
      
      timett2 <- time[peakLoc[1]:yp2x]
      tt2x <- timett2[tt2]
      
      halfwidth <- (tt2x - tt1x)*1000
      
      AP_Rise_Half <- Max_Rise/halfwidth
      
      Max_Rise_OG <- Max_Rise;
      
      ##ISIRatio
      end <- length(peakLoc)
      sumisi <- (peakLoc[end]-peakLoc[end-1])+ (peakLoc[end-1]-peakLoc[end-2])+ (peakLoc[end-2]-peakLoc[end-3]);
      ISIRatio <- (peakLoc[2]-peakLoc[1])/( sumisi/3);
      
      n3<- 3
      incProgress(1/n3, detail = paste("Starting Analysis"))
      
      ##Max Adaptation
      tnl<- length(tn)/2
      
      sumadapt <- 1/((peakLoc[end]-peakLoc[end-1])/tnl)+ 1/((peakLoc[end-1]-peakLoc[end-2])/tnl)+ 1/((peakLoc[end-2]-peakLoc[end-3])/tnl);
      
      
      MaxAdaptation <- (1/((peakLoc[2]-peakLoc[1])/tnl)- (sumadapt/3) );
      
      ##Cv2 and CV1
      
      nCV <- length(peakLoc);
      nCV1 <- nCV-1
      CV= vector();
      
      for (i in 1:nCV1){
        
        CV<- append(CV, (peakLoc[i+1]-peakLoc[i])/50)
        
      }
      
      CVn <- vector();
      n2CV <- length(CV);
      n2CV1 <- n2CV-1;
      
      for (i in 1:n2CV1){
        
        CVn <- append(CVn, 2*( abs((CV[i]-CV[i+1]))/((CV[i]+CV[i+1]))))
        
      }
      
      n3CV <- length(CVn);
      CVall <- mean(CVn[1:n3CV]);
      CV1 <- mean(CVn[2:n3CV]);
      CV2 <- mean(CVn[3:n3CV]);
      
      #Rm, RMP, Ct
      abftest$epoch$fEpochInitLevel
      abftest$epoch$fEpochLevelInc
      c2 = vector();
      for (i in sweep.max){
        
        n <-  i;
        n1L <- {n-1}*lengthofepisode+1
        n2L <- lengthofepisode*n
        tn <- abftest$traces[1,n1L:n2L];
        peak1list <- peakdet(tn, 20)
        example1  <- peak1list[[1]]
        if (length(example1[,1]) == 1){
          c2 <- i
          
        } else if (length(example1[,1]) > 1){
          c2<- i
          break
        }
        
        
      }
      creq <- c2-1;
      iv.v <-  vector();
      iv.i <-  vector();
      withProgress(message = 'Analyzing traces', value = 0,{
        n2 <- creq
        for (i in 1:creq){
          n <-  i;
          n1L <- {n-1}*lengthofepisode+1
          n2L <- lengthofepisode*n
          tn <- abftest$traces[1,n1L:n2L];
          lt <- length(tn);
          peakrmp <- tn[lt/2];
          peaki <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(i-1)
          iv.v<- append(iv.v,peakrmp)
          iv.i<- append(iv.i,peaki)
          
          incProgress(1/n2, detail = paste("Doing part", i))
          
        }})
      
      P <- polyfit(iv.i,iv.v,1)
      Rm <- P[1]*1000
      RMP <- P[2]
      
      Ct <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(c2-1)
      
      
      
      n3<- 5
      incProgress(1/n3, detail = paste("Halfway there!"))
      
      # Ctx2 Calculation
      
      Ctx2 <- Ct*2
      
      c2 <- vector();
      for (i in sweep.max){
        n <-  i;
        n1L <- {n-1}*lengthofepisode+1
        n2L <- lengthofepisode*n
        tn <- abftest$traces[1,n1L:n2L];
        peak1 <- tn
        #peak1list <- peakdet(as.vector(peak1, mode = "any"), 20)
        #example1  <- peak1list[[1]]
        CtM <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(i-1)
        if (CtM == Ctx2){
          c2 <- i
        } else if (CtM > Ctx2){
          c2<- i
          break
        }
      }
      
      n <-  c2-1;
      n1L <- {n-1}*lengthofepisode+1
      n2L <- lengthofepisode*n
      tn <- abftest$traces[1,n1L:n2L];
      peak2 <- tn
      peak2list <- peakdet(as.vector(peak2, mode = "any"), 20)
      apnumberx <- peak2list[[1]]
      apnumber <- length(apnumberx[,1]) 
      peakLoc <- apnumberx[,1]
      
      ##ISIRatio
      end <- length(peakLoc)
      sumisiCtx2 <- (peakLoc[end]-peakLoc[end-1])+ (peakLoc[end-1]-peakLoc[end-2])+ (peakLoc[end-2]-peakLoc[end-3]);
      ISIRatioCtx2 <- (peakLoc[2]-peakLoc[1])/( sumisiCtx2/3);
      
      
      
      ##Max Adaptation
      
      tnl<- length(tn)/2
      sumadaptCtx2 <- 1/((peakLoc[end]-peakLoc[end-1])/tnl)+ 1/((peakLoc[end-1]-peakLoc[end-2])/tnl)+ 1/((peakLoc[end-2]-peakLoc[end-3])/tnl);
      
      
      MaxAdaptationCtx2 <- (1/((peakLoc[2]-peakLoc[1])/tnl)- (sumadaptCtx2/3) );
      
      ##Adding in new IN PArams
      
      ##Cv2 and CV1
      
      nCV <- length(peakLoc);
      nCV1 <- nCV-1
      CV= vector();
      
      for (i in 1:nCV1){
        
        CV<- append(CV, (peakLoc[i+1]-peakLoc[i])/50)
        
      }
      
      CVn <- vector();
      n2CV <- length(CV);
      n2CV1 <- n2CV-1;
      
      for (i in 1:n2CV1){
        
        CVn <- append(CVn, 2*( abs((CV[i]-CV[i+1]))/((CV[i]+CV[i+1]))))
        
      }
      
      n3CV <- length(CVn);
      CVall_Ctx2 <- mean(CVn[1:n3CV]);
      CV1_Ctx2 <- mean(CVn[2:n3CV]);
      CV2_Ctx2 <- mean(CVn[3:n3CV]);
      
      ## AP Amp 
      
      
      peak1list <- peakdet(tn, 10)
      example1  <- peak1list[[1]]
      
      
      APAmp <- example1[1,2] - apthres
      APAmp2 <- example1[2,2] - apthres
      APAmp3 <- example1[3,2] - apthres
      APAmpEnd <- example1[end,2] - apthres
      
      AP_1_2_R_Ctx2 <- APAmp/APAmp2;
      AP_2_3_R_Ctx2 <- APAmp2/APAmp3;
      AP_end_3_R_Ctx2 <- APAmp/APAmpEnd;
      
      
      
      
      ###Max Ap
      
      apnumberlist <- vector();
      iv.i.entire <- vector();
      for (i in sweep.max){
        
        n <-  i;
        n1L <- {n-1}*lengthofepisode+1
        n2L <- lengthofepisode*n
        tn <- abftest$traces[1,n1L:n2L];
        peak1 <- tn
        peak1list <- peakdet(peak1, 20)
        example1  <- peak1list[[1]]
        apnumber2 <- length(example1[,1]) 
        apnumberlist<- append(apnumberlist,apnumber2)
        peaki2 <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(i-1)
        iv.i.entire<- append(iv.i.entire,peaki2)
      }  
      
      n3<- 7
      incProgress(1/n3, detail = paste("Almost!"))
      
      apnumberlist  
      iv.i.entire  
      
      #dfap <- data.frame(iv.i.entire,apnumberlist)
      #ggplot(dfap,aes(x=iv.i.entire, y=apnumberlist)) + geom_line()
      
      dyap <- diff(apnumberlist)/diff(iv.i.entire)
      
      maxaploc <- which.max(apnumberlist)
      maxapsweep <- maxaploc
      
      n <-  maxapsweep;
      n1L <- {n-1}*lengthofepisode+1
      n2L <- lengthofepisode*n
      tn <- abftest$traces[1,n1L:n2L];
      peak1 <- tn
      peak1list <- peakdet(peak1, 20)
      example1  <- peak1list[[1]]
      apnumbermax <- length(example1[,1]) 
      peakLoc <- example1[,1]
      #ggplot(dfap,aes(x=iv.i.entire, y=apnumberlist)) + geom_line()
      
      ##ISIRatio
      end <- length(peakLoc)
      sumisimax <- (peakLoc[end]-peakLoc[end-1])+ (peakLoc[end-1]-peakLoc[end-2])+ (peakLoc[end-2]-peakLoc[end-3]);
      ISIRatiomax <- (peakLoc[2]-peakLoc[1])/( sumisimax/3);
      
      
      
      ##Max Adaptation
      
      tnl<- length(tn)/2
      sumadaptmax <- 1/((peakLoc[end]-peakLoc[end-1])/tnl)+ 1/((peakLoc[end-1]-peakLoc[end-2])/tnl)+ 1/((peakLoc[end-2]-peakLoc[end-3])/tnl);
      
      
      MaxAdaptationmax <- (1/((peakLoc[2]-peakLoc[1])/tnl)- (sumadaptmax/3) );
      
      ##New Params MaxAPNo
      
      ##Cv2 and CV1
      
      nCV <- length(peakLoc);
      nCV1 <- nCV-1
      CV= vector();
      
      for (i in 1:nCV1){
        
        CV<- append(CV, (peakLoc[i+1]-peakLoc[i])/50)
        
      }
      
      CVn <- vector();
      n2CV <- length(CV);
      n2CV1 <- n2CV-1;
      
      for (i in 1:n2CV1){
        
        CVn <- append(CVn, 2*( abs((CV[i]-CV[i+1]))/((CV[i]+CV[i+1]))))
        
      }
      
      n3CV <- length(CVn);
      CVall_Max <- mean(CVn[1:n3CV]);
      CV1_Max <- mean(CVn[2:n3CV]);
      CV2_Max <- mean(CVn[3:n3CV]);
      
      ## AP Amp 
      
      
      peak1list <- peakdet(tn, 10)
      example1  <- peak1list[[1]]
      
      
      APAmp <- example1[1,2] - apthres
      APAmp2 <- example1[2,2] - apthres
      APAmp3 <- example1[3,2] - apthres
      APAmpEnd <- example1[end,2] - apthres
      
      AP_1_2_R_Max <- APAmp/APAmp2;
      AP_2_3_R_Max <- APAmp2/APAmp3;
      AP_end_3_R_Max <- APAmp/APAmpEnd;
      
      Value <- c(RMP, Rm, Ct, apnumber, Max_Fall, AHPAmp, ISIRatioCtx2, MaxAdaptationCtx2,CVall_Ctx2,CV1_Ctx2,CV2_Ctx2,AP_1_2_R_Ctx2,AP_2_3_R_Ctx2,AP_end_3_R_Ctx2)
      Property.Name <- c("RMP","Rm","Ct","AP Number at Ctx2","Max Fall Rate","AHP Amp","ISIRatio at Ctx2","Max Adaptation at Ctx2","Ctx2 CVAll","Ctx2 CV1","Ctx2 CV2","Ctx2 1/2AP","Ctx2 2/3AP","Ctx2 1/end AP")
      
      data1 <- data.frame(Property.Name, Value);
      
      n3<- 10
      incProgress(1/n3, detail = paste("Done!"))
    })
    Value
    }
  })
  
  output$extractedproperties <- DT::renderDataTable({
    req(input$file1)
    #req(ctx2check())
    if(length(tracenumber()) == 0){
      Value_Table <- 0
      Property.Name <- c("Please use a file with at least 11 action potentials")
      dataT <- data.frame(Property.Name, Value_Table);
      
      DT::datatable(dataT, extensions = 'Buttons', 
                    options = list(dom = 'Bfrtip',
                                   buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), pageLength = 25))  
      
      
    }
    
    else if (is.null(ctx2check())){
      
      Value_Table <- apcalc()
      Property.Name <- c("RMP (mV)","Rm (MOhm)","Ct (pA)","Max Fall Rate (mV/ms)","AHP Amp (mV)","ISI Ratio","Max Adapt (Hz)","CV2","CV2- 1st/2nd AP","AP Rise to Decay Ratio","CV2 of 1st AP","AP Threshold/Amplitude of AP","AP Rise/HalfWidth","1st/2nd AP","2nd/3rd AP","Max AP Number","Max Adaptation at Max AP No","ISIRatio at Max AP No","HalfWidth (ms)")
      dataT <- data.frame(Property.Name, Value_Table);
      
      DT::datatable(dataT, extensions = 'Buttons', 
                    options = list(dom = 'Bfrtip',
                                   buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), pageLength = 25))  
      
    }else {
    
    
    Value_Table <- apcalc()
    Property.Name <- c("RMP","Rm","Ct","AP Number at Ctx2","Max Fall Rate","AHP Amp","ISIRatio at Ctx2","Max Adaptation at Ctx2","Ctx2 CVAll","Ctx2 CV1","Ctx2 CV2","Ctx2 1/2AP","Ctx2 2/3AP","Ctx2 1/end AP")
    dataT <- data.frame(Property.Name, Value_Table);
    
    DT::datatable(dataT, extensions = 'Buttons', 
                  options = list(dom = 'Bfrtip',
                                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), pageLength = 50))
  }
    })
  
# Plot FI curve 
  output$plotfi <- renderPlotly({
    
    req(input$file1)
    
    apnumberlist <- vector();
    iv.i.entire <- vector();
    
    abftest <- abfloadv2(input$file1$datapath)
    numberoftraces <- abftest$episodes
    lengthofepisode <- length(abftest$traces[1,])/numberoftraces
    
    sweep.max <- seq(1,numberoftraces,1)
    withProgress(message = 'Ploting FI Curve', value = 0,{
      
    for (i in sweep.max){
      
      n <-  i;
      n1L <- {n-1}*lengthofepisode+1
      n2L <- lengthofepisode*n
      tn <- abftest$traces[1,n1L:n2L];
      peak1 <- tn
      peak1list <- peakdet(peak1, 20)
      example1  <- peak1list[[1]]
      apnumber2 <- length(example1[,1]) 
      apnumberlist<- append(apnumberlist,apnumber2)
      peaki2 <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(i-1)
      iv.i.entire<- append(iv.i.entire,peaki2)
      incProgress(1/n, detail = paste("Doing part", i))
      
      
    } })  
    apnumberlist  
    iv.i.entire  
    
    dfap <- data.frame(iv.i.entire,apnumberlist)
    plotfi <- ggplot(dfap,aes(x=iv.i.entire, y=apnumberlist)) + geom_line() + labs(x = "Current Injection",y = "AP Number")
    plotfi + labs(x = "Current Injection",y = "AP Number",title = "AP Number with Current Injection")
    plotfi <- plotly_build(plotfi, registerFrames = TRUE)
    plotfi$elementId <- NULL
    plotfi
    
  })
  
  
  output$plotif <- renderPlotly({
    
    req(input$file1)
    
    apnumberlist <- vector();
    iv.i.entire <- vector();
    
    abftest <- abfloadv2(input$file1$datapath)
    numberoftraces <- abftest$episodes
    lengthofepisode <- length(abftest$traces[1,])/numberoftraces
    sweep.max <- seq(1,numberoftraces,1)
    
    
    n <-  input$obs;
    if (n > numberoftraces){
      
      text = paste("\n   There is no trace to display.\n",
                   "       Please try a step protocol with more current injections\n")
      ggplot() + 
        annotate("text", x = 4, y = 25, size=4, label = text) + 
        theme_bw() +
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.position="none",
              panel.background=element_blank(),
              panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank()) + labs(x=NULL, y=NULL)
      
    } else {
    n1L <- {n-1}*lengthofepisode+1
    n2L <- lengthofepisode*n
    tn <- abftest$traces[1,n1L:n2L];
    timeT<- abftest$s
    timeTn <- timeT[n1L:n2L]
    #df.abf <- data.frame(timeTn,tn)
    #df <- data.frame(abftest1$dataX,abftest1$dataY)
    #ggplot(df,aes(x=abftest$dataX, y=abftest$dataY)) + geom_line()
    peak1 <- tn;
    time1 <-  timeTn*1000;
    dy <- diff(peak1)/diff(time1);
    #df1 <- data.frame(time1[2:100000],dy);
    #ggplot(df1,aes(x=time1[1:99999], y=dy)) + geom_line()
    peakMag2 <- peakdet(dy,10)
    
    peak1list <- peakdet(peak1, 10)
    example1  <- peak1list[[1]]
    peakMag2x <- peakMag2[[1]]
    
    peakLoc <- example1[[1]]
    
    tnl <- length(tn)/2
    listif_1 <- c(1/(peakLoc[2]-peakLoc[1])*tnl , 1/(peakLoc[3]-peakLoc[2])*tnl,1/(peakLoc[4]-peakLoc[3])*tnl,1/(peakLoc[5]-peakLoc[4])*tnl,1/(peakLoc[6]-peakLoc[5])*tnl,1/(peakLoc[7]-peakLoc[6])*tnl,1/(peakLoc[8]-peakLoc[7])*tnl,1/(peakLoc[9]-peakLoc[8])*tnl,1/(peakLoc[10]-peakLoc[9])*tnl)
    listif_2 <- c((peakLoc[2]-peakLoc[1])/tnl,(peakLoc[3]-peakLoc[1])/tnl,(peakLoc[4]-peakLoc[1])/tnl,(peakLoc[5]-peakLoc[1])/tnl,(peakLoc[6]-peakLoc[1])/tnl,(peakLoc[7]-peakLoc[1])/tnl,(peakLoc[8]-peakLoc[1])/tnl,(peakLoc[9]-peakLoc[1])/tnl,(peakLoc[10]-peakLoc[1])/tnl)
    
    
    
    dfap <- data.frame(listif_2,listif_1)
    plotif <- ggplot(dfap,aes(x=listif_2, y=listif_1)) + geom_line() + stat_smooth(se= FALSE) + geom_point()+ labs(x = "Time (sec)",y = "Frequency (Hz)")
    plotif + labs(x = "Time (sec)",y = "Frequency (Hz)", title = "Instantaneous Frequency curve")
    plotif <- plotly_build(plotif, registerFrames = TRUE)
    plotif$elementId <- NULL
    plotif
    }
    
  })
  

  


  
# Get user input data from table in manual classifier
  
  data <- reactive({
    
    if (is.null(input$hot)) {
      Property.Name <- as.character(proplabels[2:15,1])
      Value <- data.martin.NN.Level3[1,2:15];
      
      DF1 = data.frame(data.frame(Property.Name, t(Value)))
      row.names(DF1) <- c(1:14)
      DF1
    } else {
      DF1 = hot_to_r(input$hot)
    }
    DF1
  
  })
  
# Manual classification functions
  
  output$hot <- renderRHandsontable({
    DF1 <- data()
    rhandsontable(DF1, useTypes = FALSE, selectCallback = TRUE)
  })
 
  output$hot3 <- renderRHandsontable({
    DF2 <- dataL2()
    rhandsontable(DF2, useTypes = FALSE, selectCallback = TRUE)
  })
  
  # Level 1 - Auto
  
  output$classification <- renderText({ 
    req(input$file1)
    #NewDat <- data()
    #NewDat <- t(data.martin.NN.Level1[1,2:39])
    
    ValueC <- as.numeric(apcalc());
    Property.Name <- as.character(proplabels[2:15,1])
    
    data2.scale <- ValueC
    #data2.scale <- NewDat
    trial <- t(as.data.frame(data2.scale))
    colnames(trial) <-  as.character(proplabels[2:15,1])
 
    pr.nn_ <- predict(newnetPCALevel1, trial, type = "prob")
    
    cla_result1 <- which.max(pr.nn_)
  
    if (names(cla_result1) == "PN"){
      
      pr.nn_ <- predict(newnetPCALevel2, trial, type = "prob")
      cla_result2 <- which.max(pr.nn_)
      if(names(cla_result2) == "CS"){
        
        paste("The cell is a putative claustro-subcortical non-adapting projection neuron with confidence",as.character(signif(max(pr.nn_),3)))
        
        }else if(names(cla_result2) == "CC1"){
          
          paste("The cell is a putative type 1 claustro-cortical strongly adapting projection neuron with confidence",as.character(signif(max(pr.nn_),3)))
          
          
        }
        else if(names(cla_result2) == "CC2"){
          
          paste("The cell is a putative type 2 claustro-cortical strongly adapting projection neuron with confidence",as.character(signif(max(pr.nn_),3)))
          
          
        }
        else if(names(cla_result2) == "CC3"){
          
          paste("The cell is a putative type 3 claustro-cortical strongly adapting projection neuron with confidence",as.character(signif(max(pr.nn_),3)))
          
          
        }
        else if(names(cla_result2) == "CC4"){
          
          paste("The cell is a putative type 4 claustro-cortical strongly adapting projection neuron with confidence",as.character(signif(max(pr.nn_),3)))
          
          
        }
        
        
        
      }
      
      
     else if(names(cla_result1) == "IN"){
      
       pr.nn_ <- predict(newnetPCALevel3, trial, type = "prob")
       cla_result3 <- which.max(pr.nn_)
       
       if(names(cla_result3) == "PV"){
         
         paste("The cell is a putative parvalbumin expressing (PV) interneuron with confidence",as.character(signif(max(pr.nn_),3)))
         
       }else if(names(cla_result3) == "SST"){
         
         paste("The cell is a putative somatostatin expressing (SST) interneuron with confidence", as.character(signif(max(pr.nn_),3)))
         
         
       }
       else if(names(cla_result3) == "VIP"){
         
         paste("The cell is a putative vasointestinal peptide expressing (VIP) interneuron with confidence", as.character(signif(max(pr.nn_),3)))
         
         
       }
       
      
    } 
    
  })
  
  output$classficiationaccuracy <- DT::renderDataTable({
    req(input$file1)
    ValueC <- apcalc();
    Property.Name <- as.character(proplabels[2:15,1])
    Value_2 <- ValueC;
    data2.scale <- Value_2
    trial <- t(as.data.frame(data2.scale))
    colnames(trial) <- as.character(proplabels[2:15,1])
    
    pr.nn_ <- predict(newnetPCALevel1, trial, type = "prob")
    
    cla_result1 <- which.max(pr.nn_)
    
    if (names(cla_result1) == "PN"){
      
      pr.nn_ <- predict(newnetPCALevel2, trial, type = "prob")
      cla_result2 <- which.max(pr.nn_)
      if(names(cla_result2) == "CS"){
        
        
        ConfidenceInterval <- c(pr.nn_[1],pr.nn_[2],pr.nn_[3],pr.nn_[4],pr.nn_[5])
        
        dataT2 <- as.data.frame(ConfidenceInterval)
        
        DT::datatable(dataT2, extensions = 'Buttons', 
                      options = list(dom = 'Bfrtip',
                                     buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), pageLength = 25))
        
      }else if(names(cla_result2) == "CC1"){
        
        
        ConfidenceInterval <- c(pr.nn_[1],pr.nn_[2],pr.nn_[3],pr.nn_[4],pr.nn_[5])
        
        dataT2 <- as.data.frame(ConfidenceInterval)
        
        DT::datatable(dataT2, extensions = 'Buttons', 
                      options = list(dom = 'Bfrtip',
                                     buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), pageLength = 25))
        
        
      }
      else if(names(cla_result2) == "CC2"){
        
        
        ConfidenceInterval <- c(pr.nn_[1],pr.nn_[2],pr.nn_[3],pr.nn_[4],pr.nn_[5])
        
        dataT2 <- as.data.frame(ConfidenceInterval)
        
        DT::datatable(dataT2, extensions = 'Buttons', 
                      options = list(dom = 'Bfrtip',
                                     buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), pageLength = 25))
        
        
      }
      else if(names(cla_result2) == "CC3"){
        
        
        ConfidenceInterval <- c(pr.nn_[1],pr.nn_[2],pr.nn_[3],pr.nn_[4],pr.nn_[5])
        
        dataT2 <- as.data.frame(ConfidenceInterval)
        
        DT::datatable(dataT2, extensions = 'Buttons', 
                      options = list(dom = 'Bfrtip',
                                     buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), pageLength = 25))
        
        
      }
      else if(names(cla_result2) == "CC4"){
        
        
        ConfidenceInterval <- c(pr.nn_[1],pr.nn_[2],pr.nn_[3],pr.nn_[4],pr.nn_[5])
        
        dataT2 <- as.data.frame(ConfidenceInterval)
        
        DT::datatable(dataT2, extensions = 'Buttons', 
                      options = list(dom = 'Bfrtip',
                                     buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), pageLength = 25))
        
        
      }
      
      
      
    }
    
    
    else if(names(cla_result1) == "IN"){
      
      pr.nn_ <- predict(newnetPCALevel3, trial, type = "prob")
      cla_result3 <- which.max(pr.nn_)
      
      if(names(cla_result3) == "PV"){
        
        
        ConfidenceInterval <- c(pr.nn_[1],pr.nn_[2],pr.nn_[3])
        
        dataT2 <- as.data.frame(ConfidenceInterval)
        
        DT::datatable(dataT2, extensions = 'Buttons', 
                      options = list(dom = 'Bfrtip',
                                     buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), pageLength = 25))
        
      }else if(names(cla_result3) == "SST"){
        
        
        ConfidenceInterval <- c(pr.nn_[1],pr.nn_[2],pr.nn_[3])
        
        dataT2 <- as.data.frame(ConfidenceInterval)
        
        DT::datatable(dataT2, extensions = 'Buttons', 
                      options = list(dom = 'Bfrtip',
                                     buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), pageLength = 25))
        
        
      }
      else if(names(cla_result3) == "VIP"){
        
        
        ConfidenceInterval <- c(pr.nn_[1],pr.nn_[2],pr.nn_[3])
        
        dataT2 <- as.data.frame(ConfidenceInterval)
        
        DT::datatable(dataT2, extensions = 'Buttons', 
                      options = list(dom = 'Bfrtip',
                                     buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), pageLength = 25))
        
        
      }
      
      
    } 
    
    
    

    
    
  })
  
  output$plotLIME <- renderPlot({
    
    req(input$file1)
    ValueC <- apcalc();
    
  
    prediction.lime <- ValueC
    trial <- t(as.data.frame(prediction.lime))
    colnames(trial) <- as.character(proplabels[2:15,1])
    
    pr.nn_ <- predict(newnetPCALevel1, trial, type = "prob")

    cla_result1 <- which.max(pr.nn_)

    
    if (names(cla_result1) == "PN"){
      
      explanationtest <- explain(as.data.frame(trial), explainerL2, n_labels = 1, n_features = 5)
      plot_features(explanationtest, ncol = 2)
      
    } else if(names(cla_result1) == "IN"){
      
      explanationtest <- explain(as.data.frame(trial), explainerL3, n_labels = 1, n_features = 5)
      plot_features(explanationtest, ncol = 2)
      
      
    } 
     
    
  })
  
  # Level 2 - Manual

  output$classificationmanual <- renderText({ 
    req(input$hot)
    NewDat <- data()
    ValueC <- NewDat[,2];
    
   
    Property.Name <- as.character(proplabels[2:15,1])
    
    data2.scale <- ValueC
    #data2.scale <- NewDat
    trial <- t(as.data.frame(data2.scale))
    colnames(trial) <-  as.character(proplabels[2:15,1])
    
    pr.nn_ <- predict(newnetPCALevel1, trial, type = "prob")
    
    cla_result1 <- which.max(pr.nn_)
    
    if (names(cla_result1) == "PN"){
      
      pr.nn_ <- predict(newnetPCALevel2, trial, type = "prob")
      cla_result2 <- which.max(pr.nn_)
      if(names(cla_result2) == "CS"){
        
        paste("The cell is a putative claustro-subcortical non-adapting projection neuron with confidence",as.character(signif(max(pr.nn_),3)))
        
      }else if(names(cla_result2) == "CC1"){
        
        paste("The cell is a putative type 1 claustro-cortical strongly adapting projection neuron with confidence",as.character(signif(max(pr.nn_),3)))
        
        
      }
      else if(names(cla_result2) == "CC2"){
        
        paste("The cell is a putative type 2 claustro-cortical strongly adapting projection neuron with confidence",as.character(signif(max(pr.nn_),3)))
        
        
      }
      else if(names(cla_result2) == "CC3"){
        
        paste("The cell is a putative type 3 claustro-cortical strongly adapting projection neuron with confidence",as.character(signif(max(pr.nn_),3)))
        
        
      }
      else if(names(cla_result2) == "CC4"){
        
        paste("The cell is a putative type 4 claustro-cortical strongly adapting projection neuron with confidence",as.character(signif(max(pr.nn_),3)))
        
        
      }
      
      
      
    }
    
    
    else if(names(cla_result1) == "IN"){
      
      pr.nn_ <- predict(newnetPCALevel3, trial, type = "prob")
      cla_result3 <- which.max(pr.nn_)
      
      if(names(cla_result3) == "PV"){
        
        paste("The cell is a putative parvalbumin expressing (PV) interneuron with confidence",as.character(signif(max(pr.nn_),3)))
        
      }else if(names(cla_result3) == "SST"){
        
        paste("The cell is a putative somatostatin expressing (SST) interneuron with confidence", as.character(signif(max(pr.nn_),3)))
        
        
      }
      else if(names(cla_result3) == "VIP"){
        
        paste("The cell is a putative vasointestinal peptide expressing (VIP) interneuron with confidence", as.character(signif(max(pr.nn_),3)))
        
        
      }
      
      
    } 
    
  })
  
  output$plotLIMEmanual <- renderPlot({
    
    req(input$hot)
    NewDat <- data()
    ValueC <- NewDat[,2];
    
    
    prediction.lime <- ValueC
    trial <- t(as.data.frame(prediction.lime))
    colnames(trial) <- as.character(proplabels[2:15,1])
    
    pr.nn_ <- predict(newnetPCALevel1, trial, type = "prob")
    
    cla_result1 <- which.max(pr.nn_)
    
    
    if (names(cla_result1) == "PN"){
      
      explanationtest <- explain(as.data.frame(trial), explainerL2, n_labels = 1, n_features = 5)
      plot_features(explanationtest, ncol = 2)
      
    } else if(names(cla_result1) == "IN"){
      
      explanationtest <- explain(as.data.frame(trial), explainerL3, n_labels = 1, n_features = 5)
      plot_features(explanationtest, ncol = 2)
      
      
    } 
    
    
  })
  
  
# Batch property classification functions
  
  batch <- reactive({
    
    req(input$filebatch)
    lengthoffiles <- length(input$filebatch$datapath)
    #testdrataframe = data.frame(matrix(vector(), 0, 3, dimnames=list(c(), c("Date", "File", "User"))),stringsAsFactors=F)
    n<- lengthoffiles
    testdrataframe = vector()
    
    for (im in 1:lengthoffiles){
      
      abftest <- abfloadv2(input$filebatch$datapath[im])
      numberoftraces <- abftest$episodes
      lengthofepisode <- length(abftest$traces[1,])/numberoftraces
      sweep.max <- seq(1,numberoftraces,1)
      
      if (input$tracesel == "ap11"){
        c1 = vector()
        abftest <- abfloadv2(input$filebatch$datapath[im])
        numberoftraces <- abftest$episodes
        lengthofepisode <- length(abftest$traces[1,])/numberoftraces
        sweep.max <- seq(1,numberoftraces,1)
        
        withProgress(message = 'Finding the right trace to use', value = 0,{
          n <- length(sweep.max)
          for (i in sweep.max){
            n <-  i;
            n1L <- {n-1}*lengthofepisode+1
            n2L <- lengthofepisode*n
            tn <- abftest$traces[1,n1L:n2L];
            peak1list <- peakdet(tn, 20)
            example1  <- peak1list[[1]]
            if (length(example1[,1]) == 11){
              c1 <- i
              break
            } else if (length(example1[,1]) > 11){
              c1<- i
              break
            }
            incProgress(1/n, detail = paste("Doing part", i))
            
          }
          
        })
        
        tracereq <- c1
        
      } else if (input$tracesel == "ct2"){
        
        abftest <- abfloadv2(input$filebatch$datapath[im])
        numberoftraces <- abftest$episodes
        lengthofepisode <- length(abftest$traces[1,])/numberoftraces
        sweep.max <- seq(1,numberoftraces,1)
        
        abftest$epoch$fEpochInitLevel
        abftest$epoch$fEpochLevelInc
        c2 = vector();
        for (i in sweep.max){
          
          n <-  i;
          n1L <- {n-1}*lengthofepisode+1
          n2L <- lengthofepisode*n
          tn <- abftest$traces[1,n1L:n2L];
          peak1list <- peakdet(tn, 20)
          example1  <- peak1list[[1]]
          if (length(example1[,1]) == 1){
            c2 <- i
            
          } else if (length(example1[,1]) > 1){
            c2<- i
            break
          }
          
          
        }
        creq <- c2-1;
        
        
        Ct <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(c2-1)
        
        # Ctx2
        Ctx2 <- Ct*2;
        Ctx2
        c2 <- vector();
        
        ct_times_two_c = (Ctx2 - abftest$epoch$fEpochInitLeve)/abftest$epoch$fEpochLevelInc
        tracereq <- ct_times_two_c +1
        
        
      } else if (input$tracesel == "ct1"){
        abftest <- abfloadv2(input$filebatch$datapath[im])
        numberoftraces <- abftest$episodes
        lengthofepisode <- length(abftest$traces[1,])/numberoftraces
        sweep.max <- seq(1,numberoftraces,1)
        
        abftest$epoch$fEpochInitLevel
        abftest$epoch$fEpochLevelInc
        c2 = vector();
        for (i in sweep.max){
          
          n <-  i;
          n1L <- {n-1}*lengthofepisode+1
          n2L <- lengthofepisode*n
          tn <- abftest$traces[1,n1L:n2L];
          peak1list <- peakdet(tn, 20)
          example1  <- peak1list[[1]]
          if (length(example1[,1]) == 1){
            c2 <- i
            
          } else if (length(example1[,1]) > 1){
            c2<- i
            break
          }
          
          
        }
        creq <- c2-1;
        ct_plus_one <- c2+1
        
        Ct <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(c2-1)
        
        # Ctx2
        Ctx2 <- Ct*2;
        Ct_1 <- Ct + abftest$epoch$fEpochLevelInc;
        Ct_1
        c2 <- vector();
        
        tracereq <- ct_plus_one
        
        
        
        
        
      }
      
      abftest <- abfloadv2(input$filebatch$datapath[im])
      numberoftraces <- abftest$episodes
      lengthofepisode <- length(abftest$traces[1,])/numberoftraces
      sweep.max <- seq(1,numberoftraces,1)
        
      abftest$epoch$fEpochInitLevel
      abftest$epoch$fEpochLevelInc
      c2 = vector();
      for (i in sweep.max){
          
          n <-  i;
          n1L <- {n-1}*lengthofepisode+1
          n2L <- lengthofepisode*n
          tn <- abftest$traces[1,n1L:n2L];
          peak1list <- peakdet(tn, 20)
          example1  <- peak1list[[1]]
          if (length(example1[,1]) == 1){
            c2 <- i
            
          } else if (length(example1[,1]) > 1){
            c2<- i
            break
          }
          
          
        }
        
        creq <- c2-1;
        Ct <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(c2-1)
        
        # Ctx2
        Ctx2 <- Ct*2;
        Ctx2
        c2 <- vector();
        
        ct_times_two_c = (Ctx2 - abftest$epoch$fEpochInitLeve)/abftest$epoch$fEpochLevelInc
        ct_times_two = ct_times_two_c +1
        
        if (ct_times_two > numberoftraces ){
          
          
          withProgress(message = 'Trace Analysis Begins', value = 0,{
            
            c1 = vector()
            abftest <- abfloadv2(input$filebatch$datapath[im])
            numberoftraces <- abftest$episodes
            lengthofepisode <- length(abftest$traces[1,])/numberoftraces
            sweep.max <- seq(1,numberoftraces,1)
            
            
            n <-  tracereq;
            n1L <- {n-1}*lengthofepisode+1
            n2L <- lengthofepisode*n
            tn <- abftest$traces[1,n1L:n2L];
            timeT<- abftest$s
            timeTn <- timeT[n1L:n2L]
            #df.abf <- data.frame(timeTn,tn)
            #df <- data.frame(abftest1$dataX,abftest1$dataY)
            #ggplot(df,aes(x=abftest$dataX, y=abftest$dataY)) + geom_line()
            peak1 <- tn;
            time1 <-  timeTn*1000;
            dy <- diff(peak1)/diff(time1);
            #df1 <- data.frame(time1[2:100000],dy);
            #ggplot(df1,aes(x=time1[1:99999], y=dy)) + geom_line()
            peakMag2 <- peakdet(dy,10)
            
            peakMag2x <- peakMag2[[1]]
            
            locfx <- which.max(peakMag2x[,2])
            locfx2 <- peakMag2x[locfx,1]
            
            ds <- dy[1:locfx2]
            dsf <- rev(ds)
            
            tn2<- vector();
            
            for (n in 1:locfx2){
              if (dsf[n] < 10){
                tn2 =n
                break
              }
              
            }
            
            
            
            ds2 <- peak1[1:locfx2];
            ds3 <- rev(ds2);
            apthres <- ds3[tn2]
            
            n <-  tracereq;
            
            n1L <- {n-1}*lengthofepisode+1
            n2L <- lengthofepisode*n
            tn <- abftest$traces[1,n1L:n2L];
            peak1list <- peakdet(tn, 10)
            example1  <- peak1list[[1]]
            
            df.abf <- data.frame(timeTn,tn)
            #df <- data.frame(abftest1$dataX,abftest1$dataY)
            #ggplot(df.abf,aes(x=timeTn, y=tn)) + geom_line()
            
            
            
            APAmp <- example1[1,2] - apthres
            APAmp2 <- example1[2,2] - apthres
            APAmp3 <- example1[3,2] - apthres
            
            AP_1_2_R <- APAmp/APAmp2;
            AP_2_3_R <- APAmp2/APAmp3;
            
            APThres_APAmp <- apthres/APAmp;
            
            #APRise and Fall
            #Rise
            tnMax <- tn[1:example1[1,1]];
            timeTnMax <- timeTn[1:example1[1,1]]
            peaktnMax <- tnMax;
            time1tnMax <-  timeTnMax*1000;
            dytnMax <- diff(peaktnMax)/diff(time1tnMax);
            peakMag2Max <- peakdet(dytnMax,10)
            #df <- data.frame(timeTnMax,tnMax)
            #ggplot(df,aes(x=timeTnMax, y=tnMax)) + geom_line()
            
            peakMag2x <- peakMag2Max[[1]]
            peakMag3x <- peakMag2x[[2]]
            Max_Rise <- (max(peakMag3x[]))
            
            
            #Fall
            tnFall <- tn[example1[1,1]:example1[2,1]];
            timeTnFall <- timeTn[example1[1,1]:example1[2,1]]
            peaktnFall <- tnFall;
            time1tnFall <-  timeTnFall*1000;
            dytnFall <- diff(peaktnFall)/diff(time1tnFall);
            peakMag2Fall <- peakdet(dytnFall,10)
            #df <- data.frame(timeTnFall,tnFall)
            #ggplot(df,aes(x=timeTnFall, y=tnFall)) + geom_line()
            
            peakMag2xF <- peakMag2Fall[[2]]
            peakMag3xF <- peakMag2xF[[2]]
            Max_Fall <- (min(peakMag3xF[]))
            
            
            peakLoc <- example1[[1]]
            
            
            
            #AHP
            peakAHP <-peakdet(tnFall,10)
            peakAHP2 <- peakAHP[[2]]
            peakAHP3 <-  peakAHP2[[2]]
            yAHP <- min(peakAHP3[])
            
            AHPAmp <- yAHP - apthres;
            
            AP_Rise_Decay_R <- -Max_Rise/Max_Fall;
            
            #HalfWidth
            
            half1 <- apthres + APAmp/2;
            yp1x <- peakLoc[1]-200
            yp1 <- peak1[yp1x:peakLoc[1]]
            tt1 <- which.min(abs(yp1 - half1)) 
            time <- timeTn
            
            timett1 <- time[yp1x:peakLoc[1]]
            
            tt1x <- timett1[tt1]
            yp2x <- peakLoc[1]+200
            yp2 <- peak1[peakLoc[1]:yp2x]
            
            tt2 <- which.min(abs(yp2 - half1)) 
            
            timett2 <- time[peakLoc[1]:yp2x]
            tt2x <- timett2[tt2]
            
            halfwidth <- (tt2x - tt1x)*1000
            
            AP_Rise_Half <- Max_Rise/halfwidth
            
            
            ##ISIRatio
            end <- length(peakLoc)
            sumisi <- (peakLoc[end]-peakLoc[end-1])+ (peakLoc[end-1]-peakLoc[end-2])+ (peakLoc[end-2]-peakLoc[end-3]);
            ISIRatio <- (peakLoc[2]-peakLoc[1])/( sumisi/3);
            
            n3<- 3
            incProgress(1/n3, detail = paste("Starting Analysis"))
            
            ##Max Adaptation
            tnl<- length(tn)/2
            
            sumadapt <- 1/((peakLoc[end]-peakLoc[end-1])/tnl)+ 1/((peakLoc[end-1]-peakLoc[end-2])/tnl)+ 1/((peakLoc[end-2]-peakLoc[end-3])/tnl);
            
            
            MaxAdaptation <- (1/((peakLoc[2]-peakLoc[1])/tnl)- (sumadapt/3) );
            
            ##Cv2 and CV1
            
            nCV <- length(peakLoc);
            nCV1 <- nCV-1
            CV= vector();
            
            for (i in 1:nCV1){
              
              CV<- append(CV, (peakLoc[i+1]-peakLoc[i])/50)
              
            }
            
            CVn <- vector();
            n2CV <- length(CV);
            n2CV1 <- n2CV-1;
            
            for (i in 1:n2CV1){
              
              CVn <- append(CVn, 2*( abs((CV[i]-CV[i+1]))/((CV[i]+CV[i+1]))))
              
            }
            
            n3CV <- length(CVn);
            CVall <- mean(CVn[1:n3CV]);
            CV1 <- mean(CVn[2:n3CV]);
            CV2 <- mean(CVn[3:n3CV]);
            
            #Rm, RMP, Ct
            abftest$epoch$fEpochInitLevel
            abftest$epoch$fEpochLevelInc
            c2 = vector();
            for (i in sweep.max){
              
              n <-  i;
              n1L <- {n-1}*lengthofepisode+1
              n2L <- lengthofepisode*n
              tn <- abftest$traces[1,n1L:n2L];
              peak1list <- peakdet(tn, 20)
              example1  <- peak1list[[1]]
              if (length(example1[,1]) == 1){
                c2 <- i
                
              } else if (length(example1[,1]) > 1){
                c2<- i
                break
              }
              
              
            }
            creq <- c2-1;
            iv.v <-  vector();
            iv.i <-  vector();
            withProgress(message = 'Analyzing traces', value = 0,{
              n2 <- creq
              for (i in 1:creq){
                n <-  i;
                n1L <- {n-1}*lengthofepisode+1
                n2L <- lengthofepisode*n
                tn <- abftest$traces[1,n1L:n2L];
                lt <- length(tn);
                peakrmp <- tn[lt/2];
                peaki <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(i-1)
                iv.v<- append(iv.v,peakrmp)
                iv.i<- append(iv.i,peaki)
                
                incProgress(1/n2, detail = paste("Doing part", i))
                
              }})
            
            P <- polyfit(iv.i,iv.v,1)
            Rm <- P[1]*1000
            RMP <- P[2]
            
            Ct <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(c2-1)
            
            
            
            n3<- 5
            incProgress(1/n3, detail = paste("Halfway there!"))
            
            
            ###Max Ap
            
            apnumberlist <- vector();
            iv.i.entire <- vector();
            for (i in sweep.max){
              
              n <-  i;
              n1L <- {n-1}*lengthofepisode+1
              n2L <- lengthofepisode*n
              tn <- abftest$traces[1,n1L:n2L];
              peak1 <- tn
              peak1list <- peakdet(peak1, 20)
              example1  <- peak1list[[1]]
              apnumber2 <- length(example1[,1]) 
              apnumberlist<- append(apnumberlist,apnumber2)
              peaki2 <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(i-1)
              iv.i.entire<- append(iv.i.entire,peaki2)
            }  
            
            n3<- 7
            incProgress(1/n3, detail = paste("Almost!"))
            
            apnumberlist  
            iv.i.entire  
            
            #dfap <- data.frame(iv.i.entire,apnumberlist)
            #ggplot(dfap,aes(x=iv.i.entire, y=apnumberlist)) + geom_line()
            
            dyap <- diff(apnumberlist)/diff(iv.i.entire)
            
            maxaploc <- which.max(apnumberlist)
            maxapsweep <- maxaploc
            
            n <-  maxapsweep;
            n1L <- {n-1}*lengthofepisode+1
            n2L <- lengthofepisode*n
            tn <- abftest$traces[1,n1L:n2L];
            peak1 <- tn
            peak1list <- peakdet(peak1, 20)
            example1  <- peak1list[[1]]
            apnumbermax <- length(example1[,1]) 
            peakLoc <- example1[,1]
            #ggplot(dfap,aes(x=iv.i.entire, y=apnumberlist)) + geom_line()
            
            ##ISIRatio
            end <- length(peakLoc)
            sumisimax <- (peakLoc[end]-peakLoc[end-1])+ (peakLoc[end-1]-peakLoc[end-2])+ (peakLoc[end-2]-peakLoc[end-3]);
            ISIRatiomax <- (peakLoc[2]-peakLoc[1])/( sumisimax/3);
            
            
            
            ##Max Adaptation
            
            tnl<- length(tn)/2
            sumadaptmax <- 1/((peakLoc[end]-peakLoc[end-1])/tnl)+ 1/((peakLoc[end-1]-peakLoc[end-2])/tnl)+ 1/((peakLoc[end-2]-peakLoc[end-3])/tnl);
            
            
            MaxAdaptationmax <- (1/((peakLoc[2]-peakLoc[1])/tnl)- (sumadaptmax/3) );
            
            
            
            
            Value <- c(as.character(input$filebatch$name[im]),RMP, Rm, Ct, Max_Fall, AHPAmp,ISIRatio,MaxAdaptation,CVall,CV2,AP_Rise_Decay_R,CV1,APThres_APAmp,AP_Rise_Half,AP_1_2_R,AP_2_3_R, apnumbermax,MaxAdaptationmax,ISIRatiomax,halfwidth)
            #Property.Name <- c("RMP","Rm","Ct","Max Fall Rate","AHP Amp","ISI Ratio","Max Adapt","CV2","CV2- 1st/2nd AP","AP Rise Decay R","CV2 1st AP","AP Thres/Amp","AP Rise/HalfWidth","1st/2nd AP","2nd/3rd AP","AP Number at Ctx2","Max Adaptation at Ctx2","ISIRatio at Ctx2","Max AP Number","Max Adaptation at Max AP No","ISIRatio at Max AP No", "Halfwidth")
            
            #data1 <- data.frame(Property.Name, Value);
            
            n3<- 10
            incProgress(1/n3, detail = paste("Done!"))
          })
          Value
          testdrataframe = Cbind(testdrataframe,Value)
           
          
          
          
          
          
        } else{
          
          
          
          withProgress(message = 'Trace Analysis Begins', value = 0,{
            
            c1 = vector()
            abftest <- abfloadv2(input$filebatch$datapath[im])
            numberoftraces <- abftest$episodes
            lengthofepisode <- length(abftest$traces[1,])/numberoftraces
            sweep.max <- seq(1,numberoftraces,1)
            
            
            n <-  tracereq;
            n1L <- {n-1}*lengthofepisode+1
            n2L <- lengthofepisode*n
            tn <- abftest$traces[1,n1L:n2L];
            timeT<- abftest$s
            timeTn <- timeT[n1L:n2L]
            #df.abf <- data.frame(timeTn,tn)
            #df <- data.frame(abftest1$dataX,abftest1$dataY)
            #ggplot(df,aes(x=abftest$dataX, y=abftest$dataY)) + geom_line()
            peak1 <- tn;
            time1 <-  timeTn*1000;
            dy <- diff(peak1)/diff(time1);
            #df1 <- data.frame(time1[2:100000],dy);
            #ggplot(df1,aes(x=time1[1:99999], y=dy)) + geom_line()
            peakMag2 <- peakdet(dy,10)
            
            peakMag2x <- peakMag2[[1]]
            
            locfx <- which.max(peakMag2x[,2])
            locfx2 <- peakMag2x[locfx,1]
            
            ds <- dy[1:locfx2]
            dsf <- rev(ds)
            
            tn2<- vector();
            
            for (n in 1:locfx2){
              if (dsf[n] < 10){
                tn2 =n
                break
              }
              
            }
            
            
            
            ds2 <- peak1[1:locfx2];
            ds3 <- rev(ds2);
            apthres <- ds3[tn2]
            
            n <-  tracereq;
            
            n1L <- {n-1}*lengthofepisode+1
            n2L <- lengthofepisode*n
            tn <- abftest$traces[1,n1L:n2L];
            peak1list <- peakdet(tn, 10)
            example1  <- peak1list[[1]]
            
            df.abf <- data.frame(timeTn,tn)
            #df <- data.frame(abftest1$dataX,abftest1$dataY)
            #ggplot(df.abf,aes(x=timeTn, y=tn)) + geom_line()
            
            
            
            APAmp <- example1[1,2] - apthres
            APAmp2 <- example1[2,2] - apthres
            APAmp3 <- example1[3,2] - apthres
            
            AP_1_2_R <- APAmp/APAmp2;
            AP_2_3_R <- APAmp2/APAmp3;
            
            APThres_APAmp <- apthres/APAmp;
            
            #APRise and Fall
            #Rise
            tnMax <- tn[1:example1[1,1]];
            timeTnMax <- timeTn[1:example1[1,1]]
            peaktnMax <- tnMax;
            time1tnMax <-  timeTnMax*1000;
            dytnMax <- diff(peaktnMax)/diff(time1tnMax);
            peakMag2Max <- peakdet(dytnMax,10)
            #df <- data.frame(timeTnMax,tnMax)
            #ggplot(df,aes(x=timeTnMax, y=tnMax)) + geom_line()
            
            peakMag2x <- peakMag2Max[[1]]
            peakMag3x <- peakMag2x[[2]]
            Max_Rise <- (max(peakMag3x[]))
            
            
            #Fall
            tnFall <- tn[example1[1,1]:example1[2,1]];
            timeTnFall <- timeTn[example1[1,1]:example1[2,1]]
            peaktnFall <- tnFall;
            time1tnFall <-  timeTnFall*1000;
            dytnFall <- diff(peaktnFall)/diff(time1tnFall);
            peakMag2Fall <- peakdet(dytnFall,10)
            #df <- data.frame(timeTnFall,tnFall)
            #ggplot(df,aes(x=timeTnFall, y=tnFall)) + geom_line()
            
            peakMag2xF <- peakMag2Fall[[2]]
            peakMag3xF <- peakMag2xF[[2]]
            Max_Fall <- (min(peakMag3xF[]))
            
            
            peakLoc <- example1[[1]]
            
            
            
            #AHP
            peakAHP <-peakdet(tnFall,10)
            peakAHP2 <- peakAHP[[2]]
            peakAHP3 <-  peakAHP2[[2]]
            yAHP <- min(peakAHP3[])
            
            AHPAmp <- yAHP - apthres;
            
            AP_Rise_Decay_R <- -Max_Rise/Max_Fall;
            
            #HalfWidth
            
            half1 <- apthres + APAmp/2;
            yp1x <- peakLoc[1]-200
            yp1 <- peak1[yp1x:peakLoc[1]]
            tt1 <- which.min(abs(yp1 - half1)) 
            time <- timeTn
            
            timett1 <- time[yp1x:peakLoc[1]]
            
            tt1x <- timett1[tt1]
            yp2x <- peakLoc[1]+200
            yp2 <- peak1[peakLoc[1]:yp2x]
            
            tt2 <- which.min(abs(yp2 - half1)) 
            
            timett2 <- time[peakLoc[1]:yp2x]
            tt2x <- timett2[tt2]
            
            halfwidth <- (tt2x - tt1x)*1000
            
            AP_Rise_Half <- Max_Rise/halfwidth
            
            
            ##ISIRatio
            end <- length(peakLoc)
            sumisi <- (peakLoc[end]-peakLoc[end-1])+ (peakLoc[end-1]-peakLoc[end-2])+ (peakLoc[end-2]-peakLoc[end-3]);
            ISIRatio <- (peakLoc[2]-peakLoc[1])/( sumisi/3);
            
            n3<- 3
            incProgress(1/n3, detail = paste("Starting Analysis"))
            
            ##Max Adaptation
            tnl<- length(tn)/2
            
            sumadapt <- 1/((peakLoc[end]-peakLoc[end-1])/tnl)+ 1/((peakLoc[end-1]-peakLoc[end-2])/tnl)+ 1/((peakLoc[end-2]-peakLoc[end-3])/tnl);
            
            
            MaxAdaptation <- (1/((peakLoc[2]-peakLoc[1])/tnl)- (sumadapt/3) );
            
            ##Cv2 and CV1
            
            nCV <- length(peakLoc);
            nCV1 <- nCV-1
            CV= vector();
            
            for (i in 1:nCV1){
              
              CV<- append(CV, (peakLoc[i+1]-peakLoc[i])/50)
              
            }
            
            CVn <- vector();
            n2CV <- length(CV);
            n2CV1 <- n2CV-1;
            
            for (i in 1:n2CV1){
              
              CVn <- append(CVn, 2*( abs((CV[i]-CV[i+1]))/((CV[i]+CV[i+1]))))
              
            }
            
            n3CV <- length(CVn);
            CVall <- mean(CVn[1:n3CV]);
            CV1 <- mean(CVn[2:n3CV]);
            CV2 <- mean(CVn[3:n3CV]);
            
            #Rm, RMP, Ct
            abftest$epoch$fEpochInitLevel
            abftest$epoch$fEpochLevelInc
            c2 = vector();
            for (i in sweep.max){
              
              n <-  i;
              n1L <- {n-1}*lengthofepisode+1
              n2L <- lengthofepisode*n
              tn <- abftest$traces[1,n1L:n2L];
              peak1list <- peakdet(tn, 20)
              example1  <- peak1list[[1]]
              if (length(example1[,1]) == 1){
                c2 <- i
                
              } else if (length(example1[,1]) > 1){
                c2<- i
                break
              }
              
              
            }
            creq <- c2-1;
            iv.v <-  vector();
            iv.i <-  vector();
            withProgress(message = 'Analyzing traces', value = 0,{
              n2 <- creq
              for (i in 1:creq){
                n <-  i;
                n1L <- {n-1}*lengthofepisode+1
                n2L <- lengthofepisode*n
                tn <- abftest$traces[1,n1L:n2L];
                lt <- length(tn);
                peakrmp <- tn[lt/2];
                peaki <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(i-1)
                iv.v<- append(iv.v,peakrmp)
                iv.i<- append(iv.i,peaki)
                
                incProgress(1/n2, detail = paste("Doing part", i))
                
              }})
            
            P <- polyfit(iv.i,iv.v,1)
            Rm <- P[1]*1000
            RMP <- P[2]
            
            Ct <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(c2-1)
            
            
            
            n3<- 5
            incProgress(1/n3, detail = paste("Halfway there!"))
            # Ctx2
            Ctx2 <- Ct*2
            
            c2 <- vector();
            for (i in sweep.max){
              n <-  i;
              n1L <- {n-1}*lengthofepisode+1
              n2L <- lengthofepisode*n
              tn <- abftest$traces[1,n1L:n2L];
              peak1 <- tn
              #peak1list <- peakdet(as.vector(peak1, mode = "any"), 20)
              #example1  <- peak1list[[1]]
              CtM <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(i-1)
              if (CtM == Ctx2){
                c2 <- i
              } else if (CtM > Ctx2){
                c2<- i
                break
              }
            }
            
            n <-  c2-1;
            n1L <- {n-1}*lengthofepisode+1
            n2L <- lengthofepisode*n
            tn <- abftest$traces[1,n1L:n2L];
            peak2 <- tn
            peak2list <- peakdet(as.vector(peak2, mode = "any"), 20)
            apnumberx <- peak2list[[1]]
            apnumber <- length(apnumberx[,1]) 
            peakLoc <- apnumberx[,1]
            
            ##ISIRatio
            end <- length(peakLoc)
            sumisiCtx2 <- (peakLoc[end]-peakLoc[end-1])+ (peakLoc[end-1]-peakLoc[end-2])+ (peakLoc[end-2]-peakLoc[end-3]);
            ISIRatioCtx2 <- (peakLoc[2]-peakLoc[1])/( sumisiCtx2/3);
            
            
            
            ##Max Adaptation
            
            tnl<- length(tn)/2
            sumadaptCtx2 <- 1/((peakLoc[end]-peakLoc[end-1])/tnl)+ 1/((peakLoc[end-1]-peakLoc[end-2])/tnl)+ 1/((peakLoc[end-2]-peakLoc[end-3])/tnl);
            
            
            MaxAdaptationCtx2 <- (1/((peakLoc[2]-peakLoc[1])/tnl)- (sumadaptCtx2/3) );
            
            
            ###Max Ap
            
            apnumberlist <- vector();
            iv.i.entire <- vector();
            for (i in sweep.max){
              
              n <-  i;
              n1L <- {n-1}*lengthofepisode+1
              n2L <- lengthofepisode*n
              tn <- abftest$traces[1,n1L:n2L];
              peak1 <- tn
              peak1list <- peakdet(peak1, 20)
              example1  <- peak1list[[1]]
              apnumber2 <- length(example1[,1]) 
              apnumberlist<- append(apnumberlist,apnumber2)
              peaki2 <- abftest$epoch$fEpochInitLevel + abftest$epoch$fEpochLevelInc*(i-1)
              iv.i.entire<- append(iv.i.entire,peaki2)
            }  
            
            n3<- 7
            incProgress(1/n3, detail = paste("Almost!"))
            
            apnumberlist  
            iv.i.entire  
            
            #dfap <- data.frame(iv.i.entire,apnumberlist)
            #ggplot(dfap,aes(x=iv.i.entire, y=apnumberlist)) + geom_line()
            
            dyap <- diff(apnumberlist)/diff(iv.i.entire)
            
            maxaploc <- which.max(apnumberlist)
            maxapsweep <- maxaploc
            
            n <-  maxapsweep;
            n1L <- {n-1}*lengthofepisode+1
            n2L <- lengthofepisode*n
            tn <- abftest$traces[1,n1L:n2L];
            peak1 <- tn
            peak1list <- peakdet(peak1, 20)
            example1  <- peak1list[[1]]
            apnumbermax <- length(example1[,1]) 
            peakLoc <- example1[,1]
            #ggplot(dfap,aes(x=iv.i.entire, y=apnumberlist)) + geom_line()
            
            ##ISIRatio
            end <- length(peakLoc)
            sumisimax <- (peakLoc[end]-peakLoc[end-1])+ (peakLoc[end-1]-peakLoc[end-2])+ (peakLoc[end-2]-peakLoc[end-3]);
            ISIRatiomax <- (peakLoc[2]-peakLoc[1])/( sumisimax/3);
            
            
            
            ##Max Adaptation
            
            tnl<- length(tn)/2
            sumadaptmax <- 1/((peakLoc[end]-peakLoc[end-1])/tnl)+ 1/((peakLoc[end-1]-peakLoc[end-2])/tnl)+ 1/((peakLoc[end-2]-peakLoc[end-3])/tnl);
            
            
            MaxAdaptationmax <- (1/((peakLoc[2]-peakLoc[1])/tnl)- (sumadaptmax/3) );
            
            
            
            
            Value <- c(as.character(input$filebatch$name[im]),RMP, Rm, Ct, Max_Fall, AHPAmp,ISIRatio,MaxAdaptation,CVall,CV2,AP_Rise_Decay_R,CV1,APThres_APAmp,AP_Rise_Half,AP_1_2_R,AP_2_3_R, apnumber, MaxAdaptationCtx2, ISIRatioCtx2, apnumbermax,MaxAdaptationmax,ISIRatiomax,halfwidth)
            #Property.Name <- c("RMP","Rm","Ct","Max Fall Rate","AHP Amp","ISI Ratio","Max Adapt","CV2","CV2- 1st/2nd AP","AP Rise Decay R","CV2 1st AP","AP Thres/Amp","AP Rise/HalfWidth","1st/2nd AP","2nd/3rd AP","AP Number at Ctx2","Max Adaptation at Ctx2","ISIRatio at Ctx2","Max AP Number","Max Adaptation at Max AP No","ISIRatio at Max AP No", "Halfwidth")
            
           #data1 <- data.frame(Property.Name, Value);
            
            n3<- 10
            incProgress(1/n3, detail = paste("Done!"))
          })
          Value
          testdrataframe = Cbind(testdrataframe,Value) 
          
          
          
          
          
          
          
        }
        
        
       
      
        
      
      
      
      
    }
    
    
    
    Property.Name1 <- c("FileName","RMP","Rm","Ct","Max Fall Rate","AHP Amp","ISI Ratio","Max Adapt","CV2","CV2- 1st/2nd AP","AP Rise Decay R","CV2 1st AP","AP Thres/Amp","AP Rise/HalfWidth","1st/2nd AP","2nd/3rd AP","AP Number at Ctx2","Max Adaptation at Ctx2","ISIRatio at Ctx2","Max AP Number","Max Adaptation at Max AP No","ISIRatio at Max AP No", "Halfwidth")
    #Property.Name2 <- c("FileName","RMP","Rm","Ct","Max Fall Rate","AHP Amp","ISI Ratio","Max Adapt","CV2","CV2- 1st/2nd AP","AP Rise Decay R","CV2 1st AP","AP Thres/Amp","AP Rise/HalfWidth","1st/2nd AP","2nd/3rd AP","Max AP Number","Max Adaptation at Max AP No","ISIRatio at Max AP No", "Halfwidth")
    
    testdrataframe = Cbind(testdrataframe,Property.Name1) 
    testdrataframe = as.data.frame(testdrataframe[-1])
    colnames(testdrataframe) <- unlist(testdrataframe[1,])
    data <- testdrataframe[-1,]
  })
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("Cell_Properties", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(batch(), file, row.names = FALSE)
    }
  )
  
  
}



# Run the app ----
shinyApp(ui = ui, server = server)