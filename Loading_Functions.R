# The following contain peices from abfload, a package in R used to loading ABF Files. ---------
# It currently only works for ABFV1 files and fucntions from the package are tweaked below to work for both ABFV1 and ABFV2 files

BLOCK_SIZE <- 512
STRING_JUNK_PREFIX <- 44

abfloadv2 <- function ( filename=NULL )
{
  if ( is.null(filename) )
  {
    filename <- file.choose()
  }
  
  result <- list()
  fp <- file(filename, open="rb")
  result$signature <- readChar(fp, nchars=4, useBytes=TRUE)
  result$filename <- filename
  
  if ( result$signature != "ABF2" )
  {
    close(fp)
    stop("Only ABF2 files supported at present")
  }
  
  class(result) <- "abf2"
  
  vb <- readBin(fp, what="integer", n=4, size=1)
  result$version <- vb[4] + vb[3] * 0.1 + vb[2] * 0.01 + vb[1] * 0.001
  
  result <- defRead(fp, header2Def, result)
  
  result$startTime <- result$startTime * 0.001
  
  # read the section index
  result$sections <- list()
  for ( sect.name in abfSections )
  {
    result$sections[[sect.name]] <- defRead ( fp, sectionInfo, repos=FALSE )
  }
  
  # read the strings section
  seek(fp, where=result$sections$StringSection$blockIndex * BLOCK_SIZE + STRING_JUNK_PREFIX)
  result$strings <- c(readBin(fp, "character", n=result$sections$StringSection$numEntries))
  
  # read the epoch
  
  seek(fp, where=result$section$EpochPerDACSection$blockIndex * BLOCK_SIZE + result$section$EpochPerDACSection$bytes)
  result$epoch <- defRead(fp, abfEpochDacInfo, repos=FALSE)
  
  # read the ADC defs
  result$ADC <- list()
  for ( ii in 1:result$sections$ADCSection$numEntries )
  {
    seek(fp, where=result$sections$ADCSection$blockIndex * BLOCK_SIZE + result$sections$ADCSection$bytes * (ii-1))
    result$ADC[[ii]] <- defRead(fp, abfADCInfo, repos=FALSE)
    result$ADC[[ii]]$name <- result$strings[result$ADC[[ii]]$ADCChannelNameIndex]
    result$ADC[[ii]]$units <- result$strings[result$ADC[[ii]]$ADCUnitsIndex]
  }
  
  # read the protocol block
  seek(fp, where=result$section$ProtocolSection$blockIndex * BLOCK_SIZE)
  result$protocol <- defRead(fp, abfProtocolInfo, repos=FALSE)
  
  ##read the epochperdac
  
  #seek(fp, where=result$section$EpochPerDacSection$blockIndex * BLOCK_SIZE)
  #result$epoch <- defRead(fp, abfEpochDacInfo, repos=FALSE)
  
  # and the tags, if any
  result$rawTags <- list()
  vLevel <- NA
  vUnits <- "" 
  result$tags <- data.frame(time=0, level=vLevel, units=vUnits, comment="", vChange=TRUE)
  if ( result$sections$TagSection$numEntries > 0 )
  {
    for ( ii in 1:result$sections$TagSection$numEntries )
    {
      seek(fp, where=result$sections$TagSection$blockIndex *
             BLOCK_SIZE + result$sections$TagSection$bytes * (ii-1))
      result$rawTags[[ii]] <- defRead(fp, abfTagInfo, repos=FALSE)
      
      # this is purely empirical based on what I've seen in our own files
      # I have no idea whether there's any more official scheme for this stuff
      rmx <- regexec("[^']*'([^']+)' => ([^ ]+) ([^ ]+)", result$rawTags[[ii]]$comment)
      result$rawTags[[ii]]$is.cmd <- rmx[[1]][1] != -1
      
      tagTime <- result$rawTags[[ii]]$tagTime *
        result$protocol$ADCSequenceInterval *
        1e-6 / length(result$ADC)
      
      if ( result$rawTags[[ii]]$is.cmd )
      {
        mxs <- regmatches(result$rawTags[[ii]]$comment, rmx)[[1]]
        result$rawTags[[ii]]$cmd <- mxs[2]
        result$rawTags[[ii]]$level <- as.numeric(mxs[3])
        result$rawTags[[ii]]$units <- mxs[4]
        
        vLevel <- result$rawTags[[ii]]$level
        vUnits <- result$rawTags[[ii]]$units
        
        result$tags <- rbind(result$tags,
                             data.frame(time=tagTime,
                                        level=vLevel,
                                        units=vUnits,
                                        comment="",
                                        vChange=TRUE,
                                        stringsAsFactors=FALSE))
      }
      else
      {
        result$comments <- rbind(result$comments,
                                 data.frame(time=tagTime,
                                            comment=paste(unlist(strsplit(result$rawTags[[ii]]$comment,
                                                                          " +")),
                                                          collapse=" "),
                                            stringsAsFactors=FALSE))
        result$tags <- rbind(result$tags,
                             data.frame(time=tagTime,
                                        level=vLevel,
                                        units=vUnits,
                                        comment=paste(unlist(strsplit(result$rawTags[[ii]]$comment,
                                                                      " +")),
                                                      collapse=" "),
                                        vChange=FALSE,
                                        stringsAsFactors=FALSE))
      }
    }
  }
  
  # OK, now we need to read the data
  if ( result$dataFormat == 0 )
  {
    dataSz <- 2
    dataType <- "integer"
  }
  else if ( result$dataFormat == 1 )
  {
    dataSz <- 4
    dataType <- "numeric"
  }
  else
  {
    warning(paste("unknown data format", result$dataFormat, "- data section not read"))
    close(fp)
    invisible(result)
  }
  
  headOffset <- result$sections$DataSection$blockIndex * BLOCK_SIZE
  #si <- result$protocol$ADCSequenceInterval
  
  if ( result$protocol$operationMode == 5 )
  {
    # gap free mode - this is our primary goal
    #result$dataPtsPerChannel <- result$sections$DataSection$numEntries / length(result$ADC)
    seek(fp, where=headOffset)
    result$rawdata <- readBin(fp, what=dataType, size=dataSz, n=result$sections$DataSection$numEntries)
    result$traces <- matrix(result$rawdata, nrow=length(result$ADC))
    result$s <- 0:(dim(result$traces)[2]-1) * (result$protocol$ADCSequenceInterval * 1e-6)
    
    # scale int data by ADC settings
    if ( result$dataFormat == 0 )
    {
      for ( ii in 1:length(result$ADC) )
      {
        adc <- result$ADC[[ii]]
        result$ADC[[ii]]$scaleFactor <- result$protocol$ADCRange /
          (result$protocol$ADCResolution *
             adc$instScaleFactor *
             adc$signalGain *
             adc$ADCProgGain *
             adc$teleAddGain)
        
        result$ADC[[ii]]$offset <- adc$instOffset - adc$signalOffset
        
        result$traces[ii,] <- result$traces[ii,] * result$ADC[[ii]]$scaleFactor + result$ADC[[ii]]$offset
      }
    }
  }
  else
  {
    warning(paste("unsupported data mode", result$protocol$operationMode, "- data section not read"))
  }
  
  close(fp)
  invisible(result)
}

defRead <- function ( fp, def, result=NULL, origin=0, repos=TRUE )
{
  if ( is.null(result) )
    result <- list()
  
  for ( ii in 1:length(def$field) )
  {
    if (repos) { seek(fp, where=def$offset[ii]+origin) }
    if ( def$type[ii]=="uint32" )
    {
      # R does not support unsigned long, so we have to fudge this to double
      lo <- readBin(fp, what="integer", size=2, signed=FALSE)
      hi <- readBin(fp, what="integer", size=2, signed=FALSE)
      result[[def$field[ii]]] <- 65536.0 * hi + lo
    }
    else if ( def$type[ii]=="string" )
    {
      result[[def$field[ii]]] <- readChar(fp, nchars=def$bytes[ii], useBytes=TRUE)
    }
    else if ( def$type[ii]=="skip" )
    {
      seek(fp, where=def$bytes[ii], origin="current")
    }
    else
    {
      result[[def$field[ii]]] <- readBin(fp, what=def$type[ii], size=def$bytes[ii])
    }
  }
  
  invisible(result)
}

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




