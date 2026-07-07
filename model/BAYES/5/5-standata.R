make_standata <- function(.dir) {
  ## get data file
  xdata <- readr::read_csv(file.path(.dir, "..", "..", "..", "..", "data", "derived", "RawData.csv"))

  #' Format data for Stan
  nt   <- nrow(xdata)
  iObsC <- with(xdata, (1:nrow(xdata))[EVID == 0 & CMT==1])
  iObsT<- with(xdata, (1:nrow(xdata))[EVID == 0 & CMT==3])
  nObsC <- length(iObsC)
  nObsT <- length(iObsT)
  xsub <- subset(xdata, !duplicated(ID))
  nSubjects <- length(xsub$ID)
  start <- (1:nt)[!duplicated(xdata$ID)]
  end <- c(start[-1] - 1, nt)
  cmtC = xdata$CMT
  cmtC[cmtC==1]=2
  cmtC[cmtC==3]=1
  cmtT = xdata$CMT
  cmtT[cmtT==1]=1
  cmtT[cmtT==3]=2
  
  data <- with(xdata,
                    list(nt = nt,
                         nObsC = nObsC,
                         nObsT = nObsT,
                         nSubjects = nSubjects,
                         nIIV = 8,
                         iObsC = iObsC,
                         iObsT = iObsT,
                         start = start,
                         end = end,
                         time = TIME,
                         cObsC = DV[iObsC],
                         cObsT = DV[iObsT],
                         amt =  AMT,
                         rate = RATE,
                         cmt = CMT,
                         cmtC = cmtC,
                         cmtT = cmtT,
                         evid = EVID,
                         ii = II,
                         addl = ADDL,
                         ss = ADDL*0,
                         runestimation=1))
}