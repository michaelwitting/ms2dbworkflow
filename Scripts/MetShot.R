#### code starts here ###
#some sanity checks
#check for folder structure
if (file.exists("MetShotDir") & file.exists("MassBankDir")){
  print("all folders in place...")
} else {
  print("creating folders...")
  dir.create("MetShotDir")
  dir.create("MassBankDir")
  dir.create("MassBankDir/RMassBankOptions")
}

### Part 1: MetShot ###
#load required libraries
library(MetShot)

#read input file
#format
#id \t name \t smiles \t adduct \t position
#e.g. 11  L-alanine [M+Na]+ C[C@H](N)C(O)=O 1-A,1
standards<-data.frame(read.table("clipboard", header=T, sep="\t", comment.char = "", stringsAsFactors=FALSE))

#change to MetShotDir
setwd("MetShotDir")

#collision energies and other parameters
ionMode<-"negative"
colEnergies = cbind(c(-10,-20,-30))


template =paste(getwd(), "templates/MS_neg_MetShotTemplate.m/microTOFQMaxAcquisition.method", sep="/")

#result data frame
HyStarSequence<-data.frame(position=character(0),
                     sampleName=character(0),
                     methodName=character(0))

RMassBankCmpList<-data.frame(ID=character(0),
                                   Name=character(0),
                                   SMILES=character(0),
                                   RT=numeric(0),
                                   CAS=character(0))

#create and change to directory for todays data
dir.create(as.character(Sys.Date()))
setwd(as.character(Sys.Date()))

#iterate through standards list
for(i in seq(1:nrow(standards))) {

  #calculate mass for isolation
  #can be also done with RCDK
  neutralMass<-as.numeric(gsub("1\t", "", gsub(",", ".", shell(paste("cxcalc exactmass", standards$smiles[i]), intern=T, wait=T)[2])))
  
  #currently only [M+H]+ and [M-H]- supported
  if(ionMode=="positive") {
    adductMass<-neutralMass + 1.007276
  } else if(ionMode=="negative") {
    adductMass<-neutralMass - 1.007276
  } else {
    stop("Unknown ion mode! Enter \"positive\" or \"negative\"")
  }

  #round (works better with Bruker software)
  isolationMass<-round(adductMass, digits=2)
  
  #convert to peaks of interest
  peaksOfInterest <- cbind(mzmin=isolationMass,
                           mzmed=isolationMass,
                           mzmax=isolationMass,
                           rtmin=60,
                           rtmed=65,
                           rtmax=120)
  
  #convert to picklists
  picklists <- peaklist2picklist(peaksOfInterest, widthFactor=1, fillGaps=F)  
  
  sampleName<-paste("MRM", standards$id[i], ionMode, sep="_")
  methodName<-paste("MRM", standards$id[i], ionMode, sep="_")
    
  #new method for maxis
  picklist2method(picklists[[1]],
                  methodPrefix = methodName,
                  MSmode = ionMode,
                  template = template,
                  MSMSManual_ListCollisionEnergy = colEnergies,
                  MSMSManual_ListIsolationWidth = 8,
                  instrumentprefix = "qtofmaxacq",
                  templateSegmentNr=3)
    
  #append vial position, sample name and method name to dataframe
  HyStarSequence<-rbind.data.frame(HyStarSequence, cbind.data.frame(position=standards$position[i],
                                                        sampleName=sampleName,
                                                        methodName=paste(methodName, "-", paste(colEnergies, collapse=","),"eV.m", sep="")))
  
  #create list later needed for RMassBank
  RMassBankCmpList<-rbind.data.frame(RMassBankCmpList, cbind.data.frame(ID=standards$id[i],
                                                                        Name=standards$name[i],
                                                                        SMILES=standards$smiles[i],
                                                                        RT="",
                                                                        CAS=""))
}

#write sequence table
write.table(HyStarSequence, file="sequence.tsv", sep="\t")

#change to RMassBankDir and write compound list
setwd("../../MassBankDir")

#create and change to directory for todays data
dir.create(as.character(Sys.Date()))
setwd(as.character(Sys.Date()))
dir.create("rawDataFiles")
write.table(RMassBankCmpList, file="compoundlist.csv", sep=",", row.names=F, quote=F)
setwd("..")

#proceed with part 2 RMassBank.R