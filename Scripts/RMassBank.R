###### Part 2 #######
#load required libraries
library(RMassBank)

#load RMassBank setting (differnt setting files for different methods are possible)
loadRmbSettings("RMassBankOptions/RMB_options.ini")

#create new MS/MS workspace
w <- newMsmsWorkspace()
w@files <- list.files("2016-06-25/rawDataFiles", full.names=T)

#load list from previous script with additional retention times
loadList("2016-06-25/compoundlist.csv")

#mzR workflow of RMassBank is used
#important check the ion mode
#"pH", "pNa", "pM", "pNH4", "mH", "mM", "mFA" for different ions ([M+H]+, [M+Na]+, [M]+, [M+NH4]+, [M-H]-, [M]-, [M+FA]-).
#Accession: XXYYYYZZ, XX = Letter, YYYY = ID, ZZ = Adductshifts
w <- msmsRead(w, mode="pH", files=w@files, readMethod="mzR")
w <- msmsWorkflow(w, mode="pH", steps=2:8)

#create new MassBank workspace
mb <- newMbWorkspace(w)

#execute the first steps, incl. getting synonyms etc...
mb <- mbWorkflow(mb, 1:2)

#after checking infolist.csv load new version
mb <- loadInfolist(mb,"./infolist.csv")
mb <- mbWorkflow(mb, 3:8)

