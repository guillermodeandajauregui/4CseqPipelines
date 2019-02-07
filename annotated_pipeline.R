###################################
# FourCSeq annotated pipeline
#
# an RScript with commented code
# of the FourCSeq vignette
# and reference manuals
#
# Guillermo de Anda - Jáuregui
#
###################################

library(FourCSeq)
library(NOISeq)

#example files

#################################################################
#remember, the bam files are preprocessed using the 
#HTSeq pythonpackage 
#(http://www-huber.embl.de/users/anders/HTSeq/doc/install.html
#
#python pathToScriptFile/demultiplex.py --fastq YourFASTQFile --barcode YourBarcodeFile
#
# where pathToScriptFile is [home]/R/[version]/FourCSeq/extdata/python 
# or something to that effect
#
#################################################################

#I copied the example files to a data folder for ease of use
#names are self-explanatory
referenceGenomeFile = "data/dm3_chr2L_1-6900.fa"
bamFilePath = "data/bam"
primerFile = "data/primer.fa"


#we may want to make our projectPath and fragmentDir

#will keep commented now because I have these created in my working dir. 

#dir.create(path = "exampleData/")
#dir.create(path = "exampleData/re_fragments/")

#here, we manually add metadata 
#which is experimental data information
#for the fourC object

metadata <- list(projectPath = "exampleData",  # directory to store results
                 fragmentDir = "re_fragments", # directory to store restriction fragment info
                                               # subdir of projectPath
                 referenceGenomeFile = referenceGenomeFile, #path to the reference genome or a BSgenome object.
                 reSequence1 = "GATC", # restriction enzyme recognition pattern of the first restriction enzyme
                 reSequence2 = "CATG", # reSequence2, restriction enzyme recognition pattern of the second restriction enzyme
                 primerFile = primerFile, # path to file with primer sequences for 4C library prep
                 bamFilePath = bamFilePath #path to dir where bam files are stored
                 )


#this colData is also needed for the fourC object
#it is a data frame with info for each library

colData <- DataFrame(viewpoint = "testdata", #name of the viewpoint
                     condition = factor(rep(c("WE_68h",   # id of exp condition
                                              "MESO_68h", # idem
                                              "WE_34h"),  # idem
                                            each=2), #for this example, 
                                                     #they use rep to generate a 
                                                     #3 level factor vector w/ duplicates
                                        levels = c("WE_68h", 
                                                   "MESO_68h", 
                                                   "WE_34h")
                                        ),
                     replicate = rep(c(1, 2),3), #replicate id for each exp condition
                                                 #again, using rep to fill the vector
                     bamFile = c("CRM_ap_ApME680_WE_6-8h_1_testdata.bam",    #paths to file bam files
                                 "CRM_ap_ApME680_WE_6-8h_2_testdata.bam",    #paths to file bam files
                                 "CRM_ap_ApME680_MESO_6-8h_1_testdata.bam",  #paths to file bam files
                                 "CRM_ap_ApME680_MESO_6-8h_2_testdata.bam",  #paths to file bam files
                                 "CRM_ap_ApME680_WE_3-4h_1_testdata.bam",    #paths to file bam files
                                 "CRM_ap_ApME680_WE_3-4h_2_testdata.bam"),   #paths to file bam files
                     sequencingPrimer="first" #or "second", depending on the sequencing start site 
                     )


#make the FourC object

fc <- FourC(colData = colData, 
            metadata = metadata
            )

#explore the FourC object
fc

isS4(fc) #FourC objects are S4 objects. 

#########################
#Fragment reference
#########################

fc <- addFragments(object  = fc, 
                   save    = TRUE, #if we add this, we get the fragments saved
                   minSize =  20,  #default = 20, min size of res.frag. end 
                   filter  = TRUE  #if true, frags smaller than minSize OR
                                   #not containing cutting site of 2nd res enzyme
                                   #are removed
                  )

fc #note that rowData now has 4 elements, whereas before it had none
rowRanges(fc)

#########################
#viewpoint information
#########################

#first, we use findViewpointFragments
#this will NOT return data to the environment
#rather,
#it will write out things to the project directory
#this may be time consuming 
findViewpointFragments(fc)
#now we find primerFragments.rda and .txt in the fragmentDir

#we may add the viepoint Frags to the fc object
fc <- addViewpointFrags(object = fc, 
                        primerFragFile = "primerFragments.rda" #this is the default path
                          )

#########################
#counting reads at 
#fragment ends
#########################
fc@assays #empty

fc <- countFragmentOverlaps(object = fc, 
                            trim=4, #number of bases to trim at read start,
                                    #default 0
                                    #note: the vignette uses 4, so examples weren't pretrimmed?
                            minMapq=30, #min mapping quality. Default 0. Negatives skip filtering
                            shift  =0 #max diff. in starts or ends between read and fragment
                            )

#fc is updated with two new assays in the assay slot
fc@assays #non-empty now!

#now the counts are added

fc <- combineFragEnds(object = fc, 
                      multFactor = 1, #mult.factor if only one valid end
                      filter = FALSE  #only reads from valid fragment ends
                      )

assays(fc) #now its longer
assay(fc, "counts")

############################################
# Now we load data fc, which is 
# the "ap" viewpoint that was created
# for the whole chromosmes 2L and 2R 
# of the dm3 reference genome
############################################
data(fc) #bye bye to our previous fc 

metadata(fc)$projectPath <- "exampleData"

#write out 
writeTrackFiles(object = fc, 
                assay = "counts",  #assay to be saved as track file 
                folder = "tracks", #path relative to project folder
                format = "bw",     #or bedGraph
                removeZeros = TRUE #if true, fragments with zero counts are removed
                )

writeTrackFiles(fc, format='bedGraph')

################################################
#smooth PCR arctifact spikes
################################################

fc <- smoothCounts(object = fc, 
                   assay = "counts", 
                   binWidth = 5 #integer vector of odd  numbers
                   )
fc@assays # a new assays is added, counts_5

#####
#plot count values 
#for replicate reproducibility
#####

#this will plot for 2 replicates
plotScatter(fc[,c("ap_WE_68h_1", "ap_WE_68h_2")], 
            assay = "counts",
            xlab="Replicate1", 
            ylab="Replicate2", 
            asp=1
            )

#if the conditions are not defined, the function plots for all conditions
#without markings
plotScatter(fc,
            xlab="Replicate1", 
            ylab="Replicate2", 
            asp=1
)

### this function ain't great, 
### it'd be better to extract the data by ourselves...

##not done, but something along the lines of 
#my_data = assays(fc)[["counts"]]
# ggplot(my_data, aes(condition1, condition2)) +
#  geom_density_2d() + 
#  stat_density_2d()
## but haven't made it work yet

################################################
#Detecting interactions 
################################################

fcf <- getZScores(object = fc, 
                  removeZeros = TRUE,
                  minCount    = 40, 
                  minDist     = NULL,
                  fitFun      = "distFitMonotoneSymmetric", 
                  sdFun       = mad
                  )

fcf #a new class fourC object 
    #only fragments kept fitting the criteria from getZscore

zScore <- assay(fcf, "zScore")
zScore %>% head

#the suggestion is plotting with hist 
hist(zScore[,"ap_MESO_68h_1"], breaks = 100)
hist(zScore[,"ap_MESO_68h_2"], breaks = 100)

#but we can access this through assay 
qplot(zScore[,"ap_MESO_68h_1"], geom="histogram", binwidth = 0.05) 
qplot(zScore[,"ap_MESO_68h_2"], geom="histogram", binwidth = 0.05) 


####
#check, whether normal assumption for calculating the p-values is justified
####

qqnorm(zScore[,"ap_MESO_68h_1"],
       main="Normal Q-Q Plot - ap_MESO_68h_1")
abline(a=0, b=1)

qqnorm(zScore[,"ap_MESO_68h_2"],
       main="Normal Q-Q Plot - ap_MESO_68h_2")
abline(a=0, b=1)

###
#Interacting regions
###

#contains booleans 
#indicating whether an interaction has been called for a fragment ornot.
fcf <- addPeaks(object = fcf, 
                zScoreThresh = 3, 
                fdrThresh=0.01
                )

###
#Visualize the fits used to calculate the z-scores
###
plotFits(fcf[,1], main="")
?plotFits

###
#
###
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
plotZScores(fcf[,c("ap_WE_68h_1", "ap_WE_68h_2")],
            txdb=TxDb.Dmelanogaster.UCSC.dm3.ensGene)