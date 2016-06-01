#load kopls library and data (consensus input data and class file)
library(kopls)
vanc<-read.csv('/home/louic/Desktop/mouse_pellets/16s_data/consensus_opls/t2ctrlvst3vanc/vanct2t3koplsinput', header = FALSE, sep = ' ')
vancclss <- read.csv('/home/louic/Desktop/mouse_pellets/16s_data/consensus_opls/t2ctrlvst3vanc/clss', header = FALSE, sep = ' ')
#make sure data is a matrix not dataframe
v<- data.matrix(vanc)
vancclss <- data.matrix(vancclss)

#creates kernel for kopls cv
Ktr<-koplsKernel(v,NA, 'g', sigma)
#do cross validation fo model 
modelCV<-koplsCV(Ktr,vancclss,1,3,nrcv=7,cvType='mccv',preProcK='no',
+                  preProcY='no', cvFrac = 0.33,modelType='da')


koplsPlotScores(model, x = NA, xsub = "p", y  = NA, ysub = 'o')