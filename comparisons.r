library("ggplot2")
library("ggrepel")
library("reshape2")
#getting p values
test.mouse<- function(ctrl,treat) {
  p.values <- matrix(nrow = dim(ctrl)[2])
  for (i in 1:dim(ctrl)[2] ) {
    
   p.values[i]<-t.test(treat[,i], ctrl[,i])$p.value
   }
   return(p.values)
}
#getting fold changes
fold.change.mouse<- function(ctrl,treat){
  fold.change<-matrix(nrow =dim(ctrl)[2])
  for (i in 1:dim(ctrl)[2]) {
    fold.change[i] <- mean(treat[,i])/mean(ctrl[,i])
    
      
  }
  return(fold.change)
}
#function for volcano plots
volcano.plot<-function(data) {
p<- ggplot(data, (aes(log2(fold.change), -log(p.adj) )))
p+
geom_point(aes(col = condition), size = rel(4))+
geom_hline(yintercept = -log(0.05)) + 
geom_vline(xintercept = 1.5) +
geom_vline(xintercept = -1.5) +
labs(title = "Altered Metabolites", x = "log2(Fold Change)", y = "-log(p Value) FDR Adjusted") +
geom_text_repel( data = filter(data, (p.adj<0.05 | log2(fold.change) > 1.5 | log2(fold.change) < -1.5)), 
		  aes(label = metabolite))+
theme(plot.title = element_text(size=(rel(2))), 
      axis.text = element_text(size = rel(1.4)),
      axis.title.x = element_text(size = rel(1.4)), 
      axis.title.y = element_text(size = rel(1.4)),
      legend.text = element_text(size = rel(1.4)), 
      legend.title = element_text(size = rel(1.4))
  ) 
}

#loading data and separating it into groups for comparisons
#there's probably a much better way to do this.
mouse<-read.csv("/mnt/disk2/tokuwa_metabolomics_data/mouse_pellet/mouse_pellets/assigned_metabolites/norm_data_r.csv", header = TRUE)
rownames(mouse)<- mouse[,1]
mouse<-mouse[,2:51]

t2ctrl<-mouse[grep("T2G1", rownames(mouse)),] 
t2cipro<-mouse[grep("T2G2", rownames(mouse)),] 
t2metro<-mouse[grep("T2G4", rownames(mouse)),] 
t2vanco<-mouse[grep("T2G3", rownames(mouse)),] 
t3ctrl<-mouse[grep("T3G1", rownames(mouse)),] 
t3cipro<-mouse[grep("T3G2", rownames(mouse)),]
t3metro<-mouse[grep("T3G4", rownames(mouse)),]
t3vanco<-mouse[grep("T3G3", rownames(mouse)),]
t4ctrl<-mouse[grep("T4G1", rownames(mouse)),]
t4cipro<-mouse[grep("T4G2", rownames(mouse)),]
t4metro<-mouse[grep("T4G4", rownames(mouse)),]
t4vanco<-mouse[grep("T4G3", rownames(mouse)),]
    
mouse.p.values<-matrix(nrow = 50, ncol = 9)
mouse.p.values[,1]<- test.mouse(t2ctrl,t2cipro)
mouse.p.values[,2]<- test.mouse(t2ctrl,t2metro)
mouse.p.values[,3]<- test.mouse(t2ctrl,t2vanco)
mouse.p.values[,4]<- test.mouse(t3ctrl,t3cipro)
mouse.p.values[,5]<- test.mouse(t3ctrl,t3metro)
mouse.p.values[,6]<- test.mouse(t3ctrl,t3vanco)
mouse.p.values[,7]<- test.mouse(t4ctrl,t4cipro)
mouse.p.values[,8]<- test.mouse(t4ctrl,t4metro)
mouse.p.values[,9]<- test.mouse(t4ctrl,t4vanco)

mouse.p.adj.values <- as.matrix(p.adjust(mouse.p.values, "fdr"))
dim(mouse.p.adj.values)<- dim(mouse.p.values)
mouse.p.adj.values<- data.frame(mouse.p.adj.values)
names(mouse.p.adj.values)<- c("t2cipro", "t2metro", "t2vanco", "t3cipro", "t3metro", "t3vanco", "t4cipro", "t4metro", "t4vanco")
rownames(mouse.p.adj.values)<- names(mouse)

mouse.fold.changes <- matrix(nrow = 50, ncol = 9)

mouse.fold.changes[,1]<- fold.change.mouse(t2ctrl,t2cipro)
mouse.fold.changes[,2]<- fold.change.mouse(t2ctrl,t2metro)
mouse.fold.changes[,3]<- fold.change.mouse(t2ctrl,t2vanco)
mouse.fold.changes[,4]<- fold.change.mouse(t3ctrl,t3cipro)
mouse.fold.changes[,5]<- fold.change.mouse(t3ctrl,t3metro)
mouse.fold.changes[,6]<- fold.change.mouse(t3ctrl,t3vanco)
mouse.fold.changes[,7]<- fold.change.mouse(t4ctrl,t4cipro)
mouse.fold.changes[,8]<- fold.change.mouse(t4ctrl,t4metro)
mouse.fold.changes[,9]<- fold.change.mouse(t4ctrl,t4vanco)
mouse.fold.changes<-data.frame(mouse.fold.changes)
names(mouse.fold.changes)<- c("t2cipro", "t2metro", "t2vanco", "t3cipro", "t3metro", "t3vanco", "t4cipro", "t4metro", "t4vanco")

rownames(mouse.fold.changes)<- names(mouse)
rownames(mouse.p.adj.values)<- names(mouse)
mouse.fold.changes<-as.matrix(mouse.fold.changes)
mouse.p.adj.values<-as.matrix(mouse.p.adj.values)

fm<- melt(mouse.fold.changes)
pm<-melt(mouse.p.adj.values)

mouse.met<- data.frame(c(fm,pm))
names(mouse.met)<-c("metabolite", "condition", "fold.change", "nothing", "nothing2", "p.adj")


#volcano plot function
volcano.plot(mouse.met[grep("t2", mouse.met$condition),])
volcano.plot(mouse.met[grep("t3", mouse.met$condition),])
volcano.plot(mouse.met[grep("t4", mouse.met$condition),])
volcano.plot(mouse.met)
