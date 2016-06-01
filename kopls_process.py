import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/louic/Desktop/metabolomics-gui/tokuwa_scripts/kopls')
import nmr_data_preprocess as nmr
import rv_modified as rv
import center_kernel_matrix as ckm
#import nmr data pre-proccesss script with normalisation and scaling

nmr = np.loadtxt('/home/louic/Desktop/mouse_pellets/all_cpmg/bin_shift_data')
cls = np.loadtxt('/home/louic/Desktop/mouse_pellets/all_cpmg/clss')
micro = np.loadtxt('/home/louic/Desktop/mouse_pellets/16s_data/data_rm.csv', delimiter = ',')
chosen_classes = [5,7]
ppm = nmr[0]
miclabel  = np.array(range(1,35))  #number of taxons of microbiome data
nmr = nmr[1:122]
#vanc = [nmr[i] for i in range(len(cls)) if cls[i] == 5 or cls[i] ==7] #this is for van t2 vs t3
vanc = [nmr[i] for i in range(len(cls)) if cls[i] == chosen_classes[0] or cls[i] ==chosen_classes[1]] #this is for ctrl t4 vanc t4
vanc = np.array(vanc)

#normalise and scale nmr data
vanc = scalepqn(vanc)
vanc = pareto_scale(vanc)

vmicro = [micro[i] for i in range(len(cls)) if cls[i] == chosen_classes[0] or cls[i] ==chosen_classes[1]]
vmicro = np.array(vmicro)

vcls = [cls[i] for i in range(len(cls)) if cls[i] == chosen_classes[0] or cls[i] ==chosen_classes[1]]
vcls = np.array(vcls)

new_cls = np.unique(vcls)
for i in range(len(vcls)):
  if vcls[i] == new_cls[0]:
    vcls[i] = 0
  else:
    vcls[i] = 1
n = len(vcls)
vcls = vcls.reshape(n,1)

vn = vanc.dot(vanc.T)
vm = vmicro.dot(vmicro.T)

rvn = rv_modified(vanc, vcls)
rvm = rv_modified(vmicro, vcls)

block1 = rvn*vn

block2 = rvm*vm

blocksum = block1+block2
#center matrix before exporting to R
blocksum = center_kernel_matrix(blocksum)

np.savetxt('/home/louic/Desktop/mouse_pellets/16s_data/vanct3_nmr_micro_uncentered_kernel', blocksum)
np.savetxt('/home/louic/Desktop/mouse_pellets/16s_data/vanct3_nmr_micro_uncentered_kernel_Y', vcls)
#now go into R and do KOPLS, kernel isn't centered but is already scaled

###R enviornment

data <- read.csv('/home/louic/Desktop/mouse_pellets/16s_data/vanct3_nmr_micro_uncentered_kernel', header = FALSE, sep = ' ')
cls <- read.csv('/home/louic/Desktop/mouse_pellets/16s_data/vanct3_nmr_micro_uncentered_kernel_Y', header = FALSE, sep = ' ')

data<- data.matrix(data)
cls <- data.matrix(cls)
cls<- cls+1
dataCenter<-koplsCenterKTrTr(data)

#do and plot PCA of block
svd.res<-svd(data, nu=2, nv=2)  # add colors to plot if neccessary for presentation
#plot PCA of original data
plot(svd.res$u[,1], svd.res$u[,2], col=cls, pch=16, xlab="PC1", ylab="PC2", main="PCA of original data set")
#create kopls model
model <- koplsModel(data,cls, 1,3, preProcK = 'no', preProcY = 'no')
#plot scores of kopls model
 koplsPlotScores(model, col=cls, pch=pch.vec)
#export Tp[,1] and to[,1] for graphing in python
tp<- as.data.frame(model$Tp[1])[,1]
to<- as.data.frame(model$To[,1])

write.table(tp, '/home/louic/Desktop/mouse_pellets/16s_data/t3vanc_Tp1', sep = ',')
write.table(to, '/home/louic/Desktop/mouse_pellets/16s_data/t3vanc_T01', sep = ',')
#plot with ggplot?


#modify data so numpy can read it, i'm lazy right now and just doing excel
#back to python
tpto = np.loadtxt('/home/louic/Desktop/mouse_pellets/kopls_t3vanc/tpto.csv', delimiter =',')

tp =tpto[:,0]
to = tpto[:,1]
#back calculate loadings
def calc_loadings(Xa, t, RV):
  x = (RV*(Xa.T).dot(t)) / ((t.T).dot(t))
  return x

tp_load_nmr = calc_loadings(vanc, tp, rvn)
to_load_nmr = calc_loadings(vanc, to, rvn)
tp_load_mic = calc_loadings(vmicro, tp,rvm)
to_load_mic = calc_loadings(vmicro,to,rvm)

tmax = tp_load_nmr.max()
tmin = tp_load_nmr.min()
kmax = to_load_nmr.max()
kmin = to_load_nmr.min()

threshold = 0.9
s = 80
fig,ax = plt.subplots()
ax.scatter(tp_load_nmr, to_load_nmr, s=s, color= 'blue', alpha = 0.3)
ax.scatter(tp_load_mic, to_load_mic,s=s,color =  'red', alpha = 0.75)
for i in range(len(miclabel)):
  ax.annotate(miclabel[i], xy = (tp_load_mic[i], to_load_mic[i]), xytext = (tp_load_mic[i], to_load_mic[i]))
for i in range(len(ppm)):
  if tp_load_nmr[i] > tmax*threshold or tp_load_nmr[i] < tmin*threshold or to_load_nmr[i] < kmin*threshold or to_load_nmr[i]>kmax*threshold:
    print i, ppm[i]
    ax.annotate(' %.3f'%ppm[i], (tp_load_nmr[i],to_load_nmr[i]))

ax.set_xlabel('Predictive component (Tp,1)')
ax.set_ylabel('Orthogonal Component (To,1)')
ax.set_title('KOPLS Loadings (Microbiome + NMR): T3 ctrl vs T3 vanco-imi')

plt.show()
  