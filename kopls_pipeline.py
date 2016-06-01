import numpy as np
import nmr_data_preprocess as nmr
import rv_modified as rv
import center_kernel_matrix as ckm
# get data
nmr_data = np.loadtxt('/home/louic/Desktop/mouse_pellets/all_cpmg/bin_shift_data')
cls = np.loadtxt('/home/louic/Desktop/mouse_pellets/all_cpmg/clss')
micro = np.loadtxt('/home/louic/Desktop/mouse_pellets/16s_data/data_rm.csv', delimiter = ',')

ppm = nmr_data[0]
miclabel  = np.array(range(1,35))  #number of taxons of microbiome data
nmr_data = nmr_data[1:122]
vanc = [nmr_data[i] for i in range(len(cls)) if cls[i] == 5 or cls[i] ==7]
vanc = np.array(vanc)

#normalise and scale nmr data
vanc = nmr.scalepqn(vanc)
vanc = nmr.pareto_scale(vanc)

vmicro = [micro[i] for i in range(len(cls)) if cls[i] == 5 or cls[i] ==7]
vmicro = np.array(vmicro)

vcls = [cls[i] for i in range(len(cls)) if cls[i] == 5 or cls[i] ==7]
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

rvn = rv.rv_modified(vanc, vcls)
rvm = rv.rv_modified(vmicro, vcls)

block1 = rvn*vn

block2 = rvm*vm

blocksum = block1+block2
#assuming data is already normalised and scaled

#create kernel -

#calcluate rv

#create addition matrix
#center matrix

#repeat cv
  #leave sample out
  
  #create training and test sets (%samples in training vs. test). needs randomisation with constraints
    #test number of components (by F1)
      #create model
	#test model
	#get F1
	#choose components and use that model
        #get r2x, dq2 save data

#extract loadings, do backscale loadings?

plot data (scores, loadings from best model, press, predictions)