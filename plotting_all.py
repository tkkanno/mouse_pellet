import numpy as np
import matplotlib.pyplot as plt
metabolites ='/home/louic/Desktop/mouse_pellets/metabolite_list'
metabolites = np.loadtxt(metabolites)
direct = '/home/louic/Desktop/mouse_pellets/cpmg_mouse_pellet/align'

nmr_data = np.loadtxt('/home/louic/Desktop/mouse_pellets/all_cpmg/bin_shift_data')
spect = np.loadtxt('/home/louic/Desktop/mouse_pellets/all_cpmg/shifted_data')
ppm = data[0]
data = data[1:,:]

def collect_annotated_metabolites(ppm, data, metabolites):
  annotated_metabolites = []
  for i in metabolites:
    idx = np.argmin(np.abs(ppm - i))
    print "%4f - %4f = %4f" %(i, ppm[idx], (i-ppm[idx]))
    temp = data[:,idx]
    annotated_metabolites.append(temp)
  annotated_metabolites = np.vstack(annotated_metabolites)
  return annotated_metabolites

def plot_all(ppm,data, ppms, spect):
  data = (1000/np.max(data))*data
  spect = (1000/np.max(spect))*spect
  fig,ax = plt.subplots(2, sharex=True, sharey = True)
  m = np.mean(data,0)
  n = np.mean(spect,0)
  for i in range(len(data)):
    ax[0].plot(ppm, data[i])
    ax[1].plot(ppms, spect[i])
  ax[0].plot(ppm, m+50, 'r')
  ax[1].plot(ppms, n+100, 'r')
  ax[0].invert_xaxis()
  plt.show()
  #plots the aligned spectra and the binned spectra
  #along with the averaged spectra, helpful for visualisation when doing assignments
  
plot_all(ppm,data,ppms, spect)