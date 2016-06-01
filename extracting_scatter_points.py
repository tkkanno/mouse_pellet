import numpy as np
import matplotlib.pyplot as plt
def scaletointegral(ydata):
        integral = ydata.sum(1) / 100
        newydata = np.transpose(ydata.T / integral)
        return np.array(newydata)
        
def scalepqn(ydata):
      # 1: scale to integral
      newydata = scaletointegral(ydata)
      # 2: reference spectrum is median of all samples
      yref = np.median(ydata, 0)
      # 3: quotient of test spectra with ref spectra
      yquot = newydata / yref
      # 4: median of those quotients
      ymed = np.median(yquot, 0)
      # 5: divide by this median
      newydata = newydata / ymed
      return np.array(newydata)
      
def mean_center(data):
  data = data - np.mean(data,0)
  data = np.nan_to_num(data)
  return data
  
def scale_pareto(data):
    data = mean_center(data) / np.sqrt(np.std(data,0))
    data = np.nan_to_num(data)
    return data
    
directory = '/mnt/disk2/tokuwa_metabolomics_data/mouse_pellet/t2vst3/metro/'
metabolite_list = 'ppm_heatmap'
datafile = 'shifted_data'
classfile = 'class'
labelfile  = 'text_metabolites'

datlist = np.loadtxt(directory+metabolite_list)
data = np.loadtxt(directory+datafile)
clss = np.loadtxt(directory+classfile)
labels = np.genfromtxt(directory+labelfile, dtype = 'str')

ppm = data[0,:]
data = data[1:,:]
norm_dat = scalepqn(data)
scale_dat = scale_pareto(data)

x,y,z =[],[],[]
for i in datlist:
  idx = np.argmin(np.abs(ppm-i))
  x.append(data[:,idx])
  y.append(norm_dat[:,idx])
  z.append(scale_dat[:,idx])

x = np.vstack(x)
y = np.vstack(y)
z = np.vstack(z)

n= clss.shape[0]
clss = np.reshape(clss, (n,1))
x, y, z = x.T, y.T, z.T
x = np.hstack((clss,x))
y = np.hstack((clss,y))
z = np.hstack((clss,z))

x = x[x[:,0].argsort()]
y = y[y[:,0].argsort()]
z = z[z[:,0].argsort()]


np.savetxt(directory+'raw_metabolite_data',  x, delimiter = ',')
np.savetxt(directory+'norm_metabolite_data', y, delimiter = ',')
np.savetxt(directory+'scale_metabolite_data',z, delimiter = ',')

a = x[:,1:]
b = y[:,1:]
c = z[:,1:]

#fig,ax = plt.subplots()

#for i in range(len(data)):
  #ax.plot(ppm, data[i])
  
#for i in range(len(datlist)):
  #idx = np.argmin(np.abs(ppm-datlist[i]))
  #print i, datlist[i], idx, data[:,idx].max()
  #ax.annotate(("%s, %.4f")%(i, data[:,idx].max()), xy = (datlist[i], data[:,idx].max()), \
		  #xycoords = 'data', xytext =(datlist[i], data[:,idx].max() + 10 ),\
		  #arrowprops = dict(arrowstyle = "->"),rotation = 270)
#plt.show()
####generating 50 boxplots if you so desire from this data. is impossible to read but ok
##m = 10
##n = 5
#count = 0
#fig = plt.subplots(m,n)
#for i in range(m):
  #for j in range(n):

    #fig[1][i][j].boxplot((c[0:8,count], c[8:16,count]))
    #fig[1][i][j].text(.5,.75,'%s, %.3f ppm'%(labels[count], datlist[count]), \
				 #horizontalalignment = 'center', \
				 #transform = fig[1][i][j].transAxes)
    #count +=1
    #print i, j, count


#plt.show()




