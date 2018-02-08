import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.metrics import r2_score
import pickle
import pdb
import sklearn
import scipy
import scipy.io
from sklearn import metrics
import os,glob
from matplotlib.pyplot import cm
import matplotlib
matplotlib.rcParams.update({'font.size': 22})
import pandas as pd
import itertools
import seaborn as sns
from scipy import stats
from scipy.stats import norm
from math import exp, sqrt
from usefulfns import getcorr
import pickle
from sklearn.neighbors import KernelDensity

sns.set(font_scale = 1.5)
custom = {'axes.linewidth':5,'font.family':'sans-serif','font.sans-serif':['STHeiti']}
sns.set_style('white',custom)

# this is the one where we're going to take GLM classifier
pickle_in = open("/Volumes/norman/amennen/PythonMot5/evidencebystim_glmclassifier_alpha100_intercept_motion.pickle","rb")
evbystim = pickle.load(pickle_in)
# specify now which computer you're using!
motpath = '/Volumes/norman/amennen/motStudy05_transferred/'
behavioral_data = '/Volumes/norman/amennen/motStudy05_transferred/BehavioralData/'
savepath = '/Volumes/norman/amennen/PythonMot5/'
flatui = ["#DB5461", "#593C8F"]
all_sub = np.array([1,3,4,5,6,8,10,11,12,13,14,16,17,19,20,21,23,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40])
subtouse = all_sub
#subtouse =np.array([1,3,4,5,6,8,10,12,16,17,19,20,21,23,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40])
nSub= np.int(len(subtouse))
npairs = np.int(nSub/2)
RT_sub = np.array([1, 3, 4,5,6,8,10,12,13,14,19,21,23,26,29,32])
YC_sub = np.array([20,40,36,34,11,27,17,16,35,30,39,25,37,31,38,33])
RT_ind = np.searchsorted(all_sub,RT_sub)
YC_ind = np.searchsorted(all_sub,YC_sub)
subarray = np.zeros(nSub)
subarray[YC_ind] = 1

bw = .05

# load subjects data and plot each-histogram and smoothing
windowsize = 0.05
min = -1.5
max = -1*min + windowsize # to go one over
catrange = np.arange(min,max,windowsize)
cr2 = np.reshape(catrange,(len(catrange),1))

s=0
s_ind = subtouse[s]
sub = "Subject%01d" % s_ind
st=3      #thissep = allSep[st,:,s]
thissep = evbystim[sub][:,st]
x2 = np.reshape(thissep, (len(thissep), 1))
kde = KernelDensity(kernel='gaussian',bandwidth=bw).fit(x2)
allvals = np.exp(kde.score_samples(cr2))

# now plot original and gaussian
bins = np.arange(min,max+windowsize,windowsize)
z = plt.hist(x2,bins)
norm_counts = z[0]/np.sum(z[0])

plt.figure()
plt.plot(cr2,norm_counts)
plt.plot(cr2,allvals, 'r')
plt.show()

print(np.std(evbystim[sub]))
print(np.shape(evbystim[sub]))
print(1.06 * np.std(evbystim[sub]) * (np.shape(evbystim[sub])[0]*np.shape(evbystim[sub])[1])**-.2)
largeev = np.zeros((nSub,10*12))
for i in np.arange(nSub):
    s_ind = subtouse[i]
    sub = "Subject%01d" % s_ind
    thisev = np.reshape(evbystim[sub],(1,10*12))
    largeev[i,:] = thisev
large_vec = np.reshape(largeev,(1,10*12*nSub))
np.mean(largeev)
np.std(largeev)
formula = 1.06 * np.std(evbystim[sub]) * 120**-.2
# but I guess the kernel should be applied to each individual thing?
# applied on subject basis or no?