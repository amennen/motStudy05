import numpy as np
import scipy
from scipy import stats
from scipy.stats import norm
def getcorr(subjectMatrix,TRmatrix):
    nwin = np.shape(TRmatrix)[1]
    nsub = np.shape(subjectMatrix)[1]
    allr = np.zeros((nwin, nsub))
    allp = np.zeros((nwin, nsub))
    for s in np.arange(nsub):
        thissim = subjectMatrix[:, s]
        for w in np.arange(nwin):
            thisTR = TRmatrix[:, w, s]
            allr[w, s], allp[w, s] = scipy.stats.pearsonr(thissim, thisTR)
    return allr