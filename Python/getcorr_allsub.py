import numpy as np
import scipy
import numpy
from scipy import stats
from scipy.stats import norm

def getcorr_vector(subjectMatrix,TRmatrix):
    nwin = np.shape(TRmatrix)[1]
    nsub = np.shape(subjectMatrix)[1]
    allr = np.zeros((nwin))
    allp = np.zeros((nwin))
    megadatamatrix = np.empty((0, nwin))
    for s in np.arange(nsub):
        subdata = TRmatrix[:, :, s]
        megadatamatrix = np.concatenate((megadatamatrix, subdata), axis=0)
    thissim = subjectMatrix.flatten(order="F")
    for w in np.arange(nwin):
        thisTR = megadatamatrix[:, w]
        nas = np.logical_or(np.isnan(thissim), np.isnan(thisTR))
        if len(np.argwhere(~nas)[:,0]) > 2:
            allr[w], allp[w] = scipy.stats.pearsonr(thissim[~nas],thisTR[~nas])
        else: # not enough samples
            allr[w] = np.nan
            allp[w] = np.nan
    return allr
# taken from: http://scipy-cookbook.readthedocs.io/items/SignalSmooth.html


