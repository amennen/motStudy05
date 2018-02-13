import numpy as np
import scipy
import numpy
from scipy import stats
from scipy.stats import norm

def getcorr(subjectMatrix,TRmatrix, subjname,TRname):
    nwin = np.shape(TRmatrix)[1]
    nsub = np.shape(subjectMatrix)[1]
    allr = np.zeros((nwin, nsub))
    allp = np.zeros((nwin, nsub))
    for s in np.arange(nsub):
        thissim = subjectMatrix[:, s]
        for w in np.arange(nwin):
            thisTR = TRmatrix[:, w, s]
            nas = np.logical_or(np.isnan(thissim), np.isnan(thisTR))
            if len(np.argwhere(~nas)[:,0]) > 2:
                allr[w, s], allp[w, s] = scipy.stats.pearsonr(thissim[~nas],thisTR[~nas])
               # print('TEST')
            else: # not enough samples
                allr[w,s] = np.nan
                allp[w,s] = np.nan

    pvals = np.zeros((nwin))
    print('pvalues for subjectmatrix: ' + subjname + ' and TR matrix: ' + TRname)
    for w in np.arange(nwin):
        d = allr.T[:, w]
        nas = np.isnan(d)
        pvals[w] = stats.ttest_1samp(d[~nas], 0).pvalue
    print(pvals)

    return allr
# taken from: http://scipy-cookbook.readthedocs.io/items/SignalSmooth.html


