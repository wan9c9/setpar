import numpy as np
import statsmodels.formula.api as smf
import statsmodels.tsa.stattools as stattools


def sim_tpar(r,para,n):
    ny = int(n*1.5)
    y = np.zeros(ny)
    lmbd = 1.0
    for i in xrange(1,ny):
        if y[i-1] <= r:
            lmbd=para[0]+para[1]*lmbd+para[2]*y[i-1]
        else: 
            lmbd=para[3]+para[4]*lmbd+para[4]*y[i-1]
        x=np.random.poisson(lmbd,1)
        y[i]=x[0]
    return y[ny-n+1:]

r=6
para=[0.5,0.8,0.7,0.2,0.2,0.1]
obs = sim_tpar(6,para,10000)
ac = stattools.acf(obs)
print ac



