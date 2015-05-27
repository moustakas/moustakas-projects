"""Build figures for the talks I'm giving at the DESI collaboration
meeting at FNAL in May 2015.

jm15may25siena

To build the ELG spectrum do:
~/repos/git/desihub/desisim/bin/simulate-templates.py --nmodel 5 --objtype ELG --maxwave 25000
"""

import os
import numpy as np
import seaborn as sns
from astropy.io import fits
from astropy.table import Table
from desispec.io.util import header2wave
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from astropy.modeling import models, fitting
from desisim.templates import read_base_templates

sns.set(style='ticks', context='talk', font_scale=1.5, palette='deep')

talkdir = os.path.join(os.getenv('IM_RESEARCH_DIR'),'talks','desi','15may_desi')
print(talkdir)

# Velocity width

cflux, cwave, cmeta = read_base_templates(objtype='elg',observed=True)
sigma = cmeta['SIGMA_KMS']
sigminmax = np.log10([8.0,500.0])
binsz = 0.05
nbin = long((sigminmax[1]-sigminmax[0])/binsz)
ngal, bins = np.histogram(np.log10(sigma), bins=nbin, range=sigminmax)
cbins = (bins[:-1] + bins[1:]) / 2.0
   
ginit = models.Gaussian1D(ngal.max(), mean=70.0, stddev=20.0)
gfit = fitting.LevMarLSQFitter()
gauss = gfit(ginit,cbins,ngal)
    
xgauss = np.linspace(sigminmax[0],sigminmax[1],nbin*10)

fig = plt.figure(figsize=(8,6))
plt.bar(10**cbins, ngal, align='center', width=binsz*np.log(10)*10**cbins)
sns.despine()
plt.plot(10**xgauss, gauss(xgauss), '-', lw=3, color='firebrick')
plt.xlim(10**sigminmax)
plt.xlabel('$log_{10}\ Emission-line\ Velocity\ Width\ \sigma\ (km\ s^{-1})$')
plt.ylabel('$Number\ of\ Galaxies$')
plt.tick_params(axis='both', which='major', labelsize=14)
fig.subplots_adjust(bottom=0.15,left=0.15)
fig.savefig(talkdir+'/linesigma.pdf')


#plt.text(0.95,0.9,'$<log_{10}$'+' $\sigma>$ = '+
#         '{:.3f}$\pm${:.3f} km/s'.format(gauss.mean.value,gauss.stddev.value),
#         horizontalalignment='right',color='black',
#         transform=plt.gca().transAxes, fontsize=18)



# Example spectrum
tempfile = talkdir+'/elg-templates.fits'
flux, hdr = fits.getdata(tempfile, 0, header=True)
meta = Table(fits.getdata(tempfile, 1))
wave = header2wave(hdr)
print(flux.shape)

fig = plt.figure(figsize=(8,6))
ax = fig.gca()
plt.loglog(wave/1E4, 1E17*flux[2,:])
plt.axis('tight')
sns.despine()
plt.xlabel('$Wavelength\ (\mu$'+'$m)$')
plt.ylabel('$F_{\lambda}\ (10^{-17} erg\,s^{-1}\,cm^{-2}\,\AA^{-1}$)')
plt.xlim([0.3,2.5])
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.set_xticks([0.4,0.7,1.0,2.0])
fig.subplots_adjust(bottom=0.17,left=0.19)
fig.savefig(talkdir+'/elg.pdf')
