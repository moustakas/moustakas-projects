pro compute_mz_limits, doplot=doplot
; jm09mar25nyu - compute the limiting absolute magnitude and stellar
;   mass for the emission-line galaxy sample

    mzpath = ages_path(/projects)+'mz/'
    kk = read_mz_emline_sample(/mzhiiplus_ancillary)

; if I_lim=19.95 is the apparent magnitude limit of the survey, I is
; the apparent magnitude of the source, and M_g is the absolute
; magnitude (computed using K-correct), then the limiting absolute
; g-band magnitude is given simply by: M_g_lim = M_g - (I-I_lim)
;
; the limiting stellar mass, log_M_lim, is similarly given by:
; log_M_lim = log_M + 0.4(I-I_lim)
;
; compute the 50% and 95% limiting absolute magnitudes and stellar
; masses in ZBIN wide bins of redshift

    zbin = 0.05
    I_lim = 19.95
    zaxis = im_array(0.05,0.75,0.01)

; absolute magnitude
    mglim = kk.ugriz_absmag[1] - (kk.phot_I-I_lim)
    gg = im_medxbin(kk.z,mglim,zbin,weight=kk.spec_weight)
    mglim_50 = smooth(interpol(gg.medy,gg.binctr,zaxis),10)
    mglim_75 = smooth(interpol(gg.sigy75,gg.binctr,zaxis),10)
    mglim_95 = interpol(gg.sigy95,gg.binctr,zaxis)

; mass    
    masslim = kk.isedfit_mass + 0.4*(kk.phot_I-I_lim)
    mm = im_medxbin(kk.z,masslim,zbin,weight=kk.spec_weight,maxx=0.65)
    mmlim_50 = smooth(interpol(mm.medy,mm.binctr,zaxis),10)
    mmlim_75 = smooth(interpol(mm.sigy25,mm.binctr,zaxis),10)
    mmlim_95 = interpol(mm.sigy05,mm.binctr,zaxis)

; write out    
    limits = {zaxis: zaxis, $
      mglim_50: mglim_50, mglim_75: mglim_75, mglim_95: mglim_95, $
      mmlim_50: mmlim_50, mmlim_75: mmlim_75, mmlim_95: mmlim_95}
    im_mwrfits, limits, mzpath+'mz_limits.fits'
    
; make a plot    
    im_plotconfig, 0, psfile=mzpath+'mz_limits.ps'
    djs_plot, kk.z, kk.ugriz_absmag[1], psym=2, yrange=[-15,-24], $
      xtitle='Redshift', ytitle='M_{0.1g}', sym=0.3, xsty=1, ysty=1, $
      xrange=[0,0.8]
    im_legend, ['50%','75%','95%'], /right, /bottom, $
      box=0, color=['blue','orange','dark green'], line=[0,0,0], pspacing=1.2
    djs_oplot, gg.binctr, gg.medy, psym=4, color='red', sym=2.0
    djs_oplot, zaxis, mglim_50, color='blue', line=0
    djs_oplot, zaxis, mglim_75, color='orange', line=0
    djs_oplot, zaxis, mglim_95, color='dark green', line=0

    djs_plot, kk.z, kk.isedfit_mass, psym=2, yrange=[8,12], $
      xtitle='Redshift', ytitle='log (M_{*}/M'+sunsymbol()+')', sym=0.3, $
      xsty=1, ysty=1, xrange=[0,0.8]
    im_legend, ['50%','75%','95%'], /right, /bottom, $
      box=0, color=['blue','orange','dark green'], line=[0,0,0], pspacing=1.2
    djs_oplot, mm.binctr, mm.medy, psym=4, color='red', sym=2.0
    djs_oplot, zaxis, mmlim_50, color='blue', line=0
    djs_oplot, zaxis, mmlim_75, color='orange', line=0
    djs_oplot, zaxis, mmlim_95, color='dark green', line=0

    im_plotconfig, /psclose
    
stop    
    
return
end
    
