pro atlas_varykl_specfit
; jm05oct20uofa - the referee wants to know the effect of using a
;                 different k(lambda) on the continuum reddening
;                 results; so redo the fitting with the O'Donnell and
;                 the SMC curves, and compare them against the Charlot
;                 & Fall results

    datapath = atlas_path(/atlas1d)
    specfitpath = atlas_path(/specfit)+'klambda_tests/'
    atlas = atlas_read_info()

    int = where(atlas.drift and (atlas.drift_agnflag eq 0L))
    speclist = strtrim(atlas[int].drift_file,2)

; O'Donnell (1994)    
    
    suffix = 'integrated_atlas_odonnell'

    specdata = ispeclinefit(speclist,specres=8.0,snrcut=1.0,dustmodel=0,$
      datapath=datapath,linepath=specfitpath,suffix=suffix,Zmulti=Zmulti,$
      /odonnell,/zcrosscor,/postscript,/write,vmaxshift=1000.0,/nologfile,$
      starvdisp=100.0)

    suffix = 'integrated_atlas_smc'

    specdata = ispeclinefit(speclist,specres=8.0,snrcut=1.0,dustmodel=0,$
      datapath=datapath,linepath=specfitpath,suffix=suffix,Zmulti=Zmulti,$
      /smc,/zcrosscor,/postscript,/write,vmaxshift=1000.0,/nologfile,$
      starvdisp=100.0)

; analysis

    smc = mrdfits('53664_integrated_atlas_smc_specdata.fits.gz',1)
    mw = mrdfits('53664_integrated_atlas_odonnell_specdata.fits.gz',1)
    cf = read_integrated()
    cf = cf[sort(cf.ra)]
;   niceprint, cf.galaxy, mw.galaxy

; chi2    
    
    jj = im_stats(cf.continuum_chi2/smc.continuum_chi2,/verbose)
    jj = im_stats(cf.continuum_chi2/mw.continuum_chi2,/verbose)
    jj = im_stats(smc.continuum_chi2/mw.continuum_chi2,/verbose)
    jj = im_stats(cf.continuum_ebv[0]-smc.continuum_ebv[0],/verbose)
    jj = im_stats(cf.continuum_ebv[0]-mw.continuum_ebv[0],/verbose)
    jj = im_stats(smc.continuum_ebv[0]-mw.continuum_ebv[0],/verbose)

    g = where((smc.h_alpha[1] gt 0.0) and (cf.h_alpha[1] gt 0.0))
    jj = im_stats(alog10(smc[g].h_alpha[0]/cf[g].h_alpha[0]),/verbose)
    g = where((smc.h_alpha_ew[1] gt 0.0) and (cf.h_alpha_ew[1] gt 0.0))
    jj = im_stats(alog10(smc[g].h_alpha_ew[0]/cf[g].h_alpha_ew[0]),/verbose)

    g = where((smc.h_beta[1] gt 0.0) and (cf.h_beta[1] gt 0.0))
    jj = im_stats(alog10(smc[g].h_beta[0]/cf[g].h_beta[0]),/verbose)
    g = where((smc.h_beta_ew[1] gt 0.0) and (cf.h_beta_ew[1] gt 0.0))
    jj = im_stats(alog10(smc[g].h_beta_ew[0]/cf[g].h_beta_ew[0]),/verbose)

    g = where((smc.oiii_5007[1] gt 0.0) and (cf.oiii_5007[1] gt 0.0))
    jj = im_stats(alog10(smc[g].oiii_5007[0]/cf[g].oiii_5007[0]),/verbose)
    g = where((smc.oiii_5007_ew[1] gt 0.0) and (cf.oiii_5007_ew[1] gt 0.0))
    jj = im_stats(alog10(smc[g].oiii_5007_ew[0]/cf[g].oiii_5007_ew[0]),/verbose)

    g = where((smc.oii_3727[1] gt 0.0) and (cf.oii_3727[1] gt 0.0))
    jj = im_stats(alog10(smc[g].oii_3727[0]/cf[g].oii_3727[0]),/verbose)
    g = where((smc.oii_3727_ew[1] gt 0.0) and (cf.oii_3727_ew[1] gt 0.0))
    jj = im_stats(alog10(smc[g].oii_3727_ew[0]/cf[g].oii_3727_ew[0]),/verbose)

    plothist, cf.continuum_ebv[0]-smc.continuum_ebv[0], bin=0.01
    plot, cf.continuum_ebv[0], smc.continuum_ebv[0], ps=4, xsty=3, ysty=3
    plot, cf.continuum_ebv[0], smc.continuum_ebv[0]-cf.continuum_ebv[0], ps=4, xsty=3, ysty=3
    plot, cf.h_alpha_ew[0], alog10(smc.h_alpha_ew[0]/cf.h_alpha_ew[0]), ps=4, xsty=3, ysty=3, /xlog
    plot, cf.h_alpha_ew[0], alog10(smc.h_alpha[0]/cf.h_alpha[0]), ps=4, xsty=3, ysty=3, /xlog

    plothist, cf.continuum_ebv[0]-mw.continuum_ebv[0], bin=0.05
    plot, cf.continuum_ebv[0], mw.continuum_ebv[0], ps=4, xsty=3, ysty=3
    plot, cf.continuum_ebv[0], mw.continuum_ebv[0]-cf.continuum_ebv[0], ps=4, xsty=3, ysty=3
    plot, cf.h_alpha_ew[0], alog10(mw.h_alpha_ew[0]/cf.h_alpha_ew[0]), ps=4, xsty=3, ysty=3, /xlog
    plot, cf.h_alpha_ew[0], alog10(mw.h_alpha[0]/cf.h_alpha[0]), ps=4, xsty=3, ysty=3, /xlog
    
    plothist, smc.continuum_ebv[0]-mw.continuum_ebv[0], bin=0.05
    plot, smc.continuum_ebv[0], mw.continuum_ebv[0], ps=4, xsty=3, ysty=3
    plot, smc.continuum_ebv[0], mw.continuum_ebv[0]-smc.continuum_ebv[0], ps=4, xsty=3, ysty=3
    plot, smc.h_alpha_ew[0], alog10(mw.h_alpha_ew[0]/smc.h_alpha_ew[0]), ps=4, xsty=3, ysty=3, /xlog
    plot, smc.h_alpha_ew[0], alog10(mw.h_alpha[0]/smc.h_alpha[0]), ps=4, xsty=3, ysty=3, /xlog
    
return
end    
