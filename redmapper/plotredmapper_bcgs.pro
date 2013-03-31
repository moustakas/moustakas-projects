pro plotredmapper_bcgs
; jm13mar29siena - make some plots of the BCG properties
    
    prefix = 'redmapper'
    ver = 'v5.2'
    
    catalogs_dir = getenv('REDMAPPER_DATA')+'/catalogs/'
    isedfit_dir = getenv('REDMAPPER_DATA')+'/'

    bcgs = mrdfits(catalogs_dir+'dr8_run_redmapper_'+ver+'_lgt20_catalog.fits.gz',1)
    ised = mrdfits(isedfit_dir+'redmapper_fsps_chab_charlot_sfhgrid01.fits.gz',1)
    phot = mrdfits(catalogs_dir+'redmapper_'+ver+'_photometry.fits.gz',1)
    phot = phot[where(phot.isbcg)]

; hack!
    spherematch, phot.ra, phot.dec, bcgs.ra, bcgs.dec, 1D/3600.0, m1, m2
    srt = sort(m1)
    bcgs = bcgs[m2[srt]]

; --------------------------------------------------
; SFR vs stellar mass
    djs_plot, ised.mass_50, ised.sfr100, psym=8

    ww = where(total(ised.maggies[0:1] gt 0,1) gt 0 and ised.zobj lt 0.3)
    djs_plot, ised[ww].mass_50, ised[ww].sfr100_50-ised[ww].mass_50+9, psym=8, yr=[-10,1]

    ww = where(total(ised.maggies[0:1] gt 0,1) gt 0 and ised.zobj gt 0.3)
    djs_plot, ised[ww].mass_50, ised[ww].sfr100_50-ised[ww].mass_50+9, psym=8, symsize=0.5, yr=[-10,1]


    djs_plot, ised.mass_50, bcgs.lambda_chisq, psym=8, symsize=0.5, /ylog
    
stop    
    
return
end
    
