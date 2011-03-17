pro bcgs_plots, ps=ps
; jm10jul22ucsd - basic plots for the BCGS project with Anthony 

    bcgspath = ages_path(/projects)+'bcgs/'
    paperpath = bcgspath
;   paperpath = ages_path(/papers)+'bcgs/'

;; read the 37 kpc K-band photometry and compute the aperture
;; correction
;    sample = rsex('bcgs_sample_v3.sex')
;    phot = mrdfits(bcgspath+'bcgs_photometry_v3.fits.gz',1)
;    all = djs_readlines(bcgspath+'all.Ks.out')
;    ngal = n_elements(all)
;    ktot = dblarr(ngal)
;    for ii = 0, ngal-1 do ktot[ii] = (strsplit(all[ii],' ',/extract))[38]
;
;    im_plothist, 10^(-0.4*(phot.ks_mag_aper_08-ktot))
    
; read the stellar masses    
    grid1 = mrdfits(bcgspath+'BwRIJHKsirac_bc03_chab_calzetti_sfhgrid01.fits.gz',1)
    grid2 = mrdfits(bcgspath+'BwRIJHKsirac_bc03_chab_calzetti_sfhgrid02.fits.gz',1)
    grid4 = mrdfits(bcgspath+'BwRIJHKsirac_bc03_chab_calzetti_sfhgrid04.fits.gz',1)
;   grid3 = mrdfits('BwRIJHKsirac_bc03_chab_sfhgrid03.fits.gz',1)
    sdss_grid2 = mrdfits(bcgspath+'ugrizJHKs_bc03_chab_calzetti_sfhgrid02.fits.gz',1)

    ps = 1
    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

; --------------------------------------------------
; age vs redshift
    zrange = [0.0,1.85]
    agerange = [0.0,12.0]
    ztitle = 'Redshift'
    agetitle = 'Age (Gyr)'

    zaxis = range(0.0,2.0,100)
    
    psfile = paperpath+'z_vs_age'+suffix
    im_plotconfig, 0, pos, psfile=psfile, height=5.0

    djs_plot, [0], [0], /nodata, position=pos, $
      xrange=zrange, yrange=agerange, xsty=1, ysty=1, $
      xtitle=ztitle, ytitle=agetitle
    oploterror, grid2.zobj, grid2.age_avg, grid2.age_err, $
      psym=symcat(16), color=fsc_color('dodger blue',100), $
      errcolor=fsc_color('dodger blue',100), symsize=1.3, $
      errthick=4
;   oploterror, grid4.zobj, grid4.age_avg, grid4.age_err, $
;     psym=symcat(15), color=fsc_color('firebrick',100), $
;     errcolor=fsc_color('firebrick',100), symsize=1.3, $
;     errthick=4
    oploterror, sdss_grid2.zobj, sdss_grid2.age_avg, sdss_grid2.age_err, $
      psym=symcat(15), color=fsc_color('forest green',100), $
      errcolor=fsc_color('forest green',100), symsize=1.3, $
      errthick=4
;   djs_oplot, sdss_grid2.zobj, sdss_grid2.age_avg, $
;     psym=symcat(14), symsize=1.5
;   djs_oplot, grid1.zobj, grid1.age_avg, psym=7
    djs_oplot, zaxis, getage(zaxis), line=0
    djs_oplot, zaxis, getage(zaxis)-getage(3.0), line=5
    im_legend, ['Bootes','Stott+08'], /right, /top, box=0, $
      color=['dodger blue','forest green'], psym=[16,15]
;   im_legend, ['Z=0.02','Z=0.05'], /right, /top, box=0, $
;     color=['dodger blue','firebrick'], psym=[16,15]

; --------------------------------------------------
; mass vs redshift
    zrange = [0.0,1.85]
    massrange = [10.5,12.3]
    ztitle = 'Redshift'
    masstitle = 'log (M/M_{'+sunsymbol()+'})'

    psfile = paperpath+'z_vs_mass'+suffix
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, xmargin=[1.2,0.3]

    djs_plot, [0], [0], /nodata, position=pos, $
      xrange=zrange, yrange=massrange, xsty=1, ysty=1, $
      xtitle=ztitle, ytitle=masstitle
    oploterror, grid2.zobj, grid2.mass_avg, grid2.mass_err, $
      psym=symcat(16), color=fsc_color('dodger blue',100), $
      errcolor=fsc_color('dodger blue',100), symsize=1.3, $
      errthick=4
;   oploterror, grid4.zobj, grid4.mass_avg, grid4.mass_err, $
;     psym=symcat(15), color=fsc_color('firebrick',100), $
;     errcolor=fsc_color('firebrick',100), symsize=1.3, $
;     errthick=4
    oploterror, sdss_grid2.zobj, sdss_grid2.mass_avg, sdss_grid2.mass_err, $
      psym=symcat(15), color=fsc_color('forest green',100), $
      errcolor=fsc_color('forest green',100), symsize=1.3, $
      errthick=4
;   djs_oplot, sdss_grid2.zobj, sdss_grid2.age_avg, $
;     psym=symcat(14), symsize=1.5
;   im_legend, ['Z=0.02','Z=0.05'], /right, /bottom, box=0, $
;     color=['dodger blue','firebrick'], psym=[16,15]
    im_legend, ['Bootes','Stott+08'], /right, /top, box=0, $
      color=['dodger blue','forest green'], psym=[16,15]

    im_plotconfig, /psclose


;   plot, grid4.chi2, grid2.chi2, psym=6, /xlog, /ylog, xr=[0.1,300], yr=[0.1,300], xsty=1, ysty=1 & oplot, 10^!x.crange, 10^!y.crange
    
    
return
end
    
