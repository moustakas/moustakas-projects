pro plotbcgmstar_sersic, pdf=pdf, build_models=build_models, clobber=clobber
; jm14aug10siena - make plots of our Sersic fitting results

    sample = read_bcgmstar_sample()
    ncl = n_elements(sample)
    
    paperpath = bcgmstar_path(/paper)
    ellpath = bcgmstar_path(/ellipse)
    sersicpath = bcgmstar_path(/sersic)

    filt = bcgmstar_filterlist(weff=weff)
    
; ---------------------------------------------------------------------------
; fundamental plane in I814W
    re = fltarr(ncl)
    Ie = fltarr(ncl)
    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       ser = read_bcgmstar_sersic(cluster,band='f814w')
       re[ic] = ser.sersic_all_n
       Ie[ic] = 10^(-0.4*ser.sersic_all_sbe)
    endfor

stop    
    
    psfile = paperpath+'bcg_fp.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.5

    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=[0.3,2.2], yrange=[0,11], /xlog, $
      xtitle='Wavelength (\mu'+'m)', ytitle='Sersic n'
    im_plotconfig, psfile=psfile, /psclose, /pdf


; ---------------------------------------------------------------------------
; Sersic n vs wavelength (showing the effect of alpha)
    wave = range(min(weff),max(weff),100)
    
    psfile = paperpath+'bcg_sersic.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.5

    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=[0.3,2.2], yrange=[0,11], /xlog, $
      xtitle='Wavelength (\mu'+'m)', ytitle='Sersic n'
    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       sersic = mrdfits(sersicpath+cluster+'-allsersic-results.fits.gz',1,/silent)
       print, cluster, sersic.sersic_chi2/sersic.sersic_dof, format='(A20,F8.3)'
       djs_oplot, weff/1D4, sersic.sersic_n_ref*(weff/sersic.sersic_wave_ref)^sersic.sersic_alpha1, $
         psym=-symcat(16), color=cgcolor('dark grey')
    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf

return
end
    
    
