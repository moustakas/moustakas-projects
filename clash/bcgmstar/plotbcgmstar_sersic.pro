function read_kormendy
; read and parse the Kormendy+09 data

    nrr = 1000
    modelrr = [0D,range(1E-4,300,nrr,/log)] ; equivalent radius [kpc]

    korm = rsex(bcgmstar_path(/propath)+'kormendy09.dat')
    ngal = n_elements(korm)
    korm = struct_addtags(korm,replicate({sbe_mean: 0.0},ngal))

    arcsec2kpc = dangular(70*korm.dist/im_light(/kms),/kpc)/206265D
    re_arcsec = 10D^korm.logre/arcsec2kpc
    
    korm.sbe_mean = korm.mue+2.5*alog10(2.0*!pi*re_arcsec^2)
    
return, korm    
end

pro plotbcgmstar_sersic, pdf=pdf, build_models=build_models, clobber=clobber
; jm14aug10siena - make plots of our Sersic fitting results

    sample = read_bcgmstar_sample()
    ncl = n_elements(sample)
    
    paperpath = bcgmstar_path(/paper)
    sersicpath = bcgmstar_path(/sersic)

    filt = bcgmstar_filterlist(weff=weff)

    korm = read_kormendy()

stop    
    
; ---------------------------------------------------------------------------
; fundamental plane
    re = fltarr(ncl)
    mue_f814w = fltarr(ncl)
    mue_f160w = fltarr(ncl)
    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       ser = read_bcgmstar_sersic(cluster,phot=phot,band=['f814w','f160w'])
       re[ic] = phot[0].re_mean
       mue_f814w[ic] = -2.5*alog10(phot[0].sbe_mean)
       mue_f160w[ic] = -2.5*alog10(phot[1].sbe_mean)
    endfor
    niceprint, sample.shortname, re, mue_f814w

    xrange = [18,28]
    yrange = [0,3]
    
    psfile = paperpath+'bcg_fp.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0

    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, $
      ytitle='log_{10} (r_{e} / h^{-1} kpc)', xtitle='<\mu_{e}> (mag arcsec^{-2})'
    djs_oplot, mue_f814w, alog10(re), psym=symcat(16), $
      symsize=2, color=cgcolor('orange')
;   djs_oplot, mue_f160w, re, psym=symcat(15), symsize=2, color=cgcolor('forest green')

; overplot the FP from the literature; Table 4 of Bernardi+03 has a
; nice summary of previous results; generally everyone agrees that
; b=0.8, so assume that here
    mm = range(xrange[0],xrange[1],500)
    cc = im_linefit(alog10(10^(-0.4*mue_f814w)),alog10(re),$
      coeff_guess=[1.0,-0.8],coeff_fixed=[0,1])
    djs_oplot, mm, poly(alog10(10^(-0.4*mm)),cc), line=0
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

stop
    

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
    
    
