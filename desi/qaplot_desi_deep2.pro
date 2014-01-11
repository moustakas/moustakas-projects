pro qaplot_desi_deep2
; jm14jan09siena - build some QAplots for the DESI/DEEP2 sample

    light = im_light(/Ang)
    fluxscale = 1D18

    prefix = 'desi_deep2'
    isedfit_dir = getenv('IM_PROJECTS_DIR')+'/desi/deep2/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

; read the full DEEP2 photometric catalog
    junk = read_deep2_zcat(photo=pcat)
;   pcat = mrdfits(deep2_path(/cat)+'pcat_ext.fits.gz',1)
    deep2_to_maggies, pcat, pmm
    pgood = where(total(pmm[0:2,*] gt 0,1) eq 3)

; read the output of BUILD_DESI_DEEP2_SAMPLE    
    zcat = mrdfits(isedfit_dir+'deep2_zcat.fits.gz',1)
    ispec = mrdfits(isedfit_dir+'deep2_ispec.fits.gz',1)
    ised = read_isedfit(isedfit_paramfile,isedfit_dir=isedfit_dir)
    ngal = n_elements(zcat)

; read my simulation outputs
    sim = mrdfits(isedfit_dir+'simulate_desi_deep2.fits.gz',1)
    wave = mrdfits(isedfit_dir+'simulate_desi_deep2.fits.gz',2)

    weff = k_lambda_eff(filterlist=deep2_filterlist())
    nband = n_elements(weff)
    
; --------------------------------------------------
    psfile = isedfit_dir+'qa_desi_deep2.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.8, $
      xmargin=[1.3,0.4]

; B-R vs R-I    
    im_hogg_scatterplot, -2.5*alog10(pmm[1,pgood]/pmm[2,pgood]), $
      -2.5*alog10(pmm[0,pgood]/pmm[1,pgood]), position=pos, $
      xsty=1, ysty=1, xrange=[-0.3,1.8], yrange=[-0.5,5], /internal, $
      xtitle='R - I', ytitle='B - R', /nogrey, /nooutlier;, xnpix=151, ynpix=151
    im_hogg_scatterplot, -2.5*alog10(ised.bestmaggies[1]/ised.bestmaggies[2]), $
      -2.5*alog10(ised.bestmaggies[0]/ised.bestmaggies[1]), /overplot, /nogrey, $
      contour_color=cgcolor('firebrick'), cthick=8
    im_legend, ['Full DEEP2/DR4 Sample (z=0-2)','[OII] Sample (z=0.7-1.6, S/N>3)'], $
      /right, /top, box=0, color=[cgcolor('black'),cgcolor('firebrick')], $
      line=0, thick=8, pspacing=1.8, charsize=1.5, margin=0

;; redshift vs line-flux
;    im_hogg_scatterplot, ispec.z, ispec.oii_3727_1[0]+ispec.oii_3727_2[0], $
;      position=pos, xsty=1, ysty=1, xrange=[8,12], yrange=[0.5,3], /internal, $
;      xtitle=textoidl('log (M_{*} / M_{\odot})'), outpsym=symcat(16), /outlier, $
;      outsymsize=0.5, outcolor=cgcolor('black'), $
;      ytitle=textoidl('log EW([O II]) (\AA)'), /nogrey
    
; mass vs EW([OII])
    im_hogg_scatterplot, ised.mstar_50, $
      alog10(ispec.oii_3727_1_ew[0]+ispec.oii_3727_2_ew[0]), $
      position=pos, xsty=1, ysty=1, xrange=[8,12], yrange=[0.5,3], /internal, $
      xtitle=textoidl('log (M_{*} / M_{\odot})'), outpsym=symcat(16), /outlier, $
      outsymsize=0.5, outcolor=cgcolor('black'), $
      ytitle=textoidl('log EW([O II]) (\AA)'), /nogrey

; a couple pages of spectra    
    im_plotconfig, 0, pos2, height=3.0

    indx = [0,1,2]
    for ii = 0, n_elements(indx)-1 do begin
       mag = maggies2mag(sim[indx[ii]].maggies,magerr=magerr,$
         ivarmaggies=sim[indx[ii]].ivarmaggies)
       modelmag = maggies2mag(sim[indx[ii]].bestmaggies)

       abflux = -2.5*alog10((sim[indx[ii]].flux*wave.wave^2.0/light)>1D-50)-48.6
       yrange = [max(abflux),min(abflux)]

       djs_plot, [0], [0], /nodata, position=pos2, xsty=1, ysty=1, $
         xrange=[0.2,1.3], yrange=yrange, xtitle='Wavelength (\mu'+'m)', $
         ytitle='AB Magnitude'
       djs_oplot, wave.wave/1D4, abflux, color=cgcolor('grey')
            
       oploterror, weff/1D4, mag, magerr, psym=symcat(16), $
         color=cgcolor('dodger blue'), errcolor=cgcolor('dodger blue'), $
         symsize=1.5
       djs_oplot, weff/1D4, modelmag, psym=symcat(6,thick=4), symsize=2.0

       im_legend, 'z = '+strtrim(string(sim[ii].z,format='(F12.4)'),2), /left, $
         /top, box=0, margin=0

;; make an inset centered on [OII]
;       xr = 3727.0*(1+sim[indx[ii]].z)*[0.99,1.01]
;       get_element, wave.wave, xr, kk
;       yr = minmax(sim[indx[ii]].flux[kk[0]:kk[1]]*fluxscale)
;       
;       djs_plot, [0], [0], /nodata, /noerase, position=[0.5,0.9,0.3,0.8], $
;         xsty=1, ysty=3, xrange=xr/1D4, yrange=yr, /norm
;       djs_oplot, wave.wave/1D4, sim[indx[ii]].flux*fluxscale
    endfor
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    

return
end
    

    
