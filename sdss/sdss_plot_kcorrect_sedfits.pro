pro sdss_plot_kcorrect_sedfits, sdss1
; jm06jun19uofa - written
; jm06sep18nyu - updated

    vname = 'default.nolines'
    
    light = 2.99792458D18       ; speed of light [A/s]

    path = sdss_path()
    if (n_elements(sdss1) eq 0L) then sdss1 = read_sdss()

    srt = sort(sdss1.kcorr_chi2) & sdss = sdss1[srt]
;   sdss = sdss1

;   w = where(sdss1.kcorr_chi2 gt 100.0)
;   struct_print, struct_trimtags(sdss1,select=['sdss_id','ra','dec','kcorr_mass','kcorr_chi2']);, file='/global/data/scr/ioannis/sdss_crummy_photometry.dat'

    ngalaxy = n_elements(sdss)

    filters = [$
      'sdss_u0',$
      'sdss_g0',$
      'sdss_r0',$
      'sdss_i0',$
      'sdss_z0' $ ; low S/N
    ]+'.par'

    filtinfo = im_filterspecs(filterlist=filters,/verbose)
    nfilter = n_elements(filtinfo)

    k_load_vmatrix, vmatrix, lambda, vfile=vfile, vpath=vpath, lfile=lfile, vname=vname
 
    plotsym, 0, 1.5, /fill

    for i = 0L, ngalaxy-1L do begin

       if (sdss[i].kcorr_mass gt 0.0) then begin

          wave = lambda*(1.0+sdss[i].z)
          restflux = vmatrix#sdss[i].kcorr_coeffs ; f_lambda
          flux = restflux/(1.0+sdss[i].z)
          flux_fnu = flux*wave^2.0/light*10.0^(0.4*48.6) ; f_nu

          good = where((sdss[i].kcorr_mobs_ab gt -900.0),ngood)
;         good = where((sdss[i].kcorr_maggies gt 0.0) and (sdss[i].kcorr_maggies_err lt 1E15),ngood)

          if (ngood ne 0L) then begin

             maggies = 10.0^(-0.4*sdss[i].kcorr_mobs_ab[good])
             maggies_err = 0.4*maggies*1.0/sqrt(sdss[i].kcorr_mobs_ab_ivar[good])
             
;            maggies_flam = maggies*10^(-0.4*48.6)*light/rebin(filtinfo.weff,nfilter,ngalaxy)^2
             
             xrange = [min(filtinfo[good].weff-1.5*filtinfo[good].fwhm),$
               max(filtinfo[good].weff+1.5*filtinfo[good].fwhm)] ;/(1.0+zobj)
             xrange[0] = xrange[0]>90.0
             get_element, wave, xrange, xx
          
             yrange = [min(maggies-maggies_err),max(maggies+maggies_err)]
             yrange[0] = yrange[0]<min(flux_fnu[xx[0]:xx[1]])
             yrange[1] = yrange[1]>max(flux_fnu[xx[0]:xx[1]])
             
             plot, [0], [0], /nodata, xsty=3, ysty=3, charsize=1.5, charthick=2.0, xthick=2.0, $
               ythick=2.0, xtitle='Observed Wavelength (\AA)]', ytitle=textoidl('f_{\nu}'), $
               yrange=yrange, xrange=xrange, title=string(sdss[i].sdss_id,format='(I6.6)')+' z = '+$
               strtrim(string(sdss[i].z,format='(F12.3)'),2)
             oplot, wave, flux_fnu, line=0.1
;            oploterror, filtinfo[good].weff, sdss[i].kcorr_maggies[good], filtinfo[good].fwhm, $
;              sdss[i].kcorr_maggies_err[good], ps=8, color=djs_icolor('red'), errcolor=djs_icolor('red'), $
;              thick=2.0, errthick=2.0

             sdssfilt = where(strmatch(filtinfo[good].filter,'*sdss*',/fold),nsdss)
             if (nsdss ne 0L) then $
               oploterror, filtinfo[good[sdssfilt]].weff, maggies[sdssfilt], filtinfo[good[sdssfilt]].fwhm, $
               maggies_err[sdssfilt], ps=8, color=djs_icolor('blue'), errcolor=djs_icolor('blue'), $
               thick=2.0, errthick=2.0

             label = ['M = '+strtrim(string(sdss[i].kcorr_mass,format='(F12.2)'),2)+' M'+sunsymbol(),$
               '\chi^{2}_{\nu} = '+strtrim(string(sdss[i].kcorr_chi2,format='(F12.2)'),2)]
             legend, textoidl(label), /left, /top, box=0, charsize=1.6, charthick=2.0
             cc = get_kbrd(1)
          
          endif

       endif
          
    endfor 

return
end
    
