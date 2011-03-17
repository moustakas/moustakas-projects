pro deep2_plot_kcorrect_sedfits, deep1
; jm08may25nyu - based on earlier code

    vname = 'default.nolines'
    light = 2.99792458D18       ; speed of light [A/s]
    scale = 1D18
    
    if (n_elements(deep1) eq 0L) then $
      deep = read_deep2(/kcorr) else $
      deep = deep1
    ngalaxy = n_elements(deep)

    filters = [$
      'deep_B',$
      'deep_R',$
      'deep_I' $
      ]+'.par'
    filtinfo = im_filterspecs(filterlist=filters,/verbose)
    nfilter = n_elements(filtinfo)

    k_load_vmatrix, vmatrix, lambda, vfile=vfile, $
      vpath=vpath, lfile=lfile, vname=vname

    im_window, 0, xratio=0.7, yratio=0.5
    
    plotsym, 0, 1.5, /fill
    for ii = 0L, ngalaxy-1L do begin

       if (deep[ii].mass gt 0.0) then begin

          wave = lambda*(1.0+deep[ii].z)
          restflux = scale*vmatrix#deep[ii].coeffs ; f_lambda
          flux_flam = restflux/(1.0+deep[ii].z)
          flux_fnu = flux_flam*wave^2.0/light*10.0^(0.4*48.6) ; f_nu
          flux_ab = -2.5*alog10(flux_fnu>1D-60)

          good = where((deep[ii].abmaggies gt 0.0),ngood)
          if (ngood ne 0L) then begin

             weff = filtinfo[good].weff
             fnu2flam = 10^(-0.4*48.6)*light/weff^2.0

             abmaggies = deep[ii].abmaggies[good]
             abmaggies_ivar = deep[ii].abmaggies_ivar[good]
             maggies_flam = scale*abmaggies*fnu2flam
             maggies_flam_ivar = abmaggies_ivar/fnu2flam^2.0/scale^2.0

             xrange = [min(filtinfo[good].weff-1.5*filtinfo[good].fwhm),$
               max(filtinfo[good].weff+1.5*filtinfo[good].fwhm)] ;/(1.0+zobj)
             xrange[0] = xrange[0]>90.0
             get_element, wave, xrange, xx

             yrange = minmax(flux_flam[xx[0]:xx[1]])
;            yrange = [min(mab-1.5*mab_err),max(mab+1.5*mab_err)]
;            yrange[0] = yrange[0]<min(flux_ab[xx[0]:xx[1]])
;            yrange[1] = yrange[1]>max(flux_ab[xx[0]:xx[1]])
;            yrange = reverse(yrange)
             
             djs_plot, [0], [0], /nodata, xsty=1, ysty=1, charsize=2.0, charthick=2.0, xthick=2.0, $
               ythick=2.0, xtitle='Observed Wavelength (\AA)', $
               ytitle=textoidl('f_{\lambda} (10^{-18} '+flam_units()+')'), $
               yrange=yrange, xrange=xrange, title='DEEP2/'+$
               string(deep[ii].objno,format='(I0)')+', z = '+$
               strtrim(string(deep[ii].z,format='(F12.3)'),2)
             oplot, wave, flux_flam, line=0.1
             oploterror, filtinfo[good].weff, maggies_flam, filtinfo[good].fwhm, $
               1.0/sqrt(maggies_flam_ivar), ps=8, color=djs_icolor('dark green'), $
               errcolor=djs_icolor('dark green'), thick=2.0, errthick=2.0

;            label = ['M = '+strtrim(string(alog10(deep[ii].mass),format='(F12.2)'),2)+' M'+sunsymbol(),$
;              '\chi^{2}_{\nu} = '+strtrim(string(deep[ii].chi2,format='(F12.2)'),2)]
;            legend, textoidl(label), /left, /top, box=0, charsize=1.6, charthick=2.0
             cc = get_kbrd(1)
          
          endif

       endif
          
    endfor 

return
end
    
