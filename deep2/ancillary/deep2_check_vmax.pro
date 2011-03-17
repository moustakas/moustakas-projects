pro deep2_check_vmax, deep, vmax, postscript=postscript
; jm07nov20nyu - check the completeness predicted by the VMAX code

    path = deep2_path(/analysis)
    pspath = '~/deep/'

    if (n_elements(deep) eq 0L) then deep = mrdfits(path+'deep2.dr3.kcorr.fits.gz',1,/silent)
    if (n_elements(vmax) eq 0L) then $
      vmax = mrdfits(path+'deep2_vmax_0.7_0.9.fits',1,/silent)

; pick redshift and M_B ranges
    
    zlo = 0.7 & zhi = 0.9
;   mb1lo = -18.8 & mb1hi = -19.5
;   mb1lo = -18.5 & mb1hi = -19.0
    mb1lo = -19.0 & mb1hi = -19.3

    mbstr = string(abs(mb1hi),format='(F4.1)')+'_'+string(abs(mb1lo),format='(F4.1)')
    
    if keyword_set(postscript) then begin
       dfpsplot, pspath+'check_vmax_'+mbstr+'.ps', /square, /color
       postthick1 = 4.0
       postthick2 = 3.0
    endif else begin
       im_window, 0, xr=0.4, /square
       postthick1 = 2.0
       postthick2 = 2.0
    endelse

    djs_plot, deep.z, deep.m_b, ps=3, xsty=1, ysty=1, xrange=[0.6,1.5], yrange=[-18,-24], $
      xthick=postthick1, ythick=postthick1, charsize=2.0, charthick=postthick2, xtitle='Redshift', $
      ytitle='M_{B}'
    djs_oplot, zlo*[1,1], [mb1lo,mb1hi], line=0, thick=postthick1, color='red'
    djs_oplot, zhi*[1,1], [mb1lo,mb1hi], line=0, thick=postthick1, color='red'
    djs_oplot, [zlo,zhi], mb1lo*[1,1],   line=0, thick=postthick1, color='red'
    djs_oplot, [zlo,zhi], mb1hi*[1,1],   line=0, thick=postthick1, color='red'
 
    if (not keyword_set(postscript)) then cc = get_kbrd(1)

    indx = where((deep.z gt zlo) and (deep.z lt zhi) and (deep.m_b lt mb1lo) and $
      (deep.m_b gt mb1hi) and (vmax.vmax gt 0.0))

    z = deep[indx].z
    mb = deep[indx].m_b
    bmr = vmax[indx].magb-vmax[indx].magr
    rmi = vmax[indx].magr-vmax[indx].magi
    rmag = vmax[indx].magr

    bmr_model = vmax[indx].bmr_stats[3]
    rmi_model = vmax[indx].rmi_stats[3]
    rmag_model = vmax[indx].rmag_stats[3]
    
    bmr_model_err  = ( vmax[indx].bmr_stats[6]- vmax[indx].bmr_stats[0])/2.0
    rmi_model_err  = ( vmax[indx].rmi_stats[6]- vmax[indx].rmi_stats[0])/2.0
    rmag_model_err = (vmax[indx].rmag_stats[6]-vmax[indx].rmag_stats[0])/2.0
    
; check the model vs galaxy colors

    rmicut = findgen((2.0-(-2.0))/0.001)*0.001+(-2.0)
    bmrcut = poly(rmicut,[-0.25,2.35])

    plotsym, 0, 0.5, /fill
    
    djs_plot, rmi, bmr, xrange=[-0.4,1.2], yrange=[-0.3,2.0], xsty=1, ysty=1, $
      xthick=postthick1, ythick=postthick1, charsize=2.0, charthick=postthick2, $
      xtitle='R - I', ytitle='B - R', psym=8, color='dark red'
;   hogg_scatterplot, rmi, bmr, /nogrey, levels=errorf(0.5*[1.0,2.0,3.0]), $
;     xrange=[-0.3,1.2], yrange=[-0.3,1.9], xsty=1, ysty=1, xthick=postthick1, ythick=postthick1, $
;     charsize=2.0, charthick=postthick2, xtitle='R - I', ytitle='B - R'
    djs_oplot, rmi_model, bmr_model, psym=8, color=fsc_color('dodger blue',10), sym=0.7
;   oploterror, rmi_model, bmr_model, rmi_model_err, bmr_model_err, psym=8, nskip=5, $
;     errcolor=fsc_color('powder blue',10), errthick=1.5, color=fsc_color('powder blue',10)
    djs_oplot, [!x.crange[0],interpol(rmicut,bmrcut,0.5)], 0.5*[1,1], line=0, thick=postthick1
;   djs_oplot, 1.15*[1,1], !y.crange, line=0, thick=postthick1
    djs_oplot, [interpol(rmicut,bmrcut,0.5),max(rmicut)], [0.5,max(bmrcut)], line=0, thick=postthick1
    
;   djs_plot, [0], [0], /nodata, xrange=[-0.5,1.5], yrange=[-0.5,2.5], $
;     xsty=1, ysty=1, xthick=postthick1, ythick=postthick1, $
;     charsize=2.0, charthick=postthick2, xtitle='R - I', ytitle='B - R'
;   djs_oplot, rmi, bmr, psym=8
;   djs_oplot, rmi_model, bmr_model, psym=8, color='green'

    legend, textoidl([string(zlo,format='(F4.2)')+' < z < '+string(zhi,format='(F4.2)'),$
      string(mb1hi,format='(F5.1)')+' < M_{B} < '+string(mb1lo,format='(F5.1)')]), $
      /left, /top, box=0, charsize=2.0, charthick=postthick2

    if (not keyword_set(postscript)) then cc = get_kbrd(1)

; check the completeness model    
    
    zbinsize = 0.01 ; from IM_DEEP2_VMAX
    nzbins = fix((zhi-zlo)/zbinsize)+1L

    zarray = findgen(nzbins)*zbinsize+zlo
    model = fltarr(nzbins)
    model_err = fltarr(nzbins)
    for iz = 0L, nzbins-1L do begin
       good = where(vmax[indx].completeness[iz] gt -1.0,ngood)
       if (ngood ne 0L) then begin
          resistant_mean, vmax[indx[good]].completeness[iz], 3.0, mn, sig
          model[iz] = mn
          model_err[iz] = robust_sigma(vmax[indx[good]].completeness[iz]) ; sig
       endif
    endfor
    
    zbinsize2 = 0.0025

    im_plothist, z, binsize=zbinsize2, xhist, zhist, histmin=zlo, histmax=zhi, binedge=0, /noplot
    yrange = [-2,max(zhist)] ; minmax(zhist/total(zhist))
    djs_plot, [0], [0], /nodata, yrange=yrange, xrange=[zlo,zhi], xsty=3, ysty=1, xthick=postthick1, ythick=postthick1, $
      charsize=2.0, charthick=postthick2, xtitle='Redshift', ytitle='Number of Galaxies'
    djs_oplot, xhist, zhist, ps=10, color=djs_icolor('red'), thick=postthick1, line=1
    factor = mean(zhist[where((xhist gt 0.70) and (xhist lt 0.76))])/max(model)
    oploterror, zarray, model*factor, model_err*factor, ps=10, thick=3, line=0, $
      errthick=1.5

    legend, textoidl([string(zlo,format='(F4.2)')+' < z < '+string(zhi,format='(F4.2)'),$
      string(mb1hi,format='(F5.1)')+' < M_{B} < '+string(mb1lo,format='(F5.1)')]), $
      /right, /top, box=0, charsize=2.0, charthick=postthick2
    
    if keyword_set(postscript) then dfpsclose

stop    

return
end
    
