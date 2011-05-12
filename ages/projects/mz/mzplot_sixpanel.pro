pro mzplot_sixpanel, zobj, xall, yall, weightall, oh_err=oh_err, $
  xtitle=xtitle1, ytitle=ytitle1, title=title1, psfile=psfile, pos=pos, $
  _extra=extra, mzlocal=mzlocal, lzlocal=lzlocal, mzevol=mzevol, $
  lzevol=lzevol, localline=localline, localcolor=localcolor, $
  evolline=evolline, evolcolor=evolcolor, postscript=postscript, $
  mztest=mztest
; jm10oct14ucsd - simple wrapper to make a 6-panel plot for each 
; redshift bin that gets used repeatedly 

    massaxis = range(8.4,11.4,500)
    im_plotconfig, 7, pos, psfile=psfile, charsize=1.6, $
      height=2.7*[1,1,1]

    zbins = mz_zbins(nzbins)
    for iz = 0, nzbins-1 do begin
       these = where((zobj ge zbins[iz].zlo) and $
         (zobj lt zbins[iz].zup),nthese)
       xx = xall[these]
       yy = yall[these]
       ww = weightall[these]
       ee = oh_err[these]
       
       if odd(iz) then begin
          ytitle = '' & ytickname = replicate(' ',10)
       endif else begin
          ytitle = ytitle1 & delvarx, ytickname
       endelse
       if (iz lt 4) then begin
          xtitle = '' & xtickname = replicate(' ',10)
       endif else delvarx, xtickname 

       if (nthese gt 20) then begin
          mzplot_scatterplot, xx, yy, weight=ww, noerase=(iz gt 0), $
            position=pos[*,iz], xrange=xrange, yrange=yrange, xsty=1, ysty=1, $
            xtitle=xtitle, ytitle=ytitle, xtickname=xtickname, $
            ytickname=ytickname, _extra=extra, ccolor=djs_icolor('grey'), /nogrey
;         abin = im_medxbin(xx,yy,0.15,weight=ww/ee^2,/verbose,minpts=3)
;         oploterror, abin.xbin, abin.meany, abin.sigymean, psym=symcat(9,thick=5), $
;           symsize=1.1, thick=6, color=fsc_color('blue',101), $
;           errcolor=fsc_color('blue',101)
;         djs_oplot, mm.xbin, mm.medy, psym=6, symsize=0.5
       endif else begin
          djs_plot, [0], [0], /nodata, noerase=(iz gt 0), $
            position=pos[*,iz], xrange=xrange, yrange=yrange, xsty=1, ysty=1, $
            xtitle=xtitle, ytitle=ytitle, xtickname=xtickname, $
            ytickname=ytickname, _extra=extra
       endelse
       legend, 'z='+strtrim(string(zbins[iz].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[iz].zup,format='(F12.2)'),2), $
         /left, /top, box=0, margin=0, charsize=1.4

; overplot various relations, if desired

; ----------
; MZ relation
       if (n_elements(mzlocal) ne 0) then begin
          ohlocal = mz_closedbox(massaxis,mzlocal.coeff)
          good = where(ohlocal gt !y.crange[0]+0.07)
          djs_oplot, massaxis[good], ohlocal[good], $
            line=localline, color=localcolor, thick=8
; overplot the evolutionary model derived in FIT_MZLZEVOL (see also
; MZPLOT_OHEVOL)
          if (iz gt 0) then begin
             ohmodel = ohlocal + (zbins[iz].zbin-mzevol.qz0)*$
               poly(massaxis-mzevol.dlogohdz_normmass,mzevol.dlogohdz_coeff)
             keep = where((ohmodel gt !y.crange[0]+0.07))
             djs_oplot, massaxis[keep], ohmodel[keep], line=evolline, $
               color=fsc_color(evolcolor,101), thick=8
          endif
       endif
       
;       if (n_elements(mzevol) ne 0) then begin
;          if (iz gt 0) then begin
;;; free R0
;;             ohmodel = mzevol_func(massaxis,mzevol.coeffs,$
;;               z=zbins[iz].zbin,qz0=mzevol.qz0)
;;             good = where(ohmodel gt !y.crange[0]+0.07)
;;             djs_oplot, massaxis[good], ohmodel[good], line=3, $
;;               color=fsc_color('navy',101), thick=8
;; R0=0
;             ohmodel = mzevol_func(massaxis,mzevol.coeffs_r0zero,$
;               z=zbins[iz].zbin,qz0=mzevol.qz0)
;             good = where(ohmodel gt !y.crange[0]+0.07)
;             djs_oplot, massaxis[good], ohmodel[good], line=5, $
;               color=fsc_color(evolcolor,101), thick=8
;          endif 
;      endif
; ----------
; LZ relation
;      if (n_elements(lzlocal) ne 0) then begin
;         oplot_lzfit, lzlocal.coeff, band='B', linestyle=localline, $
;           linecolor=localcolor
;      endif
       if (n_elements(lzevol) ne 0) then begin
          oplot_lzfit, lzevol.coeffs[*,0]*[1,1,0,1], band='B', linestyle=localline, $
            linecolor=localcolor, z=zbins[iz].zbin, qz0=lzevol.qz0, /evol 
          if (iz gt 0) then begin
; q0=0.0 mag/z
             oplot_lzfit, lzevol.coeffs[*,0], band='B', linestyle=evolline, $
               linecolor=evolcolor, z=zbins[iz].zbin, qz0=lzevol.qz0, /evol
; q0=1.6 mag/z
             oplot_lzfit, lzevol.coeffs[*,2]*[1,1,0,1], band='B', linestyle=1, $
               linecolor=evolcolor, z=zbins[iz].zbin, qz0=lzevol.qz0, /evol
          endif
       endif
    endfor
;   if (n_elements(ytitle1) ne 0) then xyouts, pos[0,0]-0.1, $
;     pos[1,0], align=0.5, /norm, orientation=90, ytitle1
;   if (n_elements(ytitle1) ne 0) then xyouts, pos[0,2]-0.1, $
;     pos[1,2], align=0.5, /norm, orientation=90, ytitle1

    if (n_elements(ytitle1) ne 0) then xyouts, pos[0,5], $
      pos[1,5]-0.07, align=0.5, /norm, xtitle1
    if (n_elements(title1) ne 0) then xyouts, pos[2,0], $
      pos[3,0]+0.02, align=0.5, /norm, title1
    im_plotconfig, /psclose, psfile=psfile, gzip=keyword_set(postscript)
    
return
end
    
