pro mzplot_sixpanel, zobj, xall, yall, weightall, xtitle=xtitle1, $
  ytitle=ytitle1, title=title1, psfile=psfile, pos=pos, _extra=extra, $
  mzlocal=mzlocal, lzlocal=lzlocal, mzevol=mzevol, lzevol=lzevol, $
  localline=localline, localcolor=localcolor, evolline=evolline, $
  evolcolor=evolcolor
; jm10oct14ucsd - simple wrapper to make a 6-panel plot for each 
; redshift bin that gets used repeatedly 

    massaxis = range(8.4,11.4,500)
    im_plotconfig, 7, pos, psfile=psfile, charsize=1.6, $
      height=2.7*[1,1,1]

    zbins = mz_zbins(nzbins)
    for iz = 0, nzbins-1 do begin
       these = where((zobj gt zbins[iz].zlo) and $
         (zobj lt zbins[iz].zup),nthese)
       xx = xall[these]
       yy = yall[these]
       ww = weightall[these]
       
       if odd(iz) then begin
          ytitle = '' & ytickname = replicate(' ',10)
       endif else begin
          ytitle = ytitle1 & delvarx, ytickname
       endelse
       if (iz lt 4) then begin
          xtitle = '' & xtickname = replicate(' ',10)
       endif else delvarx, xtickname 

       if (nthese gt 40L) then begin
          mzplot_scatterplot, xx, yy, weight=ww, noerase=(iz gt 0), $
            position=pos[*,iz], xrange=xrange, yrange=yrange, xsty=1, ysty=1, $
            xtitle=xtitle, ytitle=ytitle, xtickname=xtickname, $
            ytickname=ytickname, _extra=extra
          mm = im_medxbin(xx,yy,0.1,weight=ww)
          djs_oplot, mm.xbin, mm.medy, psym=6, symsize=0.5
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
;      if (n_elements(mzlocal) ne 0) then begin
;         djs_oplot, massaxis, mz_brokenpl(massaxis,mzlocal.coeff_bin), $
;           line=localline, color=localcolor, thick=8
;      endif
       if (n_elements(mzevol) ne 0) then begin
          djs_oplot, massaxis, mzevol_func(massaxis,mzevol.coeffs[*,0]*[1,1,1,0,0],$
            z=zbins[iz].zbin,qz0=mzevol.qz0), line=localline, $
            color=fsc_color(localcolor,101), thick=8
          if (iz gt 0) then begin
; r0=0.0 dex/z
             djs_oplot, massaxis, mzevol_func(massaxis,mzevol.coeffs[*,0],$
               z=zbins[iz].zbin,qz0=mzevol.qz0), line=evolline, $
               color=fsc_color(evolcolor,101), thick=8
; r0=-0.5 dex/z
             djs_oplot, massaxis, mzevol_func(massaxis,mzevol.coeffs[*,2],$;*[1,1,1,0,1],$
               z=zbins[iz].zbin,qz0=mzevol.qz0), line=1, $
               color=fsc_color(evolcolor,101), thick=8
          endif
       endif
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
    im_plotconfig, /psclose
    
return
end
    
