;+
; NAME:
;   PLOT_COSMICIMF_LF24
;
; PURPOSE:
;   Plot the mass function results from COSMICIMF_LF24.
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Feb 05, UCSD
;
; Copyright (C) 2010, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

pro plot_cosmicimf_lf24, ps=ps, keynote=keynote
; jm10mar25ucsd - plot the L(24) LFs in AGES

    cosmicimfpath = ages_path(/projects)+'cosmicimf/'
    paperpath = ages_path(/papers)+'cosmicimf/'
    if keyword_set(keynote) then paperpath = $
      getenv('RESEARCHPATH')+'/talks/2010/10apr_florida/keynote/'
    zbins = cosmicimf_zbins(nzbins,/lf24)

    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

    lf24_data = mrdfits(cosmicimfpath+'lf24_data.fits.gz',1)
    lf24_fit = mrdfits(cosmicimfpath+'lf24_fit.fits.gz',1)

    l24axis = im_array(6.0,15.0,0.02)

; --------------------------------------------------
; LFs in four redshift bins
    psfile = paperpath+'lf24'+suffix
    im_plotconfig, 5, pos, psfile=psfile, charsize=1.6, $
      height=3.0*[1,1], keynote=keynote

    if keyword_set(keynote) then keycolor = djs_icolor('white')

    shupe = lf24_shupe()
    wiphu = lf24_wiphu()
    wlimit = [9.0,9.7,10.2,10.55]
    for jj = 0, nzbins-1 do begin
       if odd(jj) then ytickname = replicate(' ',10) else $
         delvarx, ytickname 
       if (jj lt 2) then xtickname = replicate(' ',10) else $
         delvarx, xtickname 
       djs_plot, [0], [0], /nodata, noerase=(jj gt 0), $
         position=pos[*,jj], xsty=1, ysty=1, /ylog, $
         yrange=[1E-6,0.009], xrange=[8.8,11.99], xtickname=xtickname, $
         ytickname=ytickname, xtitle='', ytitle='', $
         xtickinterval=1.0
       legend, 'z='+strtrim(string(zbins[jj].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[jj].zup,format='(F12.2)'),2), $
         /right, /top, box=0, margin=0, charsize=1.5

; plot the LF_data, distinguishing between objects that are above and
; below the completeness limit
       these = where((lf24_data[jj].phi gt 0.0),nthese)
       fullbin = lf24_data[jj].fullbin[these]
       l24 = lf24_data[jj].l24[these]
       phi = lf24_data[jj].phi[these]
       phierr = lf24_data[jj].phierr[these]
       phierr_cv = lf24_data[jj].phierr_cv[these]
       
       above = where(fullbin eq 1)
       below = where((fullbin eq 0) and lf24_data[jj].number[these] gt 1,nbelow)
;      oploterror, l24[above], phi[above], phierr[above], $
;        psym=symcat(9,thick=!p.thick), symsize=1.2
;      oploterror, l24[below], phi[below], phierr[below], $
;        psym=symcat(9,thick=!p.thick), symsize=1.2, color=djs_icolor('grey'), $
;        errcolor=djs_icolor('grey')

; overplot the local LF
;      djs_oplot, l24axis, lf_double_powerlaw(l24axis,wiphu[0].phistar,$
;        wiphu[0].lstar,wiphu[0].alpha,wiphu[0].beta), $
;        line=5, thick=8, color=fsc_color('dodger blue',102) ; color=fsc_color(wiphu[0].color2,101)
       djs_oplot, l24axis, lf_double_powerlaw(l24axis,shupe), $
         line=5, color=fsc_color('dodger blue',102), thick=8
;        line=5, color=fsc_color(shupe[0].color2,102), thick=6

; overplot my best-fit and Wiphu's best fit
       djs_oplot, l24axis, lf_double_powerlaw(l24axis,lf24_fit[jj]), $
         line=0, color=fsc_color('tomato',101), thick=8
;      djs_oplot, l24axis, lf_double_powerlaw(l24axis,wiphu[jj+1].phistar,$ ; offset from local
;        wiphu[jj+1].lstar,wiphu[jj+1].alpha,wiphu[jj+1].beta), $
;        line=5, color=fsc_color(wiphu[0].color1,101)

       if (jj eq 2) then begin
;         im_legend, 'z=0', /left, /bottom, box=0, line=5, $
;           thick=8, color='dodger blue', margin=0, charsize=1.4, $
;           pspacing=1.4, textcolor=keycolor
         im_legend, ['Shupe+98 (z=0)','Rujopakarn+10'], /left, /bottom, $
           box=0, line=[5,0], thick=8, color=['dodger blue','tomato'], $
           margin=0, charsize=1.2, pspacing=1.4, textcolor=keycolor, $
           charthick=3.0
       endif
       
; overplot Wiphu's data
       if (jj eq 0) then begin
;         wlum = im_array(8.1,10.7,0.2)
;         wphi = [-2.44,-2.44,-2.49,-2.34,-2.47,-2.45,$ ; with AGN
;           -2.58,-2.69,-2.92,-3.08,-3.35,-3.72,-4.15,-4.33]
          wlum = [7.9,8.1,8.3,8.5,8.7,8.9,9.1,9.3,9.5,9.7,$
            9.9,10.1,10.3,10.5,10.7]
          wphierr = [-2.95,-3.07,-3.23,-3.40,-3.44,-3.63,$
            -3.68,-3.76,-3.80,-3.93,-4.02,-4.18,-4.39,-4.58,-4.78]
          wphi = [-2.55,-2.44,-2.45,-2.49,-2.34,-2.48,-2.46,$
            -2.59,-2.71,-2.94,-3.11,-3.44,-3.85,-4.23,-4.63]
       endif
       if (jj eq 1) then begin
;         wlum = im_array(9.1,11.3,0.2)
;         wphi = [-2.88,-2.68,-2.79,-2.89,-2.97,-3.14,$
;           -3.42,-3.70,-4.04,-4.40,-4.51,-5.21]
          wlum = [9.1,9.3,9.5,9.7,9.9,10.1,10.3,10.5,10.7,10.9,11.1]
          wphierr = [-3.65,-3.80,-4.02,-4.12,-4.22,-4.34,-4.48,-4.63,-4.85,-5.12,-5.21]
          wphi = [-2.88,-2.69,-2.80,-2.89,-2.99,-3.17,$
            -3.45,-3.74,-4.19,-4.74,-4.91]
       endif
       if (jj eq 2) then begin
;         wlum = im_array(9.7,11.3,0.2)
;         wphi = [ -3.03,-2.88,-3.09,-3.15,-3.45,-3.81,$
;           -4.15,-4.64,-4.80]
          wlum = [9.7,9.9,10.1,10.3,10.5,10.7,10.9,11.1,11.3]
          wphierr = [-3.79,-3.94,-4.23,-3.89,-4.37,-4.84,-5.07,-5.34,-5.46]
          wphi = [-3.03,-2.89,-3.10,-3.19,-3.52,-3.89,-4.33,-4.86,-5.12]
       endif
       if (jj eq 3) then begin
;         wlum = im_array(10.3,11.9,0.2)
;         wphi = [-3.60,-3.49,-3.59,-3.94,-4.23,-4.48,-4.97,-4.94,-5.53]
          wlum = [10.3,10.5,10.7,10.9,11.1,11.3,11.5,11.7]
          wphierr = [-4.50,-4.60,-4.53,-4.87,-5.08,-5.15,-5.58,-5.77]
          wphi = [-3.60,-3.49,-3.63,-4.00,-4.33,-4.53,-5.16,-5.53]
       endif
       wphi = 10^wphi
       wphierr = 10^wphierr
       niceprint, wlum, wphi, wphi & print
       above = where(wlum gt wlimit[jj],comp=below)
       oploterror, wlum[above], wphi[above], wphierr[above], $
         psym=symcat(6,thick=8), symsize=1.2, color=keycolor, $
         errcolor=keycolor, errthick=8
       oploterror, wlum[below], wphi[below], wphierr[below], $
         psym=symcat(6,thick=5), symsize=1.2, color=djs_icolor('grey'), $
         errcolor=djs_icolor('grey'), errthick=5
    endfor
    xyouts, pos[0,0]-0.1, pos[1,0], textoidl('\Phi(L_{24}/h_{70}^3 Mpc^{-3} dex^{-1})'), $
      align=0.5, orientation=90, /normal
    xyouts, pos[0,3], pos[1,3]-0.1, textoidl('log (L_{24}/L_{\odot})'), $
      align=0.5, /normal
    
    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.eps','.pdf'), /sh
       rmfile, psfile
    endif

return
end
