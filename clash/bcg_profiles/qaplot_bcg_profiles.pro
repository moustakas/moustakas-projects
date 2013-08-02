pro qaplot_bcg_profiles
; jm13may22siena - build a QAplot of all the BCG surface brightness
; profiles 
    
    qapath = getenv('CLASH_DATA')+'/bcg_profiles/'

    filt = clash_filterlist(short=short,weff=weff)
    f814w = where(short eq 'f814w')
    f110w = where(short eq 'f110w')
    f160w = where(short eq 'f160w')

; plot all the clusters; F110W is missing from MACS1423
    clash = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    ncl = n_elements(clash)

    psfile = qapath+'qa_bcg_profiles.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.2
    djs_plot, [0], [0], /nodata, position=pos, $
      xsty=1, ysty=1, xtitle='Semimajor Axis (kpc)', $
      ytitle='F110W - F160W (AB mag)', xrange=[0.5,150], /xlog, $
      yrange=[0,0.7]
;   im_legend, allcl, /right, /top, box=0
    
    sberrmax = 1.0
    for ii = 0, ncl-1 do begin
       data = read_bcg_profiles(clash[ii].shortname)
       if size(data,/type) eq 8 then begin ; good photometry
          ww = where(data[f110w].sma gt -90.0 and data[f160w].sma gt -90.0 and $
            data[f110w].mu_err lt sberrmax and data[f160w].mu_err lt sberrmax,nww)
;         splog, clash[ii].cluster_short, nww
          if nww eq 0 then begin
             splog, 'Insufficient photometry for cluster '+clash[ii].cluster
          endif else begin
             if clash[ii].starforming then begin
                color = 'dodger blue'
                line = 5
             endif else begin
                color = 'firebrick'
                line = 0
             endelse
             ircolor = data[f110w].mu[ww]-data[f160w].mu[ww]
             ircolor_err = sqrt(data[f110w].mu_err[ww]^2+data[f160w].mu_err[ww]^2)
             djs_oplot, data[f160w].sma[ww], ircolor, psym=-8, $
               color=im_color(color), line=line
;            oploterror, data[f160w].sma[ww], ircolor, ircolor_err, psym=8

;; pack the results into a handy structure for writing out
;             out = replicate({sma: 0.0, f110w: 0.0, f110w_err: 0.0, $
;               f160w: 0.0, f160w_err: 0.0},nww)
;;            out = replicate({sma: 0.0, f814w: 0.0, f110w: 0.0, f160w: 0.0},nww)
;             out.sma = data[f160w].sma[ww]
;;            out.f814w = data[f814w].mu[ww]
;             out.f110w = data[f110w].mu[ww]
;             out.f110w_err = data[f110w].mu_err[ww]
;             out.f160w = data[f160w].mu[ww]
;             out.f160w_err = data[f160w].mu_err[ww]
;             wsex, out, outfile=qapath+clash[ii].shortname+'_sbprofiles.txt'
          endelse
       endif
    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf

return
end
    
;    psfile = qapath+'qa_bcg_profiles.ps'
;    im_plotconfig, 0, pos, psfile=psfile, height=5.2
;    for ii = 0, 3 do begin
;;   for ii = 0, ncl-1 do begin
;       splog, allcl[ii]
;       
;       data = read_bcg_profiles(allcl[ii])
;       f160w = where(short eq 'f160w')
;       f110w = where(short eq 'f110w')
;       ww = where(data[f110w].sma gt -90.0 and data[f160w].sma gt -90.0)
;       col = data[f110w].mu[ww]-data[f160w].mu[ww]
;       col_err = sqrt(data[f110w].mu_err[ww]^2+data[f160w].mu_err[ww]^2)
;       yrange = [-2,2]
;       xrange = [0,max(data[f110w].sma)*1.2]
;
;       djs_plot, [0], [0], /nodata, position=pos, $
;         xsty=1, ysty=1, xtitle='Semi-major Axis (arcsec)', $
;         ytitle='\mu', xrange=xrange, yrange=yrange
;       oploterror, data[f160w].sma[ww], col, col_err, psym=8
;       im_legend, allcl[ii], /right, /top, box=0
;    endfor
;    im_plotconfig, psfile=psfile, /psclose, /pdf
