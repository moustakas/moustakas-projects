pro plotsings_log12oh_ohmethod, ps=ps
; jm10mar11ucsd - demonstration of how we measure abundances 

    pspath = sings_path(/papers)+'log12oh/FIG_LOG12OH/'
    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

    psfile = pspath+'r23_vs_12oh'+suffix
    im_plotconfig, 7, pos, psfile=psfile, charsize=1.4, $
      xmargin=[0.8,0.8], xspace=0.1, width=[3.4,3.4], $
      height=[3.0,3.0,3.0], thick=6.0

    pos1 = pos[*,2*lindgen(3)]
    pos2 = pos[*,2*lindgen(3)+1]

    ohrange1 = [7.2,9.5]
    ohrange2 = [6.6,9.7]

    logq = alog10([5D6,4E7,1.5D8]) ; alog10(4E7)
    logr23 = findgen((1.1-(-0.5))/0.001+1)*0.001-0.5
    linestyle = [0,1,2]
    
; define the R23 and O32 values    
    r23obs = [4.0,10.0,15.0]
    r23obs_err = [1.0,2.0,4.0]
    o32obs = [0.5,0.6,0.7]
    o32obs_err = [0.05,0.05,0.1]

    logr23obs = alog10(r23obs)
    logo32obs = alog10(o32obs)
    logr23obs_err = r23obs_err/r23obs/alog(10.0)
    logo32obs_err = o32obs_err/o32obs/alog(10.0)

    line = replicate({oiii: [0.0,0.0], oiii: [0.0,0.0], hbeta: [1.0,0.0]},3)
fix this
    line.oiii[0] = r23
    line.oiii[0] = r23
    kk04 = monte_log12oh_kk04(r23obs,r23obs_err,o32obs,o32obs_err,$
      nmonte=1000,result_monte=kk04_monte,seed=100.0)
    help, kk04, /str

    for jj = 0L, 2L do begin ; loop on each value
    
; plot the R23-O/H relation and the KK04 theoretical calibration

;      djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
;        position=pos1[*,jj], xsty=4, ysty=4
;      polyfill, [logr23obs-logr23obs_err,logr23obs+logr23obs_err,logr23obs+logr23obs_err,logr23obs-logr23obs_err], $
;        [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], $
;        noclip=0, clip=[!x.crange[0],!y.crange[0],!x.crange[1],!y.crange[1]], $
;        /fill, color=fsc_color('light grey',20)

       if (jj eq 2L) then begin
          delvarx, xtickname
          xtitle1 = 'log (R_{23})'
          xtitle2 = '12 + log (O/H)'
       endif else begin
          xtickname = replicate(' ',10)
          delvarx, xtitle1, xtitle2
       endelse

       xrange = [-0.3,1.3]
       yrange = ohrange1
       
       djs_plot, [0], [0], /nodata, noerase=(jj gt 0), xrange=xrange, yrange=yrange, $
         xtitle=xtitle1, ytitle='12 + log (O/H)', $
         position=pos1[*,jj], xsty=3, ysty=3, xtickname=xtickname
;      im_legend, '(a)', /left, /top, box=0, charsize=charsize_4, charthick=postthick2, margin=0

       case jj of
          0L: label = '(O/H)_{upper}>(O/H)_{lower}'
          1L: label = '(O/H)_{upper}<(O/H)_{lower}'
          2L: label = '(O/H)_{upper}<<(O/H)_{lower}'
       endcase
       im_legend, label, /right, /top, box=0, charsize=1.1, margin=0
         

       for iq = 0L, n_elements(logq)-1L do begin
          logoh_upper = 9.72D - 0.777*logr23 - 0.951*logr23^2 - 0.072*logr23^3 - 0.811*logr23^4 - $
            logq[iq]*(0.0737 - 0.0713*logr23 - 0.141*logr23^2 + 0.0373*logr23^3 - 0.058*logr23^4)
          logoh_lower = 9.40D + 4.65D*logr23 - 3.17D*logr23^2 - logq[iq]*(0.272D + 0.547D*logr23 - 0.513D*logr23^2)
          good1 = where((logoh_upper gt logoh_lower))
          good2 = where((logoh_lower[good1] gt 7.5))

          djs_oplot, logr23[good1], logoh_upper[good1], linestyle=linestyle[iq]
          djs_oplot, logr23[good1], logoh_lower[good1], linestyle=linestyle[iq]
;         djs_oplot, logr23[good1[good2]], logoh_lower[good1[good2]], line=0
       endfor

; plot the best abundance
       plotsym, 0, 1.6, fill=1, color=im_color('red',10)
       oploterror, logr23obs[jj], kk04[jj].log12oh_upper, logr23obs_err[jj], $
         kk04[jj].log12oh_upper_err, psym=8, errthick=errthick1, errcolor=im_color('red',10)
       plotsym, 4, 1.8, fill=1, color=im_color('blue',11)
       oploterror, logr23obs[jj], kk04[jj].log12oh_lower, logr23obs_err[jj], $
         kk04[jj].log12oh_lower_err, psym=8, errthick=errthick1, errcolor=im_color('blue',11)

; now plot the distribution of O/H values
       
       binsize = 0.05
       
       im_plothist, kk04_monte[jj].log12oh_lower, bin=binsize, xbin_lower, ybin_lower, /noplot
       im_plothist, kk04_monte[jj].log12oh_upper, bin=binsize, xbin_upper, ybin_upper, /noplot
       nbin_lower = n_elements(xbin_lower)
       nbin_upper = n_elements(xbin_upper)
       
       yrange = [0.0,max(ybin_lower)>max(ybin_upper)]*1.05
       xrange = ohrange2

       djs_plot, [0], [0], /nodata, /noerase, yrange=yrange, $
         xrange=xrange, position=pos2[*,jj], xsty=5, ysty=5

; shade where the two histograms overlap

       xar = xbin_lower
       yar = fltarr(nbin_lower)
       for i = 0L, nbin_lower-1L do begin
          ix = where((xbin_upper gt xbin_lower[i]-binsize/2.0) and (xbin_upper lt xbin_lower[i]+binsize/2.0))
          if (ix[0] ne -1L) then yar[i] = min([ybin_lower[i],ybin_upper[ix]])
       endfor

       zix = where(yar eq 0.0)  ; now remove leading/trailing zeros
       shx = zix-shift(zix,-1)
       thx = where(shx ne -1L)
;      xar = xar[zix[thx[0]]+1L:nbin_lower-1L]
;      yar = yar[zix[thx[0]]+1L:nbin_lower-1L]
;      oplot, xar, yar, psym=4, symsize=2 ; do a sanity-check overplot...

       Xfill = transpose([[Xar-binsize/2.0],[Xar+binsize/2.0]])
       Xfill = reform(Xfill, n_elements(Xfill))
       Xfill = [Xfill[0], Xfill, Xfill[n_elements(Xfill)-1]]
       Yfill = transpose([[Yar],[Yar]])
       Yfill = reform(Yfill, n_elements(Yfill))
       Yfill = [0, Yfill, 0]
       Xfill = Xfill > !X.CRANGE[0] < !X.CRANGE[1]
       Yfill = Yfill > !Y.CRANGE[0] < !Y.CRANGE[1]

       polyfill, xfill, yfill, color=fsc_color('light grey',5)

; now plot the data and the axes
       im_plothist, kk04_monte[jj].log12oh_lower, bin=binsize, /overplot, $
         linestyle=2
       im_plothist, kk04_monte[jj].log12oh_upper, bin=binsize, /overplot, $
         linestyle=0

       djs_plot, [0], [0], /nodata, /noerase, yrange=yrange, xrange=xrange, $
         xtitle=xtitle2, ytitle='', ytickname=replicate(' ',10), $
         position=pos2[*,jj], xsty=1, ysty=1, xtickname=xtickname
       axis, /yaxis, ysty=1, yrange=yrange, ytitle='Number'

; overplot a Gaussian fit
       case jj of
          0L: begin
             x1 = kk04[jj].log12oh_upper-5.0*kk04[jj].log12oh_upper_err
             x2 = kk04[jj].log12oh_upper+5.0*kk04[jj].log12oh_upper_err
             xaxis = findgen((x2-x1)/0.005+1)*0.005+x1
             djs_oplot, xaxis, gauss1(xaxis,[kk04[jj].log12oh_upper,$
               kk04[jj].log12oh_upper_err,max(ybin_upper)],/peak), $
               color=djs_icolor('dark green')
          end
          1L: begin
             x1 = kk04[jj].log12oh_avg-5.0*kk04[jj].log12oh_avg_err
             x2 = kk04[jj].log12oh_avg+5.0*kk04[jj].log12oh_avg_err
             xaxis = findgen((x2-x1)/0.005+1)*0.005+x1
             djs_oplot, xaxis, gauss1(xaxis,[kk04[jj].log12oh_avg,$
               kk04[jj].log12oh_avg_err,max(yfill)],/peak), $
               color=djs_icolor('dark green')
          end
          else: 
       endcase

; label the best-fit values       

       case jj of
          0L: final = ['(O/H)='+strtrim(string(kk04[jj].log12oh_upper,format='(F12.2)'),2)+$
            '\pm'+strtrim(string(kk04[jj].log12oh_upper_err,format='(F12.2)'),2),$
            'log(U)='+strtrim(string(kk04[jj].logu_upper,format='(F12.2)'),2)+$
            '\pm'+strtrim(string(kk04[jj].logu_upper_err,format='(F12.2)'),2)]
          1L: final = ['(O/H)='+strtrim(string(kk04[jj].log12oh_avg,format='(F12.2)'),2)+$
            '\pm'+strtrim(string(kk04[jj].log12oh_avg_err,format='(F12.2)'),2),$
            'log(U)='+strtrim(string(kk04[jj].logu_avg,format='(F12.2)'),2)+$
            '\pm'+strtrim(string(kk04[jj].logu_avg_err,format='(F12.2)'),2)]
          2L: final = ['No (O/H) possible','No log(U) possible']
       endcase
       
       im_legend, final, /left, /top, box=0, charsize=1.1, margin=0
         

;      legend, '(b)', /left, /top, box=0, charsize=charsize_4, charthick=postthick2, margin=0
;      legend, ['Upper','Lower'], /left, /top, box=0, charsize=1.2, $
;        charthick=postthick2, thick=postthick1, line=[0,2], margin=0
       
; add the y-titles

;      xyouts, pos1[0,0]-0.08, pos[1,0], textoidl('12 + log (O/H)'), charsize=charsize_4, $
;        align=0.5, orientation=90, /normal, charthick=postthick2
;      xyouts, pos1[0,1]-0.08, pos1[1,1], textoidl('12 + log (O/H)'), charsize=charsize_4, $
;        align=0.5, orientation=90, /normal, charthick=postthick2

;      xyouts, pos2[2,0]+0.08, pos[1,0], textoidl('Number of Realizations'), charsize=charsize_4, $
;        align=0.5, orientation=90, /normal, charthick=postthick2
;      xyouts, pos2[2,1]+0.08, pos1[1,1], textoidl('Number of Realizations'), charsize=charsize_4, $
;        align=0.5, orientation=90, /normal, charthick=postthick2

    endfor
       
    im_plotconfig, /psclose

return
end
