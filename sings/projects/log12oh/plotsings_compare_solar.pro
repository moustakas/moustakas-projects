pro doplot_compare, data1, data2, title=title
; driving plotting routine; main routine below

    ngalaxy = n_elements(data2)
    snrcut = 0.0 ; 5.0

    pagemaker, nx=4, ny=4, position=pos1, /normal, xmargin=[1.2,0.2], $ 
      ymargin=[0.4,1.2], yspace=0.0, xspace=0.0, xpage=xpage, ypage=ypage
    pagemaker, nx=3, ny=3, position=pos2, /normal, xmargin=[1.1,0.2], $ 
      ymargin=[0.5,1.1], xpage=xpage, ypage=ypage, xspace=0.0, yspace=0.0
    pagemaker, nx=1, ny=3, position=pos3, /normal, xmargin=[1.1,0.2], $
      ymargin=[1.0,1.1], xpage=xpage, ypage=ypage, xspace=0.0, yspace=0.0

    if (ngalaxy gt 300L) then psym = 3L else begin
       plotsym, 0, 0.7, /fill
       psym = 8L
    endelse

    line = ['oii_3727','h_gamma','h_beta','oiii_5007','h_alpha','nii_6584','sii_6716','sii_6731']+'_ew'
    linelabel = 'EW('+['[O II]','H\gamma','H\beta','[O III] \lambda5007',$
      'H\alpha','[N II] \lambda6584','[S II] \lambda6716','[S II] \lambda6731']+')'

; ---------------------------------------------------------------------------
; compare emission lines EWs: value vs value + residuals
    xrange = [0.1,500]
    xrange2 = xrange
    yrange = xrange
    yrange2 = 125*[-1,1]

    posindx1 = [0,1,2,3,8,9,10,11]
    posindx2 = [4,5,6,7,12,13,14,15]
    
    for itag = 0L, n_elements(line)-1L do begin
       data11 = struct_trimtags(data1,select=line[itag])
       data21 = struct_trimtags(data2,select=line[itag])

       x = reform((data11.(0))[0,*])
       y = reform((data21.(0))[0,*])
       xerr = reform((data11.(0))[1,*])
       yerr = reform((data21.(0))[1,*])

       good1 = where((x/xerr gt 1.0) and (y/yerr gt 1.0) and $
         (xerr gt 0.0) and (abs(x) lt 1E4) and (abs(y) lt 1E4),ngood1)
       good2 = where((x/xerr gt snrcut) and (y/yerr gt snrcut) and $
         (xerr gt 0.0) and (abs(x) lt 1E4) and (abs(y) lt 1E4),ngood2)
;      splog, linelabel[itag], ngood1, ngood2

; some statistics       
       if (ngood1 ne 0) then begin
          xresid1 = y[good1]
          yresid1 = 100.0*(y[good1]-x[good1])/x[good1]
       endif
       if (ngood2 ne 0) then begin
          xresid2 = y[good2]
          yresid2 = 100.0*(y[good2]-x[good2])/x[good2]
       endif

       if (ngood2 gt 3) then begin
          stats = im_stats(yresid2,verbose=0,/baremin,sigrej=snrcut)
          xstr = strtrim(string(stats.mean_rej,format='(F12.1)'),2)+$
            '\pm'+strtrim(string(stats.sigma_rej,format='(F12.1)'),2)+$
            ' ('+strtrim(string(stats.median_rej,format='(F12.1)'),2)+')'
       endif

; make the plot       
       if ((itag mod 4) eq 0) then begin ; with [SII]
          delvarx, ytickname
          ytitle = 'EW [Zsolar]'
          ytitle2 = 'Residuals [%]'
       endif else begin
          ytickname = replicate(' ',10)
          ytitle = ''
          ytitle2 = ''
       endelse

       xtickname = replicate(' ',10)
       xtitle = ''

       if (itag eq 1) then thistitle = title else delvarx, thistitle
       
       djs_plot, x[good1], y[good1], noerase=(itag ne 0), xtitle=xtitle, ytitle=ytitle, $
         psym=symcat(16), xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, /xlog, /ylog, $
         position=pos1[*,posindx1[itag]], xtickname=xtickname, ytickname=ytickname, $
         charsize=1.4
       if (ngood2 ne 0) then djs_oplot, x[good2], y[good2], psym=symcat(16), color='red'
       djs_oplot, 10^!x.crange, 10^!y.crange, line=0
       legend, textoidl(linelabel[itag]), /left, /top, box=0, charsize=1.0, margin=0

       if (itag eq 0) then xyouts, (pos1[2,15]-pos1[0,0])/2.0+pos1[0,0], pos1[3,0]*1.01, $ ; with [SII]
         title, align=0.5, charsize=1.3, /norm

       if (posindx2[itag] ge 12) then begin ; with [SII]
          delvarx, xtickname
          xtitle = 'EW [Zgrid]'
       endif
       
       djs_plot, xresid1, yresid1, xtitle=xtitle, ytitle=ytitle2, psym=symcat(16), /noerase, $
         xrange=xrange2, yrange=yrange2, xstyle=3, ystyle=3, position=pos1[*,posindx2[itag]], $
         /xlog, xtickname=xtickname, ytickname=ytickname, charsize=1.4
       if (ngood2 ne 0) then djs_oplot, xresid2, yresid2, psym=symcat(16), color='red'
       djs_oplot, 10^!x.crange, [0,0], line=0
       if (ngood2 gt 3) then legend, textoidl(xstr), /right, /bottom, box=0, $
         charsize=1.0, margin=0

    endfor 
       
return
end
    
pro plotsings_compare_solar
; jm10jul01ucsd - compare the emission-line strengths based on the
; multi-metallicity and solar-only models

    nuclear = read_sings_gandalf(/nuclear)
    drift20 = read_sings_gandalf(/drift20)
    drift56 = read_sings_gandalf(/drift56)

    nuclear_solar = read_sings_gandalf(/nuclear,/solar)
    drift20_solar = read_sings_gandalf(/drift20,/solar)
    drift56_solar = read_sings_gandalf(/drift56,/solar)

    psfile = sings_path(/projects)+'log12oh/compare_solar.ps'
    im_plotconfig, 0, psfile=psfile
    doplot_compare, nuclear, nuclear_solar, title='Nuclear'
    doplot_compare, drift20, drift20_solar, title='Circumnuclear'
    doplot_compare, drift56, drift56_solar, title='Radial-Strip'
    im_plotconfig, psfile=psfile, /psclose, /gzip
    
return
end
