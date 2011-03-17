pro doplot_compare, ancillary, ispec, unfluxed, title=title, abscor=abscor, $
  ylabel=ylabel, histograms=histograms, postscript=postscript, silent=silent
; driving plotting routine; main routine below

    ngalaxy = n_elements(unfluxed)
    snrcut = 5.0

    pagemaker, nx=3, ny=4, position=pos1, /normal, xmargin=[1.2,0.2], $
      ymargin=[0.4,1.2], yspace=0.0, xspace=0.0, xpage=xpage, ypage=ypage
;   pagemaker, nx=4, ny=4, position=pos1, /normal, xmargin=[1.2,0.2], $ ; with [SII]
;     ymargin=[0.4,1.2], yspace=0.0, xspace=0.0, xpage=xpage, ypage=ypage
    pagemaker, nx=3, ny=2, position=pos2, /normal, xmargin=[1.1,0.2], $
      ymargin=[0.5,1.1], xpage=xpage, ypage=ypage, xspace=0.0, yspace=0.0
;   pagemaker, nx=3, ny=3, position=pos2, /normal, xmargin=[1.1,0.2], $ ; with [SII]
;     ymargin=[0.5,1.1], xpage=xpage, ypage=ypage, xspace=0.0, yspace=0.0
    pagemaker, nx=1, ny=3, position=pos3, /normal, xmargin=[1.1,0.2], $
      ymargin=[1.0,1.1], xpage=xpage, ypage=ypage, xspace=0.0, yspace=0.0

    if keyword_set(postscript) then postthick = 4.0 else postthick = 1.5

    if (ngalaxy gt 300L) then psym = 3L else begin
       plotsym, 0, 0.7, /fill
       psym = 8L
    endelse

    line = ['oii_3727','h_gamma','h_beta','oiii_5007','h_alpha','nii_6584']+'_ew'
    linelabel = 'EW('+['[O II]','H\gamma','H\beta','[O III] \lambda5007','H\alpha','[N II] \lambda6584']+')'
;   line = ['oii_3727','h_gamma','h_beta','oiii_5007','h_alpha','nii_6584','sii_6716','sii_6731']+'_ew'
;   linelabel = 'EW('+['[O II]','H\gamma','H\beta','[O III] \lambda5007',$
;     'H\alpha','[N II] \lambda6584','[S II] \lambda6716','[S II] \lambda6731']+')'

    ewhacor = 0.66 ; 0.62 ; 1.382
    ewhbcor = 1.80 ; 1.78 ; 2.170
    ewhgcor = 1.71 ; 1.70 ; 1.725
;   ewhacor = 0.0 & ewhbcor = 0.0 & ewhgcor = 0.0 
    
;;; ---------------------------------------------------------------------------
;;; plot the distribution of Balmer absorption-line EWs
;;; ---------------------------------------------------------------------------
;;
;;    snrcut2 = 3.0 ; 4.0
;;    binsize = 0.1
;;    xrange = [-1,6]
;;    
;;; H-alpha    
;;
;;    good = where((ispec.babs_h_alpha_ew[1] gt 0.0),ngood)
;;    high = where((ispec[good].babs_h_alpha_ew[0]/ispec[good].babs_h_alpha_ew[1]) gt snrcut2,nhigh)
;;    if (not keyword_set(silent)) then splog, 'H-alpha', ngood, nhigh
;;    im_plothist, ispec[good].babs_h_alpha_ew[0], bin=binsize, xha, yha, /noplot
;;
;;    yrange = [0,max(yha)*1.2]
;;
;;    djs_plot, [0], [0], /nodata, xtitle='', ytitle='Number', xrange=xrange, $
;;      yrange=yrange, xstyle=1, ystyle=1, xthick=postthick, ythick=postthick, $
;;      charsize=1.4, charthick=postthick, xtickname=replicate(' ',10), $
;;      position=pos3[*,0], title='Balmer Absorption-Line EWs !c'+title
;;    im_plothist, ispec[good].babs_h_alpha_ew[0], bin=binsize, xha, yha, thick=postthick, /overplot
;;    if (nhigh gt 3L) then $
;;      im_plothist, ispec[good[high]].babs_h_alpha_ew[0], bin=binsize, xha_hi, yha_hi, $
;;      /fill, fcolor=fsc_color('light grey',3), thick=postthick, /overplot
;;
;;    stats = im_stats(ispec[good].babs_h_alpha_ew[0],sigrej=3.0)
;;;   stats = im_stats(ispec[good[high]].babs_h_alpha_ew[0],sigrej=3.0)
;;    legend, textoidl(['H\alpha '+strtrim(string(stats.mean_rej,format='(F12.3)'),2)+'\pm'+$
;;      strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+' \AA ('+$
;;      strtrim(string(stats.median_rej,format='(F12.3)'),2)+' \AA)']), /left, /top, $
;;      box=0, charsize=1.4, charthick=postthick;, /clear
;;
;;; H-beta    
;;
;;    good = where((ispec.babs_h_beta_ew[1] gt 0.0),ngood)
;;    high = where((ispec[good].babs_h_beta_ew[0]/ispec[good].babs_h_beta_ew[1]) gt snrcut2,nhigh)
;;    if (not keyword_set(silent)) then splog, 'H-beta', ngood, nhigh
;;    im_plothist, ispec[good].babs_h_beta_ew[0], bin=binsize, xhb, yhb, /noplot
;;
;;    djs_plot, [0], [0], /nodata, /noerase, xtitle='', ytitle='Number', xrange=xrange, $
;;      yrange=yrange, xstyle=1, ystyle=1, xthick=postthick, ythick=postthick, $
;;      charsize=1.4, charthick=postthick, xtickname=replicate(' ',10), $
;;      position=pos3[*,1]
;;    im_plothist, ispec[good].babs_h_beta_ew[0], bin=binsize, xhb, yhb, thick=postthick, /overplot
;;    if (nhigh gt 3L) then $
;;      im_plothist, ispec[good[high]].babs_h_beta_ew[0], bin=binsize, xhb_hi, yhb_hi, $
;;      /fill, fcolor=fsc_color('light grey',3), thick=postthick, /overplot
;;
;;    stats = im_stats(ispec[good].babs_h_beta_ew[0],sigrej=3.0)
;;;   stats = im_stats(ispec[good[high]].babs_h_beta_ew[0],sigrej=3.0)
;;    legend, textoidl(['H\beta '+strtrim(string(stats.mean_rej,format='(F12.3)'),2)+'\pm'+$
;;      strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+ '\AA ('+$
;;      strtrim(string(stats.median_rej,format='(F12.3)'),2)+' \AA)']), /left, /top, $
;;      box=0, charsize=1.4, charthick=postthick;, /clear
;;
;;; H-gamma
;;
;;    good = where((ispec.babs_h_gamma_ew[1] gt 0.0),ngood)
;;    high = where((ispec[good].babs_h_gamma_ew[0]/ispec[good].babs_h_gamma_ew[1]) gt snrcut2,nhigh)
;;    if (not keyword_set(silent)) then splog, 'H-gamma', ngood, nhigh
;;    im_plothist, ispec[good].babs_h_gamma_ew[0], bin=binsize, xhg, yhg, /noplot
;;
;;    djs_plot, [0], [0], /nodata, /noerase, xtitle='EW (\AA)', ytitle='Number', $
;;      xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, xthick=postthick, $
;;      ythick=postthick, charsize=1.4, charthick=postthick, position=pos3[*,2]
;;    im_plothist, ispec[good].babs_h_gamma_ew[0], bin=binsize, xhg, yhg, thick=postthick, /overplot
;;    if (nhigh gt 3L) then $
;;      im_plothist, ispec[good[high]].babs_h_gamma_ew[0], bin=binsize, xhg_hi, yhg_hi, $
;;      /fill, fcolor=fsc_color('light grey',3), thick=postthick, /overplot
;;
;;    stats = im_stats(ispec[good].babs_h_gamma_ew[0],sigrej=3.0)
;;;   stats = im_stats(ispec[good[high]].babs_h_gamma_ew[0],sigrej=3.0)
;;    legend, textoidl(['H\gamma '+strtrim(string(stats.mean_rej,format='(F12.3)'),2)+'\pm'+$
;;      strtrim(string(stats.sigma_rej,format='(F12.3)'),2)+ '\AA ('+$
;;      strtrim(string(stats.median_rej,format='(F12.3)'),2)+' \AA)']), /left, /top, $
;;      box=0, charsize=1.4, charthick=postthick;, /clear
;;
;;    if (not keyword_set(postscript)) then begin
;;       splog, 'Press any key to continue.'
;;       cc = get_kbrd(1)
;;    endif
       
; ---------------------------------------------------------------------------    
; compare emission lines using histogram distributions
; ---------------------------------------------------------------------------    

    if keyword_set(histograms) then begin
       
       binsize = 0.4
       xrange = 5.0*[-1,1]
       xgauss = findgen((xrange[1]-xrange[0])/0.01+1)*0.01+xrange[0]
       ygauss = gauss1(xgauss,[0.0,1.0,1.0])
       ygauss = ygauss/max(ygauss)

       for itag = 0L, n_elements(line)-1L do begin

          ispec1 = struct_trimtags(ispec,select=line[itag])
          unfluxed1 = struct_trimtags(unfluxed,select=line[itag])

          x = reform((ispec1.(0))[0,*])
          y = reform((unfluxed1.(0))[0,*])
          xerr = reform((ispec1.(0))[1,*])
          yerr = reform((unfluxed1.(0))[1,*])
          
          good1 = where((x/xerr gt 1.0) and (y/yerr gt 1.0) and $
            (xerr gt 0.0) and (abs(x) lt 1E4) and (abs(y) lt 1E4),ngood1)
          good2 = where((x/xerr gt snrcut) and (y/yerr gt snrcut) and $
            (xerr gt 0.0) and (abs(x) lt 1E4) and (abs(y) lt 1E4),ngood2)
          if (not keyword_set(silent)) then splog, linelabel[itag], ngood1, ngood2

; correct for stellar absorption          

          if keyword_set(abscor) then begin
             
             hgmatch = strmatch(linelabel[itag],'*gamma*',/fold)
             hbmatch = strmatch(linelabel[itag],'*beta*',/fold)
             hamatch = strmatch(linelabel[itag],'*alpha*',/fold)
             match = hamatch or hbmatch or hgmatch

             if hamatch then y = y + ewhacor
             if hbmatch then y = y + ewhbcor
             if hgmatch then y = y + ewhgcor

          endif

; some statistics       
          
          if (ngood2 gt 3L) then begin
             ratio2 = (y[good2]-x[good2])/xerr[good2]
             im_plothist, ratio2, bin=binsize, xhist2, yhist2, /noplot, /nan ;, /peak
             stats2 = im_stats(ratio2,sigrej=3.0) ;,/verbose)
          endif

          if (ngood1 gt 3L) then begin
             ratio1 = (y[good1]-x[good1])/xerr[good1]
             im_plothist, ratio1, bin=binsize, xhist1, yhist1, /noplot, /nan ;, /peak
             stats1 = im_stats(ratio1,sigrej=3.0) ;,/verbose)
          endif

; make the plot          
          
          if (itag lt 3L) then begin
;         if (itag lt 5L) then begin ; with [SII]
             xtickname = replicate(' ',10)
             xtitle = ''
          endif else begin
             delvarx, xtickname
             xtitle = '\Delta(EW)/\sigma(EW)'
          endelse
          
          if ((itag mod 3L) eq 0L) then begin
             delvarx, ytickname
             ytitle = 'Relative Fraction'
          endif else begin
             ytickname = replicate(' ',10)
             ytitle = ''
          endelse
          
          djs_plot, [0], [0], /nodata, noerase=(itag gt 0L), xtitle=xtitle, $
            ytitle=ytitle, xrange=xrange, yrange=[0.0,1.3], xstyle=1, ystyle=1, $
            xthick=postthick, ythick=postthick, charsize=1.4, charthick=postthick, $
            position=pos2[*,itag], xtickname=xtickname, ytickname=ytickname
          if (itag eq 0L) then xyouts, (pos2[2,2]-pos2[0,0])/2.0+pos2[0,0], pos2[3,0]*1.01, $
            title, align=0.5, charthick=postthick, charsize=1.4, /norm
          if (ngood2 gt 3L) then im_plothist, ratio2, bin=binsize, /nan, /fill, $ ; shaded S/N=5 histogram
            fcolor=fsc_color('light grey',3), normfactor=max(yhist1), /overplot
          if (ngood1 gt 3L) then im_plothist, ratio1, bin=binsize, /nan, $ ; S/N=1 histogram
            normfactor=max(yhist1), /overplot              
; plot a unit Gaussian
          djs_oplot, xgauss, ygauss, line=1, thick=postthick, color='red'
;         djs_oplot, [0,0], !y.crange, line=0, thick=1.0
; legend
          if (ngood2 gt 3L) then begin
             djs_oplot, stats2.mean_rej*[1,1], !y.crange, line=2, thick=postthick
             label = [linelabel[itag],strtrim(string(stats2.mean_rej,format='(F12.2)'),2)+'\pm'+$
               strtrim(string(stats2.sigma_rej,format='(F12.2)'),2)+' ('+$
               strtrim(string(stats2.median_rej,format='(F12.2)'),2)+')']
             legend, textoidl(label), /left, /top, box=0, charsize=1.0, charthick=postthick, $
               clear=keyword_set(postscript)
          endif

       endfor
       
       if (not keyword_set(postscript)) then begin
          splog, 'Press any key to continue.'
          cc = get_kbrd(1)
       endif

    endif
       
; ---------------------------------------------------------------------------
; compare emission lines EWs: value vs value + residuals
; ---------------------------------------------------------------------------
    
    xrange = [0.1,500]
    xrange2 = xrange
    yrange = xrange
    yrange2 = 125*[-1,1]

    posindx1 = [0,1,2,6,7,8]
    posindx2 = [3,4,5,9,10,11]
;   posindx1 = [0,1,2,3,8,9,10,11] ; with [SII]
;   posindx2 = [4,5,6,7,12,13,14,15] ; with [SII]
    
    for itag = 0L, n_elements(line)-1L do begin

       ispec1 = struct_trimtags(ispec,select=line[itag])
       unfluxed1 = struct_trimtags(unfluxed,select=line[itag])

       x = reform((ispec1.(0))[0,*])
       y = reform((unfluxed1.(0))[0,*])
       xerr = reform((ispec1.(0))[1,*])
       yerr = reform((unfluxed1.(0))[1,*])

       good1 = where((x/xerr gt 1.0) and (y/yerr gt 1.0) and $
         (xerr gt 0.0) and (abs(x) lt 1E4) and (abs(y) lt 1E4),ngood1)
       good2 = where((x/xerr gt snrcut) and (y/yerr gt snrcut) and $
         (xerr gt 0.0) and (abs(x) lt 1E4) and (abs(y) lt 1E4),ngood2)
;      splog, linelabel[itag], ngood1, ngood2

; correct for stellar absorption          
       
       if keyword_set(abscor) then begin

          hgmatch = strmatch(linelabel[itag],'*gamma*',/fold)
          hbmatch = strmatch(linelabel[itag],'*beta*',/fold)
          hamatch = strmatch(linelabel[itag],'*alpha*',/fold)
          match = hamatch or hbmatch or hgmatch

          if hamatch then y = y + ewhacor
          if hbmatch then y = y + ewhbcor
          if hgmatch then y = y + ewhgcor

       endif

; some statistics       
       
       if (ngood1 ne 0L) then begin
          xresid1 = y[good1]
          yresid1 = 100.0*(y[good1]-x[good1])/x[good1]
       endif

       if (ngood2 ne 0L) then begin
          xresid2 = y[good2]
          yresid2 = 100.0*(y[good2]-x[good2])/x[good2]
       endif

       if (ngood2 gt 3L) then begin
          stats = im_stats(yresid2,verbose=0,/baremin,sigrej=snrcut)
          xstr = strtrim(string(stats.mean_rej,format='(F12.1)'),2)+$
            '\pm'+strtrim(string(stats.sigma_rej,format='(F12.1)'),2)+$
            ' ('+strtrim(string(stats.median_rej,format='(F12.1)'),2)+')'
       endif

;      splog, 'Objects with large residuals after S/N cuts:'
;      srt = reverse(sort(abs(yresid2))) & srt = srt[0:25<(ngood2-1L)]
;      niceprint, unfluxed[good2[srt]].pass, unfluxed[good2[srt]].aper, ispec[good2[srt]].fit_id*2, yresid2[srt], $
;        x[good2[srt]], y[good2[srt]], xerr[good2[srt]], yerr[good2[srt]], unfluxed[good2[srt]].z, srt

; make the plot       
       
       if ((itag mod 3L) eq 0L) then begin
;      if ((itag mod 4L) eq 0L) then begin ; with [SII]
          delvarx, ytickname
          ytitle = 'EW ['+ylabel+']'
          ytitle2 = 'Residuals [%]'
       endif else begin
          ytickname = replicate(' ',10)
          ytitle = ''
          ytitle2 = ''
       endelse

       xtickname = replicate(' ',10)
       xtitle = ''

;      xrange = [(min(x[good1])<min(y[good1]))>0.1,((max(x[good1])>max(y[good1]))*1.2)<500.0]

       if (itag eq 1L) then thistitle = title else delvarx, thistitle
       
       djs_plot, x[good1], y[good1], noerase=(itag ne 0), xtitle=xtitle, ytitle=ytitle, $
         psym=psym, xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, /xlog, /ylog, $
         xthick=postthick, ythick=postthick, charsize=1.4, charthick=postthick, position=pos1[*,posindx1[itag]], $
         xtickname=xtickname, ytickname=ytickname, title=thistitle
       if (ngood2 ne 0L) then djs_oplot, x[good2], y[good2], psym=psym, color='red'
       djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick
       legend, textoidl(linelabel[itag]), /left, /top, box=0, charsize=1.0, charthick=postthick, $
         margin=0, clear=keyword_set(postscript)

;      if (itag eq 0L) then xyouts, (pos1[2,15]-pos1[0,0])/2.0+pos1[0,0], pos1[3,0]*1.01, $ ; with [SII]
;        title, align=0.5, charthick=postthick, charsize=1.3, /norm

       if (posindx2[itag] ge 8L) then begin
;      if (posindx2[itag] ge 12L) then begin ; with [SII]
          delvarx, xtickname
          xtitle = 'EW [ispec]'
       endif
       
       djs_plot, xresid1, yresid1, xtitle=xtitle, ytitle=ytitle2, psym=psym, /noerase, $
         xrange=xrange2, yrange=yrange2, xstyle=3, ystyle=3, position=pos1[*,posindx2[itag]], $
         xthick=postthick, ythick=postthick, charsize=1.4, charthick=postthick, /xlog, $
         xtickname=xtickname, ytickname=ytickname
       if (ngood2 ne 0L) then djs_oplot, xresid2, yresid2, psym=psym, color='red'
       djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick
       if (ngood2 gt 3L) then legend, textoidl(xstr), /right, /bottom, box=0, $
         charsize=1.0, charthick=postthick, margin=0, clear=keyword_set(postscript)

    endfor 
    
    if (not keyword_set(postscript)) then begin
       splog, 'Press any key to continue.'
       cc = get_kbrd(1)
    endif
       
;;; ---------------------------------------------------------------------------
;;; compare emission-line S/N
;;; ---------------------------------------------------------------------------
;;    
;;    xrange = [0.5,500]
;;    xrange2 = xrange
;;    yrange = xrange
;;    yrange2 = 125*[-1,1]
;;
;;    for itag = 0L, n_elements(line)-1L do begin
;;
;;       ispec1 = struct_trimtags(ispec,select=line[itag])
;;       unfluxed1 = struct_trimtags(unfluxed,select=line[itag])
;;
;;       x = reform((ispec1.(0))[0,*])
;;       y = reform((unfluxed1.(0))[0,*])
;;       xerr = reform((ispec1.(0))[1,*])
;;       yerr = reform((unfluxed1.(0))[1,*])
;;
;;       good1 = where((x/xerr gt 1.0) and (y/yerr gt 1.0) and $
;;         (xerr gt 0.0) and (abs(x) lt 1E4) and (abs(y) lt 1E4),ngood1)
;;       good2 = where((x/xerr gt snrcut) and (y/yerr gt snrcut) and $
;;         (xerr gt 0.0) and (abs(x) lt 1E4) and (abs(y) lt 1E4),ngood2)
;;       if (not keyword_set(silent)) then splog, linelabel[itag], ngood1, ngood2
;;
;;; correct for stellar absorption          
;;       
;;       if keyword_set(abscor) then begin
;;
;;          hgmatch = strmatch(linelabel[itag],'*gamma*',/fold)
;;          hbmatch = strmatch(linelabel[itag],'*beta*',/fold)
;;          hamatch = strmatch(linelabel[itag],'*alpha*',/fold)
;;          match = hamatch or hbmatch or hgmatch
;;
;;          if hamatch then y = y + ewhacor
;;          if hbmatch then y = y + ewhbcor
;;          if hgmatch then y = y + ewhgcor
;;
;;       endif
;;          
;;       if (ngood2 gt 3L) then begin
;;          ratio2 = ((y[good2]/yerr[good2])-(x[good2]/xerr[good2]))/(x[good2]/xerr[good2])
;;          im_plothist, ratio2, bin=binsize, xhist2, yhist2, /noplot, /nan ;, /peak
;;          stats2 = im_stats(100*ratio2,sigrej=3.0) ;,/verbose)
;;          xstr = strtrim(string(stats2.mean_rej,format='(F12.1)'),2)+$
;;            '\pm'+strtrim(string(stats2.sigma_rej,format='(F12.1)'),2)+$
;;            '% ('+strtrim(string(stats2.median_rej,format='(F12.1)'),2)+'%)'
;;       endif
;;
;;       if (ngood1 gt 3L) then begin
;;          ratio1 = ((y[good1]/yerr[good1])-(x[good1]/xerr[good1]))/(x[good1]/xerr[good1])
;;          im_plothist, ratio1, bin=binsize, xhist1, yhist1, /noplot, /nan ;, /peak
;;          stats1 = im_stats(ratio1,sigrej=3.0) ;,/verbose)
;;       endif
;;
;;; make the plot       
;;       
;;       if (itag lt 5L) then begin
;;          xtickname = replicate(' ',10)
;;          xtitle = ''
;;       endif else begin
;;          delvarx, xtickname
;;          xtitle = '\delta(EW)/EW [ispec]'
;;       endelse
;;       
;;       if ((itag mod 3L) eq 0L) then begin
;;          delvarx, ytickname
;;          ytitle = '\delta(EW)/EW [b-spline]'
;;       endif else begin
;;          ytickname = replicate(' ',10)
;;          ytitle = ''
;;       endelse
;;
;;;      xrange = [(min(x[good1])<min(y[good1]))>0.1,((max(x[good1])>max(y[good1]))*1.2)<500.0]
;;
;;       djs_plot, x[good1]/xerr[good1], y[good1]/yerr[good1], noerase=(itag ne 0), xtitle=xtitle, ytitle=ytitle, $
;;         psym=psym, xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, /xlog, /ylog, $
;;         xthick=postthick, ythick=postthick, charsize=1.4, charthick=postthick, position=pos2[*,itag], $
;;         xtickname=xtickname, ytickname=ytickname
;;       if (ngood2 ne 0L) then djs_oplot, x[good2]/xerr[good2], y[good2]/yerr[good2], psym=psym, color='red'
;;       djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick
;;       legend, textoidl(linelabel[itag]), /left, /top, box=0, charsize=1.0, charthick=postthick, $
;;         margin=0, clear=keyword_set(postscript)
;;       if (ngood2 gt 3L) then legend, textoidl(xstr), /right, /bottom, box=0, $
;;         charsize=1.0, charthick=postthick, margin=0, clear=keyword_set(postscript)
;;       
;;    endfor 
;;    
;;    if (not keyword_set(postscript)) then begin
;;       splog, 'Press any key to continue.'
;;       cc = get_kbrd(1)
;;    endif

;;; ---------------------------------------------------------------------------
;;; compare emission lines EW *errors*; value vs value + residuals
;;; ---------------------------------------------------------------------------
;;    
;;    xrange = [0.1,200]
;;    xrange2 = xrange
;;    yrange = xrange
;;    yrange2 = 125*[-1,1]
;;
;;    posindx1 = [0,1,2,3,8,9,10,11]
;;    posindx2 = [4,5,6,7,12,13,14,15]
;;    
;;    for itag = 0L, n_elements(line)-1L do begin
;;
;;       ispec1 = struct_trimtags(ispec,select=line[itag])
;;       unfluxed1 = struct_trimtags(unfluxed,select=line[itag])
;;
;;       x = reform((ispec1.(0))[0,*])
;;       y = reform((unfluxed1.(0))[0,*])
;;       xerr = reform((ispec1.(0))[1,*])
;;       yerr = reform((unfluxed1.(0))[1,*])
;;
;;       good1 = where((x/xerr gt 1.0) and (y/yerr gt 1.0) and $
;;         (xerr gt 0.0) and (abs(x) lt 1E4) and (abs(y) lt 1E4),ngood1)
;;       good2 = where((x/xerr gt snrcut) and (y/yerr gt snrcut) and $
;;         (xerr gt 0.0) and (abs(x) lt 1E4) and (abs(y) lt 1E4),ngood2)
;;;      splog, linelabel[itag], ngood1, ngood2
;;
;;       if (ngood1 ne 0L) then begin
;;          xresid1 = yerr[good1] ; NOTE I'M USING ERRORS HERE!
;;          yresid1 = 100.0*(yerr[good1]-xerr[good1])/xerr[good1]
;;       endif
;;
;;       if (ngood2 ne 0L) then begin
;;          xresid2 = yerr[good2] ; NOTE I'M USING ERRORS HERE!
;;          yresid2 = 100.0*(yerr[good2]-xerr[good2])/xerr[good2]
;;       endif
;;
;;       if (ngood2 gt 3L) then begin
;;          stats = im_stats(yresid2,verbose=0,/baremin,sigrej=snrcut)
;;          xstr = strtrim(string(stats.mean_rej,format='(F12.1)'),2)+$
;;            '\pm'+strtrim(string(stats.sigma_rej,format='(F12.1)'),2)+$
;;            ' ('+strtrim(string(stats.median_rej,format='(F12.1)'),2)+')'
;;       endif
;;
;;; make the plot       
;;       
;;       if ((itag mod 4L) eq 0L) then begin
;;          delvarx, ytickname
;;          ytitle = '\delta (EW) [b-spline]'
;;          ytitle2 = 'Residuals (%)'
;;       endif else begin
;;          ytickname = replicate(' ',10)
;;          ytitle = ''
;;          ytitle2 = ''
;;       endelse
;;
;;       xtickname = replicate(' ',10)
;;       xtitle = ''
;;
;;;      xrange = [(min(x[good1])<min(y[good1]))>0.1,((max(x[good1])>max(y[good1]))*1.2)<500.0]
;;
;;       djs_plot, xerr[good1], yerr[good1], noerase=(itag ne 0), xtitle=xtitle, ytitle=ytitle, $ ; ERRORS!
;;         psym=psym, xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, /xlog, /ylog, $
;;         xthick=postthick, ythick=postthick, charsize=1.4, charthick=postthick, position=pos1[*,posindx1[itag]], $
;;         xtickname=xtickname, ytickname=ytickname
;;       if (ngood2 ne 0L) then djs_oplot, xerr[good2], yerr[good2], psym=psym, color='red'
;;       djs_oplot, 10^!x.crange, 10^!y.crange, line=0, thick=postthick
;;       legend, textoidl(linelabel[itag]), /left, /top, box=0, charsize=1.0, charthick=postthick, $
;;         margin=0, clear=keyword_set(postscript)
;;       
;;       if (itag eq 0L) then xyouts, (pos1[2,15]-pos1[0,0])/2.0+pos1[0,0], pos1[3,0]*1.01, $
;;         title, align=0.5, charthick=postthick, charsize=1.3, /norm
;;
;;       if (posindx2[itag] ge 12L) then begin
;;          delvarx, xtickname
;;          xtitle = '\delta (EW) [ispec]'
;;       endif
;;       
;;       djs_plot, xresid1, yresid1, xtitle=xtitle, ytitle=ytitle2, psym=psym, /noerase, $
;;         xrange=xrange2, yrange=yrange2, xstyle=3, ystyle=3, position=pos1[*,posindx2[itag]], $
;;         xthick=postthick, ythick=postthick, charsize=1.4, charthick=postthick, /xlog, $
;;         xtickname=xtickname, ytickname=ytickname
;;       if (ngood2 ne 0L) then djs_oplot, xresid2, yresid2, psym=psym, color='red'
;;       djs_oplot, 10^!x.crange, [0,0], line=0, thick=postthick
;;       if (ngood2 gt 3L) then legend, textoidl(xstr), /left, /bottom, box=0, $
;;         charsize=1.0, charthick=postthick, margin=0, clear=keyword_set(postscript)
;;
;;    endfor 
;;    
;;    if (not keyword_set(postscript)) then begin
;;       splog, 'Press any key to continue.'
;;       cc = get_kbrd(1)
;;    endif
       
return
end
    
pro ages_compare_specdata, postscript=postscript
; jm07apr13nyu - compare the EWs measured from the fluxed & unfluxed plates
; jm08nov01nyu - rewritten a bit
    
    specfitpath = ages_path(/specfit)
    analysis_path = ages_path(/analysis)

    ancillary1 = read_ages(/ancillary)
    ispec1 = read_ages(/ispec)
    unfluxed1 = read_ages(/unfluxed)
;   parent = where((ancillary1.main_flag eq 1L) and (ancillary1.z gt 0.01) and (ancillary1.z lt 0.8),nparent)
    parent = where((ancillary1.main_flag eq 1) and (ancillary1.z gt 0.01) and (ancillary1.z lt 0.8) and $
      (ancillary1.pass ne 106) and (ancillary1.pass ne 110) and (ancillary1.pass ne 209) and $
      (ancillary1.pass ne 310) and (ancillary1.pass ne 311),nparent)
    ancillary = (temporary(ancillary1))[parent]
    unfluxed = (temporary(unfluxed1))[parent]
    ispec = (temporary(ispec1))[parent]

; compare all the plates together, and then each plate individually 

    allpass = ancillary.pass
    upass = allpass[uniq(allpass,sort(allpass))]
    npass = n_elements(upass)

    if (not keyword_set(postscript)) then im_window, 0, /square
    
; ###########################################################################
; ispec vs unfluxed - do not compare unfluxed plates
; ###########################################################################
    
    if keyword_set(postscript) then begin
       psname = analysis_path+'ages_compare_specdata_ispec_vs_unfluxed.ps'
       splog, 'Writing '+psname
       dfpsplot, psname, /color
    endif

    doplot_compare, ancillary, ispec, unfluxed, postscript=postscript, $
      /abscor, /histograms, ylabel='b-spline', title='All Plates: ispec vs b-spline'

    for ipass = 0L, npass-1L do begin
       print, format='("Plate ",I0,"/",I0,".",A10,$)', ipass+1, npass, string(13b)
       title = 'Plate '+string(upass[ipass],format='(I0)')+': ispec vs b-spline'
       passindx = where(upass[ipass] eq allpass)
       doplot_compare, ancillary[passindx], ispec[passindx], unfluxed[passindx], $
         title=title, ylabel='b-spline', /abscor, /histograms, postscript=postscript, /silent
    endfor
    
    if keyword_set(postscript) then begin
       dfpsclose
       spawn, 'gzip -f '+psname, /sh
    endif

stop    
    
; ###########################################################################
; ispec vs ispec_restrict_old
; ###########################################################################
    
    if keyword_set(postscript) then begin
       psname = analysis_path+'ages_compare_specdata_ispec_vs_ispec_restrict_old.ps'
       splog, 'Writing '+psname
       dfpsplot, psname, /color
    endif

    doplot_compare, ancillary, ispec, ispec_restrict_old, postscript=postscript, $
      ylabel='ispec-restrict-old', title='All Plates: ispec vs ispec-restrict-old'

    for ipass = 0L, npass-1L do begin
       print, format='("Plate ",I0,"/",I0,".",A10,$)', ipass+1, npass, string(13b)
       title = 'Plate '+string(upass[ipass],format='(I0)')+': ispec vs ispec-restrict-old'
       passindx = where(upass[ipass] eq allpass)
       doplot_compare, ancillary[passindx], ispec[passindx], ispec_restrict_old[passindx], $
         title=title, ylabel='ispec-restrict-old', postscript=postscript, /silent
    endfor
    
    if keyword_set(postscript) then begin
       dfpsclose
       spawn, 'gzip -f '+psname, /sh
    endif

stop    
    
; ###########################################################################
; ispec vs ispec_Zmulti    
; ###########################################################################
    
    if keyword_set(postscript) then begin
       psname = analysis_path+'ages_compare_specdata_ispec_vs_ispec_Zmulti.ps'
       splog, 'Writing '+psname
       dfpsplot, psname, /color
    endif

    doplot_compare, ancillary, ispec, ispec_Zmulti, postscript=postscript, $
      ylabel='ispec-Zmulti', title='All Plates: ispec vs ispec-Zmulti'

    for ipass = 0L, npass-1L do begin
       print, format='("Plate ",I0,"/",I0,".",A10,$)', ipass+1, npass, string(13b)
       title = 'Plate '+string(upass[ipass],format='(I0)')+': ispec vs ispec-Zmulti'
       passindx = where(upass[ipass] eq allpass)
       doplot_compare, ancillary[passindx], ispec[passindx], ispec_Zmulti[passindx], $
         title=title, ylabel='ispec-Zmulti', postscript=postscript, /silent
    endfor
    
    if keyword_set(postscript) then begin
       dfpsclose
       spawn, 'gzip -f '+psname, /sh
    endif

stop       

return
end
