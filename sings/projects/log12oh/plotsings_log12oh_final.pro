pro plotsings_log12oh_final, ps=ps
; jm10mar10ucsd - build the plots containing the final metallicities

    metpath = sings_path(/projects)+'log12oh/'
    pspath = sings_path(/papers)+'log12oh/FIG_LOG12OH/'
    version = sings_log12oh_version()

    finalfile = metpath+'sings_log12oh_final_'+version+'.fits.gz'
    final = mrdfits(finalfile,1)
    final_kk04 = mrdfits(finalfile,2)
    final_pt05 = mrdfits(finalfile,3)
    lzfit = mrdfits(metpath+'lzfit_'+version+'.fits.gz',1)
    ngal = n_elements(final)
    
    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

    mbpivot = lzfit[0].pivot
    mbaxis = findgen(((-13.2)-(-22.0))/0.01+1)*0.01-22.0

    earlysym1 = 16 & earlypsize1 = 1.7 & earlycolor1 = 'navy' ; 'dodger blue'
    midsym1   = 15 & midpsize1   = 1.5 & midcolor1   = 'tan' ; 'salmon' ; 'black'
    latesym1  = 14 & latepsize1  = 1.8 & latecolor1  = 'firebrick' ; tan

;   earlysym1 = 16 & earlypsize1 = 1.8 & earlycolor1 = 'dodger blue'
;   midsym1   = 15 & midpsize1   = 1.7 & midcolor1   = 'black'
;   latesym1  = 4 & latepsize1  = 2.0 & latecolor1  = 'firebrick'

    nuclearsym = 34 & nuclearpsize = 1.7 & nuclearcolor = 'blue'
    circumsym  = 17 & circumpsize  = 1.5 & circumcolor = 'magenta'
    stripsym   = 14 & strippsize   = 1.4 & stripcolor = 'orange'
    charsym    = 15 & charpsize    = 1.0 & charcolor = 'forest green'
    centsym    = 46 & centpsize    = 1.8 & centcolor = 'dodger blue'
    avgsym     = 16 & avgpsize     = 1.1 & avgcolor = 'firebrick'

; ---------------------------------------------------------------------------    
; all the abundances vs R/R25

    xrange1 = [0.0025,1.0]
    yrange1 = [7.6,9.6]
    yrange2 = [6.9,9.1]

    xtitle1 = '\rho/\rho_{25}'
    ytitle1 = '12 + log (O/H) [KK04]'
    ytitle2 = '12 + log (O/H) [PT05]'

    psfile = pspath+'sings_rr25_vs_12oh_all'+suffix
    im_plotconfig, 6, pos, psfile=psfile, $
      charsize=1.6, height=3.5*[1,1], width=8

; #########################
; KK04    
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      xtitle='', ytitle=ytitle1, xrange=xrange1, yrange=yrange1, $
      xtickname=replicate(' ',10), /xlog
    djs_oplot, 0.1*[1,1], !y.crange, line=5

    nuclear = where(strmatch(final_kk04.comment,'*nuclear*',/fold),nnuclear)
    circum = where(strmatch(final_kk04.comment,'*circum*',/fold),ncircum)
    strip = where(strmatch(final_kk04.comment,'*strip*',/fold),nstrip)
    char = where(strmatch(final_kk04.comment,'*char*',/fold),nchar)
    cent = where(strmatch(final_kk04.comment,'*cent*',/fold),ncent)
    avg = where(strmatch(final_kk04.comment,'*avg*',/fold),navg)

    rchar = (1+randomu(seed,nchar)*0.2)*0.0+1.0
    rcent = (1+randomu(seed,ncent)*0.15)*0.0+1.0
    ravg = (1+randomu(seed,navg)*0.2)*0.0+1.0

    if (nnuclear ne 0) then oploterror, final_kk04[nuclear].rr25, final_kk04[nuclear].oh, $
      final_kk04[nuclear].oh_err, color=fsc_color(nuclearcolor,100), $
      errcolor=fsc_color(nuclearcolor,100), errthick=!p.thick, psym=symcat(nuclearsym), $
      symsize=nuclearpsize
    if (ncircum ne 0) then oploterror, final_kk04[circum].rr25, final_kk04[circum].oh, $
      final_kk04[circum].oh_err, color=fsc_color(circumcolor,100), $
      errcolor=fsc_color(circumcolor,100), errthick=!p.thick, psym=symcat(circumsym), $
      symsize=circumpsize
    if (nstrip ne 0) then oploterror, final_kk04[strip].rr25, final_kk04[strip].oh, $
      final_kk04[strip].oh_err, color=fsc_color(stripcolor,100), $
      errcolor=fsc_color(stripcolor,100), errthick=!p.thick, psym=symcat(stripsym), $
      symsize=strippsize
    if (nchar ne 0) then oploterror, final_kk04[char].rr25*rchar, final_kk04[char].oh, $
      final_kk04[char].oh_err, color=fsc_color(charcolor,100), $
      errcolor=fsc_color(charcolor,100), errthick=!p.thick, psym=symcat(charsym), $
      symsize=charpsize
    if (ncent ne 0) then oploterror, (final_kk04[cent].rr25+0.003)*rcent, final_kk04[cent].oh, $
      final_kk04[cent].oh_err, color=fsc_color(centcolor,100), $
      errcolor=fsc_color(centcolor,100), errthick=!p.thick, psym=symcat(centsym), $
      symsize=centpsize
    if (navg ne 0) then oploterror, final_kk04[avg].rr25*ravg, final_kk04[avg].oh, $
      final_kk04[avg].oh_err, color=fsc_color(avgcolor,100), $
      errcolor=fsc_color(avgcolor,100), errthick=!p.thick, psym=symcat(avgsym), $
      symsize=avgpsize

; #########################
; PT05    
    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      xtitle=xtitle1, ytitle=ytitle2, xrange=xrange1, yrange=yrange2, /xlog
    djs_oplot, 0.1*[1,1], !y.crange, line=5

    nuclear = where(strmatch(final_pt05.comment,'*nuclear*',/fold),nnuclear)
    circum = where(strmatch(final_pt05.comment,'*circum*',/fold),ncircum)
    strip = where(strmatch(final_pt05.comment,'*strip*',/fold),nstrip)
    char = where(strmatch(final_pt05.comment,'*char*',/fold),nchar)
    cent = where(strmatch(final_pt05.comment,'*cent*',/fold),ncent)
    avg = where(strmatch(final_pt05.comment,'*avg*',/fold),navg)

    rchar = (1+randomu(seed,nchar)*0.2)*0.0+1.0
    rcent = (1+randomu(seed,ncent)*0.15)*0.0+1.0
    ravg = (1+randomu(seed,navg)*0.2)*0.0+1.0
    
    if (nnuclear ne 0) then oploterror, final_pt05[nuclear].rr25, final_pt05[nuclear].oh, $
      final_pt05[nuclear].oh_err, color=fsc_color(nuclearcolor,100), $
      errcolor=fsc_color(nuclearcolor,100), errthick=!p.thick, psym=symcat(nuclearsym), $
      symsize=nuclearpsize
    if (ncircum ne 0) then oploterror, final_pt05[circum].rr25, final_pt05[circum].oh, $
      final_pt05[circum].oh_err, color=fsc_color(circumcolor,100), $
      errcolor=fsc_color(circumcolor,100), errthick=!p.thick, psym=symcat(circumsym), $
      symsize=circumpsize
    if (nstrip ne 0) then oploterror, final_pt05[strip].rr25, final_pt05[strip].oh, $
      final_pt05[strip].oh_err, color=fsc_color(stripcolor,100), $
      errcolor=fsc_color(stripcolor,100), errthick=!p.thick, psym=symcat(stripsym), $
      symsize=strippsize
    if (nchar ne 0) then oploterror, final_pt05[char].rr25*rchar, final_pt05[char].oh, $
      final_pt05[char].oh_err, color=fsc_color(charcolor,100), $
      errcolor=fsc_color(charcolor,100), errthick=!p.thick, psym=symcat(charsym), $
      symsize=charpsize
    if (ncent ne 0) then oploterror, (final_pt05[cent].rr25+0.003)*rcent, final_pt05[cent].oh, $
      final_pt05[cent].oh_err, color=fsc_color(centcolor,100), $
      errcolor=fsc_color(centcolor,100), errthick=!p.thick, psym=symcat(centsym), $
      symsize=centpsize
    if (navg ne 0) then oploterror, final_pt05[avg].rr25*ravg, final_pt05[avg].oh, $
      final_pt05[avg].oh_err, color=fsc_color(avgcolor,100), $
      errcolor=fsc_color(avgcolor,100), errthick=!p.thick, psym=symcat(avgsym), $
      symsize=avgpsize

;   im_legend, ['Nuclear','Circumnuclear','Radial Strip','HII: Average',$
;     'HII: \rho=0.4\rho_{25}'], /left, /bottom, box=0, $
;     margin=0, symsize=1.5, spacing=1.8, charsize=1.3, $
;     psym=[nuclearsym,circumsym,stripsym,avgsym,charsym], $
;     color=[nuclearcolor,circumcolor,stripcolor,avgcolor,charcolor]

    im_legend, ['Nuclear','Circumnuclear','Radial Strip','HII: Average',$
      'HII: \rho=0.4\rho_{25}','HII: \rho=0'], /left, /bottom, box=0, $
      margin=0, spacing=1.8, charsize=1.3, $
      psym=[nuclearsym,circumsym,stripsym,avgsym,charsym,centsym], $
      symsize=[nuclearpsize,circumpsize,strippsize,avgpsize,charpsize,centpsize], $
      color=[nuclearcolor,circumcolor,stripcolor,avgcolor,charcolor,centcolor]

    im_plotconfig, /psclose

; ---------------------------------------------------------------------------
; plot the B-band LZ relation    

    psfile = pspath+'sings_lz_ohfinal'+suffix
    im_plotconfig, 12, pos, psfile=psfile, width=4.2*[1,1], $
      height=4.2, charsize=1.8, charthick=3.0

    xtitle = 'M_{B}'
    ytitle = '12 + log (O/H)_{final}'

    xrange = [-12.8,-22.5]
    yrange = [7.3,9.6]

; KK04    
    indx = where((final.mb gt -900.0) and (final.log12oh_kk04_char[0] gt -900.0),nindx)
    x = final[indx].mb
    xerr = final[indx].mb_err
    y = final[indx].log12oh_kk04_char[0]
    yerr = final[indx].log12oh_kk04_char[1]

    early = where((final[indx].t le 1.0),nearly) ; E-->Sa
    mid   = where((final[indx].t ge 2.0) and (final[indx].t le 4.0),nmid) ; Sab-->Sbc
    late  = where((final[indx].t ge 5.0),nlate) ; Sc-->Im

    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, position=pos[*,0]
    lzfit_kk04 = poly(mbaxis-mbpivot,lzfit[0].coeff)
    djs_oplot, mbaxis, lzfit_kk04, line=0
    djs_oplot, mbaxis, lzfit_kk04+lzfit[0].scatter, line=1
    djs_oplot, mbaxis, lzfit_kk04-lzfit[0].scatter, line=1
    
    oploterror, x[mid], y[mid], xerr[mid], yerr[mid], $
      color=fsc_color(midcolor1,100), errcolor=fsc_color(midcolor1,100), $
      errthick=errthick1, psym=symcat(midsym1,thick=symthick1), symsize=midpsize1
    oploterror, x[late], y[late], xerr[late], yerr[late], $
      color=fsc_color(latecolor1,100), errcolor=fsc_color(latecolor1,100), $
      errthick=errthick1, psym=symcat(latesym1,thick=symthick1), symsize=latepsize1
    oploterror, x[early], y[early], xerr[early], yerr[early], $
      color=fsc_color(earlycolor1,100), errcolor=fsc_color(earlycolor1,100), $
      errthick=errthick1, psym=symcat(earlysym1,thick=symthick1), symsize=earlypsize1

    legend, 'KK04', /left, /top, box=0, charsize=1.5, margin=0
    
    im_legend, ['E/S0','Sa-Sbc','Sc-Im'], psym=[earlysym1,midsym1,latesym1], $
      /right, /bottom, box=0, fill=[1,1,1], spacing=1.8, symsize=[1.8,1.5,2.1], $
      color=[earlycolor1,midcolor1,latecolor1], symthick=symthick1

    xyouts, -14.5, 9.1, 'Ho IX', align=0.5, /data, charsize=1.1
    
; PT05    
    indx = where((final.mb gt -900.0) and (final.log12oh_pt05_char[0] gt -900.0),nindx)
    x = final[indx].mb
    xerr = final[indx].mb_err
    y = final[indx].log12oh_pt05_char[0]
    yerr = final[indx].log12oh_pt05_char[1]

    early = where((final[indx].t le 1.0),nearly) ; E-->Sa
    mid   = where((final[indx].t ge 2.0) and (final[indx].t le 4.0),nmid) ; Sab-->Sbc
    late  = where((final[indx].t ge 5.0),nlate) ; Sc-->Im

    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle, ytitle='', xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, position=pos[*,1], $
      ytickname=replicate(' ',10)
    legend, 'PT05', /left, /top, box=0, charsize=1.5, margin=0

    lzfit_pt05 = poly(mbaxis-mbpivot,lzfit[1].coeff)
    djs_oplot, mbaxis, lzfit_pt05, line=0
    djs_oplot, mbaxis, lzfit_pt05+lzfit[1].scatter, line=1
    djs_oplot, mbaxis, lzfit_pt05-lzfit[1].scatter, line=1

    oploterror, x[mid], y[mid], xerr[mid], yerr[mid], $
      color=fsc_color(midcolor1,100), errcolor=fsc_color(midcolor1,100), $
      errthick=errthick1, psym=symcat(midsym1,thick=symthick1), symsize=midpsize1
    oploterror, x[late], y[late], xerr[late], yerr[late], $
      color=fsc_color(latecolor1,100), errcolor=fsc_color(latecolor1,100), $
      errthick=errthick1, psym=symcat(latesym1,thick=symthick1), symsize=latepsize1
    oploterror, x[early], y[early], xerr[early], yerr[early], $
      color=fsc_color(earlycolor1,100), errcolor=fsc_color(earlycolor1,100), $
      errthick=errthick1, psym=symcat(earlysym1,thick=symthick1), symsize=earlypsize1

    im_plotconfig, /psclose

return
end
    
