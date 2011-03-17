pro plotsings_log12oh_hiiregions, ps=ps
; jm10mar10ucsd - build the HII-region plots
    
    metpath = sings_path(/projects)+'log12oh/'
    pspath = sings_path(/papers)+'log12oh/FIG_LOG12OH/'

    version = sings_log12oh_version()
    galinfo = mrdfits(metpath+'sings_log12oh_'+version+'.fits.gz',1)
    hiifile = metpath+'sings_log12oh_hiiregions_'+version+'.fits.gz'
    splog, 'Reading '+hiifile
    hiiinfo = mrdfits(hiifile,1)
    hii = mrdfits(hiifile,2)
    
    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

    ohaxis = findgen((10.0-0.0)/0.01)*0.01+0.0
    rr25axis = findgen((1.05-0.0)/0.01)*0.01+0.0
    rr25range = [-0.15,1.2]
    sloperange = [-1.3,0.3]
;   sloperange = [-1.3,0.3]
    yresidrange = 0.46*[-1,1]
    mbrange = [-21.9,-18.1]

    kk04psize1 = 0.7 & kk04psize2 = 1.1 & kk04psize3 = 2.2 & kk04psize4 = 1.7
    pt05psize1 = 0.7 & pt05psize2 = 1.1 & pt05psize3 = 2.2 & pt05psize4 = 1.7
    kk04sym1 = 15 & kk04color1 = 'firebrick' ; 'dark green' ; 
    kk04sym2 =  6 & kk04color2 = 'firebrick' ; 'dark green' ; 
    pt05sym1 = 16 & pt05color1 = 'blue' ; 'orchid' ; 
    pt05sym2 =  9 & pt05color2 = 'blue' ; 'orchid' ; 
;   pt05sym1 = 16 & pt05color1 = 'royal blue' ; 'orchid' ; 
;   pt05sym2 =  9 & pt05color2 = 'royal blue' ; 'orchid' ; 

    earlysym = 16 & earlypsize = 1.8 & earlycolor = 'navy' ; 'dodger blue'
    midsym   = 15 & midpsize   = 2.0 & midcolor   = 'tan' ; 'salmon' ; 'black'
    latesym  = 14 & latepsize  = 2.6 & latecolor  = 'firebrick' ; tan
    earlysym2 = 9 & midsym2 = 6 & latesym2 = 4

    symthick1 = 5
    errthick1 = 5
    
; ---------------------------------------------------------------------------
; Slope [KK04] vs Slope [PT05]

    psfile = pspath+'slope_kk04_vs_slope_pt05'+suffix
    im_plotconfig, 12, pos, psfile=psfile, width=4.2*[1,1], height=4.2, $
      charsize=1.6, xspace=1.2

    plotsym, 0, 1.3, /fill

    indx = where((hiiinfo.hii_kk04_slope[0] gt -900.0) and $
      (hiiinfo.hii_pt05_slope[0] gt -900.0),nindx)

    early = where((hiiinfo[indx].t ge 2.0) and (hiiinfo[indx].t le 3.0),nearly)
    mid   = where((hiiinfo[indx].t ge 4.0) and (hiiinfo[indx].t le 5.0),nmid)
    late  = where((hiiinfo[indx].t ge 6.0) and (hiiinfo[indx].t le 7.0),nlate)

; compare the characteristic abundances
    x = hiiinfo[indx].hii_kk04_log12oh_char[0]
    xerr = hiiinfo[indx].hii_kk04_log12oh_char[1]
    y = hiiinfo[indx].hii_pt05_log12oh_char[0]
    yerr = hiiinfo[indx].hii_pt05_log12oh_char[1]
    galaxy = strtrim(hiiinfo[indx].galaxy,2)
    type = strtrim(hiiinfo[indx].type,2)

    resid = abs(y-x) & srt = sort(resid)
    splog, 'Residual statistics: characteristic 12+log(O/H)'
    resid_stats = im_stats(resid,/verbose)
;   hi = where(x gt 8.97,comp=lo)
;   resid_stats = im_stats(resid[hi],/verbose)
;   resid_stats = im_stats(resid[lo],/verbose)
    splog, 'Residuals: characteristic 12+log(O/H)'
    for ii = 0, nindx-1 do print, galaxy[srt[ii]], type[srt[ii]], y[srt[ii]], $
      yerr[srt[ii]], x[srt[ii]], xerr[srt[ii]], resid[srt[ii]], format='(A10,A10,5F12.5)'

    xrange = [8.75,9.25]
    yrange = xrange - mean(resid)
    splog, median(resid), mean(resid)

    xtitle = '12+log(O/H)_{char} [KK04]'
    ytitle = '12+log(O/H)_{char} [PT05]'

    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,0]
    oplot, xrange, yrange, line=2
; early
    bar = where((strmatch(hiiinfo[indx[early]].type,'*B*') eq 1L),nbar,comp=nobar,ncomp=nnobar)
    if (nbar ne 0L) then begin
       oploterror, x[early[bar]], y[early[bar]], xerr[early[bar]], yerr[early[bar]], $
         color=fsc_color(earlycolor,100), errcolor=fsc_color(earlycolor,100), errthick=errthick1, $
         psym=symcat(earlysym2,thick=symthick1), symsize=earlypsize
    endif
    if (nnobar ne 0L) then begin
       oploterror, x[early[nobar]], y[early[nobar]], xerr[early[nobar]], yerr[early[nobar]], $
         color=fsc_color(earlycolor,100), errcolor=fsc_color(earlycolor,100), errthick=errthick1, $
         psym=symcat(earlysym), symsize=earlypsize
    endif
; mid
    bar = where((strmatch(hiiinfo[indx[mid]].type,'*B*') eq 1L),nbar,comp=nobar,ncomp=nnobar)
    if (nbar ne 0L) then begin
       oploterror, x[mid[bar]], y[mid[bar]], xerr[mid[bar]], yerr[mid[bar]], $
         color=fsc_color(midcolor,100), errcolor=fsc_color(midcolor,100), errthick=errthick1, $
         psym=symcat(midsym2,thick=symthick1), symsize=midpsize
    endif
    if (nnobar ne 0L) then begin
       oploterror, x[mid[nobar]], y[mid[nobar]], xerr[mid[nobar]], yerr[mid[nobar]], $
         color=fsc_color(midcolor,100), errcolor=fsc_color(midcolor,100), errthick=errthick1, $
         psym=symcat(midsym,thick=symthick1), symsize=midpsize
    endif
; late
    bar = where((strmatch(hiiinfo[indx[late]].type,'*B*') eq 1L),nbar,comp=nobar,ncomp=nnobar)
    if (nbar ne 0L) then begin
       oploterror, x[late[bar]], y[late[bar]], xerr[late[bar]], yerr[late[bar]], $
         color=fsc_color(latecolor,100), errcolor=fsc_color(latecolor,100), errthick=errthick1, $
         psym=symcat(latesym2,thick=symthick1), symsize=latepsize
    endif
    if (nnobar ne 0L) then begin
       oploterror, x[late[nobar]], y[late[nobar]], xerr[late[nobar]], yerr[late[nobar]], $
         color=fsc_color(latecolor,100), errcolor=fsc_color(latecolor,100), errthick=errthick1, $
         psym=symcat(latesym,thick=symthick1), symsize=latepsize
    endif

    im_legend, ['Sab-Sb','Sbc-Sc','Scd-Sd'], psym=[16,15,14], $
      /right, /bottom, box=0, fill=[1,1,1], charsize=1.6, $
      spacing=1.8, color=[earlycolor,midcolor,latecolor], $
      symsize=0.85*[earlypsize,midpsize,latepsize], margin=0

; compare the slopes    
    x = hiiinfo[indx].hii_kk04_slope[0]
    xerr = hiiinfo[indx].hii_kk04_slope[1]
    y = hiiinfo[indx].hii_pt05_slope[0]
    yerr = hiiinfo[indx].hii_pt05_slope[1]
    galaxy = strtrim(hiiinfo[indx].galaxy,2)
    type = strtrim(hiiinfo[indx].type,2)

    resid = y-x & srt = sort(abs(resid))
;   resid = abs(y-x) & srt = sort(resid)
    splog, 'Residual statistics: slope vs slope'
    good = where((galaxy ne 'NGC3621'))
;   good = where((galaxy ne 'NGC3621') and (galaxy ne 'NGC7331') and (galaxy ne 'NGC3521'))
    resid_stats = im_stats(resid[good],/verbose)
;   resid_stats = im_stats(resid,/verbose,sigrej=1.5,mask=mask)
    splog, 'Residuals: slope vs slope'
    for ii = 0, nindx-1 do print, galaxy[srt[ii]], type[srt[ii]], y[srt[ii]], $
      yerr[srt[ii]], x[srt[ii]], xerr[srt[ii]], resid[srt[ii]], format='(A10,A10,5F12.5)'

    rcor = r_correlate(x[good],y[good],zd=zd,probd=probd)
    splog, 'Spearman rank test: ', rcor, probd, zd

; some stats:    
    good = lindgen(nindx)
    mn = sings_weighted_mean(y[good],yerr[good],wmean=mn_err,wquant=med,quant=0.5)
    splog, 'PT05 slope: min, max, mean, median, sigma, weighted median, mean, and error: '
    splog, minmax(y[good]), djs_mean(y[good]), djs_median(y[good]), djsig(y[good]), med, mn, mn_err

    mn = sings_weighted_mean(x[good],xerr[good],wmean=mn_err,wquant=med,quant=0.5)
    splog, 'KK04 slope: min, max, mean, median, sigma, weighted median, mean, and error: '
    splog, minmax(x[good]), djs_mean(x[good]), djs_median(x[good]), djsig(x[good]), med, mn, mn_err
    
    good = where($
      (strmatch(galaxy,'*3521*') eq 0) and $
      (strmatch(galaxy,'*7331*') eq 0) and $
      (strmatch(galaxy,'*2841*') eq 0),ngood)
;     (strmatch(galaxy,'*3621*') eq 0),ngood)
    resid = y[good]-x[good]
    resid_err = sqrt(yerr[good]^2+xerr[good]^2)
    mn = sings_weighted_mean(resid,resid_err,wmean_err=mn_err,wquant=med,quant=0.5)
    splog, 'Slope difference: mean, median, sigma, weighted median, mean, and error:'
    splog, djs_mean(resid), djs_median(resid), djsig(resid), $
      djsig(resid)/sqrt(ngood), med, mn, mn_err
    
    xrange = sloperange
    yrange = xrange

    xtitle = 'Gradient Slope (dex \rho_{25}^{-1}) [KK04]'
    ytitle = 'Gradient Slope (dex \rho_{25}^{-1}) [PT05]'

    djs_plot, [0], [0], /nodata, /noerase, $
      xtitle=xtitle, ytitle=ytitle, $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,1]
    oplot, !x.crange, !y.crange, line=0
; early
    bar = where((strmatch(hiiinfo[indx[early]].type,'*B*') eq 1L),nbar,comp=nobar,ncomp=nnobar)
    if (nbar ne 0L) then begin
       oploterror, x[early[bar]], y[early[bar]], xerr[early[bar]], yerr[early[bar]], $
         color=fsc_color(earlycolor,100), errcolor=fsc_color(earlycolor,100), errthick=errthick1, $
         psym=symcat(earlysym2,thick=symthick1), symsize=earlypsize
;      n4725 = speclinefit_locate(hiiinfo[indx[early[bar]]],'ngc4725')
;      xyouts, x[early[bar[n4725]]]+0.03, y[early[[bar[n4725]]]]+0.05, 'NGC4725', $
;        align=0.0, charsize=1.2
    endif
    if (nnobar ne 0L) then begin
       oploterror, x[early[nobar]], y[early[nobar]], xerr[early[nobar]], yerr[early[nobar]], $
         color=fsc_color(earlycolor,100), errcolor=fsc_color(earlycolor,100), errthick=errthick1, $
         psym=symcat(earlysym), symsize=earlypsize
;      n7331 = speclinefit_locate(hiiinfo[indx[early[nobar]]],'ngc7331')
;      xyouts, x[early[nobar[n7331]]]+0.05, y[early[[nobar[n7331]]]]+0.05, 'NGC7331', $
;        align=1.0, orientation=-45, charsize=1.2
    endif
; mid
    bar = where((strmatch(hiiinfo[indx[mid]].type,'*B*') eq 1L),nbar,comp=nobar,ncomp=nnobar)
    if (nbar ne 0L) then begin
       oploterror, x[mid[bar]], y[mid[bar]], xerr[mid[bar]], yerr[mid[bar]], $
         color=fsc_color(midcolor,100), errcolor=fsc_color(midcolor,100), errthick=errthick1, $
         psym=symcat(midsym2,thick=symthick1), symsize=midpsize
;      n3521 = speclinefit_locate(hiiinfo[indx[mid[bar]]],'ngc3521')
;      xyouts, x[mid[bar[n3521]]]-0.05, y[mid[[bar[n3521]]]]-0.05, 'NGC3521', $
;        align=1.0, orientation=45, charsize=1.2
    endif
    if (nnobar ne 0L) then begin
       oploterror, x[mid[nobar]], y[mid[nobar]], xerr[mid[nobar]], yerr[mid[nobar]], $
         color=fsc_color(midcolor,100), errcolor=fsc_color(midcolor,100), errthick=errthick1, $
         psym=symcat(midsym,thick=symthick1), symsize=midpsize
    endif
; late
    bar = where((strmatch(hiiinfo[indx[late]].type,'*B*') eq 1L),nbar,comp=nobar,ncomp=nnobar)
    if (nbar ne 0L) then begin
       oploterror, x[late[bar]], y[late[bar]], xerr[late[bar]], yerr[late[bar]], $
         color=fsc_color(latecolor,100), errcolor=fsc_color(latecolor,100), errthick=errthick1, $
         psym=symcat(latesym2,thick=symthick1), symsize=latepsize
    endif
    if (nnobar ne 0L) then begin
       oploterror, x[late[nobar]], y[late[nobar]], xerr[late[nobar]], yerr[late[nobar]], $
         color=fsc_color(latecolor,100), errcolor=fsc_color(latecolor,100), errthick=errthick1, $
         psym=symcat(latesym,thick=symthick1), symsize=latepsize
       n3621 = speclinefit_locate(hiiinfo[indx[late[nobar]]],'ngc3621')
       xyouts, x[late[nobar[n3621]]]-0.01, y[late[[nobar[n3621]]]]-0.1, 'NGC3621', $
         align=1.0, charsize=1.2
    endif

    im_plotconfig, /psclose

; ---------------------------------------------------------------------------        
; multi-panel abundance gradient plot

    grad = where(hiiinfo.gradient_flag,ngrad)

    psfile = pspath+'sings_gradients'+suffix
    im_plotconfig, 21, pos, psfile=psfile, charsize=1.2

    ncols = 4
    nrows = 6
    
    xtitle = '\rho / \rho_{25}'       
    ytitle = '12 + log (O/H)'

    xrange = rr25range
    yrange = [7.6,9.6]

    for ii = 0, ngrad-1 do begin

       galaxy = strtrim(hiiinfo[grad[ii]].galaxy,2)
       type = strtrim(repstr(hiiinfo[grad[ii]].type,'~pec',''),2)
       match = where(galaxy eq strtrim(hii.galaxy,2),nmatch)
       thesehii = hii[match]
       niceprint, galaxy, nmatch;, thesehii.reference

; now make the plot       
       if (ii mod ncols) eq 0L then begin
          delvarx, ytickname
          ytitle1 = ytitle
       endif else begin
          ytickname = replicate(' ',10)
          ytitle1 = ''
       endelse

       if (ii ge 17) then begin
;      if (ii ge 18) then begin
          delvarx, xtickname
          xtitle1 = xtitle
       endif else begin
          xtickname = replicate(' ',10)
          xtitle1 = ''
       endelse
       
       djs_plot, [0], [0], /nodata, xtitle=xtitle1, ytitle=ytitle1, $
         xstyle=1, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,ii], $
         ytickname=ytickname, xtickname=xtickname, noerase=(ii gt 0L), xtickinterval=0.5
       legend, galaxy+' ('+type+')', /left, /bottom, box=0, charsize=1.0, $
         spacing=-1.0, margin=-0.5;, /clear

; #########################       
; PT05 - shaded region showing the 1-sigma range of lines (see
; SINGS_LOG12OH_HIIREGIONS)
       covar = hiiinfo[grad[ii]].hii_pt05_gradient_covar
       coeff = [hiiinfo[grad[ii]].hii_pt05_log12oh_nuclear[0],hiiinfo[grad[ii]].hii_pt05_slope[0]]
       rand = mrandomn(seed,covar,1000)
       rand[*,0] = rand[*,0] + coeff[0]
       rand[*,1] = rand[*,1] + coeff[1]
       
       ell = covar2ellipse(covar,nsigma=1.0) ; get models within 1-sigma
       indx = get_ellipse_indices(rand[*,0],rand[*,1],$
         major=ell.major,minor=ell.minor,angle=ell.angle, $
         xcenter=coeff[0],ycenter=coeff[1])
       nindx = n_elements(indx)

       allgrad = fltarr(n_elements(rr25axis),nindx)
       for jj = 0, nindx-1 do allgrad[*,jj] = poly(rr25axis,rand[indx[jj],*])
       mingrad = min(allgrad,dim=2)
       maxgrad = max(allgrad,dim=2)
       polyfill, [rr25axis,reverse(rr25axis)],[mingrad,reverse(maxgrad)], $
         /data, /fill, color=fsc_color('powder blue',100), noclip=0
       
; KK04 - plot a shaded region showing the 1-sigma range of lines (see
; SINGS_LOG12OH_HIIREGIONS)
       covar = hiiinfo[grad[ii]].hii_kk04_gradient_covar
       coeff = [hiiinfo[grad[ii]].hii_kk04_log12oh_nuclear[0],hiiinfo[grad[ii]].hii_kk04_slope[0]]
       rand = mrandomn(seed,covar,1000)
       rand[*,0] = rand[*,0] + coeff[0]
       rand[*,1] = rand[*,1] + coeff[1]
       
       ell = covar2ellipse(covar,nsigma=1.0) ; get models within 1-sigma
       indx = get_ellipse_indices(rand[*,0],rand[*,1],$
         major=ell.major,minor=ell.minor,angle=ell.angle, $
         xcenter=coeff[0],ycenter=coeff[1])
       nindx = n_elements(indx)

       allgrad = fltarr(n_elements(rr25axis),nindx)
       for jj = 0, nindx-1 do allgrad[*,jj] = poly(rr25axis,rand[indx[jj],*])
       mingrad = min(allgrad,dim=2)
       maxgrad = max(allgrad,dim=2)
       polyfill, [rr25axis,reverse(rr25axis)],[mingrad,reverse(maxgrad)], $
         /data, /fill, color=fsc_color('light salmon',100), noclip=0

; KK04 - overplot the data and the best-fitting line
       good = where((thesehii.log12oh_kk04[0] gt -900.0) and $
         (strtrim(thesehii.r23_branch_kk04,2) ne 'A') and $
         (strtrim(thesehii.r23_branch_kk04,2) ne 'X'))
       ambig = where((thesehii.log12oh_kk04[0] gt -900.0) and $
         ((strtrim(thesehii.r23_branch_kk04,2) eq 'A') or $
         (strtrim(thesehii.r23_branch_kk04,2) eq 'X')),nambig)
       oploterror, thesehii[good].rc3_rr25, thesehii[good].log12oh_kk04[0], $
         thesehii[good].log12oh_kk04[1], errthick=errthick1, color=fsc_color(kk04color1,100), $
         errcolor=fsc_color(kk04color1,100), psym=symcat(kk04sym1), symsize=kk04psize1
       if (nambig eq 1) then plots, thesehii[ambig].rc3_rr25, thesehii[ambig].log12oh_kk04[0], $
         color=fsc_color(kk04color2,100), psym=symcat(kk04sym2,thick=4), symsize=kk04psize1
       if (nambig gt 1) then djs_oplot, thesehii[ambig].rc3_rr25, thesehii[ambig].log12oh_kk04[0], $
         color=fsc_color(kk04color2,100), psym=symcat(kk04sym2,thick=4), symsize=kk04psize1
       coeff = [hiiinfo[grad[ii]].hii_kk04_log12oh_nuclear[0],hiiinfo[grad[ii]].hii_kk04_slope[0]]
       djs_oplot, rr25axis, poly(rr25axis,coeff), line=0

; PT05 - overplot the data and the best-fitting line
       good = where((thesehii.log12oh_pt05[0] gt -900.0) and $
         (strtrim(thesehii.r23_branch_pt05,2) ne 'A') and $
         (strtrim(thesehii.r23_branch_pt05,2) ne 'X'))
       ambig = where((thesehii.log12oh_pt05[0] gt -900.0) and $
         ((strtrim(thesehii.r23_branch_pt05,2) eq 'A') or $
         (strtrim(thesehii.r23_branch_pt05,2) eq 'X')),nambig)
       oploterror, thesehii[good].rc3_rr25, thesehii[good].log12oh_pt05[0], $
         thesehii[good].log12oh_pt05[1], errthick=errthick1, color=fsc_color(pt05color1,100), $
         errcolor=fsc_color(pt05color1,100), psym=symcat(pt05sym1), symsize=pt05psize1
       if (nambig eq 1) then plots, thesehii[ambig].rc3_rr25, thesehii[ambig].log12oh_pt05[0], $
         color=fsc_color(pt05color2,100), psym=symcat(pt05sym2,thick=4), symsize=pt05psize1
       if (nambig gt 1) then djs_oplot, thesehii[ambig].rc3_rr25, thesehii[ambig].log12oh_pt05[0], $
         color=fsc_color(pt05color2,100), psym=symcat(pt05sym2,thick=4), symsize=pt05psize1
       coeff = [hiiinfo[grad[ii]].hii_pt05_log12oh_nuclear[0],hiiinfo[grad[ii]].hii_pt05_slope[0]]
       djs_oplot, rr25axis, poly(rr25axis,coeff), line=2
    endfor
    
    im_plotconfig, /psclose

; ---------------------------------------------------------------------------
; Slope vs type

    psfile = pspath+'slope_vs_type'+suffix
    im_plotconfig, 12, pos, psfile=psfile, width=4.2*[1,1], $
      height=4.2, charsize=1.8, charthick=3.0, xmargin=[1.3,0.4]

    indx = where((hiiinfo.hii_kk04_slope[0] gt -900.0) and (hiiinfo.hii_pt05_slope[0] gt -900.0),nindx)
    early = where((hiiinfo[indx].t ge 2.0) and (hiiinfo[indx].t le 3.0),nearly)
    mid   = where((hiiinfo[indx].t ge 4.0) and (hiiinfo[indx].t le 5.0),nmid)
    late  = where((hiiinfo[indx].t ge 6.0) and (hiiinfo[indx].t le 7.0),nlate)

    x = galinfo[indx].t
    xerr = x*0.0
    galaxy = strtrim(hiiinfo[indx].galaxy,2)
    type = strtrim(hiiinfo[indx].type,2)
    ut = x[uniq(x,sort(x))]

    xrange = [1.5,7.5]
    yrange = sloperange

    xtitle = 'Hubble Type'
    ytitle = 'Gradient Slope (dex \rho_{25}^{-1})'

; -------------------------    
; KK04
    y = hiiinfo[indx].hii_kk04_slope[0]
    yerr = hiiinfo[indx].hii_kk04_slope[1]

    djs_plot, [0], [0], /nodata, position=pos[*,0], $
      xtitle=xtitle, ytitle=ytitle, $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange
    legend, 'KK04', /left, /top, box=0, margin=0, charsize=1.5
    
; early
    bar = where((strmatch(hiiinfo[indx[early]].type,'*B*') eq 1L),nbar,comp=nobar,ncomp=nnobar)
    if (nbar ne 0L) then begin
       oploterror, x[early[bar]], y[early[bar]], xerr[early[bar]], yerr[early[bar]], $
         color=fsc_color(earlycolor,100), errcolor=fsc_color(earlycolor,100), errthick=errthick1, $
         psym=symcat(earlysym2,thick=symthick1), symsize=earlypsize
    endif
    if (nnobar ne 0L) then begin
       oploterror, x[early[nobar]], y[early[nobar]], xerr[early[nobar]], yerr[early[nobar]], $
         color=fsc_color(earlycolor,100), errcolor=fsc_color(earlycolor,100), errthick=errthick1, $
         psym=symcat(earlysym), symsize=earlypsize
    endif
; mid
    bar = where((strmatch(hiiinfo[indx[mid]].type,'*B*') eq 1L),nbar,comp=nobar,ncomp=nnobar)
    if (nbar ne 0L) then begin
       oploterror, x[mid[bar]], y[mid[bar]], xerr[mid[bar]], yerr[mid[bar]], $
         color=fsc_color(midcolor,100), errcolor=fsc_color(midcolor,100), errthick=errthick1, $
         psym=symcat(midsym2,thick=symthick1), symsize=midpsize
    endif
    if (nnobar ne 0L) then begin
       oploterror, x[mid[nobar]], y[mid[nobar]], xerr[mid[nobar]], yerr[mid[nobar]], $
         color=fsc_color(midcolor,100), errcolor=fsc_color(midcolor,100), errthick=errthick1, $
         psym=symcat(midsym,thick=symthick1), symsize=midpsize
    endif
; late
    bar = where((strmatch(hiiinfo[indx[late]].type,'*B*') eq 1L),nbar,comp=nobar,ncomp=nnobar)
    if (nbar ne 0L) then begin
       oploterror, x[late[bar]], y[late[bar]], xerr[late[bar]], yerr[late[bar]], $
         color=fsc_color(latecolor,100), errcolor=fsc_color(latecolor,100), errthick=errthick1, $
         psym=symcat(latesym2,thick=symthick1), symsize=latepsize
    endif
    if (nnobar ne 0L) then begin
       oploterror, x[late[nobar]], y[late[nobar]], xerr[late[nobar]], yerr[late[nobar]], $
         color=fsc_color(latecolor,100), errcolor=fsc_color(latecolor,100), errthick=errthick1, $
         psym=symcat(latesym,thick=symthick1), symsize=latepsize
    endif

; compute the average in bins of T    
    mn = fltarr(n_elements(ut))
    mn_err = fltarr(n_elements(ut))
    for ii = 0, n_elements(ut)-1 do begin
       these = where(ut[ii] eq x)
       mn[ii] = sings_weighted_mean(y[these],yerr[these],wmean_err=mn_err1)
       mn_err[ii] = mn_err1
    endfor
    oploterror, ut-0.2, mn, mn_err, psym=-symcat(6,thick=4), symsize=2.0
    
; -------------------------    
; pt05
    y = hiiinfo[indx].hii_pt05_slope[0]
    yerr = hiiinfo[indx].hii_pt05_slope[1]

    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], $
      xtitle=xtitle, ytitle='', ytickname=replicate(' ',10), $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange
    legend, 'pt05', /left, /top, box=0, margin=0, charsize=1.5

    im_legend, ['sab-sb','sbc-sc','scd-sd'], psym=[16,15,14], $
      /right, /bottom, box=0, fill=[1,1,1], charsize=1.6, $
      spacing=1.8, color=[earlycolor,midcolor,latecolor], $
      symsize=0.85*[earlypsize,midpsize,latepsize], margin=0

; early
    bar = where((strmatch(hiiinfo[indx[early]].type,'*b*') eq 1l),nbar,comp=nobar,ncomp=nnobar)
    if (nbar ne 0l) then begin
       oploterror, x[early[bar]], y[early[bar]], xerr[early[bar]], yerr[early[bar]], $
         color=fsc_color(earlycolor,100), errcolor=fsc_color(earlycolor,100), errthick=errthick1, $
         psym=symcat(earlysym2,thick=symthick1), symsize=earlypsize
    endif
    if (nnobar ne 0l) then begin
       oploterror, x[early[nobar]], y[early[nobar]], xerr[early[nobar]], yerr[early[nobar]], $
         color=fsc_color(earlycolor,100), errcolor=fsc_color(earlycolor,100), errthick=errthick1, $
         psym=symcat(earlysym), symsize=earlypsize
    endif
; mid
    bar = where((strmatch(hiiinfo[indx[mid]].type,'*b*') eq 1l),nbar,comp=nobar,ncomp=nnobar)
    if (nbar ne 0l) then begin
       oploterror, x[mid[bar]], y[mid[bar]], xerr[mid[bar]], yerr[mid[bar]], $
         color=fsc_color(midcolor,100), errcolor=fsc_color(midcolor,100), errthick=errthick1, $
         psym=symcat(midsym2,thick=symthick1), symsize=midpsize
    endif
    if (nnobar ne 0l) then begin
       oploterror, x[mid[nobar]], y[mid[nobar]], xerr[mid[nobar]], yerr[mid[nobar]], $
         color=fsc_color(midcolor,100), errcolor=fsc_color(midcolor,100), errthick=errthick1, $
         psym=symcat(midsym,thick=symthick1), symsize=midpsize
    endif
; late
    bar = where((strmatch(hiiinfo[indx[late]].type,'*b*') eq 1l),nbar,comp=nobar,ncomp=nnobar)
    if (nbar ne 0l) then begin
       oploterror, x[late[bar]], y[late[bar]], xerr[late[bar]], yerr[late[bar]], $
         color=fsc_color(latecolor,100), errcolor=fsc_color(latecolor,100), errthick=errthick1, $
         psym=symcat(latesym2,thick=symthick1), symsize=latepsize
    endif
    if (nnobar ne 0l) then begin
       oploterror, x[late[nobar]], y[late[nobar]], xerr[late[nobar]], yerr[late[nobar]], $
         color=fsc_color(latecolor,100), errcolor=fsc_color(latecolor,100), errthick=errthick1, $
         psym=symcat(latesym,thick=symthick1), symsize=latepsize
    endif

; compute the average in bins of t    
    mn = fltarr(n_elements(ut))
    mn_err = fltarr(n_elements(ut))
    for ii = 0, n_elements(ut)-1 do begin
       these = where(ut[ii] eq x)
       mn[ii] = sings_weighted_mean(y[these],yerr[these],wmean_err=mn_err1)
       mn_err[ii] = mn_err1
    endfor
    oploterror, ut+0.2, mn, mn_err, psym=-symcat(6,thick=4), symsize=2.0
    
    im_plotconfig, /psclose

; ---------------------------------------------------------------------------
; Slope vs MB

    psfile = pspath+'slope_vs_mb'+suffix
    im_plotconfig, 12, pos, psfile=psfile, width=4.2*[1,1], $
      height=4.2, charsize=1.8, charthick=3.0, xmargin=[1.3,0.4]

    indx = where((hiiinfo.hii_kk04_slope[0] gt -900.0) and (hiiinfo.hii_pt05_slope[0] gt -900.0),nindx)
    early = where((hiiinfo[indx].t ge 2.0) and (hiiinfo[indx].t le 3.0),nearly)
    mid   = where((hiiinfo[indx].t ge 4.0) and (hiiinfo[indx].t le 5.0),nmid)
    late  = where((hiiinfo[indx].t ge 6.0) and (hiiinfo[indx].t le 7.0),nlate)

    x = galinfo[indx].mb
    xerr = galinfo[indx].mb_err
    galaxy = strtrim(hiiinfo[indx].galaxy,2)
    type = strtrim(hiiinfo[indx].type,2)

    xrange = mbrange
    yrange = sloperange

    xtitle = 'M_{B}'
    ytitle = 'Gradient Slope (dex \rho_{25}^{-1})'

; -------------------------    
; KK04
    y = hiiinfo[indx].hii_kk04_slope[0]
    yerr = hiiinfo[indx].hii_kk04_slope[1]

    djs_plot, [0], [0], /nodata, position=pos[*,0], $
      xtitle=xtitle, ytitle=ytitle, $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange
    legend, 'KK04', /left, /top, box=0, margin=0, charsize=1.5
    
; early
    bar = where((strmatch(hiiinfo[indx[early]].type,'*B*') eq 1L),nbar,comp=nobar,ncomp=nnobar)
    if (nbar ne 0L) then begin
       oploterror, x[early[bar]], y[early[bar]], xerr[early[bar]], yerr[early[bar]], $
         color=fsc_color(earlycolor,100), errcolor=fsc_color(earlycolor,100), errthick=errthick1, $
         psym=symcat(earlysym2,thick=symthick1), symsize=earlypsize
    endif
    if (nnobar ne 0L) then begin
       oploterror, x[early[nobar]], y[early[nobar]], xerr[early[nobar]], yerr[early[nobar]], $
         color=fsc_color(earlycolor,100), errcolor=fsc_color(earlycolor,100), errthick=errthick1, $
         psym=symcat(earlysym), symsize=earlypsize
    endif
; mid
    bar = where((strmatch(hiiinfo[indx[mid]].type,'*B*') eq 1L),nbar,comp=nobar,ncomp=nnobar)
    if (nbar ne 0L) then begin
       oploterror, x[mid[bar]], y[mid[bar]], xerr[mid[bar]], yerr[mid[bar]], $
         color=fsc_color(midcolor,100), errcolor=fsc_color(midcolor,100), errthick=errthick1, $
         psym=symcat(midsym2,thick=symthick1), symsize=midpsize
    endif
    if (nnobar ne 0L) then begin
       oploterror, x[mid[nobar]], y[mid[nobar]], xerr[mid[nobar]], yerr[mid[nobar]], $
         color=fsc_color(midcolor,100), errcolor=fsc_color(midcolor,100), errthick=errthick1, $
         psym=symcat(midsym,thick=symthick1), symsize=midpsize
    endif
; late
    bar = where((strmatch(hiiinfo[indx[late]].type,'*B*') eq 1L),nbar,comp=nobar,ncomp=nnobar)
    if (nbar ne 0L) then begin
       oploterror, x[late[bar]], y[late[bar]], xerr[late[bar]], yerr[late[bar]], $
         color=fsc_color(latecolor,100), errcolor=fsc_color(latecolor,100), errthick=errthick1, $
         psym=symcat(latesym2,thick=symthick1), symsize=latepsize
    endif
    if (nnobar ne 0L) then begin
       oploterror, x[late[nobar]], y[late[nobar]], xerr[late[nobar]], yerr[late[nobar]], $
         color=fsc_color(latecolor,100), errcolor=fsc_color(latecolor,100), errthick=errthick1, $
         psym=symcat(latesym,thick=symthick1), symsize=latepsize
    endif

; -------------------------    
; PT05
    y = hiiinfo[indx].hii_pt05_slope[0]
    yerr = hiiinfo[indx].hii_pt05_slope[1]

    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], $
      xtitle=xtitle, ytitle='', ytickname=replicate(' ',10), $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange
    legend, 'PT05', /left, /top, box=0, margin=0, charsize=1.5

    im_legend, ['Sab-Sb','Sbc-Sc','Scd-Sd'], psym=[16,15,14], $
      /right, /bottom, box=0, fill=[1,1,1], charsize=1.6, $
      spacing=1.8, color=[earlycolor,midcolor,latecolor], $
      symsize=0.85*[earlypsize,midpsize,latepsize], margin=0

; early
    bar = where((strmatch(hiiinfo[indx[early]].type,'*B*') eq 1L),nbar,comp=nobar,ncomp=nnobar)
    if (nbar ne 0L) then begin
       oploterror, x[early[bar]], y[early[bar]], xerr[early[bar]], yerr[early[bar]], $
         color=fsc_color(earlycolor,100), errcolor=fsc_color(earlycolor,100), errthick=errthick1, $
         psym=symcat(earlysym2,thick=symthick1), symsize=earlypsize
    endif
    if (nnobar ne 0L) then begin
       oploterror, x[early[nobar]], y[early[nobar]], xerr[early[nobar]], yerr[early[nobar]], $
         color=fsc_color(earlycolor,100), errcolor=fsc_color(earlycolor,100), errthick=errthick1, $
         psym=symcat(earlysym), symsize=earlypsize
    endif
; mid
    bar = where((strmatch(hiiinfo[indx[mid]].type,'*B*') eq 1L),nbar,comp=nobar,ncomp=nnobar)
    if (nbar ne 0L) then begin
       oploterror, x[mid[bar]], y[mid[bar]], xerr[mid[bar]], yerr[mid[bar]], $
         color=fsc_color(midcolor,100), errcolor=fsc_color(midcolor,100), errthick=errthick1, $
         psym=symcat(midsym2,thick=symthick1), symsize=midpsize
    endif
    if (nnobar ne 0L) then begin
       oploterror, x[mid[nobar]], y[mid[nobar]], xerr[mid[nobar]], yerr[mid[nobar]], $
         color=fsc_color(midcolor,100), errcolor=fsc_color(midcolor,100), errthick=errthick1, $
         psym=symcat(midsym,thick=symthick1), symsize=midpsize
    endif
; late
    bar = where((strmatch(hiiinfo[indx[late]].type,'*B*') eq 1L),nbar,comp=nobar,ncomp=nnobar)
    if (nbar ne 0L) then begin
       oploterror, x[late[bar]], y[late[bar]], xerr[late[bar]], yerr[late[bar]], $
         color=fsc_color(latecolor,100), errcolor=fsc_color(latecolor,100), errthick=errthick1, $
         psym=symcat(latesym2,thick=symthick1), symsize=latepsize
    endif
    if (nnobar ne 0L) then begin
       oploterror, x[late[nobar]], y[late[nobar]], xerr[late[nobar]], yerr[late[nobar]], $
         color=fsc_color(latecolor,100), errcolor=fsc_color(latecolor,100), errthick=errthick1, $
         psym=symcat(latesym,thick=symthick1), symsize=latepsize
    endif

    im_plotconfig, /psclose

; ---------------------------------------------------------------------------    
; Characteristic and Nuclear vs Average

    psfile = pspath+'sings_hiiavg_vs_hiichar_hiinuc'+suffix
    im_plotconfig, 5, pos, psfile=psfile, charsize=1.5, xmargin=[1.0,1.0], $
      width=3.25*[1,1], height=[2.75,1.75]

    xrange = [7.85,9.6]
    yrange = xrange

    yresidtitle = '\Delta log (O/H)'
    yresidrange = 0.56*[-1,1]

; avg vs nuclear
    
    pt05_indx = where((hiiinfo.hii_pt05_log12oh_char[0] gt -900.0) and $
      (hiiinfo.hii_pt05_log12oh_nuclear[0] gt -900.0) and $
      (hiiinfo.hii_pt05_log12oh_avg[0] gt -900.0),pt05_nindx)
    xpt05 = hiiinfo[pt05_indx].hii_pt05_log12oh_avg[0]
    ypt05 = hiiinfo[pt05_indx].hii_pt05_log12oh_nuclear[0]
    xpt05_err = hiiinfo[pt05_indx].hii_pt05_log12oh_avg[1]
    ypt05_err = hiiinfo[pt05_indx].hii_pt05_log12oh_nuclear[1]
    yresidpt05 = ypt05-xpt05
    yresidpt05_err = sqrt(ypt05_err^2+xpt05_err^2)
    
    splog, 'PT05: Average vs nuclear'
    stats = im_stats(yresidpt05)
    xstr_pt05 = strtrim(string(stats.median,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma,format='(F12.2)'),2)+')'
    splog, xstr_pt05
    
    kk04_indx = where((hiiinfo.hii_kk04_log12oh_char[0] gt -900.0) and $
      (hiiinfo.hii_kk04_log12oh_nuclear[0] gt -900.0) and $
      (hiiinfo.hii_kk04_log12oh_avg[0] gt -900.0),kk04_nindx)
    xkk04 = hiiinfo[kk04_indx].hii_kk04_log12oh_avg[0]
    ykk04 = hiiinfo[kk04_indx].hii_kk04_log12oh_nuclear[0]
    xkk04_err = hiiinfo[kk04_indx].hii_kk04_log12oh_avg[1]
    ykk04_err = hiiinfo[kk04_indx].hii_kk04_log12oh_nuclear[1]
    yresidkk04 = ykk04-xkk04
    yresidkk04_err = sqrt(ykk04_err^2+xkk04_err^2)
    
    splog, 'KK04: Average vs nuclear'
    stats = im_stats(yresidkk04)
    xstr_kk04 = strtrim(string(stats.median,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma,format='(F12.2)'),2)+')'
    splog, xstr_kk04

    stats = im_stats([yresidkk04,yresidpt05])
    xstr = '<\Delta>='+strtrim(string(stats.mean,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma,format='(F12.2)'),2)
;   xstr = '<\Delta>='+strtrim(string(stats.median,format='(F12.2)'),2)+' ('+$
;     strtrim(string(stats.mean,format='(F12.2)'),2)+$
;     '\pm'+strtrim(string(stats.sigma,format='(F12.2)'),2)+')'
    
    xtitle = '12 + log (O/H)_{HII avg}'
    ytitle = '12 + log (O/H)_{\rho=0}'
    
    djs_plot, [0], [0], /nodata, xtitle='', ytitle=ytitle, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, position=pos[*,0], xtickname=replicate(' ',10)
    djs_oplot, ohaxis, ohaxis, line=0
;   legend, '(a)', /left, /top, box=0, charsize=charsize_5, charthick=postthick2, margin=0

    oploterror, xkk04, ykk04, xkk04_err, ykk04_err, $
      errthick=errthick1, color=fsc_color(kk04color1,100), errcolor=fsc_color(kk04color1,100), $
      psym=symcat(kk04sym1), symsize=kk04psize2
    oploterror, xpt05, ypt05, xpt05_err, ypt05_err, $
      errthick=errthick1, color=fsc_color(pt05color1,100), errcolor=fsc_color(pt05color1,100), $
      psym=symcat(pt05sym1), symsize=pt05psize2

    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle, ytitle=yresidtitle, xsty=1, ysty=1, $
      xrange=xrange, yrange=yresidrange, position=pos[*,2]
    djs_oplot, ohaxis, ohaxis*0.0, line=0
    legend, textoidl(xstr), /right, /bottom, box=0, margin=0, charsize=1.2

    oplot, xkk04, yresidkk04, color=fsc_color(kk04color1,100), $
      psym=symcat(kk04sym1), symsize=kk04psize2
    oplot, xpt05, yresidpt05, color=fsc_color(pt05color1,100), $
      psym=symcat(pt05sym1), symsize=pt05psize2

; avg vs char
    pt05_indx = where((hiiinfo.hii_pt05_log12oh_char[0] gt -900.0) and $
      (hiiinfo.hii_pt05_log12oh_nuclear[0] gt -900.0) and $
      (hiiinfo.hii_pt05_log12oh_avg[0] gt -900.0),pt05_nindx)
    xpt05 = hiiinfo[pt05_indx].hii_pt05_log12oh_avg[0]
    ypt05 = hiiinfo[pt05_indx].hii_pt05_log12oh_char[0]
    xpt05_err = hiiinfo[pt05_indx].hii_pt05_log12oh_avg[1]
    ypt05_err = hiiinfo[pt05_indx].hii_pt05_log12oh_char[1]
    yresidpt05 = ypt05-xpt05
    yresidpt05_err = sqrt(ypt05_err^2+xpt05_err^2)

    splog, 'PT05: Average vs characteristic'
    stats = im_stats(yresidpt05)
    xstr_pt05 = strtrim(string(stats.median,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma,format='(F12.2)'),2)+')'
    splog, xstr_pt05

    kk04_indx = where((hiiinfo.hii_kk04_log12oh_char[0] gt -900.0) and $
      (hiiinfo.hii_kk04_log12oh_nuclear[0] gt -900.0) and $
      (hiiinfo.hii_kk04_log12oh_avg[0] gt -900.0),kk04_nindx)
    xkk04 = hiiinfo[kk04_indx].hii_kk04_log12oh_avg[0]
    ykk04 = hiiinfo[kk04_indx].hii_kk04_log12oh_char[0]
    xkk04_err = hiiinfo[kk04_indx].hii_kk04_log12oh_avg[1]
    ykk04_err = hiiinfo[kk04_indx].hii_kk04_log12oh_char[1]
    yresidkk04 = ykk04-xkk04
    yresidkk04_err = sqrt(ykk04_err^2+xkk04_err^2)

    splog, 'KK04: Average vs characteristic'
    stats = im_stats(yresidkk04)
    xstr_kk04 = strtrim(string(stats.median,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma,format='(F12.2)'),2)+')'
    splog, xstr_kk04

    stats = im_stats([yresidkk04,yresidpt05])
    xstr = '<\Delta>='+strtrim(string(stats.mean,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma,format='(F12.2)'),2)
;   xstr = '<\Delta>='+strtrim(string(stats.median,format='(F12.2)'),2)+' ('+$
;     strtrim(string(stats.mean,format='(F12.2)'),2)+$
;     '\pm'+strtrim(string(stats.sigma,format='(F12.2)'),2)+')'
    
    xtitle = '12 + log (O/H)_{HII avg}'
    ytitle = '12 + log (O/H)_{char}'
    
    djs_plot, [0], [0], /nodata, /noerase, xsty=1, ysty=1, xtitle='', ytitle='', $
      xrange=xrange, yrange=yrange, position=pos[*,1], $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    djs_oplot, ohaxis, ohaxis, line=0
    axis, /yaxis, ysty=1, yrange=yrange, ytitle=textoidl(ytitle)

    oploterror, xkk04, ykk04, xkk04_err, ykk04_err, $
      errthick=errthick1, color=fsc_color(kk04color1,100), errcolor=fsc_color(kk04color1,100), $
      psym=symcat(kk04sym1), symsize=kk04psize2
    oploterror, xpt05, ypt05, xpt05_err, ypt05_err, $
      errthick=errthick1, color=fsc_color(pt05color1,100), errcolor=fsc_color(pt05color1,100), $
      psym=symcat(pt05sym1), symsize=pt05psize2
;   legend, '(b)', /left, /top, box=0, charsize=charsize_5, charthick=postthick2, margin=0

    im_legend, ['KK04','PT05'], psym=[kk04sym1,pt05sym1], /right, /bottom, $
      box=0, fill=[1,1], spacing=1.8, symsize=1.3, color=[kk04color1,pt05color1], margin=0, charsize=1.2
    
    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle, ytitle='', xsty=1, ysty=1, $
      xrange=xrange, yrange=yresidrange, position=pos[*,3], ytickname=replicate(' ',10)
    djs_oplot, ohaxis, ohaxis*0.0, line=0
    axis, /yaxis, ysty=1, yrange=yresidrange, ytitle=textoidl(yresidtitle)
    legend, textoidl(xstr), /right, /bottom, box=0, margin=0, charsize=1.2

    oplot, xkk04, yresidkk04, color=fsc_color(kk04color1,100), $
      psym=symcat(kk04sym1), symsize=kk04psize2
    oplot, xpt05, yresidpt05, color=fsc_color(pt05color1,100), $
      psym=symcat(pt05sym1), symsize=pt05psize2
    
    im_plotconfig, /psclose

return
end
