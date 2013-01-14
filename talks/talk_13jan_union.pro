pro polyfill_legend, xc, yc, dx, dy, color=color, hatch=hatch, outcolor=outcolor
    xx = [xc-dx,xc+dx,xc+dx,xc-dx]
    yy = [yc-dy,yc-dy,yc+dy,yc+dy]
    if keyword_set(hatch) then begin
       polyfill, xx, yy, /data, color=im_color(color), noclip=0, $
         /line_fill, orientation=45, spacing=0.1
       polyfill, xx, yy, /data, color=im_color(color), noclip=0, $
         /line_fill, orientation=135, spacing=0.1
    endif else begin
       polyfill, xx, yy, /data, color=im_color(color), noclip=0, /fill
    endelse
    djs_oplot, [xc-dx,xc+dx], (yc-dy)*[1,1], color=im_color(outcolor), line=0, thick=4
    djs_oplot, [xc-dx,xc+dx], (yc+dy)*[1,1], color=im_color(outcolor), line=0, thick=4
    djs_oplot, (xc-dx)*[1,1], [yc-dy,yc+dy], color=im_color(outcolor), line=0, thick=4
    djs_oplot, (xc+dx)*[1,1], [yc-dy,yc+dy], color=im_color(outcolor), line=0, thick=4
return
end    
    
pro mfpolyfill_legend, xc, yc, dx, dy, color=color, hatch=hatch, outcolor=outcolor
    xx = [xc-dx,xc+dx,xc+dx,xc-dx]
    yy = [yc-dy,yc-dy,yc+dy,yc+dy]
    if keyword_set(hatch) then begin
       polyfill, xx, yy, /data, color=im_color(color), noclip=0, $
         /line_fill, orientation=45, spacing=0.1
       polyfill, xx, yy, /data, color=im_color(color), noclip=0, $
         /line_fill, orientation=135, spacing=0.1
    endif else begin
       polyfill, xx, yy, /data, color=im_color(color), noclip=0, /fill
    endelse
    if n_elements(outcolor) ne 0 then begin
       djs_oplot, [xc-dx,xc+dx], (yc-dy)*[1,1], color=im_color(outcolor), line=0, thick=4
       djs_oplot, [xc-dx,xc+dx], (yc+dy)*[1,1], color=im_color(outcolor), line=0, thick=4
       djs_oplot, (xc-dx)*[1,1], [yc-dy,yc+dy], color=im_color(outcolor), line=0, thick=4
       djs_oplot, (xc+dx)*[1,1], [yc-dy,yc+dy], color=im_color(outcolor), line=0, thick=4
    endif
return
end    
    
pro talk_13jan_union, keynote=keynote, noevol=noevol
; jm12jan28ucsd - miscellaneous plots for my 2012 Feb talk at Siena

    common com_siena, model, mstar, isedfit
    
    mfpath = mf_path()
    mfspath = mf_path(/mfs)
    isedpath = mf_path(/isedfit)
    mzpath = mz_path()
    
    talkpath = getenv('IM_RESEARCH_DIR')+'/talks/2013/13jan_union/'

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('black')

    evol_zbins = mf_zbins(nevolzbins,/evolved)
    
    
; --------------------------------------------------
; Figure 16: fractional quenching rate versus mass and redshift 
    result = mrdfits(mfpath+'mf_quenching_rate_supergrid01.fits.gz',1)

    scale = 100D
    yrange = scale*[-0.0,0.23]
    xrange = [-0.02,0.83]
    
;   mass = [10.25,10.75,11.25]
;   masslo = [10.0,10.5,11.0]
;   massup = [10.5,11.0,11.5]
;   masserr = [0.25,0.25,0.25]
    mass = [9.75,10.25,10.75,11.25]
    masslo = [9.5,10.0,10.5,11.0]
    massup = [10.0,10.5,11.0,11.5]
    masserr = [0.25,0.25,0.25,0.25]
    nmass = n_elements(mass)

;    masscolors = ['tan','steel blue','grey']
;    psym = [16,14,17]
;    line = [0,5,3]
;    symsize1 = [2.3,2.3,2.6]
    
   masscolors = ['tomato','tan','sky blue','grey80']
;  masscolorsout = ['tomato','tan','sky blue','grey80']
   masscolorsout = ['firebrick','brown','dodger blue','grey40']
   psym = [16,17,14,15]
   line = [4,5,3,0]
   lineout = [0,0,0,0]
;  line = [0,5,3,4]
   symsize1 = [2.3,2.6,2.5,2.3]

    zpoly = result.zpoly
    zmax = [result.zmax_95_10,result.zmax_10_105,result.zmax_105_11,result.zmax_11_115]
    allquench = scale*[[result.fquench_95_10],[result.fquench_10_105],$
      [result.fquench_105_11],[result.fquench_11_115]]
    allquenchup = scale*[[result.fquenchup_95_10],[result.fquenchup_10_105],$
      [result.fquenchup_105_11],[result.fquenchup_11_115]]
    allquenchlo = scale*[[result.fquenchlo_95_10],[result.fquenchlo_10_105],$
      [result.fquenchlo_105_11],[result.fquenchlo_11_115]]

;   zmax = [result.zmax_10_105,result.zmax_105_11,result.zmax_11_115]
;   allquench = scale*[[result.fquench_10_105],$
;     [result.fquench_105_11],[result.fquench_11_115]]
;   allquenchup = scale*[[result.fquenchup_10_105],$
;     [result.fquenchup_105_11],[result.fquenchup_11_115]]
;   allquenchlo = scale*[[result.fquenchlo_10_105],$
;     [result.fquenchlo_105_11],[result.fquenchlo_11_115]]

    psfile = talkpath+'fquench.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.4,0.4], $
      height=4.8, width=6.7, charsize=3, keynote=keynote

    plot, [0], [0], /nodata, position=pos, xsty=9, ysty=9, $
      yrange=yrange, xrange=xrange, color=keycolor
;   djs_oplot, !x.crange, [0,0], line=0, thick=6, color='grey'

    label = ['log (M/M'+sunsymbol()+') =  ',replicate('             ',4)]+$
;   label = ['\Delta'+'log(M) = ',replicate('             ',4)]+$
      strtrim(string(masslo,format='(F12.1)'),2)+'-'+strtrim(string(massup,format='(F12.1)'),2)
;   label = ['log (M/M'+sunsymbol()+') = ',replicate('             ',4)]+$
;     strtrim(string(mass,format='(F12.2)'),2)
    polyfill_legend, xrange[1]-0.125, yrange[1]-1.7, (xrange[1]-xrange[0])*0.095, (yrange[1]-yrange[0])*0.027, $
      color=masscolors[0], outcolor=masscolorsout[0]
    polyfill_legend, xrange[1]-0.125, yrange[1]-3.4, (xrange[1]-xrange[0])*0.095, (yrange[1]-yrange[0])*0.027, $
      color=masscolors[1], outcolor=masscolorsout[1], /hatch
    polyfill_legend, xrange[1]-0.125, yrange[1]-5.0, (xrange[1]-xrange[0])*0.095, (yrange[1]-yrange[0])*0.027, $
      color=masscolors[2], outcolor=masscolorsout[2]
    polyfill_legend, xrange[1]-0.125, yrange[1]-6.7, (xrange[1]-xrange[0])*0.095, (yrange[1]-yrange[0])*0.027, $
      color=masscolors[3], outcolor=masscolorsout[3], /hatch
    im_legend, label, /right, /top, box=0, $ ; color=masscolors, $
      margin=-0.2, charsize=1.6, thick=8, textcolor=keycolor, color=keycolor, $
      pspacing=1.6, spacing=2.5, symsize=symsize1*0.5;, psym=-psym; line=line, 

    off = [0.0,0.005,0.00,0.01]*0
    hatch = [0,1,0,1]
    lthick = 4*[1,1,1,1]
    for mm = nmass-1, 0, -1 do begin
       good = where(zpoly ge 0.001 and zpoly le zmax[mm],ngood)
       if mm eq 0 then good = where(zpoly le zmax[mm] and zpoly ge 0.011,ngood)
       if mm eq 1 then good = where(zpoly le zmax[mm] and zpoly ge 0.007,ngood)
       if mm eq 2 then good = where(zpoly le zmax[mm]-0.01 and zpoly ge 0.003,ngood)
       if mm eq 3 then good = where(zpoly le zmax[mm],ngood)

       if hatch[mm] then begin
          polyfill, [zpoly[good],reverse(zpoly[good])]+off[mm], [allquenchlo[good,mm],reverse(allquenchup[good,mm])], $
            /data, color=im_color(masscolors[mm]), noclip=0, /line_fill, orientation=45, spacing=0.1
          polyfill, [zpoly[good],reverse(zpoly[good])]+off[mm], [allquenchlo[good,mm],reverse(allquenchup[good,mm])], $
            /data, color=im_color(masscolors[mm]), noclip=0, /line_fill, orientation=135, spacing=0.1
       endif else begin
          polyfill, [zpoly[good],reverse(zpoly[good])]+off[mm], [allquenchlo[good,mm],reverse(allquenchup[good,mm])], $
            /data, color=im_color(masscolors[mm]), noclip=0, /fill ;/line_fill, orientation=45, spacing=spacing
       endelse
; add the boundaries
       djs_oplot, zpoly[good]+off[mm], allquenchlo[good,mm], color=im_color(masscolorsout[mm]), line=lineout[mm], thick=lthick[mm]
       djs_oplot, zpoly[good]+off[mm], allquenchup[good,mm], color=im_color(masscolorsout[mm]), line=lineout[mm], thick=lthick[mm]
       djs_oplot, zpoly[good[0]]*[1,1]+off[mm], [allquenchlo[good[0],mm],allquenchup[good[0],mm]], $
         color=im_color(masscolorsout[mm]), line=lineout[mm], thick=lthick[mm]
;      oplot, zpoly[good[ngood-1]]*[1,1]+off[mm], [allquenchlo[good[ngood-1],mm],allquenchup[good[ngood-1],mm]], $
;        color=im_color(masscolorsout[mm]), line=lineout[mm], thick=lthick[mm]
    endfor

    for mm = nmass-1, 0, -1 do begin
       good = where(zpoly ge 0.001 and zpoly le zmax[mm],ngood)
       if mm eq 0 then good = where(zpoly le zmax[mm] and zpoly ge 0.011,ngood)
       if mm eq 1 then good = where(zpoly le zmax[mm] and zpoly ge 0.007,ngood)
       if mm eq 2 then good = where(zpoly le zmax[mm]-0.01 and zpoly ge 0.003,ngood)
       if mm eq 3 then good = where(zpoly le zmax[mm],ngood)
;      oplot, zpoly[good], allquench[good,mm], line=line[mm], thick=8, color=keycolor
       gd = where(evol_zbins.zbin lt zmax[mm])
;      djs_oplot, evol_zbins[gd].zbin, interpol(allquench[good,mm],zpoly[good],evol_zbins[gd].zbin);, $
;        psym=symcat(psym[mm]), symsize=symsize1*0.6
    endfor
    
;; overlay some of the lines
;    for mm = nmass-1, 0, -1 do begin
;       if mm eq 2 then begin
;;      if mm eq 2 or mm eq 3 then begin
;          if mm eq 2 then good = where(zpoly le zmax[mm]-0.01 and zpoly ge 0.003,ngood)
;          if mm eq 3 then good = where(zpoly le zmax[mm],ngood)
;;         good = where(zpoly le zmax[mm],ngood)
;          djs_oplot, zpoly[good]+off[mm], allquenchlo[good,mm], color=im_color(masscolorsout[mm]), line=line[mm], thick=6
;          djs_oplot, zpoly[good]+off[mm], allquenchup[good,mm], color=im_color(masscolorsout[mm]), line=line[mm], thick=6
;          djs_oplot, zpoly[good[0]]*[1,1]+off[mm], [allquenchlo[good[0],mm],allquenchup[good[0],mm]], $
;            color=im_color(masscolorsout[mm]), line=line[mm], thick=6
;          djs_oplot, zpoly[good[ngood-1]]*[1,1]+off[mm], [allquenchlo[good[ngood-1],mm],allquenchup[good[ngood-1],mm]], $
;            color=im_color(masscolorsout[mm]), line=line[mm], thick=6
;       endif
;    endfor

; overlay the axes
    plot, [0], [0], /nodata, /noerase, position=pos, xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle='Redshift', $
      ytitle=textoidl('F_{quench} (% Gyr^{-1})'), color=keycolor

    im_plotconfig, /psclose, psfile=psfile, pdf=pdf, keynote=keynote

stop    
    
; --------------------------------------------------
; Figure 12: number-density evolution for all, quiescent, and
; star-forming galaxies in four bins of mass 
    all = [read_mf_vmax('sdss',/rho,/log),$
      read_mf_vmax(/rho,/avgfield,/log)]
    qq = [read_mf_vmax('sdss',/rho,/log,/quiescent),$
      read_mf_vmax(/rho,/avgfield,/log,/quiescent)]
    sf = [read_mf_vmax('sdss',/rho,/log,/active),$
      read_mf_vmax(/rho,/avgfield,/log,/active)]
;   qq = mrdfits(mfpath+'rhoden_numden_vs_redshift_quiescent.fits.gz',1)
;   sf = mrdfits(mfpath+'rhoden_numden_vs_redshift_active.fits.gz',1)

    result = mrdfits(mfpath+'fit_rhoden_numden_supergrid01.fits.gz',1)

    zbins = [mf_zbins(/sdss),mf_zbins()]
    zz = zbins.zbin & zsig = zbins.zsigma*0
    
;   drory = [[get_mf_literature(/drory,/quiescent)],$
;     [get_mf_literature(/drory,/active)]]

    nmonte = 100
    nsample = 3
;   tot = [[qq],[sf],[all]]
;   sample = ['qq','sf','all']
;   psym = [14,16,15]
;   color = ['dark red','steel blue','black']
;   symsize = [2.0,1.8,1.2]

    tot = [[all],[qq],[sf]]
    sample = ['all','qq','sf']
    psym = [15,14,16]
    color = ['white','dark red','steel blue']
    fillcolor = ['grey10','tomato','powder blue']
    symsize = [1.5,2.2,2.0]
    line = [0,3,5]

    ytickint = 1
    psym_open = [4,9,6]

    psfile = talkpath+'z_vs_numden_bymass.ps'
    im_plotconfig, 2, pos, psfile=psfile, xmargin=[1.3,1.5], $
      width=4*[1,1], height=2.9*[1,1], charsize=1.9, yspace=0.1, $
      xspace=0.2, keynote=keynote

    xrange = [-0.03,1.05]
    numrange = [-4.6,-2.3]
;   numrange = [-4.7,-2.0]

    zaxis = range(0,1,50)
    
; #########################
; 9.5-10
    plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      yrange=numrange, xrange=xrange, xtitle='', $
      ytitle='', ytickinterval=ytickint, xtickname=replicate(' ',10), color=keycolor
    im_legend, 'log (M/M'+sunsymbol()+')=9.5-10', /right, /bottom, box=0, margin=0, charsize=1.7, textcolor=keycolor

; overlay the error bands
    for ii = 1, nsample-1 do begin
;   for ii = 0, nsample-1 do begin
       zaxis = range(0.0,result[ii].zmax_lim_95_10,50)
       coeff = result[ii].numdencoeff_95_10
       covar = result[ii].numdencovar_95_10
       rand = mrandomn(seed,covar,nmonte)
       rand[*,0] = rand[*,0] + coeff[0]
       rand[*,1] = rand[*,1] + coeff[1]
       ell = covar2ellipse(covar,nsigma=1.0) ; get models within 1-sigma
       indx = get_ellipse_indices(rand[*,0],rand[*,1],$
         major=ell.major,minor=ell.minor,angle=ell.angle, $
         xcenter=coeff[0],ycenter=coeff[1])
       nindx = n_elements(indx)
       shade = fltarr(n_elements(zaxis),nindx)
       for bb = 0, nindx-1 do shade[*,bb] = poly(alog10(1+zaxis),rand[indx[bb],*])
       polyfill, [zaxis,reverse(zaxis)], [min(shade,dim=2),reverse(max(shade,dim=2))], $
         noclip=0, /fill, color=im_color(fillcolor[ii])
;      for bb = 0, nindx-1 do djs_oplot, zaxis, poly(alog10(1+zaxis),rand[indx[bb],*]), color=im_color(fillcolor[ii])
    end
       
; now plot the data    
    for ii = 1, nsample-1 do begin
;   for ii = 0, nsample-1 do begin
       gd = where(zbins.zup le result[ii].zmax_95_10,ngd)
       lim = where(zbins.zup gt result[ii].zmax_95_10 and zbins.zup le result[ii].zmax_lim_95_10,nlim)

       wmean = im_weighted_mean(tot[gd,ii].num_95_10,errors=tot[gd,ii].numerr_stat_95_10,wsigma=wsig)
       splog, 'M = 9.5-10, '+sample[ii]+', '+strtrim(result[ii].numdencoeff_95_10[1],2)+'+/-'+$
         strtrim(result[ii].numdencoeff_err_95_10[1],2)+', '+strtrim(result[ii].numdencoeff_95_10[0],2)+'+/-'+$
         strtrim(result[ii].numdencoeff_err_95_10[0],2)+', '+strtrim(wmean,2)+'+/-'+strtrim(wsig,2)

       z0 = 0.4 & fact = -result[ii].numdencoeff_95_10[1]*alog10(1+z0)
       facterr = result[ii].numdencoeff_err_95_10[1]*alog10(1+z0)
       if 10^fact gt 2.0 then splog, 'M = 9.5-10, '+sample[ii]+' to z='+strtrim(z0,2)+$
         ', factor '+strtrim(10^fact,2)+'+/-'+strtrim(alog(10)*facterr*10^fact,2) else $
           splog, 'M = 9.5-10, '+sample[ii]+' to z='+strtrim(z0,2)+$
         ', percent '+strtrim(100*(10^(fact)-1),2)+'+/-'+strtrim(100*alog(10)*facterr,2)
       if ii ne nsample-1 then print

       zaxis = range(0.0,result[ii].zmax_lim_95_10,50)
;      oplot, zaxis, poly(alog10(1+zaxis),result[ii].numdencoeff_95_10), $
;        line=line[ii], color=im_color(color[ii])

       oploterror, zz[gd], tot[gd,ii].num_95_10, zsig[gd], tot[gd,ii].numerr_stat_95_10, $
         psym=symcat(psym[ii]), symsize=symsize[ii], color=im_color(color[ii]), errcolor=im_color(color[ii]), $
         errstyle=0, errthick=8, nohat=nohat, line=line[ii], thick=6

       if nlim ne 0 then begin
          oploterror, zz[lim], tot[lim,ii].num_95_10, zsig[lim], tot[lim,ii].numerr_stat_95_10*0, $
            psym=symcat(psym[ii]), symsize=symsize[ii], color=im_color(color[ii]), errcolor=im_color(color[ii]), $
            errstyle=0, errthick=8, nohat=nohat, line=line[ii], thick=6
          im_symbols, 112, psize=symsize*2, color=im_color(color[ii]), thick=8
          oplot, [zz[lim]], [tot[lim,ii].num_95_10], psym=8, color=keycolor
;         djs_oplot, [zz[gd[ngd-1]],zz[lim[0]]], [tot[gd[ngd-1],ii].num_95_10,tot[lim[0],ii].num_95_10], $
;           line=line[ii], color=im_color(color[ii]), thick=6
       endif
    endfor 
    print, '###########################################################################'

; #########################
; 10-10.5
    plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      yrange=numrange, xrange=xrange, xtitle='', ytickname=replicate(' ',10), $
      ytitle='', ytickinterval=ytickint, xtickname=replicate(' ',10), color=keycolor
;   axis, yrange=numrange, ysty=1, /yaxis, color=keycolor
    im_legend, 'log (M/M'+sunsymbol()+')=10-10.5', /right, /bottom, box=0, margin=0, charsize=1.7, textcolor=keycolor

; overlay the error bands
    for ii = 1, nsample-1 do begin
;   for ii = 0, nsample-1 do begin
       zaxis = range(0.0,result[ii].zmax_lim_10_105,50)
       coeff = result[ii].numdencoeff_10_105
       covar = result[ii].numdencovar_10_105
       rand = mrandomn(seed,covar,nmonte)
       rand[*,0] = rand[*,0] + coeff[0]
       rand[*,1] = rand[*,1] + coeff[1]
       ell = covar2ellipse(covar,nsigma=1.0) ; get models within 1-sigma
       indx = get_ellipse_indices(rand[*,0],rand[*,1],$
         major=ell.major,minor=ell.minor,angle=ell.angle, $
         xcenter=coeff[0],ycenter=coeff[1])
       nindx = n_elements(indx)
       shade = fltarr(n_elements(zaxis),nindx)
       for bb = 0, nindx-1 do shade[*,bb] = poly(alog10(1+zaxis),rand[indx[bb],*])
       polyfill, [zaxis,reverse(zaxis)], [min(shade,dim=2),reverse(max(shade,dim=2))], $
         noclip=0, /fill, color=im_color(fillcolor[ii])
;      for bb = 0, nindx-1 do djs_oplot, zaxis, poly(alog10(1+zaxis),rand[indx[bb],*]), color=im_color(fillcolor[ii])
    end
       
; now plot the data    
    for ii = 1, nsample-1 do begin
;   for ii = 0, nsample-1 do begin
       gd = where(zbins.zup le result[ii].zmax_10_105,ngd)
       lim = where(zbins.zup gt result[ii].zmax_10_105 and zbins.zup le result[ii].zmax_lim_10_105,nlim)

       wmean = im_weighted_mean(tot[gd,ii].num_10_105,errors=tot[gd,ii].numerr_stat_10_105,wsigma=wsig)
       splog, 'M = 10-10.5, '+sample[ii]+', '+strtrim(result[ii].numdencoeff_10_105[1],2)+'+/-'+$
         strtrim(result[ii].numdencoeff_err_10_105[1],2)+', '+strtrim(result[ii].numdencoeff_10_105[0],2)+'+/-'+$
         strtrim(result[ii].numdencoeff_err_10_105[0],2)+', '+strtrim(wmean,2)+'+/-'+strtrim(wsig,2)

       z0 = 0.6 & fact = -result[ii].numdencoeff_10_105[1]*alog10(1+z0)
       facterr = result[ii].numdencoeff_err_10_105[1]*alog10(1+z0)
       if 10^fact gt 2.0 then splog, 'M = 10-10.5, '+sample[ii]+' to z='+strtrim(z0,2)+$
         ', factor '+strtrim(10^fact,2)+'+/-'+strtrim(alog(10)*facterr*10^fact,2) else $
           splog, 'M = 10-10.5, '+sample[ii]+' to z='+strtrim(z0,2)+$
         ', percent '+strtrim(100*(10^(fact)-1),2)+'+/-'+strtrim(100*alog(10)*facterr,2)
       if ii ne nsample-1 then print

       zaxis = range(0.0,result[ii].zmax_lim_10_105,50)
;      djs_oplot, zaxis, poly(alog10(1+zaxis),result[ii].numdencoeff_10_105), $
;        line=line[ii], color=im_color(color[ii])

       oploterror, zz[gd], tot[gd,ii].num_10_105, zsig[gd], tot[gd,ii].numerr_stat_10_105, $
         psym=symcat(psym[ii]), symsize=symsize[ii], color=im_color(color[ii]), errcolor=im_color(color[ii]), $
         errstyle=0, errthick=8, nohat=nohat, line=line[ii], thick=6

       if nlim ne 0 then begin
          oploterror, zz[lim], tot[lim,ii].num_10_105, zsig[lim], tot[lim,ii].numerr_stat_10_105*0, $
            psym=symcat(psym[ii]), symsize=symsize[ii], color=im_color(color[ii]), errcolor=im_color(color[ii]), $
            errstyle=0, errthick=8, nohat=nohat, line=line[ii], thick=6
          im_symbols, 112, psize=symsize*2, color=im_color(color[ii]), thick=8
          oplot, [zz[lim]], [tot[lim,ii].num_10_105], psym=8, color=keycolor
;         djs_oplot, [zz[gd[ngd-1]],zz[lim[0]]], [tot[gd[ngd-1],ii].num_10_105,tot[lim[0],ii].num_10_105], $
;           line=line[ii], color=im_color(color[ii]), thick=6
       endif
    endfor 
    print, '###########################################################################'

; #########################
; 10.5-11
    plot, [0], [0], /nodata, /noerase, position=pos[*,2], xsty=1, ysty=1, $
      yrange=numrange, xrange=xrange, xtitle='Redshift', $
      ytitle='', ytickinterval=ytickint, color=keycolor, xtickinterval=0.5
    im_legend, 'log (M/M'+sunsymbol()+')=10.5-11', /right, /bottom, box=0, margin=0, charsize=1.7, textcolor=keycolor

; overlay the error bands
    for ii = 1, nsample-1 do begin
;   for ii = 0, nsample-1 do begin
       zaxis = range(0.0,result[ii].zmax_lim_105_11,50)
       coeff = result[ii].numdencoeff_105_11
       covar = result[ii].numdencovar_105_11
       rand = mrandomn(seed,covar,nmonte)
       rand[*,0] = rand[*,0] + coeff[0]
       rand[*,1] = rand[*,1] + coeff[1]
       ell = covar2ellipse(covar,nsigma=1.0) ; get models within 1-sigma
       indx = get_ellipse_indices(rand[*,0],rand[*,1],$
         major=ell.major,minor=ell.minor,angle=ell.angle, $
         xcenter=coeff[0],ycenter=coeff[1])
       nindx = n_elements(indx)
       shade = fltarr(n_elements(zaxis),nindx)
       for bb = 0, nindx-1 do shade[*,bb] = poly(alog10(1+zaxis),rand[indx[bb],*])
       polyfill, [zaxis,reverse(zaxis)], [min(shade,dim=2),reverse(max(shade,dim=2))], $
         noclip=0, /fill, color=im_color(fillcolor[ii])
;      for bb = 0, nindx-1 do djs_oplot, zaxis, poly(alog10(1+zaxis),rand[indx[bb],*]), color=im_color(fillcolor[ii])
    end
       
; now plot the data    
    for ii = 1, nsample-1 do begin
;   for ii = 0, nsample-1 do begin
       gd = where(zbins.zup le result[ii].zmax_105_11,ngd)
       lim = where(zbins.zup gt result[ii].zmax_105_11 and zbins.zup le result[ii].zmax_lim_105_11,nlim)

       wmean = im_weighted_mean(tot[gd,ii].num_105_11,errors=tot[gd,ii].numerr_stat_105_11,wsigma=wsig)
       splog, 'M = 10.5-11, '+sample[ii]+', '+strtrim(result[ii].numdencoeff_105_11[1],2)+'+/-'+$
         strtrim(result[ii].numdencoeff_err_105_11[1],2)+', '+strtrim(result[ii].numdencoeff_105_11[0],2)+'+/-'+$
         strtrim(result[ii].numdencoeff_err_105_11[0],2)+', '+strtrim(wmean,2)+'+/-'+strtrim(wsig,2)

       z0 = 0.8 & fact = -result[ii].numdencoeff_105_11[1]*alog10(1+z0)
       facterr = result[ii].numdencoeff_err_105_11[1]*alog10(1+z0)
       if 10^fact gt 2.0 then splog, 'M = 10.5-11, '+sample[ii]+' to z='+strtrim(z0,2)+$
         ', factor '+strtrim(10^fact,2)+'+/-'+strtrim(alog(10)*facterr*10^fact,2) else $
           splog, 'M = 10.5-11, '+sample[ii]+' to z='+strtrim(z0,2)+$
         ', percent '+strtrim(100*(10^(fact)-1),2)+'+/-'+strtrim(100*alog(10)*facterr,2)
       if ii ne nsample-1 then print

       zaxis = range(0.0,result[ii].zmax_lim_105_11,50)
;      djs_oplot, zaxis, poly(alog10(1+zaxis),result[ii].numdencoeff_105_11), $
;        line=line[ii], color=im_color(color[ii])

       oploterror, zz[gd], tot[gd,ii].num_105_11, zsig[gd], tot[gd,ii].numerr_stat_105_11, $
         psym=symcat(psym[ii]), symsize=symsize[ii], color=im_color(color[ii]), errcolor=im_color(color[ii]), $
         errstyle=0, errthick=8, nohat=nohat, line=line[ii], thick=6

       if nlim ne 0 then begin
          oploterror, zz[lim], tot[lim,ii].num_105_11, zsig[lim], tot[lim,ii].numerr_stat_105_11*0, $
            psym=-symcat(psym[ii]), symsize=symsize[ii], color=im_color(color[ii]), errcolor=im_color(color[ii]), $
            errstyle=0, errthick=8, nohat=nohat, line=line[ii], thick=6
          im_symbols, 112, psize=symsize*2, color=im_color(color[ii]), thick=8
          oplot, [zz[lim]], [tot[lim,ii].num_105_11], psym=8, color=keycolor
;         djs_oplot, [zz[gd[ngd-1]],zz[lim[0]]], [tot[gd[ngd-1],ii].num_105_11,tot[lim[0],ii].num_105_11], $
;           line=line[ii], color=im_color(color[ii]), thick=6
       endif
    endfor 
    print, '###########################################################################'
    
; #########################
; 11-11.5
    plot, [0], [0], /nodata, /noerase, position=pos[*,3], xsty=1, ysty=1, $
      yrange=numrange, xrange=xrange, xtitle='Redshift', $
      ytickname=replicate(' ',10), ytitle='', ytickinterval=ytickint, color=keycolor, xtickinterval=0.5
;   axis, yrange=numrange, ysty=1, /yaxis, color=keycolor
    im_legend, 'log (M/M'+sunsymbol()+')=11-11.5', /right, /bottom, box=0, margin=0, charsize=1.7, textcolor=keycolor

;   mfpolyfill_legend, 0.9, -2.46, 0.11, 0.07, color=fillcolor[0]
;   mfpolyfill_legend, 0.9, -2.69, 0.11, 0.07, color=fillcolor[1]
;   mfpolyfill_legend, 0.9, -2.91, 0.11, 0.07, color=fillcolor[2]

;    im_legend, ['All','Quiescent','Star-Forming'], /right, /top, box=0, $
;      psym=-psym, symsize=symsize*0.9, charsize=1.4, $
;      color=color, spacing=2.05, margin=-0.5, pspacing=1.6, thick=8, line=line, textcolor=keycolor
    
;    mfpolyfill_legend, 0.9, -2.46, 0.11, 0.07, color=fillcolor[1]
;    mfpolyfill_legend, 0.9, -2.69, 0.11, 0.07, color=fillcolor[2]
;    im_legend, ['Quiescent','Star-Forming'], /right, /top, box=0, $
;      psym=-psym[1:2], symsize=symsize[1:2]*0.9, charsize=1.4, $
;      color=color[1:2], spacing=2.05, margin=-0.5, pspacing=1.6, thick=8, line=line[1:2], textcolor=keycolor

;   im_legend, ['All','Quiescent','Star-Forming'], /left, /top, box=0, $
;     psym=-psym, symsize=symsize*0.9, charsize=1.4, $
;     color=color, spacing=2, margin=-0.5, line=line, pspacing=1.6, thick=8

; overlay the error bands
    for ii = 1, nsample-1 do begin
;   for ii = 0, nsample-1 do begin
       zaxis = range(0.0,result[ii].zmax_lim_11_115,50)
       coeff = result[ii].numdencoeff_11_115
       covar = result[ii].numdencovar_11_115
       rand = mrandomn(seed,covar,nmonte)
       rand[*,0] = rand[*,0] + coeff[0]
       rand[*,1] = rand[*,1] + coeff[1]
       ell = covar2ellipse(covar,nsigma=1.0) ; get models within 1-sigma
       indx = get_ellipse_indices(rand[*,0],rand[*,1],$
         major=ell.major,minor=ell.minor,angle=ell.angle, $
         xcenter=coeff[0],ycenter=coeff[1])
       nindx = n_elements(indx)
       shade = fltarr(n_elements(zaxis),nindx)
       for bb = 0, nindx-1 do shade[*,bb] = poly(alog10(1+zaxis),rand[indx[bb],*])
       polyfill, [zaxis,reverse(zaxis)], [min(shade,dim=2),reverse(max(shade,dim=2))], $
         noclip=0, /fill, color=im_color(fillcolor[ii])
;      for bb = 0, nindx-1 do djs_oplot, zaxis, poly(alog10(1+zaxis),rand[indx[bb],*]), color=im_color(fillcolor[ii])
    end
       
    for ii = 1, nsample-1 do begin
;   for ii = 0, nsample-1 do begin
       gd = where(zbins.zup le result[ii].zmax_11_115,ngd)
       lim = where(zbins.zup gt result[ii].zmax_11_115 and zbins.zup le result[ii].zmax_lim_11_115,nlim)

       wmean = im_weighted_mean(tot[gd,ii].num_11_115,errors=tot[gd,ii].numerr_stat_11_115,wsigma=wsig)
       splog, 'M = 11-11.5, '+sample[ii]+', '+strtrim(result[ii].numdencoeff_11_115[1],2)+'+/-'+$
         strtrim(result[ii].numdencoeff_err_11_115[1],2)+', '+strtrim(result[ii].numdencoeff_11_115[0],2)+'+/-'+$
         strtrim(result[ii].numdencoeff_err_11_115[0],2)+', '+strtrim(wmean,2)+'+/-'+strtrim(wsig,2)

       z0 = 1.0 & fact = -result[ii].numdencoeff_11_115[1]*alog10(1+z0)
       facterr = result[ii].numdencoeff_err_11_115[1]*alog10(1+z0)
       if 10^fact gt 2.0 then splog, 'M = 11-11.5, '+sample[ii]+' to z='+strtrim(z0,2)+$
         ', factor '+strtrim(10^fact,2)+'+/-'+strtrim(alog(10)*facterr*10^fact,2) else $
           splog, 'M = 11-11.5, '+sample[ii]+' to z='+strtrim(z0,2)+$
         ', percent '+strtrim(100*(10^(fact)-1),2)+'+/-'+strtrim(100*alog(10)*facterr,2)
       if ii ne nsample-1 then print
       
       zaxis = range(0.0,result[ii].zmax_lim_11_115,50)
;      djs_oplot, zaxis, poly(alog10(1+zaxis),result[ii].numdencoeff_11_115), $
;        line=line[ii], color=im_color(color[ii])

       oploterror, zz[gd], tot[gd,ii].num_11_115, zsig[gd], tot[gd,ii].numerr_stat_11_115, $
         psym=symcat(psym[ii]), symsize=symsize[ii], color=im_color(color[ii]), errcolor=im_color(color[ii]), $
         errstyle=0, errthick=8, nohat=nohat, line=line[ii], thick=6

       if nlim ne 0 then begin
          oploterror, zz[lim], tot[lim,ii].num_11_115, zsig[lim], tot[lim,ii].numerr_stat_11_115*0, $
            psym=symcat(psym[ii]), symsize=symsize[ii], color=im_color(color[ii]), errcolor=im_color(color[ii]), $
            errstyle=0, errthick=8, nohat=nohat, line=line[ii], thick=6
          im_symbols, 112, psize=symsize*2, color=im_color(color[ii]), thick=8
          oplot, [zz[lim]], [tot[lim,ii].num_11_115], psym=8, color=keycolor
;         djs_oplot, [zz[gd[ngd-1]],zz[lim[0]]], [tot[gd[ngd-1],ii].num_11_115,tot[lim[0],ii].num_105_11], $
 ;          line=line[ii], color=im_color(color[ii]), thick=6
       endif
    endfor 

    xyouts, pos[0,0]-0.08, djs_mean([pos[1,0],pos[3,2]]), mfplot_numdentitle(), $
      align=0.5, /normal, orientation=90, color=keycolor
;   xyouts, pos[2,1]+0.1, djs_mean([pos[1,0],pos[3,2]]), mfplot_numdentitle(), $
;     align=0.5, /normal, orientation=270, color=keycolor

    splog, 'Factor by which QQ galaxies outnumber SF galaxies at z=1:'
    splog, poly(alog10(2),result[1].numdencoeff_11_115)-poly(alog10(2),result[2].numdencoeff_11_115) 
    
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf, keynote=keynote

stop    
    
; ---------------------------------------------------------------------------
; evolution of the qq MF in one panel
    zbins = mf_zbins(nzbins)
    xrange = [8.7,12.2]
    yrange = [-5,-2]

    ss = read_mf_vmax('sdss',noevol=noevol,/log,/quiescent)
    pp = read_mf_vmax(noevol=noevol,/avgfield,/log,/quiescent)

    if keyword_set(keynote) then col1 = 'white' else col1 = 'black'
    col = [col1,'coral','steel blue','forest green','orchid','orange','tan']
    for ii = -1, nzbins-1 do begin
       psfile = talkpath+'qq_evol_zbin'+strtrim(ii+1,2)+'.ps'
       im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, height=5.5

       plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
         ytitle=mfplot_phititle(), $
         xtickinterval=1, xminor=4, color=keycolor, ytickinterval=1

       good = where(ss.limit eq 1 and ss.number ge 3)
       phimass = ss.mass[good]
       polyfill, [phimass,reverse(phimass)],[ss.phi_lower_stat[good]-0.03,$
         reverse(ss.phi_upper_stat[good]+0.03)], $
         /data, color=im_color(col[0],10), noclip=0, /fill

       if ii ge 0 then for iz = 0, ii do begin
          good = where(pp[iz].limit eq 1 and pp[iz].number ge 3)
          phimass = pp[iz].mass[good]
          polyfill, [phimass,reverse(phimass)],[pp[iz].phi_lower_stat[good]-0.03,$
            reverse(pp[iz].phi_upper_stat[good]+0.03)], $
            /data, color=im_color(col[iz+1]), noclip=0, /fill
       endfor
       im_plotconfig, /psclose, psfile=psfile, keynote=keynote
    endfor


; ---------------------------------------------------------------------------
; sdss - compare quiescent and active samples
    psfile = talkpath+'mf_sdss.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.5,0.2], $
      width=6.8, height=6.0, charsize=3.0, keynote=keynote

    xrange = [8.8,12.2]
    yrange = [-6,-2]
    maxis1 = range(xrange[0]+0.15,12.5,100)

    subsample = ['active','quiescent']
    mfcolor = ['powder blue','tomato']

    psymsize = [0.9,1.1,0.9]*1.2
    mfpsym = [16,14,15]

    plot, [0], [0], /nodata, position=pos, xsty=9, ysty=9, $
      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
      ytitle=mfplot_phititle(), xtickinterval=1, color=keycolor
;   im_legend, ['All','Quiescent','Star-Forming'], /left, /bottom, box=0, $
;     color=reverse(mfcolor), psym=reverse(mfpsym), symsize=[1.7,2.1,1.9], $
;     symthick=8, spacing=2.4, charsize=1.8
    
    for jj = 0, n_elements(subsample)-1 do begin
       case jj of
          0: begin
             mfdata = read_mf_vmax('sdss',/active,/log)
          end
          1: begin
             mfdata = read_mf_vmax('sdss',/quiescent,/log)
          end
       endcase
       
       good = where(mfdata.limit eq 1)
       few = where(mfdata.limit eq 1 and mfdata.number le 3,nfew)

       phimass = mfdata.mass[good]
       polyfill, [phimass,reverse(phimass)],[mfdata.phi_lower_stat[good],$
         reverse(mfdata.phi_upper_stat[good])], $
         /data, color=im_color(mfcolor[jj]), noclip=0, /fill
       
;      oploterror, mfdata.mass[good], mfdata.phi[good], mfdata.phierr_upper[good], $
;        psym=symcat(mfpsym[jj],thick=5), symsize=psymsize[jj], errthick=5, $
;        color=im_color(mfcolor[jj],100+jj), errcolor=im_color(mfcolor[jj],100+jj), /hibar, /nohat
;      oploterror, mfdata.mass[good], mfdata.phi[good], mfdata.phierr_lower[good], psym=3, $
;        symsize=psymsize[jj], errthick=5, color=im_color(mfcolor[jj],100+jj), $
;        errcolor=im_color(mfcolor[jj],100+jj), /lobar, /nohat
    endfor

    im_plotconfig, /psclose, psfile=psfile, keynote=keynote


; ---------------------------------------------------------------------------
    

    
return
end
