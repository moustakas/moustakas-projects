pro oplot_salpeter_box, rhomin, rhomax, color=color, $
  outcolor=outcolor, _extra=extra
    ds = 0.03
    slope = 1.35
    polyfill, [slope-ds,slope-ds,slope+ds,slope+ds,slope-ds], $
      [rhomin,rhomax,rhomax,rhomin,rhomin], $
      /data, /fill, color=fsc_color(color,100), noclip=0
    djs_oplot, [slope-ds,slope-ds], [rhomin,rhomax], line=0, color=outcolor
    djs_oplot, [slope-ds,slope+ds], [rhomax,rhomax], line=0, color=outcolor
    djs_oplot, [slope+ds,slope+ds], [rhomax,rhomin], line=0, color=outcolor
    djs_oplot, [slope+ds,slope-ds], [rhomin,rhomin], line=0, color=outcolor
return
end
    
pro shade_rhostar, gamma, rhostarmin, rhostarmax, color=color, fill=fill
    if keyword_set(fill) then begin ; solid fill
       polyfill, [gamma,reverse(gamma)], $
         [rhostarmin,reverse(rhostarmax)], $
         /data, /fill, color=fsc_color(color,100), noclip=0
    endif else begin ; line-fill
       polyfill, [gamma,reverse(gamma)], $
         [rhostarmin,reverse(rhostarmax)], $
         /data, /line_fill, linestyle=0, orientation=45, spacing=0.1, $
         color=fsc_color(color,101), noclip=0, thick=1
       polyfill, [gamma,reverse(gamma)], $
         [rhostarmin,reverse(rhostarmax)], $
         /data, /line_fill, linestyle=0, orientation=135, spacing=0.1, $
         color=fsc_color(color,101), noclip=0, thick=1
    endelse
return
end

function get_rhostar, data, rhostar_err=rhostar_err, $
  local=local, sfh1=sfh1, sfh5=sfh5
    if keyword_set(sfh1) or keyword_set(sfh5) then begin
       if keyword_set(sfh1) then begin
          rhostar = data.sfh_rhostar1
          rhostar_err = data.sfh_rhostar1_err
       endif else begin
          rhostar = data.sfh_rhostar5
          rhostar_err = data.sfh_rhostar5_err
       endelse
    endif else begin
       if keyword_set(local) then begin
          rhostar = data.rhostar_local
          rhostar_err = data.rhostar_local_err
       endif else begin
          rhostar = data.rhostar
          rhostar_err = data.rhostar_err
       endelse
    endelse
; take the log        
    rhostar_err = rhostar_err/rhostar/alog(10)
    rhostar = alog10(rhostar)
return, rhostar
end

function get_rho_lit, data, rho_err=rho_err
    rho = data.rho
    rho_err = data.rho_err
; take the log        
    rho_err = rho_err/rho/alog(10)
    rho = alog10(rho)
return, rho
end

pro plot_cosmicimf_final, ps=ps, keynote=keynote
; jm10mar26ucsd - build the final plots

    cosmicimfpath = ages_path(/projects)+'cosmicimf/'
    paperpath = ages_path(/papers)+'cosmicimf/'
    if keyword_set(keynote) then paperpath = $
      getenv('RESEARCHPATH')+'/meetings/10apr_florida/keynote/'

    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

    final = mrdfits(cosmicimfpath+'final.fits.gz',1)
    models1 = mrdfits(cosmicimfpath+'final_models_zform1.2.fits.gz',1)
    models5 = mrdfits(cosmicimfpath+'final_models_zform5.0.fits.gz',1)
    dmsalp = mrdfits(cosmicimfpath+'dmsalp.fits.gz',1)
    zbins = cosmicimf_zbins(nzbins)
    peg = read_cosmicimf_sample(/pegase)
    nimf = n_elements(peg)

    zaxis = im_array(0,5.0,0.02)

; ---------------------------------------------------------------------------
; multipanel plot showing rho* vs redshift    
    psfile = paperpath+'final_hiz'+suffix
    im_plotconfig, 14, pos, psfile=psfile, keynote=keynote;, $
;     width=[2.4,2.4,2.4], xmargin=[1.1,0.2], height=[2.0,2.0]
    if keyword_set(keynote) then keycolor = djs_icolor('white')
    
;   xrange = [-0.08,4]
;   yrange = [7.0,9.0]
    xrange = [-0.08,1.4]
    yrange = [7.7,9.0]

; literature data
    bell = mf_bell()
    cole = mf_cole()
    borch = mf_borch()
    perez = mf_perez()
    bundy = mf_bundy()
    pozzetti07 = mf_07pozzetti()

    wilkins = rsex(getenv('PAPERSPATH')+'/literature/data/08wilkins.sex')
    wz = total([[wilkins.zmin],[wilkins.zmax]],2)/2.0
    wz_err = (wilkins.zmax-wilkins.zmin)/2.0
    wrho = 1.36D11*wilkins.omega
    wrho_err = 1.36D11*wilkins.omega_err
    
; Salpeter + a select set of slopes; we want to include Gamma=1.15
; (the best fit), so we have to interpolate
    gamma = [0.95,1.05,1.15,1.25,1.35]
    indx = [0.0,findex(peg[1:nimf-1].slope-1,gamma)+1] ; fractional index
;   print, interpolate(peg.slope-1,indx)
;   get_element, peg.slope-1, [1.0,1.1,1.2,1.3,1.4], indx
;   indx = [0,indx]
    nindx = n_elements(indx)

    for ii = 0, nindx-1 do begin
       if (ii eq 0) or (ii eq 3) then begin
          delvarx, ytickname 
       endif else begin
          ytitle = ''
          ytickname = replicate(' ',10)
       endelse
       if (ii ge 3) then begin
          xtitle = 'Redshift'
          delvarx, xtickname 
       endif else begin
          xtitle = ''
          xtickname = replicate(' ',10)
       endelse

       djs_plot, [0], [0], /nodata, noerase=(ii gt 0), position=pos[*,ii], $
         xsty=5, ysty=5, xrange=xrange, yrange=yrange

; show the models       
;      scale1 = 1.0+0.0*final[0].rhostar_local[indx[ii]]/models1[0].sfh_rhostar[0,indx[ii]]
;      for jj = 0, n_elements(models1)-1 do djs_oplot, models1[jj].zaxis, $
;        alog10(scale1*models1[jj].sfh_rhostar[*,indx[ii]]), line=0, color='green'

;      scale5 = final[0].rhostar_local[indx[ii]]/models5[0].sfh_rhostar[0,indx[ii]]
;      for jj = 0, n_elements(models5)-1 do djs_oplot, models5[jj].zaxis, $
;        alog10(scale5*models5[jj].sfh_rhostar[*,indx[ii]]), line=5, $
;        color='dark green', thick=5
       rhomin = alog10(interpolate(models5[1].sfh_rhostar,indx[ii]))
       rhomax = alog10(interpolate(models5[2].sfh_rhostar,indx[ii]))
       if keyword_set(keynote) then col = 'powder blue' else col = 'powder blue'
       polyfill, [models5[0].zaxis,reverse(models5[0].zaxis)],$
         [rhomin,reverse(rhomax)], /data, /fill, noclip=0, $
         color=fsc_color(col,101)
       
; Salpeter conversion plus data from wilkins       
       conv = 10^(-interpolate(dmsalp.dmsalp_sfhgrid03,indx[ii])) ; salpeter --> IMF[ii]
       conv_err = alog(10)*interpolate(dmsalp.dmsalp_sfhgrid03_err,indx[ii])*conv
       if keyword_set(keynote) then col = 'lime green' else col = 'lime green'
       oploterror, wz, alog10(wrho*conv), wz_err, wrho_err/wrho/alog(10), $
         psym=symcat(9), color=fsc_color(col,101), $
         errcolor=fsc_color(col,101), /nohat, symsize=0.5, $
         thick=4

; this is integrating the SFR density from z=0 to z=zform; save for later       
;      djs_oplot, models5[0].zaxis, alog10(interpolate(models5[0].sfh_rhostar_z0,indx[ii])), $
;        line=0, thick=5

;      if (ii eq 0) then begin ; from Wilkins+08 converted to Salpeter
;         djs_oplot, zaxis, alog10(1.36D11*0.0023*exp(-0.68*zaxis^1.2)/0.58), $
;           line=0, thick=8, color=keycolor
;      endif
       
; measured mass density
       splog, 'HACK!!!!!'
       if keyword_set(keynote) then col = 'orange red' else col = 'firebrick'
       for jj = 0, nzbins-1 do begin
;         rhostar = alog10(final[jj].rhostar[indx[ii]])+0.14  ; hack!
          rhostar = alog10(conv*final[jj].rhostar[0])+0.14
          rhostar_err = interpolate(final[jj].rhostar_err,indx[ii])/$
            interpolate(final[jj].rhostar,indx[ii])/alog(10)
          
          oploterror, zbins[jj].zbin, rhostar, zbins[jj].zbin_err, rhostar_err, $
            psym=-symcat(16), symsize=2.0, color=fsc_color(col,101), $
            errcolor=fsc_color(col,101)
;         rhostar = get_rhostar(final[jj],rhostar_err=rhostar_err)
;         oploterror, zbins[jj].zbin, rhostar[indx[ii]], zbins[jj].zbin_err, rhostar_err[indx[ii]], $
;           psym=-symcat(16), symsize=2.0, color=fsc_color('dodger blue',101), $
;           errcolor=fsc_color('dodger blue',101)
       endfor
       rhostar_local = get_rhostar(final[0],rhostar_err=rhostar_local_err,/local)
       rhostar_local = interpolate(rhostar_local,indx[ii])
       rhostar_local_err = interpolate(rhostar_local_err,indx[ii])
       if keyword_set(keynote) then col = 'blue' else col = 'navy'
       oploterror, [0.02], rhostar_local, rhostar_local_err, $
         psym=-symcat(6,thick=8), symsize=2.5, color=fsc_color(col,101), $
         errcolor=fsc_color(col,101)

; overplot the axes
       djs_plot, [0], [0], /nodata, /noerase, position=pos[*,ii], $
         xsty=1, ysty=1, xrange=xrange, yrange=yrange, xtickname=xtickname, $
         ytickname=ytickname, xtitle=xtitle, ytitle=ytitle, xtickinterval=0.5

       if (ii eq 0) then label = 'Salpeter' else $ ; (0.1-100 M_{\odot})' else $
         label = '\Gamma='+string(interpolate(peg.slope-1,indx[ii]),format='(F4.2)')
       im_legend, label, /right, /top, box=0, margin=0, charsize=1.5, textcolor=keycolor
    endfor

; ytitle
    xyouts, pos[0,0]-0.06, pos[1,0], textoidl('log (\rho_{*}/h_{70} M_{\odot} Mpc^{-3})'), $
      align=0.5, orientation=90, /normal
    
    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf -dEPSCrop '+psfile+' '+repstr(psfile,suffix,'.pdf'), /sh
       rmfile, psfile
    endif

stop    
    
; ---------------------------------------------------------------------------
    psfile = paperpath+'final_local'+suffix
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], $
      width=6.8, height=6.8, keynote=keynote
    if keyword_set(keynote) then keycolor = djs_icolor('white')

    gamma = peg.slope-1
    indx = lindgen(nimf-1)+1 ; exclude Salpeter
    gammaxis = im_array(min(gamma),max(gamma),0.02)
    
    djs_plot, [0], [0], /nodata, position=pos, xsty=5, ysty=5, $
      xrange=[0.5,2.0], yrange=[7.5,9.5]
;   legend, 'z=0.1', /right, /bottom, box=0, margin=0
; SFH mass density, zform=1.2
    rhomin1 = interpol(alog10(final[0].sfh_rhostar1_min[indx]),gamma[indx],gammaxis)
    rhomax1 = interpol(alog10(final[0].sfh_rhostar1_max[indx]),gamma[indx],gammaxis)
    rhobest1 = interpol(alog10(final[0].sfh_rhostar1_best[indx]),gamma[indx],gammaxis)
    rhomin1_twosigma = rhomin1+2*(rhomin1-rhobest1)
;   shade_rhostar, gammaxis, rhomin1, rhomax1, color='wheat', /fill
;   shade_rhostar, gammaxis, rhomin1-0.03, rhomin1+0.03, color='wheat', /fill
; SFH mass density, zform=5
    rhomin5 = interpol(alog10(final[0].sfh_rhostar5_min[indx]),gamma[indx],gammaxis)
    rhomax5 = interpol(alog10(final[0].sfh_rhostar5_max[indx]),gamma[indx],gammaxis)
    shade_rhostar, gammaxis, rhomin5, rhomax5, color='powder blue', /fill
; measured mass density
    rhostar = get_rhostar(final[0],rhostar_err=rhostar_err,/local)
    rhomax = interpol(rhostar[indx]+rhostar_err[indx]/2,gamma[indx],gammaxis)
    rhomin = interpol(rhostar[indx]-rhostar_err[indx]/2,gamma[indx],gammaxis)
    shade_rhostar, gammaxis, rhomin, rhomax, color='grey', /fill
; show the final results as a hatched region
; zform=5
    xx1 = interpol(gammaxis,rhomin-rhomax5,0.0)
    yy1 = interpol(rhomin,gammaxis,xx1)
    xx2 = interpol(gammaxis,rhomax-rhomax5,0.0)
    yy2 = interpol(rhomax,gammaxis,xx2)
    xx3 = interpol(gammaxis,rhomax-rhomin5,0.0)
    yy3 = interpol(rhomax,gammaxis,xx3)
    xx4 = interpol(gammaxis,rhomin-rhomin5,0.0)
    yy4 = interpol(rhomin,gammaxis,xx4)
    polyfill, [xx1,xx2,xx3,xx4], [yy1,yy2,yy3,yy4], $
      /data, /line_fill, orientation=45, color=fsc_color('navy',100), $
      spacing=0.11
    polyfill, [xx1,xx2,xx3,xx4], [yy1,yy2,yy3,yy4], $
      /data, /line_fill, orientation=135, color=fsc_color('navy',100), $
      spacing=0.11
;; zform=1.2
;    xx1 = interpol(gammaxis,rhomin-rhomax1,0.0)
;    yy1 = interpol(rhomin,gammaxis,xx1)
;    xx2 = interpol(gammaxis,rhomax-rhomax1,0.0)
;    yy2 = interpol(rhomax,gammaxis,xx2)
;    xx3 = interpol(gammaxis,rhomax-rhomin1,0.0)
;    yy3 = interpol(rhomax,gammaxis,xx3)
;    xx4 = interpol(gammaxis,rhomin-rhomin1,0.0)
;    yy4 = interpol(rhomin,gammaxis,xx4)
;    polyfill, [xx1,xx2,xx3,xx4], [yy1,yy2,yy3,yy4], $
;      /data, /line_fill, orientation=45, color=fsc_color('saddle brown',100), $
;      spacing=0.11
;    polyfill, [xx1,xx2,xx3,xx4], [yy1,yy2,yy3,yy4], $
;      /data, /line_fill, orientation=135, color=fsc_color('saddle brown',100), $
;      spacing=0.11
;; plot lines on the boundaries of the shaded regions because
;; it's perty
;    djs_oplot, gammaxis, rhomin1, line=0, thick=5, $ ; 1-sigma
;      color=fsc_color('wheat',101)
;    djs_oplot, gammaxis, rhomin1_twosigma, line=2, thick=10, $ ; 3-sigma
;      color=fsc_color('saddle brown',101)
;;   djs_oplot, gammaxis, rhobest1, line=0, thick=10, $
;;     color=fsc_color('saddle brown',101)
;;   djs_oplot, gammaxis, rhomax1, line=0, thick=5, color='saddle brown'

    djs_oplot, gammaxis, rhomin5, line=0, thick=5, color='navy'
    djs_oplot, gammaxis, rhomax5, line=0, thick=5, color='navy'
    djs_oplot, gammaxis, rhomin, line=0, thick=5, color='black'
    djs_oplot, gammaxis, rhomax, line=0, thick=5, color='black'

; overplot the axes
    djs_plot, [0], [0], /nodata, /noerase, position=pos, xsty=1, ysty=1, $
      xrange=[0.5,2.0], yrange=[7.5,9.5], $
      xtitle='IMF Slope \Gamma (0.5-120 M_{\odot})', $
      ytitle=cosmicimf_rhotitle()

; some labels
    xyouts, 0.7, 8.45, textoidl('\rho_{*} (z=0)'), $
      /data, align=0.0, charsize=1.6, orientation=8, charthick=2.5
    xyouts, 0.72, 7.7, textoidl('!MI!N\rho_{SFR} (z=0-5)'), $
      /data, align=0.0, charsize=1.6, orientation=52, charthick=2.5, $
      color=djs_icolor('black')
;   xyouts, 1.2, 8.1, textoidl('1\sigma'), /data, align=0.0, $
;     charsize=1.6, orientation=52, charthick=2.5
;   xyouts, 1.4, 8.1, textoidl('2\sigma'), /data, align=0.0, $
;     charsize=1.6, orientation=52, charthick=2.5
;   xyouts, 0.92, 7.65, textoidl('!MI!N\rho_{SFR} (z=0-1.2)'), $
;     /data, align=0.0, charsize=1.6, orientation=52, charthick=2.5

;   im_legend, ['Best fit: \Gamma=1.15\pm0.15','Upper limit: \Gamma<1.4\pm0.1'], $
    im_legend, ['\Gamma=1.15\pm0.15'], textcolor=keycolor, $
      /right, /bottom, box=0, charsize=1.8, margin=0.5
;   xyouts, 1.3, 7.7, textoidl('Upper limit: \Gamma<1.4\pm0.1'), $
;     align=0.0, charsize=1.6
;   xyouts, 1.3, 7.6, textoidl('Best fit: \Gamma=1.15\pm0.15'), $
;     align=0.0, charsize=1.6

;   oplot, 1.15*[1,1], !y.crange, line=0
;   oplot, 1.3*[1,1], !y.crange, line=5
;   oplot, 1.0*[1,1], !y.crange, line=5
    
; plot the results for Salpeter as an inset
    pos2 = [0.28,0.7,0.43,0.9]
    djs_plot, [0], [0], /nodata, /noerase, position=pos2, $
      xsty=1, ysty=1, xrange=1.35+[-0.1,0.1], yrange=[8.4,9.25], $
      charsize=1.3, xtickname=replicate(' ',10), xtitle='', $
      ytitle='log (\rho_{*})', charthick=2.5
    xyouts, (pos2[2]-pos2[0])/2.0+pos2[0], pos2[1]-0.035, 'Salpeter IMF', $
      charsize=1.3, charthick=2.5, /normal, align=0.5
;   oplot_salpeter_box, alog10(final[0].sfh_rhostar1_min[0]), $
;     alog10(final[0].sfh_rhostar1_max[0]), line=0, thick=8, $
;     color='wheat', outcolor='saddle brown'
    oplot_salpeter_box, alog10(final[0].sfh_rhostar5_min[0]), $
      alog10(final[0].sfh_rhostar5_max[0]), line=0, thick=8, $
      color='powder blue', outcolor='navy'
    rhostar = get_rhostar(final[0],rhostar_err=rhostar_err,/local)
;   plots, 1.35, rhostar[0], psym=symcat(15), symsize=2
    if keyword_set(keynote) then salpcolor = fsc_color('black',50) else $
      salpcolor = djs_icolor('default')
    oploterror, 1.35, rhostar[0], rhostar_err[0], $
      psym=symcat(16), symsize=1.6, color=salpcolor, $
      errcolor=salpcolor

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,suffix,'.pdf'), /sh
       rmfile, psfile
    endif

return
end

;
;; multipanel plot showing rho* vs redshift    
;    psfile = paperpath+'imf_rho_vs_redshift'+suffix
;    im_plotconfig, 22, pos, psfile=psfile, width=3.2*[1,1,1,1], $
;      height=3.0*[1,1,1], xmargin=[1.3,0.4]
;    
;    xrange = [-0.01,0.75]
;    yrange = [7.5,10.2]
;
;; Salpeter + every other IMF slope
;    indx = [0,lindgen(nimf/2)*2+1]
;    nindx = n_elements(indx)
;    
;    for ii = 0, nindx-1 do begin
;       if (ii eq 0) or (ii eq 4) or (ii eq 8) then begin
;          if (ii eq 4) then ytitle = textoidl('log_{10}(\rho_{*}/h_{70} M_{\odot} Mpc^{-3})')
;          delvarx, ytickname 
;       endif else begin
;          ytitle = ''
;          ytickname = replicate(' ',10)
;       endelse
;       if (ii ge 8) then begin
;          xtitle = 'Redshift'
;          delvarx, xtickname 
;       endif else begin
;          xtitle = ''
;          xtickname = replicate(' ',10)
;       endelse
;
;       djs_plot, [0], [0], /nodata, noerase=(ii gt 0), position=pos[*,ii], $
;         xsty=1, ysty=1, xrange=xrange, yrange=yrange, xtickname=xtickname, $
;         ytickname=ytickname, xtitle=xtitle, ytitle=ytitle
;         
;       rhostar = alog10(final.rhostar[indx[ii]]) ; need CV errors!
;       rhostar_err = final.rhostar_err[indx[ii]]/final.rhostar[indx[ii]]/alog(10)
;       sfh_rhostar = alog10(final.sfh_rhostar[indx[ii]]) ; need CV errors!
;       sfh_rhostar_err = final.sfh_rhostar_err[indx[ii]]/final.sfh_rhostar[indx[ii]]/alog(10)
;       oploterror, zbins.zbin, rhostar, rhostar_err, psym=-symcat(16), $
;         symsize=2.0, color=fsc_color('dodger blue',101), errcolor=fsc_color('dodger blue',101)
;       oploterror, zbins.zbin, sfh_rhostar, sfh_rhostar_err, psym=-symcat(15), $
;         symsize=2.0, color=fsc_color('firebrick',102), errcolor=fsc_color('firebrick',102)
;
;       if (ii eq 0) then label = 'Salpeter' else $ ; (0.1-100 M_{\odot})' else $
;         label = '\alpha_{2}='+string(peg[indx[ii]].slope,format='(F4.2)')
;       im_legend, label, /right, /top, box=0, margin=0, charsize=1.5
;    endfor
;       
;    im_plotconfig, /psclose
;    


;; this should be in COSMICIMF_FINAL!
;    chi2 = fltarr(nimf)
;    scale = fltarr(nimf)
;    for ii = 0, nimf-1 do begin
;       var = final.rhostar_err[ii]^2 + final.sfh_rhostar_err[ii]^2
;       scale[ii] = total(final.rhostar[ii]*final.sfh_rhostar[ii]/var,/double)/$
;         total(final.sfh_rhostar[ii]^2/var,/double)
;       chi2[ii] = total((final.rhostar[ii]-final.sfh_rhostar[ii])^2/var,/double)
;;      chi2[ii] = total((final.rhostar[ii]-scale[ii]*final.sfh_rhostar[ii])^2/var,/double)
;    endfor
;
;    findchi2min, peg[1:nimf-1].slope, chi2[1:nimf-1], $
;      chi2min, alpha2, alpha2_err
;    djs_plot, peg[1:nimf-1].slope, chi2[1:nimf-1], psym=-6, $
;      xsty=3, ysty=3, /ylog, yrange=[1,max(chi2)]
;    plots, peg[0].slope, chi2[0], psym=symcat(16), color=djs_icolor('orange'), symsize=2.5
;    djs_oplot, alpha2*[1,1], 10^!y.crange, line=0, color='red'
;    djs_oplot, !x.crange, chi2min*[1,1], line=0, color='red'

