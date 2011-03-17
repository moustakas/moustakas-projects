pro compare_dale_munoz
; jm10mar23ucsd - compare the Dale+07 and Munoz-Mateos BVRI photometry
    
    outpath = sings_path(/analysis)
    v2ab = k_vega2ab(filterlist='bessell_'+['B','V','R','I']+'.par',/kurucz,/silent)

; Dale+07    
    flux = rsex(outpath+'sings_photometry_2006dale.sex')
    ferr = rsex(outpath+'sings_photometry_2006dale.uncertainty.sex')
    dale_galaxy = repstr(repstr(repstr(repstr(repstr(flux.galaxy,'HoIX','HOLMBERGIX'),$
      'HoII','HOLMBERGII'),'HoI','HOLMBERGI'),'Mrk33','MRK0033'),'Tol89','TOLOLO89')

; Munoz-Mateos+09    
    table5 = rsex(getenv('CATALOGS_DIR')+'/09munoz/table5.sex')
    match, strtrim(strlowcase(dale_galaxy),2), strtrim(strlowcase(table5.galaxy),2), m1, m2

    table5_bvri = transpose([[table5[m2].b],[table5[m2].v],[table5[m2].r],[table5[m2].i]]) ; AB
    table5_bvri_err = transpose([[table5[m2].b_err],[table5[m2].v_err],[table5[m2].r_err],[table5[m2].i_err]])

    dale_table5_bvri = transpose([[flux[m1].b],[flux[m1].v],[flux[m1].r],[flux[m1].i]])
    dale_table5_bvri_err = transpose([[ferr[m1].b],[ferr[m1].v],[ferr[m1].r],[ferr[m1].i]])
    for ii = 0, 3 do begin
       good = where(dale_table5_bvri[ii,*] gt 0)
       dale_table5_bvri[ii,good] = -2.5*alog10(dale_table5_bvri[ii,good]*1D-23)-48.6 ; AB
       dale_table5_bvri_err[ii,good] = 2.5*dale_table5_bvri_err[ii,good]/dale_table5_bvri[ii,good]/alog(10.0)
    endfor

; also compare against the ugriz converted photometry 
    table6 = rsex(getenv('CATALOGS_DIR')+'/09munoz/table6.sex')
    match, strtrim(strlowcase(dale_galaxy),2), strtrim(strlowcase(table6.galaxy),2), m1_ugriz, m2_ugriz
    b = table6[m2_ugriz].g + 0.3915*(table6[m2_ugriz].g-table6[m2_ugriz].r-0.6102) + 0.2354 ; AB
    v = table6[m2_ugriz].g - 0.7585*(table6[m2_ugriz].g-table6[m2_ugriz].r-0.6102) - 0.3516 ; AB
;   b = table6[m2_ugriz].g + 0.3130*(table6[m2_ugriz].g-table6[m2_ugriz].r) + 0.2271 + v2ab[0] ; AB
;   v = table6[m2_ugriz].g - 0.5784*(table6[m2_ugriz].g-table6[m2_ugriz].r) - 0.0038 + v2ab[1] ; AB
    table6_bv = transpose([[b],[v]])

    dale_table6_bv = transpose([[flux[m1_ugriz].b],[flux[m1_ugriz].v]])
    for ii = 0, 1 do begin
       good = where(dale_table6_bv[ii,*] gt 0)
       dale_table6_bv[ii,good] = -2.5*alog10(dale_table6_bv[ii,good]*1D-23)-48.6 ; AB
    endfor
;   plot, table6[m2_ugriz].g, dale_table6_bv[0,*]-table6[m2_ugriz].g, psym=6
;   niceprint, table6[m2_ugriz].galaxy, table6[m2_ugriz].g, table6_bv[0,*], dale_table6_bv[0,*]
    
; make the plot
    xrange = [[7.5,15],[7,15],[6,14],[5,15]]
    yrange = xrange
    residrange = 1.4*[-1,1]
    band = ['B','V','R','I']
    psfile = outpath+'dale_vs_munoz.ps'
    im_plotconfig, 6, pos, psfile=psfile, charsize=1.5
    for ii = 0, 3 do begin
; from ugriz
       if (ii eq 0) or (ii eq 1) then begin
          good_ugriz = where((table6_bv[ii,*] gt -900.0) and (dale_table6_bv[ii,*] gt -900.0),ngood_ugriz)
          xx_ugriz = dale_table6_bv[ii,good_ugriz]
          yy_ugriz = table6_bv[ii,good_ugriz] ; - v2ab[ii]
          resid_ugriz = yy_ugriz-xx_ugriz
       endif
; BVRI       
       good = where((table5_bvri[ii,*] gt -900.0) and (dale_table5_bvri[ii,*] gt -900.0),ngood)
       xx = dale_table5_bvri[ii,good]
       yy = table5_bvri[ii,good] ; - v2ab[ii]
       resid = yy-xx
       print & splog, '### '+band[ii]
       niceprint, table5[m2[good]].galaxy, xx, yy, resid

       out = where(abs(resid) gt im_quantile(abs(resid),quant=0.65),nout) ; worst outliers

       djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
         ytitle=band[ii]+' (Munoz-Mateos+09, AB mag)', xtickname=replicate(' ',10), $
         xrange=xrange[*,ii], yrange=yrange[*,ii]
       djs_oplot, !x.crange, !y.crange, line=0, color='red'
       if (ii eq 0) or (ii eq 1) then $
         djs_oplot, xx_ugriz, yy_ugriz, psym=symcat(6,thick=5), $
         symsize=1.5, color='orange'
       djs_oplot, xx, yy, psym=symcat(16), symsize=1.5
       for jj = 0, nout-1 do xyouts, xx[out[jj]], yy[out[jj]], $
         flux[m1[good[out[jj]]]].galaxy, /data, align=0.0, charsize=1.0, $
         charthick=2.5
       

       djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
         xtitle=band[ii]+' (Dale+07, AB mag)', ytitle=band[ii]+' Residuals (mag)', $
         xrange=xrange[*,ii], yrange=residrange
       djs_oplot, !x.crange, [0,0], line=0, color='red'
       if (ii eq 0) or (ii eq 1) then $
         djs_oplot, xx_ugriz, resid_ugriz, psym=symcat(6,thick=5), $
         symsize=1.5, color='orange'
       djs_oplot, xx, resid, psym=symcat(16), symsize=1.5
       for jj = 0, nout-1 do xyouts, xx[out[jj]], resid[out[jj]], $
         flux[m1[good[out[jj]]]].galaxy, /data, align=0.0, charsize=1.0, $
         charthick=2.5
    endfor

; compare the B-V colors
    good = where((table5_bvri[0,*] gt -900.0) and (table5_bvri[1,*] gt -900.0) and $
      (dale_table5_bvri[0,*] gt -900.0) and (dale_table5_bvri[1,*] gt -900.0),ngood)
    xx = dale_table5_bvri[0,good]-dale_table5_bvri[1,good]
    yy = table5_bvri[0,good]-table5_bvri[1,good]
    resid = yy-xx
    print & splog, '### B-V'
    niceprint, table5[m2[good]].galaxy, xx, yy, resid

    good_ugriz = where((table6_bv[0,*] gt -900.0) and (dale_table6_bv[0,*] gt -900.0) and $
      (table6_bv[1,*] gt -900.0) and (dale_table6_bv[1,*] gt -900.0),ngood_ugriz)
    xx_ugriz = dale_table6_bv[0,good_ugriz]-dale_table6_bv[1,good_ugriz]
    yy_ugriz = table6_bv[0,good_ugriz]-table6_bv[1,good_ugriz]
    resid_ugriz = yy_ugriz-xx_ugriz
    
    xrange = [-0.2,1.2]
    yrange = xrange
    residrange = 0.9*[-1,1]
    
; label the worst outliers
    out = where(abs(resid) gt im_quantile(abs(resid),quant=0.65),nout)

    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      ytitle='B-V (Munoz-Mateos+09, AB mag)', xtickname=replicate(' ',10), $
      xrange=xrange, yrange=yrange
    djs_oplot, !x.crange, !y.crange, line=0, color='red'
    djs_oplot, xx_ugriz, yy_ugriz, psym=symcat(6,thick=5), $
      symsize=1.5, color='orange'
    djs_oplot, xx, yy, psym=symcat(16), symsize=1.5
    for jj = 0, nout-1 do xyouts, xx[out[jj]], yy[out[jj]], $
      flux[m1[good[out[jj]]]].galaxy, /data, align=0.0, charsize=1.0, $
      charthick=2.5

    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      xtitle='B-V (Dale+07, AB mag)', ytitle='B-V Residuals (mag)', $
      xrange=xrange, yrange=residrange
    djs_oplot, !x.crange, [0,0], line=0, color='red'
    djs_oplot, xx_ugriz, resid_ugriz, psym=symcat(6,thick=5), $
      symsize=1.5, color='orange'
    djs_oplot, xx, resid, psym=symcat(16), symsize=1.5
    for jj = 0, nout-1 do xyouts, xx[out[jj]], resid[out[jj]], $
      flux[m1[good[out[jj]]]].galaxy, /data, align=0.0, charsize=1.0, $
      charthick=2.5
    
    im_plotconfig, psfile=psfile, /psclose, /gzip

stop    
    
return
end
    
