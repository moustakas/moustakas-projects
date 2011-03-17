pro compare_dale_kara
; jm10mar23ucsd - compare the Dale+07 and Karachentsev+04 B-band photometry
    
    outpath = sings_path(/analysis)
    v2ab = k_vega2ab(filterlist='bessell_B.par',/kurucz,/silent)

; Dale+07    
    ss = sings_read_info()
    flux = rsex(outpath+'sings_photometry_2006dale.sex')
    ferr = rsex(outpath+'sings_photometry_2006dale.uncertainty.sex')

; Kara+04
    kara = read_04kara()

; match
    spherematch, 15.0*im_hms2dec(kara.raj2000), im_hms2dec(kara.dej2000), $
      15.0*im_hms2dec(ss.ra), im_hms2dec(ss.dec), 55.0/3600.0, m1, m2, d12
    niceprint, kara[m1].name, flux[m2].galaxy, d12*3600.0

    dale_b = -2.5*alog10(flux[m2].b*1D-23)-48.6 ; AB
    kara_b = kara[m1].bmag - kara[m1].abg + v2ab ; AB

; make the plot
    xrange = [6.5,17.5]
    yrange = xrange
    residrange = 1.1*[-1,1]
    psfile = outpath+'dale_vs_kara.ps'

    im_plotconfig, 6, pos, psfile=psfile, charsize=1.5

    xx = dale_b
    yy = kara_b
    resid = yy-xx
    print & splog, '### '
    niceprint, kara[m1].name, ss[m2].galaxy, xx, yy, resid

; label the worst outliers       
    out = where(abs(resid) gt im_quantile(abs(resid),quant=0.5),nout)

    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      ytitle='B (Karachentsev+04, AB mag)', xtickname=replicate(' ',10), $
      xrange=xrange, yrange=yrange
    djs_oplot, !x.crange, !y.crange, line=0, color='red'
    djs_oplot, xx, yy, psym=symcat(16), symsize=1.5
    for jj = 0, nout-1 do xyouts, xx[out[jj]], yy[out[jj]], $
      flux[m2[out[jj]]].galaxy, /data, align=0.0, charsize=1.0, $
      charthick=2.5

    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      xtitle='B (Dale+07, AB mag)', ytitle='B Residuals (mag)', $
      xrange=xrange, yrange=residrange
    djs_oplot, !x.crange, [0,0], line=0, color='red'
    djs_oplot, xx, resid, psym=symcat(16), symsize=1.5
    for jj = 0, nout-1 do xyouts, xx[out[jj]], resid[out[jj]], $
      flux[m2[out[jj]]].galaxy, /data, align=0.0, charsize=1.0, $
      charthick=2.5

    im_plotconfig, psfile=psfile, /psclose, /gzip

stop    
    
return
end
    
