pro age_metallicity, postscript=postscript
; jm03dec4uofa
; generate some diagnostic plots    

    ageZ = mrdfits('age_metallicity.fits.gz',1,/silent)
    ngalaxy = n_elements(ageZ)
    psname = 'age_metallicity.ps'
    
    ttchi2 = {$
      chi2_best: fltarr(ngalaxy), $
      chi2_diff: fltarr(2,ngalaxy)}

    for i = 0L, ngalaxy-1L do begin
       ttchi2.chi2_best[i] = ageZ[i].chi2[ageZ[i].best_indx]
       other = lindgen(3) & remove, ageZ[i].best_indx, other
       ttchi2.chi2_diff[*,i] = ageZ[i].chi2[other] - ttchi2.chi2_best[i]
    endfor

; chi2 plots

    im_openclose, psname, postscript=postscript
    
    yrange = [-0.2,1.5]
    xrange = [0.5,10.0]
    
    plotsym, 0, 1
    djs_plot, ttchi2.chi2_best, reform(ttchi2.chi2_diff[0,*]), ps=8, $
      xrange=xrange, yrange=yrange, xthick=2.0, ythick=2.0, charsize=2.0, $
      charthick=2.0, ytitle='\chi^{2} - Minimum \chi^{2}', xtitle='Minimum \chi^{2}', $
      xsty=3, ysty=3
    djs_oplot, ttchi2.chi2_best, reform(ttchi2.chi2_diff[1,*]), ps=8
    oplot, !x.crange, [0,0], line=0, thick=2.0

    if not keyword_set(postscript) then cc = get_kbrd(1)

    yrange = [0,0.6]
    xrange = [-0.2,1.5]

    im_plothist, ttchi2.chi2_diff, xrange=xrange, yrange=yrange, bin=0.1, xsty=3, $
      ysty=1, /fraction, xthick=2.0, ythick=2.0, charsize=2.0, charthick=2.0, $
      thick=4.0, ytitle='Fraction', xtitle=textoidl('\chi^{2} - Minimum \chi^{2}')

;   if not keyword_set(postscript) then cc = get_kbrd(1)
    im_openclose, psname, postscript=postscript, /close

; line flux plots
    
;    tt = {$
;      oii_best: fltarr(ngalaxy), $
;      oii_diff: fltarr(2,ngalaxy)}
;
;    for i = 0L, ngalaxy-1L do begin
;       ttchi2.chi2_best[i] = ageZ[i].chi2[ageZ[i].best_indx]
;       other = lindgen(3) & remove, ageZ[i].best_indx, other
;       ttchi2.chi2_diff[*,i] = ageZ[i].chi2[other] - ttchi2.chi2_best[i]
;    endfor

; other plots
    
;    Zline = ageZ.kd_z_combined
;    good = where(Zline gt 0.0,ngood)
;    Zline = Zline[good]
;    plotageZ = ageZ[good]
;
;    Zcont = alog10(plotageZ.Z_template/0.02)+8.69
;
;    lmc = where(plotageZ.Z_template eq 0.004,nlmc)
;    sun = where(plotageZ.Z_template eq 0.02,nsun)
;    super = where(plotageZ.Z_template eq 0.05,nsuper)

;    plot, [0], [0], xrange=[7.9,9.6], yrange=[0.0,0.7], xsty=3, ysty=3, /nodata, $
;      xthick=2.0, ythick=2.0, charsize=2.0, charthick=2.0
;    im_plothist, Zline[lmc], color=djs_icolor('red'), thick=2.0, bin=0.1, /overplot, /fraction
;    im_plothist, Zline[sun], color=djs_icolor('yellow'), thick=2.0, bin=0.1, /overplot, /fraction
;    im_plothist, Zline[super], color=djs_icolor('blue'), thick=2.0, bin=0.1, /overplot, /fraction

stop    
    
return
end    
