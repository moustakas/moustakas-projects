pro mc_abundances, nmonte=nmonte, postscript=postscript
; jm04jun28uofa
; use a Monte Carlo technique to quantify the success rate of the P01
; abundance diagnostic

    if keyword_set(postscript) then postthick = 5.0 else postthick = 2.0
    
    hii = read_hii_regions(/nolog,/nosdss)
    strong = make_linefit_structure(hii)
    keep = where((hii.ZT_log12oh gt -9000.0) and (hii.nii_6584_h_alpha gt -9000.0) and $
      (hii.nii_6584_oii gt -9000.0) and (hii.oii_h_beta gt -9000.0) and $
      (hii.oiii_5007_h_beta gt -9000.0),nhii)

    hii = hii[keep]
    strong = strong[keep]
    montestrong = strong

    linename = strong[0].linename
    nline = n_elements(linename)
    
    if (n_elements(nmonte) eq 0L) then nmonte = 100L

    Z_array = dblarr(nhii,nmonte)

    result = {$
      dZ_syst:  0.0D, $
      dZ_rand:  0.0D, $
      Z_mean:   0.0D, $
      Z_median: 0.0D, $
      Z_error:  0.0D, $
      outliers:      0.0D  $
      }
    result = replicate(result,nhii)

    splog, 'Computing Monte Carlo abundances.'
    for imonte = 0L, nmonte-1L do begin

       print, format='("Monte Carlo ",I3,"/",I3,".",A1,$)', $
         imonte+1, nmonte, string(13b)
       
       if (imonte gt 0L) then for jline = 0L, nline-1L do begin
          flux = (strong.(3+jline))[0,*]
          ferr = (strong.(3+jline))[1,*]
          newflux = flux + randomn(seed,1,nhii)*ferr
          montestrong.(3+jline) = transpose([ [ [newflux] ], [ [ferr] ] ])
       endfor

       abund = im_abundance(montestrong,snrcut=1.0,/silent)
       Z_array[*,imonte] = abund.Z_12OH_P01

    endfor

    Z_input = reform(Z_array[*,0])
    Z_Te = hii.ZT_log12oh

    for ihii = 0L, nhii-1L do begin

       array = reform(Z_array[ihii,*])
       good = where(array gt -900.0)
       
       djs_iterstat, array[good], sigrej=3.0, mean=mn, median=md, sigma=sig, mask=mask

       result[ihii].Z_mean = mn
       result[ihii].Z_median = md
       result[ihii].Z_error = sig
       result[ihii].outliers = 100.0*total(mask eq 0)/nmonte

       result[ihii].dZ_rand = mn - Z_input[ihii]
       result[ihii].dZ_syst = mn - Z_Te[ihii]

    endfor

    plotsym, 0, 0.5, /fill
    window, 0, xs=450, ys=450
    ploterror, Z_input, result.dZ_rand, result.Z_error, $
      ps=8, xsty=3, ysty=3, xrange=[7.0,9.2], yrange=[-0.7,0.7], $
      xtitle=textoidl('12 + log (O/H) [T_{e}]'), $
      ytitle='Random Error', xthick=postthick, ythick=postthick, $
      charsize=2.0, charthick=postthick
    djs_oplot, !x.crange, [0,0], line=0, thick=1.0, color='yellow'
    
    window, 2, xs=450, ys=450
    ploterror, Z_Te, result.dZ_syst, result.Z_error, $
      ps=8, xsty=3, ysty=3, xrange=[7.0,9.2], yrange=[-0.7,0.7], $
      xtitle=textoidl('12 + log (O/H) [T_{e}]'), $
      ytitle='Systematic Error', xthick=postthick, ythick=postthick, $
      charsize=2.0, charthick=postthick
    djs_oplot, !x.crange, [0,0], line=0, thick=1.0, color='yellow'

; what fraction of lower branch objects were placed on the upper
; branch and vice-versa

    lofail = where((Z_Te le 7.95) and (result.Z_mean ge 8.2),nlofail)
    hifail = where((Z_Te ge 8.2) and (result.Z_mean le 7.95),nhifail)
    
;   window, 2, xs=450, ys=450
;   plot, Z_input, result.outliers, ps=8, xsty=3, ysty=3, $
;     xrange=[7.0,9.2], yrange=[1.0,15.0], /ylog
;   djs_oplot, !x.crange, [0,0], line=0, thick=1.0, color='yellow'
    
stop    

return
end
    
