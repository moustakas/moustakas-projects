pro fit_munoz_profiles
; jm10feb28ucsd - fit the surface-brightness profiles of the galaxies
; in SINGS for which we have derived a metallicity gradient

    munozpath = getenv('CATALOGS_DIR')+'/09munoz/'
    table2 = im_read_fmr(munozpath+'table2.dat')
    table3 = im_read_fmr(munozpath+'table3.dat')

    singspath = sings_path(/projects)+'log12oh/'
    result1 = mrdfits(singspath+'sings_log12oh_v7.1.fits.gz',1)
    keep = where(result1.gradient_flag eq 1,ngal)
    result = result1[keep]

    out = {$
      galaxy:   '',$
      band:     '',$
      frac_02: 0.0,$
      frac_10: 0.0,$
      frac_40: 0.0}
    out = replicate(out,ngal)
    out.galaxy = strtrim(result.galaxy,2)

    psfile = singspath+'munoz_profiles.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.7
    for ii = 0, ngal-1 do begin
       gal = strtrim(result[ii].galaxy,2)
       r25 = result[ii].d25*60.0/2.0 ; R25 [arcsec]
       match = where(gal eq strtrim(table3.name,2),nmatch)
       if (nmatch ne 0) then begin
          rr25 = table3[match].arad/r25
          sb = table3[match].g
          sb_err = table3[match].e_g
          out[ii].band = 'g'
       endif else begin
          match = where(gal eq strtrim(table2.name,2),nmatch)
          if (nmatch ne 0) then begin
             rr25 = table2[match].arad/r25
             if (total(table2[match].b gt -900.0) ne 0.0) then begin
                sb = table2[match].b
                sb_err = table2[match].e_b
                out[ii].band = 'B'
             endif else begin
                sb = table2[match].nuv
                sb_err = table2[match].e_nuv
                out[ii].band = 'NUV'
             endelse
          endif
       endelse
; fit a model
       good = where(sb gt -900.0,ngood)
       rr25 = rr25[good]
       sb = sb[good]
       sb_err = sb_err[good]
       flux = 10^(-0.4*sb)
       ferr = sb_err*flux*alog(10.0)

       parinfo = replicate({value: 0.0D, limited: [0,0], limits: [0.0D,0.0D]},2)
       parinfo.value = [alog(1D-9),1.0]
;      parinfo.value = [1D-9,1.0]
;      parinfo[1].limited = 1
;      parinfo[1].limits = [0.0,10.0]
       range = where((rr25 gt 0.1) and (rr25 lt 0.9))

       rad = rr25[range]
       lnflux = alog(flux[range])
       lnferr = ferr[range]/flux[range]*0.0+1.0
       
       params = mpfitexpr('P[0] - X/P[1]',rad,$
         lnflux,lnferr,parinfo=parinfo)
;      params = mpfitexpr('P[0]*exp(-X/P[1])',rr25[range],$
;        flux[range],flux_err[range],parinfo=parinfo)
       rraxis = im_array(0.0,1.2,0.01)
       lnmodel = params[0] - rraxis/params[1]
       model = -2.5*alog10(exp(lnmodel))
;      djs_plot, rr25, sb, psym=6, yrange=reverse(minmax(sb[range]))
;      djs_oplot, rraxis, -2.5*alog10(model), color='red'
       
       tot = im_integral(rr25,flux,0.0,1.0)
       out[ii].frac_10 = im_integral(rr25,flux,0.0,0.1)/tot
       out[ii].frac_40 = im_integral(rr25,flux,0.0,0.4)/tot
; QAplot
       xrange = [0,1.2]
       inrange = where(rr25 lt xrange[1])
       yrange = reverse(minmax(sb[inrange]))
       ploterror, rr25, sb, sb_err, position=pos, $
         psym=6, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
         title=gal, xtitle=textoidl('R/R_{25}'), $
         ytitle=textoidl('\mu (R) ('+out[ii].band+' mag arcsec^{-2})')
       djs_oplot, rraxis, model, color='red', thick=2
;      djs_oplot, rraxis, -2.5*alog10(model), color='red', thick=2
       djs_oplot, 0.0*[1,1], !y.crange, color='grey', line=0, thick=2
       djs_oplot, 1.0*[1,1], !y.crange, color='cyan', line=0, thick=4
       djs_oplot, 0.4*[1,1], !y.crange, color='orange', line=5, thick=4
       djs_oplot, 0.1*[1,1], !y.crange, color='dark green', line=5, thick=4
    endfor
    im_plotconfig, psfile=psfile, /psclose, /gzip

    struct_print, out
    im_mwrfits, out, singspath+'munoz_profiles.fits', /clobber
    
stop    
    
return
end
    
