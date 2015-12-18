pro plot_elodieres, postscript=postscript
; jm03may21uofa
; generate a plot    

    plotsym, 0, 1, /fill

    cmrestore, 'elodie_resolution.idlsave'

    newres = resolution
    for i = 0L, nstar-1L do begin
       w = where(newres[*,i] eq 15.0,nw)
       if nw ne 0L then newres[w,i] = median(newres[*,i])
    endfor

    nwave = n_elements(wavebins)
    
    xrange = [wave1,wave2]
    yrange = minmax(sigma*fwhm2sig)

    xerr = replicate(dwave,nstar)
    yerr = replicate(0.5,nstar)
    
    meanres = fltarr(nwave)
    sigmares = fltarr(nwave)

    for j = 0L, nwave-1L do begin

       meanres[j] = djs_mean(newres[j,*])
       sigmares[j] = stddev(newres[j,*])
       
    endfor

    im_openclose, 'elodie_resolution.ps', postscript=postscript
    
    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=2.0, ythick=2.0, charsize=2.0, $
      charthick=2.0, xrange=xrange, yrange=yrange, xtitle='Wavelength ['+angstrom()+']', $
      ytitle='FWHM Resolution ['+angstrom()+']'
;   for i = 0L, nstar-1L do oplot, wavebins, newres[*,i], ps=8
;   for i = 0L, nstar-1L do oploterror, wavebins, newres[*,i], xerr, yerr, ps=8
    oploterror, wavebins, meanres, xerr, sigmares, ps=8, errthick=2.0
    oplot, !x.crange, djs_mean(meanres)*[1,1], line=2, thick=3.0

    im_openclose, postscript=postscript, /close
    
return    
end
