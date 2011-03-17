pro atlas_specres, debug=debug
; jm04may03uofa
; derive the spectral resolution for the atlas using NGC4736

; ---------------------------------------------------------------------------    
;   atlas = read_integrated(/silent)
;   srt = where(atlas.continuum_snr gt 40.0)
;   atlas = atlas[srt]
;   list = atlas_path(/atlas1d)+atlas.specfile
;   for i = 0L, n_elements(list)-1L do spawn, ['cp -fp '+list[i]+' ./'], /sh
; ---------------------------------------------------------------------------    
    
    light = 2.99792458D5        ; speed of light [km/s]
    fwhm2sig = 2.0*sqrt(2.0*alog(2.0))

    galaxy = file_search('*.fits')
;   galaxy = file_search('*.fits.gz')
;   galaxy = ['ugc_09081.fits.gz','ngc_4736.fits.gz']
    ngalaxy = n_elements(galaxy)
    
; compute the number of wavebins

    binstart = 4100.0
    binend = 6400.0

    binsize = 200.0

    bin1 = findgen(fix((binend-binstart)/binsize)+1)*binsize+binstart
    bin2 = bin1+binsize*2
    nbins = n_elements(bin1)
    
;   wavebins = [3750.0,4050.0,4350.0,4650.0,4950.0,5250.0,5550.0,5850.0,6150.0,6450.0]
;   wavebins = findgen(2400.0/binsize)*binsize+4100.0
;   nbins = n_elements(wavebins)

    result = {$
      galaxy:      '', $
      vdisp:       fltarr(nbins), $
      vdisp_err:   fltarr(nbins), $
      specres:     fltarr(nbins), $
      specres_err: fltarr(nbins), $
      specwave:    fltarr(nbins)}
    result = replicate(result,ngalaxy)
    result.galaxy = galaxy

    if keyword_set(debug) then begin
       window, 0, xs=500, ys=400
       window, 2, xs=500, ys=400
    endif

    for j = 0L, ngalaxy-1L do begin
       
       scube = rd1dspec(galaxy[j])

       objflux = scube.spec     ; object flux [erg/s/cm^2/A]
       objferr = scube.sigspec  ; object sigma
       objivar = 1.0/objferr^2.0 ; object sigma
       wave = scube.wave        ; wavelength [A]
       npix = scube.npix        ; wavelength [A]

; normalize everything    
       
       norm = djs_median(objflux)
       flux = objflux/norm
       ferr = objferr/norm
       ivar = objivar*norm^2
       
; de-redshift

       zobj = sxpar(scube.header,'Z')
       logwave = alog10(wave) - alog10(1+zobj)

       icleanup, scube
       
       for i = 0L, nbins-1L do begin

          w = where((10^logwave gt bin1[i]) and (10^logwave lt bin2[i]))
;         w = where((10^logwave gt wavebins[i]) and (10^logwave lt wavebins[i+1]))

          wlogwave = logwave[w]
          wflux = flux[w]
          wivar = ivar[w]

          result[j].specwave[i] = djs_mean(10^wlogwave)
          
; resample to constant log-lambda    
          
          coeff0 = wlogwave[0]
          coeff1 = 1D-4
          newnpix = long(1.0D + (max(wlogwave)-min(wlogwave))/coeff1)
          newlogwave = coeff0 + coeff1*findgen(newnpix)

          combine1fiber, wlogwave, wflux, wivar, newloglam=newlogwave, $
            newflux=newflux, newivar=newivar

          if keyword_set(debug) then wset, 0
          vdans = vdispfit(newflux,newivar,newlogwave,zobj=0.0,npoly=npoly,$
            eigenfile='spEigenElodie.fits',eigendir=eigendir,columns=lindgen(24),$
            yfit=dispflux,doplot=debug)

          if keyword_set(debug) then begin
             
             wset, 2
             djs_plot, 10^newlogwave, newflux, ps=10, xsty=3, ysty=3, $
               charthick=2.0, charsize=1.5
             djs_oplot, 10^newlogwave, dispflux, ps=10, color='yellow', thick=3.0
             cc = get_kbrd(1)

          endif
             
          result[j].vdisp[i] = vdans.vdisp
          result[j].vdisp_err[i] = vdans.vdisp_err
          result[j].specres[i] = result[j].specwave[i]*vdans.vdisp/light*fwhm2sig
          result[j].specres_err[i] = result[j].specwave[i]*vdans.vdisp_err/light*fwhm2sig
          
       endfor

    endfor

    window, 0
    plotsym, 8, 1.0
    indx = 0L
    ploterror, result[indx].specwave, result[indx].specres, result[indx].specres_err, $
      ps=8, xsty=3, ysty=3, xthick=2.0, ythick=2.0, charthick=2.0, $
      charsize=1.5, yrange=[5,15]
    for k = 0L, ngalaxy-1L do oploterror, result[k].specwave, $
      result[k].specres, result[k].specres_err, ps=8

    window, 2
    plot, result[0].specwave, djs_median(result.specres,2), $
      ps=8, xsty=3, ysty=3, xthick=2.0, ythick=2.0, charthick=2.0, $
      charsize=1.5, yrange=[5,15]

    stop

return
end    
