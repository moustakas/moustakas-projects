pro unpack_kenn92, doplot=doplot, wfits=wfits
; jm03may27uofa
; jm05aug02uofa - updated

; write individual FITS files for Rob's 1992 galaxy atlas which he
; gave me as a single multi-extension FITS file

    light = 2.99792458D5 ; speed of light [km/s]

    analysis_path = kenn92_path(/analysis)
    outpath = kenn92_path(/spec1d)

    rob = mrdfits(analysis_path+'rob2.fits',0,h,/silent)
    info = mrdfits(analysis_path+'kenn92_ned.fits.gz',1,/silent)
    ngalaxy = n_elements(info)
    
; the following objects were observed with the IRS 0.9-meter at very
; low spectral resolution; the remainder were observed with two
; grating systems on the 2.3-meter

    irs = [$
      'ngc4472',$
      'ngc5866',$
      'ngc3690',$
      'ngc6052',$
      'ngc3471',$
      'ngc5996',$
      'mrk0035',$
      'mrk0059',$
      'mrk0071',$
      'mrk0487',$
      'ngc7469',$
      'ngc6764',$
      'ngc3303',$
      'ngc4194']

; read the literature catalog of Kennicutt (1992a) galaxies 

    k92 = read_92kennicutt()
    
;   for i = 42, ngalaxy-1L do begin
    for i = 0L, ngalaxy-1L do begin

; cross-match with the literature data; ***only write out galaxies
; that Rob measured in Kennicutt (1992b)

       match = where(strtrim(info[i].nedgalaxy,2) eq strtrim(k92.alt_galaxy,2),nmatch)
       if (nmatch ne 1L) then continue 
;      if (nmatch eq 0L) then print, info[i].galaxy+'   '+info[i].nedgalaxy
       
; build the output spectrum       
       
       flux = reform(rob[*,i]) ; spectrum
       npix = n_elements(flux)
       wave = findgen(npix)*2.0+3650.0 ; wavelength array [Angstrom]

; repair crummy pixels

       zero = where(flux eq 0.0,nzero)
       if (nzero ne 0L) then begin
          bigindx = lindgen(npix)
          blue = where(zero lt 150L,nblue)
          if (nblue ne 0L) then remove, zero[blue], bigindx
          
          flux = flux[bigindx]
          wave = wave[bigindx]
          npix = n_elements(flux)
       endif

       zero = where(flux eq 0.0,nzero)
       if (nzero ne 0L) then begin
          bigindx = lindgen(npix)
          red = where(zero gt npix-151L,nred)
          if (nred ne 0L) then remove, zero[red], bigindx
          
          flux = flux[bigindx]
          wave = wave[bigindx]
          npix = n_elements(flux)
       endif

       neg = where(flux le 0.0,nneg)
       if nneg ne 0L then begin

          mflux = dkboxstats(flux,boxstat='mean',xwidth=30)
          flux[neg] = mflux[neg]

          newneg = where(mflux[neg] lt 0.0)
          if newneg[0] ne -1L then stop
          
       endif

; distinguish between low- and high-resolution galaxies

;      print, i, info[i].galaxy
;      if strmatch(strtrim(info[i].galaxy,2),'*0071*',/fold) then stop
       
       lowres = where(strmatch(irs,'*'+info[i].galaxy+'*',/fold) eq 1B,nlowres)
       
       if (nlowres eq 1L) then begin
          fitsfile = strcompress(strlowcase(info[i].galaxy),/remove)+'_0.9m.fits'
       endif else begin
          fitsfile = strcompress(strlowcase(info[i].galaxy),/remove)+'_2.3m.fits'
          
; degrade the spectral resolution in the blue

          degindx = where(wave lt 5050.0)
          degflux = im_gauss_broaden(wave[degindx],flux[degindx],4.5,7.0)       
          flux[degindx] = degflux

       endelse
       
; de-redden for foreground Galactic reddening

       glactc, 15.0D*im_hms2dec(info[i].ra), im_hms2dec(info[i].dec), 2000.0, gl, gb, 1, /degree
       ebv_mw = dust_getval(gl,gb,/interp)
       kl = k_lambda(wave,/odonnell,R_V=3.1)
       flux = flux*10^(0.4*ebv_mw*kl)
       
; generate the noise spectrum

       get_element, wave/(1+info[i].z), [5100,6500], xx
       djs_iterstat, flux[xx[0]:xx[1]], sigma=sig, mean=mn, sigrej=3.0
       snr = mn/sig
       
       get_element, wave, djs_mean(wave), snrindx
       ferr = flux/sqrt((flux/flux[snrindx])>0)/snr

       flux = flux*1D-15
       ferr = ferr*1D-15

; output header       
       
       mkhdr, header, float(flux), /extend
       sxdelpar, header, 'COMMENT'
       sxaddpar, header, 'OBJECT', strtrim(info[i].galaxy,2)
       sxaddpar, header, 'GALAXY', strtrim(info[i].galaxy,2), 'galaxy name as in Kennicutt 1992'

       sxaddpar, header, 'CRVAL1', min(wave), ' wavelength at CRPIX1' ; 3650
       sxaddpar, header, 'CRPIX1', 1.0, ' reference pixel number'
       sxaddpar, header, 'CD1_1', wave[1]-wave[0], ' dispersion [Angstrom/pixel]'
       sxaddpar, header, 'CDELT1', wave[1]-wave[0], ' dispersion [Angstrom/pixel]'
       sxaddpar, header, 'CTYPE1', 'LINEAR', ' projection type'
       sxaddpar, header, 'Z', info[i].z, ' NED redshift'

       sxaddpar, header, 'RA', info[i].ra, ' right ascension [HMS]'
       sxaddpar, header, 'DEC', info[i].dec, ' declination [DMS]', after='RA'
       sxaddpar, header, 'EPOCH', 2000.0, ' coordinate epoch', format='(F6.1)', after='DEC'
       sxaddpar, header, 'EQUINOX', 2000.0, ' coordinate epoch', format='(F6.1)', after='EPOCH'

       sxaddpar, header, 'SCANLEN', float(k92[match].scanlen), $
         ' scan length [arcsec]', format='(F12.2)'
          sxaddpar, header, 'POSANGLE', float(k92[match].posangle), $
                  ' slit position angle [degrees]', format='(F12.1)'
       sxaddpar, header, 'APERWID', k92[match].aperwid, $
         ' extraction aperture diameter [arcsec]', format='(F12.2)'
       
       if keyword_set(wfits) then begin

          splog, 'Writing '+outpath+fitsfile+'.'
          mwrfits, flux, outpath+fitsfile, header, /create
          mwrfits, ferr, outpath+fitsfile
;e         spawn, ['gzip -f '+outpath+fitsfile], /sh
          
; derive the galaxy redshift and update the header

;         zans = iredshift(fitsfile+'.gz',analysis_path=outpath,npoly=1,zmin=-0.01,$
;           zmax=0.01,/update,doplot=doplot)

       endif else splog, i, info[i].galaxy, snr, npix, format='(I3,A15,3x,I0,3x,I0)'

       if keyword_set(doplot) and (not keyword_set(wfits)) then begin
;         ploterror, wave, 1E15*flux, 1E15*ferr, ps=10, xsty=3, ysty=3, $
          djs_plot, wave, 1E15*flux, ps=10, xsty=3, ysty=3, $
            xthick=2.0, ythick=2.0, charsize=1.5, charthick=2.0, $ ; xrange=[6650,6900], $
            title=info[i].galaxy, xtitle='Wavelength', ytitle='Normalized Flux'
          cc = get_kbrd(1)
       endif
       
    endfor

return
end    
