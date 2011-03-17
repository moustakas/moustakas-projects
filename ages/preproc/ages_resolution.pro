pro ages_resolution

    arcwave = mrdfits('arc.fits',0,/silent)
    arcflux = mrdfits('arc.fits',1,/silent)
    arcinvvar = mrdfits('arc.fits',2,/silent)

    nfiber = 300L
    odd = lindgen(nfiber/2)*2 ; even indices
    even = odd+1              ; odd indices

    fluxodd = total(arcflux[*,odd],2)
    ferrodd = sqrt(fluxodd)
    waveodd = djs_median(arcwave[*,odd],2)
    
    minwave = min(waveodd)
    maxwave = max(waveodd)
    npix = n_elements(waveodd)
    dwave = (maxwave-minwave) / (npix-1.0) ; mean dispersion

    outwave = findgen(npix)*dwave+minwave
    noutpix = n_elements(outwave)

    outfluxodd = interpol(fluxodd,waveodd,outwave)
    outferrodd = sqrt(interpol(ferrodd^2.0,waveodd,outwave))
 
    flux = outfluxodd # (fltarr(30)+1)
    ferr = outferrodd # (fltarr(30)+1)
   
    mkhdr, header, flux, /extend
    sxaddpar, header, 'OBJECT', 'Odd Lamp'
    sxaddpar, header, 'CRVAL1', minwave, ' wavelength at CRPIX1'
    sxaddpar, header, 'CD1_1', dwave, ' dispersion [Angstrom/pixel]'
    sxaddpar, header, 'CRPIX1', 1.0, ' reference pixel number'
    sxaddpar, header, 'CTYPE1', 'LINEAR', ' projection type'
    sxaddpar, header, 'DC-FLAG', 0, ' [0 or 1] linear dispersion sampling' 
    
    mwrfits, flux, 'arc_odd.fits', header, /create
    mwrfits, ferr, 'arc_odd.fits'
    spawn, ['gzip -f arc_odd.fits'], /sh

    lampres = im_specres_lamp('arc_odd.fits.gz',psname='lampres_odd.ps',/postscript,/write)

stop
    
    
;   fluxeven = total(arcflux[*,even],2)
;   waveeven = djs_median(arcwave[*,even],2)

;   for i = 0L, nfiber-1L do begin
;
;      wave = arcwave[*,i]
;      flux = arcflux[*,i]
;      invvar = arcinvvar[*,i]
;
;      mask = invvar eq 0
;      invvar = djs_maskinterp(invvar,mask,wave,/const)
;      
;      ferr = 1.0/sqrt(invvar)
;      
;   endfor
    
return
end    
