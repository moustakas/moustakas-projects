;------------------------------------------------------------------------------
; Read the Elodie eschelle spectra and put on the SDSS pixelization.
; The spectra are originally in units of normalized F-lambda (erg/s/cm^2/Ang),
; which are also the units that we output.
;------------------------------------------------------------------------------
function im_read_elodie, filename, wave=wave, hdr=hdr

    elodie_path = filepath('',root_dir=getenv('CATALOGS_DIR'),subdirectory='elodie/LL_ELODIE')
    
    if (n_params() LT 1) then begin
       print, 'Syntax: flux = im_read_elodie(filename, [loglam=, hdr= ])'
       return, -1L
    endif

; The Elodie spectra have 13501 0.2-Ang pixels from 4100 - 6800 Ang

    objflux = mrdfits(elodie_path+filename, 0, hdr,/silent)
    if (NOT keyword_set(objflux)) then message, 'File not found: ' + filename
    objwave = sxpar(hdr,'CRVAL1') + dindgen(sxpar(hdr,'NAXIS1')) * sxpar(hdr,'CDELT1')

    igood = where(finite(objflux))
    newflux = objflux[igood]
    wave = objwave[igood]
    
;;;    dloglam = 1.d-4
;;;    nsubsamp = 8
;;;    nnew = 2220L
;;;    loglam0 = 3.6120
;;;    sigma = 1.0                 ; Gaussian sigma in pixels on the final sampled vector
;;;
;;;                                ; Convert to heliocentric wavelengths in vacuum
;;;    airtovac, objwave
;;;; It looks like the spectrum is already de-redshifted, so don't do the below
;;;;   cspeed = 2.99792458d5
;;;;   vhelio = sxpar(hdr,'VR')
;;;;   objwave = objwave * (1.d0 + vhelio / cspeed)
;;;
;;;                                ; Resample to a vector starting at bigloglam0,
;;;                                ; spaced every dloglam/nsubsamp
;;;;   bigloglam0 = loglam0 - 0.5 * (nsubsamp - 1) * dloglam / nsubsamp
;;;;   ibigpix = (alog10(objwave) - bigloglam0) / (dloglam/nsubsamp)
;;;;   newflux = fltarr(nnew*nsubsamp)
;;;;   igood = where(finite(objflux))
;;;;   populate_image, newflux, ibigpix[igood], weights=objflux[igood], $
;;;;    assign='cic'
;;;
;;;    bigloglam0 = loglam0 - 0.5 * (nsubsamp - 1) * dloglam / nsubsamp
;;;    bigloglam = bigloglam0 + lindgen(nnew*nsubsamp) * (dloglam/nsubsamp)
;;;    newflux = interpol(objflux[igood], alog10(objwave[igood]), bigloglam)
;;;
;;;                                ; Convolve with a gaussian
;;;    nkpix = long(10*sigma*nsubsamp)
;;;    kern = exp( -0.5 * (findgen(nkpix*2+1) - nkpix)^2 / (sigma*nsubsamp)^2 )
;;;    kern = kern / total(kern)
;;;    newflux = convol(newflux, kern, /center)
;;;
;;;                                ; Now rebin to the final sampling
;;;    newflux = rebin(newflux, nnew) * nsubsamp
;;;    newloglam = loglam0 + lindgen(nnew) * dloglam
;;;
;;;                                ; Set any pixels outside of the wavelength range to zero.
;;;                                ; Insist that we be at least 3*sigma from the enpoint measurements.
;;;    minloglam = alog10( min(objwave[igood]) ) + 3 * sigma * dloglam
;;;    maxloglam = alog10( max(objwave[igood]) ) - 3 * sigma * dloglam
;;;    newflux = newflux * (newloglam GT minloglam) * (newloglam LT maxloglam)

return, newflux
end
