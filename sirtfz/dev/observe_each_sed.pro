pro observe_each_sed, z, l_evol, sedcube, lf, bandcube, flux, type, indx=indx
; jm01januofa
; observe the model SEDs with noise by placing them into a unique
; luminosity function bin

; jm01may12 - added indx keyword to only allow specific LF bins    
    
    common sirtf_simulations

    if not keyword_set(indx) then indx = lindgen(lf.nlf)

    nbands = n_elements(bandcube)     ; number of filters
    nlf = n_elements(indx)            ; number of LF bins
    dlum = dluminosity(z)*3.085678D16 ; luminosity distance (m)

; perturb the central value of the LF bin by the binsize (W/Hz)

    mag_pert = lf.lum_lf[indx] + (randomu(seed,nlf)*(2D)-1D)*lf.logbinsz/2D
    lum_lf_pert = 10D^(mag_pert) * 3.826D26 * 60D / 2.99793D14

    flux = fltarr(nlf,nbands)   ; flux for each band for each LF bin
    type = intarr(nlf)          ; galaxy type (SED number)
    for k = 0L, nlf-1L do begin ; loop on the luminosity function
       
       sedindx = where(sedcube.lfbindx eq indx[k],sedcount) ; match the appropriate SEDs with the LF bin

       if sedcount eq 0L then begin ; there is no SED in this LF bin, so replace with a nearby SED

          nearindx = where(abs(sedcube.lfbindx-indx[k]) eq min(abs(sedcube.lfbindx-indx[k])),nearcount)
          if nearcount ne 0L then rsed = nearindx[floor(randomu(seed,1)*nearcount)] else $ ; random nearby SED
            message, "There was a problem finding a nearby SED."
          
       endif else rsed = sedindx[floor(randomu(seed,1)*sedcount)] ; random SED in that LF bin

       type[k] = rsed  ; SED number

       get_element, *sedcube[rsed].lambda, bandcube.lambda0, xsed ; monochromatic luminosity at each filter wavelength
       sedlum = (*sedcube[rsed].mlum)[xsed]

       for j = 0L, nbands-1L do begin ; loop on the bandpasses

          wband = *bandcube[j].wband ; bandpass wavelength
          rband = *bandcube[j].rband ; bandpass response

          if strupcase(bandcube[j].bandnames) eq '20 CM' then radio = 1B else radio = 0B
          
          restflux = band_flux(0D,*sedcube[rsed].lambda,*sedcube[rsed].mlum,wband,rband,radio=radio) ; flux at z=0
          zflux = band_flux(z,*sedcube[rsed].lambda,*sedcube[rsed].mlum,wband,rband,radio=radio)     ; flux at redshift z
       
          scale = lum_lf_pert[k]/(sedcube[rsed].lum60*60D/2.99793D14)                  ; LF-SED correction

          flux[k,j] = scale * l_evol * sedlum[j] * zflux / restflux / (4.0D*!dpi*dlum*dlum) / 1D-29 ; mJy
          
       endfor
       
    endfor
       
return
end








