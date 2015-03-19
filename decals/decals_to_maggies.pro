;+
; NAME:
;   DECALS_TO_MAGGIES
;
; PURPOSE:
;   Convert the Tractor DECam-Legacy photometry to AB maggies.  
;
; INPUTS: 
;   cat - input photometric catalog [NGAL] 
;
; KEYWORD PARAMETERS:
;
; OUTPUTS: 
;   maggies - output maggies [2,NGAL]
;   ivarmaggies - corresponding inverse variance array [8,NGAL]  
;
; COMMENTS:
;   Assume no minimum magnitude uncertainty.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2014 Sep 08, Siena
;
; Copyright (C) 2014, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

pro decals_sdss_to_maggies, cat, maggies, ivarmaggies, psf=psf

    ngal = n_elements(cat)
    dustfactor = 1D-9*10D^(0.4*cat.sdss_extinction)

    if keyword_set(psf) then begin
       maggies = cat.sdss_psfflux*dustfactor
       ivarmaggies = cat.sdss_psfflux_ivar/dustfactor^2
    endif else begin
       cmodelmaggies = cat.sdss_cmodelflux*dustfactor
       cmodelivarmaggies = cat.sdss_cmodelflux_ivar/dustfactor^2

       modelmaggies = cat.sdss_modelflux*dustfactor
       modelivarmaggies = cat.sdss_modelflux_ivar/dustfactor^2

       ratio = cmodelmaggies[2,*]/(modelmaggies[2,*]+(modelmaggies[2,*] eq 0))
       neg = where(modelmaggies[2,*] le 0)
       if (neg[0] ne -1L) then ratio[neg] = 1.0

       factor = rebin(ratio,5,ngal)
       maggies = modelmaggies*factor
       ivarmaggies = modelivarmaggies/factor^2
    endelse

    k_abfix, maggies, ivarmaggies
    k_minerror, maggies, ivarmaggies

    maggies = float(maggies)
    ivarmaggies = float(ivarmaggies)

return
end

pro decals_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist, $
  sdss=sdss, psf=psf, aperture=aperture, apmaggies=apmaggies, $
  apivarmaggies=apivarmaggies, shortwise=shortwise, decam_grz=decam_grz

    ngal = n_elements(cat)
    if (ngal eq 0L) then begin
       doc_library, 'decals_to_maggies'
       return
    endif

; reddening factors from Schlafly & Finkbeiner (2014, private
; communication) for the DECam filter curves
;   red_fac = [3.214,2.165,1.211] ; =grz
;   red_fac = [3.995,3.214,2.165,1.592,1.211,1.064] ; =ugrizY

;   euler, cat.ra, cat.dec, ll, bb, 1
;   extinction = red_fac # dust_getval(ll,bb,/interp,/noloop)

    if keyword_set(decam_grz) then these = [1,2,4] else these = [0,1,2,3,4,5]
    factor = 1D-9/cat.decam_mw_transmission[these]
    dmaggies = float(cat.decam_flux[these]*factor)
    divarmaggies = float(cat.decam_flux_ivar[these]/factor^2)
    dfilterlist = (decals_filterlist())[these]

; unWISE
    if keyword_set(shortwise) then these = [0,1] else these = [0,1,2,3]
    factor = 1D-9/cat.wise_mw_transmission[these]
    wmaggies = float(cat.wise_flux[these]*factor)
    wivarmaggies = float(cat.wise_flux_ivar[these]/factor^2)
    wfilterlist = wise_filterlist(short=keyword_set(shortwise))
    
; add SDSS
    if keyword_set(sdss) then begin
       decals_sdss_to_maggies, cat, smaggies, sivarmaggies, psf=psf
       maggies = [dmaggies,wmaggies,smaggies]
       ivarmaggies = [divarmaggies,wivarmaggies,sivarmaggies]
       filterlist = [dfilterlist,wfilterlist,sdss_filterlist()]
    endif else begin
       maggies = [dmaggies,wmaggies]
       ivarmaggies = [divarmaggies,wivarmaggies]
       filterlist = [dfilterlist,wfilterlist]
    endelse

; optionally get the DECam aperture fluxes
    if keyword_set(aperture) then begin
       naper = n_elements(cat[0].decam_apflux[*,0])

       if keyword_set(decam_grz) then these = [1,2,4] else these = [0,1,2,3,4,5]
       nband = n_elements(these)
       factor = 1D-9/cat.decam_mw_transmission[these]
       factor = rebin(reform(factor,1,nband,ngal),naper,nband,ngal)

       apmaggies = float(cat.decam_apflux[*,these]*factor)
       apivarmaggies = float(cat.decam_apflux_ivar[*,these]/factor^2)

;      these = [1,2,4]
;      nband = n_elements(these)
;      factor = 1D-9/cat.decam_mw_transmission[these]
;      factor = rebin(reform(factor,1,nband,ngal),naper,nband,ngal)
;
;      apmaggies = float(cat.decam_apflux[*,these]*factor)
;      apivarmaggies = float(cat.decam_apflux_ivar[*,these]/factor^2)
    endif
    
return   
end
