;+
; NAME:
;   LEGACYSURVEY_TO_MAGGIES
;
; PURPOSE:
;   Convert the Tractor Legacysurvey photometry to AB maggies.  
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
;   jm17sep12siena - updated to the DR4+ data model
;
; Copyright (C) 2014, 2017, John Moustakas
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

pro legacysurvey_sdss_to_maggies, cat, maggies, ivarmaggies, psf=psf

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

pro legacysurvey_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist, $
  sdss=sdss, psf=psf

    ngal = n_elements(cat)
    if (ngal eq 0L) then begin
       doc_library, 'legacysurvey_to_maggies'
       return
    endif

; reddening factors from Schlafly & Finkbeiner (2014, private
; communication) for the DECam filter curves
;   red_fac = [3.214,2.165,1.211] ; =grz
;   red_fac = [3.995,3.214,2.165,1.592,1.211,1.064] ; =ugrizY

;   euler, cat.ra, cat.dec, ll, bb, 1
;   extinction = red_fac # dust_getval(ll,bb,/interp,/noloop)

    factor = 1D-9 / transpose([ [cat.mw_transmission_g], [cat.mw_transmission_r], $
      [cat.mw_transmission_z] ])
    dmaggies = float(transpose([ [cat.flux_g], [cat.flux_r], [cat.flux_z] ]) * factor)
    divarmaggies = float(transpose([ [cat.flux_ivar_g], [cat.flux_ivar_r], $
      [cat.flux_ivar_z] ]) / factor^2)
       
    dfilterlist = legacysurvey_filterlist()

; unWISE
    factor = 1D-9 / transpose([ [cat.mw_transmission_w1], [cat.mw_transmission_w2] ])
    wmaggies = float(transpose([ [cat.flux_w1], [cat.flux_w2] ]) * factor)
    wivarmaggies = float(transpose([ [cat.flux_ivar_w1], [cat.flux_ivar_w2] ]) / factor^2)

    wfilterlist = wise_filterlist(/short)

; add SDSS (not really supported...)
    if keyword_set(sdss) then begin
       legacysurvey_sdss_to_maggies, cat, smaggies, sivarmaggies, psf=psf
       maggies = [dmaggies,wmaggies,smaggies]
       ivarmaggies = [divarmaggies,wivarmaggies,sivarmaggies]
       filterlist = [dfilterlist,wfilterlist,sdss_filterlist()]
    endif else begin
       maggies = [dmaggies,wmaggies]
       ivarmaggies = [divarmaggies,wivarmaggies]
       filterlist = [dfilterlist,wfilterlist]
    endelse

return   
end
