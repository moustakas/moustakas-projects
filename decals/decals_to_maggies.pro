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

pro decals_sdss_to_maggies, cat, maggies, ivarmaggies

    flux = 'cmodel'
    fluxname = 'sdss_'+flux+'flux'
    ivarname = 'sdss_'+flux+'flux_ivar'
    iflux = tag_indx(cat[0],fluxname)
    iivar = tag_indx(cat[0],ivarname)
    sdss_nmgy = cat.(iflux)
    sdss_ivar = cat.(iivar)

    iextinction = tag_indx(cat[0],'sdss_extinction')
    extinction = cat.(iextinction)

    sdss_nmgy = sdss_nmgy*10D^(0.4*extinction)
    sdss_ivar = sdss_ivar*10D^(-0.8*extinction)
    maggies = sdss_nmgy*1D-9
    ivarmaggies = sdss_ivar*1D18
    k_abfix, maggies, ivarmaggies
    k_minerror, maggies, ivarmaggies
    maggies = float(maggies)
    ivarmaggies = float(ivarmaggies)

return
end

pro decals_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist

    ngal = n_elements(cat)
    if (ngal eq 0L) then begin
       doc_library, 'decals_to_maggies'
       return
    endif

    decals_sdss_to_maggies, cat, smaggies, sivarmaggies

; reddening factors from Schlafly & Finkbeiner (2014, private
; communication) for the DECam filter curves
;   red_fac = [3.214,2.165,1.211] ; =grz
;   red_fac = [3.995,3.214,2.165,1.592,1.211,1.064] ; =ugrizY

;   euler, cat.ra, cat.dec, ll, bb, 1
;   extinction = red_fac # dust_getval(ll,bb,/interp,/noloop)

    these = [1,2,4]
    dmaggies = float(cat.decam_flux[these]*10D^(0.4*cat.decam_extinction[these])*1D-9)
    divarmaggies = float(cat.decam_flux_ivar[these]*10D^(-0.8*cat.decam_extinction[these])*1D18)

    maggies = [dmaggies,smaggies]
    ivarmaggies = [divarmaggies,sivarmaggies]
    filterlist = [decals_filterlist(),sdss_filterlist()]
;   filterlist = decals_filterlist()
    
return   
end
