;+
; NAME:
;   DECAMLEGACY_TO_MAGGIES
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

pro decamlegacy_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist, $
  ratag=ratag, dectag=dectag, nodust=nodust, ebv=ebv

    ngal = n_elements(cat)
    if (ngal eq 0L) then begin
       doc_library, 'decamlegacy_to_maggies'
       return
    endif

    filterlist = decamlegacy_filterlist()
;   dfilterlist = filterlist[0:2] ; just DECam/grz
    weff = k_lambda_eff(filterlist=filterlist)
    nbands = n_elements(filterlist)

;; convert the SDSS photometry to maggies
;    tags = ['modelflux','modelflux_ivar','extinction']
;    sdsscat = im_struct_trimtags(cat,select='sdss_'+tags,newtags=tags)
;    sdss_to_maggies, smaggies, sivarmaggies, calib=sdsscat, flux='model'

; now deal with the DECam photometry
    decam_flux = cat.decam_flux[[1,2,4]]
    decam_ivar = cat.decam_flux_ivar[[1,2,4]]^2
    nan = where(finite(cat.decam_flux_ivar[[1,2,4]]) eq 0)
    decam_ivar[nan] = 0

    kl = ext_ccm(weff)*3.1
;   kl = k_lambda(weff,/ccm,/silent)
    
; correct for dust
    glactc, cat.(tag_indx(cat,'dec')), cat.(tag_indx(cat,'ra')), $
      2000.0, gl, gb, 1, /deg
    ebv = dust_getval(gl,gb,/interp,/noloop)
    
; conversion factor from nanomaggies to maggies
    fact = 1D-9

; construct maggies and ivarmaggies in each band       
    dmaggies = dblarr(nbands,ngal)
    divarmaggies = dblarr(nbands,ngal)
    for ib = 0, nbands-1 do begin
       dmaggies[ib,*] = decam_flux[ib,*]*fact*10D^(0.4*(kl[ib]*ebv))
       divarmaggies[ib,*] = decam_ivar[ib,*]/(fact*10D^(0.4*(kl[ib]*ebv)))^2D
    endfor

    maggies = [dmaggies,smaggies]
    ivarmaggies = [divarmaggies,sivarmaggies]
    
return   
end
