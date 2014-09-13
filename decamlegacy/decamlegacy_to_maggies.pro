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

    filterlist = wise_filterlist(/short)
    weff = k_lambda_eff(filterlist=filterlist)
    nbands = n_elements(filterlist)

; convert the SDSS photometry to maggies
    sdsscat = im_struct_trimtags(cat,select='sdss_'+['modelflux'
    
    
    
    tags = ['sdss_'+['u','g','i','z'],'decam_'+['g','r','z']]+'_nanomaggies'
    ivartags = tags+'_ivar'

    kl = ext_ccm(weff)*3.1
;   kl = k_lambda(weff,/ccm,/silent)
    
; correct for dust
    glactc, cat.(tag_indx(cat,'dec')), cat.(tag_indx(cat,'ra')), $
      2000.0, gl, gb, 1, /deg
    ebv = dust_getval(gl,gb,/interp,/noloop)
    
; conversion factor from nanomaggies to maggies
    fact = 1D-9

; construct maggies and ivarmaggies in each band       
    maggies = dblarr(nbands,ngal)
    ivarmaggies = dblarr(nbands,ngal)
    for ib = 0, nbands-1 do begin
       ftag = tag_indx(cat[0],tags[ib])
       itag = tag_indx(cat[0],ivartags[ib])
       maggies[ib,*] = cat[*].(ftag)*fact*10D^(0.4*(kl[ib]*ebv))
       ivarmaggies[ib,*] = cat[*].(itag)/(fact*10D^(0.4*(kl[ib]*ebv)))^2D
    endfor

return   
end
