;+
; NAME:
;   SG1120_TO_MAGGIES
;
; PURPOSE:
;   Convert the SG1120 photometry to extinction-corrected AB maggies. 
;
; INPUTS: 
;   cat - input photometric catalog (see SG1120_MERGE_CATALOGS) [NGAL] 
;
; OUTPUTS: 
;   maggies - output maggies [6,NGAL]
;   ivarmaggies - corresponding inverse variance array [6,NGAL] 
;
; COMMENTS:
;   A minimum error of 0.05 mag is applied to every bandpass. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Jun 15, NYU - written
;
; Copyright (C) 2009, John Moustakas
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

pro sg1120_to_maggies, cat, maggies, ivarmaggies, sdss=sdss, $
  filterlist=filterlist, use_zpt=use_zpt

    if (n_elements(cat) eq 0) then begin
       doc_library, 'sg1120_to_maggies'
       return
    endif

; zeropoint corrections: see VIMOS_CALIBRATE (BV) and LDSS3_CALIBRATE
; (gr); the last value is for FLAMINGOS/Ks, and it is actually a
; Vega-->AB correction
;   zpt_sg1120 = [-0.298,-0.125,+0.085,+0.345,+0.235,0.0] ; v2.0 parent catalog
    zpt_sg1120 = [-0.275,-0.116,+0.097,+0.364,+0.243,+1.839]
    zpt_sdss = fltarr(5) ; no corrections necessary
    
    if keyword_set(sdss) then begin
       names = 'phot_'+['b','v','r','gprime','rprime','ks',$
         'sdssu','sdssg','sdssr','sdssi']
       zpt = [zpt_sg1120,zpt_sdss]
    endif else begin
       names = 'phot_'+['b','v','r','gprime','rprime','ks']
       zpt = zpt_sg1120
    endelse

; override the default zeropoint corrections
    nband = n_elements(names)
    if (n_elements(use_zpt) ne 0) then begin
       if (n_elements(use_zpt) ne nband) then $
         message, 'Problem here!'
       zpt = use_zpt
    endif

    filterlist = sg1120_filterlist(sdss=sdss)
    weff = k_lambda_eff(filterlist=filterlist)
    kl = k_lambda(weff,/odonnell,R_V=3.1)
    ebv = cat.ebv_mw
    
; everything is already on the AB magnitude system (see the notes
; above for FLAMINGOS)
    vega2ab = k_vega2ab(filterlist=filterlist,/kurucz,/silent)*0.0
    
    maggies = dblarr(nband,n_elements(cat))
    ivarmaggies = dblarr(nband,n_elements(cat))
    for iband = 0L, nband-1 do begin
       itag_m = tag_indx(cat[0],names[iband]+'')
       itag_msig = tag_indx(cat[0],names[iband]+'_err')
       indx = where((cat.(itag_m) gt 0.0) and (cat.(itag_m) lt 90.0) and $
         (cat.(itag_msig) ge 0.0) and (cat.(itag_msig) lt 90.0),count)
       if (count gt 0) then begin
          mag = cat[indx].(itag_m)-kl[iband]*ebv[indx]+vega2ab[iband]+zpt[iband]
          maggies[iband,indx] = 10.0^(-0.4*mag)
          sig = cat[indx].(itag_msig);> 0.001
          notzero = where((maggies[iband,indx] gt 0.0),nnotzero)
          if (nnotzero ne 0L) then ivarmaggies[iband,indx[notzero]] = $
            1.0/(0.4*alog(10.0)*(maggies[iband,indx[notzero]]*sig))^2
       endif
    endfor

; include a minimum photometric error    
    minerrors = replicate(0.05,nband)
;   minerrors = replicate(0.1,nband)

    k_minerror, maggies, ivarmaggies, minerrors

return
end
