;+
; NAME:
;   SDWFS_TO_MAGGIES
;
; PURPOSE:
;   Convert the SDWFS/IRAC photometry to AB maggies. 
;
; INPUTS: 
;   cat - input photometric catalog [NGAL] 
;
; KEYWORD PARAMETERS:
;
; OUTPUTS: 
;   maggies - output maggies [4,NGAL]
;   ivarmaggies - corresponding inverse variance array [4,NGAL]  
;
; COMMENTS:
;   A minimum error of 0.1 mag is applied to every bandpass. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Nov 14, UCSD - written
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

pro sdwfs_to_maggies, cat, maggies, ivarmaggies, $
  filterlist=filterlist

    nobj = n_elements(cat)
    if (nobj eq 0) then begin
       doc_library, 'sdwfs_to_maggies'
       return
    endif

    filterlist = irac_filterlist()
    names = 'phot_'+['ch1','ch2','ch3','ch4']
    nband = n_elements(names)
    minerrors = replicate(0.1,nband)

; convert from Vega magnitudes 
    vega2ab = k_vega2ab(filterlist=filterlist,/kurucz,/silent)

    maggies = fltarr(nband,nobj)
    ivarmaggies = fltarr(nband,nobj)
    for iband = 0L, nband-1 do begin
       ftag = tag_indx(cat[0],names[iband]+'')
       utag = tag_indx(cat[0],names[iband]+'_err')
       good = where((cat.(ftag) gt 0.0) and $
         (cat.(utag) gt 0.0) and (cat.(ftag) lt 90.0) and $
         (cat.(utag) lt 90.0),ngood)
       if (ngood ne 0) then begin
          mag = cat[good].(ftag)+vega2ab[iband]
          magerr = cat[good].(utag)
          maggies[iband,good] = 10.0^(-0.4*mag)
          notzero = where((maggies[iband,good] gt 0.0),nnotzero)
          if (nnotzero ne 0L) then ivarmaggies[iband,good[notzero]] = $
            1.0/(0.4*alog(10.0)*(maggies[iband,good[notzero]]*magerr[notzero]))^2
       endif 
    endfor

    k_minerror, maggies, ivarmaggies, minerrors

return
end
