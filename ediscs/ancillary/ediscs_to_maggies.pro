;+
; NAME:
;   EDISCS_TO_MAGGIES
;
; PURPOSE:
;   Convert the EDisCS BVRIJKs photometry to maggies.
;
; INPUTS: 
;   phot - input photometric catalog (see BUILD_EDISCS_PHOTOMETRY)
;     [NGAL]; all magnitudes are already extinction corrected and are
;     on the AB system
;
; OUTPUTS: 
;   maggies - output maggies [6,NGAL]
;   ivarmaggies - corresponding inverse variance array [6,NGAL] 
;
; COMMENTS:
;   A minimum error of 0.02 mag is applied to every bandpass. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Feb 20, NYU - written
;   jm09aug14ucsd - updated to conform to the new
;   BUILD_EDISCS_PHOTOMETRY routine
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

pro ediscs_to_maggies, phot, maggies, ivarmaggies, filterlist=filterlist

    names = 'phot_'+['b','v','r','i','j','k']
    nband = n_elements(names)
    minerrors = replicate(0.02,nband)

    filterlist = ediscs_filterlist()

    maggies = fltarr(nband,n_elements(phot))
    ivarmaggies = fltarr(nband,n_elements(phot))
    for iband = 0L, nband-1 do begin
       ftag = tag_indx(phot[0],names[iband]+'')
       utag = tag_indx(phot[0], names[iband]+'_err')
       good = where((phot.(ftag) gt 0.0) and (phot.(utag) gt 0.0),ngood)
       if (ngood ne 0) then begin
          mag = phot[good].(ftag)
          magerr = phot[good].(utag)
          maggies[iband,good] = 10.0^(-0.4*mag)
          notzero = where((maggies[iband,good] gt 0.0),nnotzero)
          if (nnotzero ne 0L) then ivarmaggies[iband,good[notzero]] = $
            1.0/(0.4*alog(10.0)*(maggies[iband,good[notzero]]*magerr[notzero]))^2
       endif
    endfor

    k_minerror, maggies, ivarmaggies, minerrors

return
end
