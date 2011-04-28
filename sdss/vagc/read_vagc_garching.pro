;+
; NAME:
;   READ_VAGC_GARCHING
;
; PURPOSE:
;   Read the line-matched output of VAGC and GARCHING files written
;   out by BUILD_VAGC_GARCHING.
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Apr 30, UCSD
;
; Copyright (C) 2010, John Moustakas
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

function read_vagc_garching, sample=sample, letter=letter, $
  poststr=poststr, postlss=postlss, mpacat=mpacat, $
  mpamassoh=mpamassoh, ispecline=ispecline, vmax_noevol=vmax_noevol, $
  vmax_evol=vmax_evol, silent=silent

    common sdss_vagc_mpa, sdss_sample, sdss_letter, sdss_poststr, $
      sdss_postlss, sdss_mpacat, sdss_mpamassoh, sdss_ispecline, $
      sdss_vmax_noevol, sdss_vmax_evol
    
    if (n_elements(sample) eq 0) then sample = 'dr72'
    if (n_elements(letter) eq 0) then letter = 'bsafe'
    if (n_elements(poststr) eq 0) then poststr = '25'
    suffix = sample+letter+poststr

    if (n_elements(sdss_sample) eq 0L) then sdss_sample = sample
    if (n_elements(sdss_letter) eq 0L) then sdss_letter = letter
    if (n_elements(sdss_poststr) eq 0L) then sdss_poststr = poststr

    vagcpath = getenv('LSS_REDUX')+'/'+sample+'/'+letter+'/'+poststr+'/'

    if keyword_set(postlss) then begin
       thisfile = vagcpath+'post_catalog.'+suffix+'.fits.gz'
       if (size(sdss_postlss,/type) ne 8) or (sample ne sdss_sample) or $
         (letter ne sdss_letter) or (poststr ne sdss_poststr) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          sdss_postlss = hogg_mrdfits(thisfile,1,silent=0,nrow=50000L)
          if (sdss_sample ne sample) then sdss_sample = sample
          if (sdss_letter ne letter) then sdss_letter = letter
          if (sdss_poststr ne poststr) then sdss_poststr = poststr
       endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
       return, sdss_postlss
    endif
    
    if keyword_set(mpacat) then begin
       thisfile = vagcpath+'mpacat.'+suffix+'.fits.gz'
       if (size(sdss_mpacat,/type) ne 8) or (sample ne sdss_sample) or $
         (letter ne sdss_letter) or (poststr ne sdss_poststr) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          sdss_mpacat = hogg_mrdfits(thisfile,1,silent=0,nrow=50000L)
          if (sdss_sample ne sample) then sdss_sample = sample
          if (sdss_letter ne letter) then sdss_letter = letter
          if (sdss_poststr ne poststr) then sdss_poststr = poststr
       endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
       return, sdss_mpacat
    endif
    
    if keyword_set(mpamassoh) then begin
       thisfile = vagcpath+'mpamassoh.'+suffix+'.fits.gz'
       if (size(sdss_mpamassoh,/type) ne 8) or (sample ne sdss_sample) or $
         (letter ne sdss_letter) or (poststr ne sdss_poststr) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          sdss_mpamassoh = hogg_mrdfits(thisfile,1,silent=0,nrow=50000L)
          if (sdss_sample ne sample) then sdss_sample = sample
          if (sdss_letter ne letter) then sdss_letter = letter
          if (sdss_poststr ne poststr) then sdss_poststr = poststr
       endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
       return, sdss_mpamassoh
    endif
    
    if keyword_set(ispecline) then begin
       thisfile = vagcpath+'ispecline.'+suffix+'.fits.gz'
       if (size(sdss_ispecline,/type) ne 8) or (sample ne sdss_sample) or $
         (letter ne sdss_letter) or (poststr ne sdss_poststr) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          sdss_ispecline = hogg_mrdfits(thisfile,1,silent=0,nrow=50000L)
          if (sdss_sample ne sample) then sdss_sample = sample
          if (sdss_letter ne letter) then sdss_letter = letter
          if (sdss_poststr ne poststr) then sdss_poststr = poststr
       endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
       return, sdss_ispecline
    endif
    
    if keyword_set(vmax_noevol) then begin
       thisfile = vagcpath+'vmax/vmax-noevol.'+suffix+'.fits.gz'
       if (size(sdss_vmax_noevol,/type) ne 8) or (sample ne sdss_sample) or $
         (letter ne sdss_letter) or (poststr ne sdss_poststr) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          sdss_vmax_noevol = hogg_mrdfits(thisfile,1,silent=0,nrow=50000L)
          if (sdss_sample ne sample) then sdss_sample = sample
          if (sdss_letter ne letter) then sdss_letter = letter
          if (sdss_poststr ne poststr) then sdss_poststr = poststr
       endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
       return, sdss_vmax_noevol
    endif

    if keyword_set(vmax_evol) then begin
       if (sample eq 'dr72') then evol_params = 'q2.00a-1.00' else message, 'Update me!'
       thisfile = vagcpath+'vmax/vmax-'+evol_params+'.'+suffix+'.fits.gz'
       if (size(sdss_vmax_evol,/type) ne 8) or (sample ne sdss_sample) or $
         (letter ne sdss_letter) or (poststr ne sdss_poststr) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          sdss_vmax_evol = hogg_mrdfits(thisfile,1,silent=0,nrow=50000L)
          if (sdss_sample ne sample) then sdss_sample = sample
          if (sdss_letter ne letter) then sdss_letter = letter
          if (sdss_poststr ne poststr) then sdss_poststr = poststr
       endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
       return, sdss_vmax_evol
    endif

return, 0    
end
    
