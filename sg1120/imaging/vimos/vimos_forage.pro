;+
; NAME:
;       VIMOS_FORAGE()
;
; PURPOSE:
;       Retrieve useful header information.
;
; CALLING SEQUENCE:
;       forage = vimos_forage(flist)
;
; INPUTS:
;       flist - list of FITS spectra
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       forage - data structure
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       SXPAR(), HEADFITS(), FILE_SEARCH()
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2007 Jan 15, NYU - based on IFORAGE()
;       jm08jul18nyu - the GAIN is actually given as CONAD
;
; Copyright (C) 2007-2008, John Moustakas
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

function vimos_forage, flist, datapath=datapath

    nflist = n_elements(flist)
    if nflist eq 0L then begin
       print, 'Syntax - forage = vimos_forage(flist)'
       return, -1
    endif

    nflist = n_elements(flist)

    forage = {$
      file:         '', $
      file_reduced: '', $
      naxis:        0L, $
      naxis1:       0L, $
      naxis2:       0L, $
      object:       '', $
      target:       '', $
      filter:       '', $
      quadrant:     '', $
      date:         '', $ ; date observed
      ra:           '', $
      dec:          '', $
      pixscale:    0.0, $
      gain:        1.0, $
      rdnoise:     0.0, $
      seeing:      1.0, $ ; average seeing
      airmass:     0.0, $
      exptime:     0.0, $
      mag_zero:    0.0, $
      magerr_zero: 0.0}

    forage = replicate(forage,nflist)

    tags = tag_names(forage[0])
    ntags = n_elements(tags)
    
    for i = 0L, nflist-1L do begin

       h = headfits(flist[i])

       for j = 0L, ntags-1L do begin
          val = sxpar(h,tags[j],count=count)
          if (count ne 0L) and (strcompress(val,/remove) ne '') then forage[i].(j) = val
          if strmatch(tags[j],'*quadrant*',/fold) then begin 
             match = where(strmatch(h,'*quadrant number*',/fold),nmatch) ; this could match twice
             if (nmatch ne 0L) then begin
                forage[i].(j) = strtrim(strmid(h[match[0]],strpos(h[match[0]],'=')+1,$
                  strpos(h[match[0]],'/')-strpos(h[match[0]],'=')-1),2)
             endif
          endif
          if strmatch(tags[j],'*filter*',/fold) then begin
             match = where(strmatch(h,'*filter? = *',/fold),nmatch) ; this could match twice
             if (nmatch ne 0L) then begin
                forage[i].(j) = strtrim(repstr(strmid(h[match[0]],strpos(h[match[0]],'=')+1,$
                  strpos(h[match[0]],'/')-strpos(h[match[0]],'=')-1),"'",''),2)
             endif
          endif
          if strmatch(tags[j],'*gain*',/fold) then begin
             match = where(strmatch(h,'*conad* = *',/fold),nmatch)
             if (nmatch ne 0L) then begin
                forage[i].(j) = strtrim(repstr(strmid(h[match],strpos(h[match],'=')+1,$
                  strpos(h[match],'/')-strpos(h[match],'=')-1),"'",''),2)
             endif
;            match = where(strmatch(h,'*gain* = *',/fold),nmatch)
;            if (nmatch ne 0L) then begin
;               forage[i].(j) = strtrim(repstr(strmid(h[match],strpos(h[match],'=')+1,$
;                 strpos(h[match],'/')-strpos(h[match],'=')-1),"'",''),2)
;            endif
          endif
          if strmatch(tags[j],'*rdnoise*',/fold) then begin
             match = where(strmatch(h,'*ron* = *',/fold),nmatch)
             if (nmatch ne 0L) then begin
                forage[i].(j) = strtrim(repstr(strmid(h[match],strpos(h[match],'=')+1,$
                  strpos(h[match],'/')-strpos(h[match],'=')-1),"'",''),2)
             endif
          endif
          if strmatch(tags[j],'*mag_zero*',/fold) then begin
             match = where(strmatch(h,'*mag zero* = *',/fold),nmatch)
             if (nmatch ne 0L) then begin
                forage[i].(j) = strtrim(repstr(strmid(h[match],strpos(h[match],'=')+1,$
                  strpos(h[match],'/')-strpos(h[match],'=')-1),"'",''),2)
             endif
          endif
          if strmatch(tags[j],'*magerr_zero*',/fold) then begin
             match = where(strmatch(h,'*magzero rms* = *',/fold),nmatch)
             if (nmatch ne 0L) then begin
                forage[i].(j) = strtrim(repstr(strmid(h[match],strpos(h[match],'=')+1,$
                  strpos(h[match],'/')-strpos(h[match],'=')-1),"'",''),2)
             endif
          endif
          if strmatch(tags[j],'*pixscale*',/fold) then begin
             match = where(strmatch(h,'*pixscale* = *',/fold),nmatch)
             if (nmatch ne 0L) then begin
                forage[i].(j) = strtrim(repstr(strmid(h[match],strpos(h[match],'=')+1,$
                  strpos(h[match],'/')-strpos(h[match],'=')-1),"'",''),2)
             endif
          endif
          if strmatch(tags[j],'*seeing*',/fold) then begin
             match1 = where(strmatch(h,'*FWHM START= *',/fold),nmatch1)
             match2 = where(strmatch(h,'*FWHM END= *',/fold),nmatch2)
             if (nmatch1 ne 0L) and (nmatch2 ne 0L) then begin
                fwhm_start = float(strtrim(repstr(strmid(h[match1],strpos(h[match1],'=')+1,$
                  strpos(h[match1],'/')-strpos(h[match1],'=')-1),"'",''),2))
                fwhm_end = float(strtrim(repstr(strmid(h[match2],strpos(h[match2],'=')+1,$
                  strpos(h[match2],'/')-strpos(h[match2],'=')-1),"'",''),2))
                if (fwhm_start gt 0.0) and (fwhm_end gt 0.0) then $
                  forage[i].(j) = (fwhm_start+fwhm_end)/2.0
                if (fwhm_start gt 0.0) and (fwhm_end le 0.0) then $
                  forage[i].(j) = fwhm_start
                if (fwhm_start le 0.0) and (fwhm_end gt 0.0) then $
                  forage[i].(j) = fwhm_end
             endif
          endif
          if strmatch(tags[j],'*target*',/fold) then begin
             match = where(strmatch(h,'*targ name*',/fold),nmatch)
             if (nmatch ne 0L) then begin
                forage[i].(j) = strtrim(repstr(strmid(h[match],strpos(h[match],'=')+1,$
                  strpos(h[match],'/')-strpos(h[match],'=')-1),"'",''),2)
             endif
          endif
          if strmatch(tags[j],'*date*',/fold) then begin
             date = sxpar(h,'DATE-OBS',count=datecount)
             if (datecount eq 1L) then forage[i].(j) = date
          endif
       endfor

       forage[i].file = file_basename(flist[i])

    endfor

return, forage
end    
