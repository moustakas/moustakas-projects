;+
; NAME:
;       LDSS3_FORAGE()
;
; PURPOSE:
;       Retrieve useful header information.
;
; CALLING SEQUENCE:
;       forage = ldss3_forage(flist)
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
;       J. Moustakas, 2006 November 6, NYU - based on IFORAGE()
;
; Copyright (C) 2006, John Moustakas
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

function ldss3_forage, flist, datapath=datapath

    nflist = n_elements(flist)
    if nflist eq 0L then begin
       print, 'Syntax - forage = ldss3_forage(flist)'
       return, -1
    endif

    nflist = n_elements(flist)

    forage = {$
      file:      '', $
      naxis:     0L, $
      naxis1:    0L, $
      naxis2:    0L, $
      object:    '', $
      filter:    '', $
      date:      '', $ ; date observed
      ra:        '', $
      dec:       '', $
      speed:     '', $
      epoch:    0.0, $
      airmass:  0.0, $
      exptime:  0.0, $  ; exposure time
      egain:    0.0, $
      enoise:   0.0, $
      biassec:   ''}

    forage = replicate(forage,nflist)

    tags = tag_names(forage[0])
    ntags = n_elements(tags)
    
    for i = 0L, nflist-1L do begin

       h = headfits(flist[i])

       for j = 0L, ntags-1L do begin
          if strmatch(tags[j],'*filter*',/fold) then begin
             match = where(strmatch(h,'*filter*',/fold),nmatch) ; sxpar can't handle e.g., r' ("r-prime")
             if (nmatch ne 0L) then begin
                str = strmid(h[match],strpos(h[match],"'"))
                forage[i].(j) = repstr(strmid(str,0,strpos(str,"'",/reverse_search)),"'",'')
             endif
          endif else begin
             val = sxpar(h,tags[j],count=count)
             if (count ne 0L) and (strcompress(val,/remove) ne '') then forage[i].(j) = val
          endelse
       endfor

       forage[i].file = file_basename(flist[i])

       date = sxpar(h,'DATE-OBS',count=ndate)
       if (ndate ne 0L) then forage[i].date = strcompress(date,/remove)

    endfor

return, forage
end    
