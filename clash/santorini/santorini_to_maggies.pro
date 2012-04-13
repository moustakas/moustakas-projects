;+
; NAME:
;   SANTORINI_TO_MAGGIES
;
; PURPOSE:
;   Convert the CLASH matched-aperture photometry to maggies. 
;
; INPUTS: 
;   clash - compatible CLASH catalog
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   maggies - 
;   ivarmaggies - 
;
; OPTIONAL OUTPUTS:
;   filterlist - 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 May 05, UCSD
;-

pro santorini_to_maggies, clash, maggies, ivar, filterlist=filterlist, $
  nominerror=nominerror, zpt=zpt

    ngal = n_elements(clash)    
    if (ngal le 0L) then begin
       doc_library, 'santorini_to_maggies'
       return
    endif

    filterlist = santorini_filterlist(short_filter=filt)
    nbands = n_elements(filterlist)

    tags = filt+'_flux'
    errtags = filt+'_fluxerr'

; conversion factor to AB mags
; 
    zpt = [$
; WFC3/UVIS
      32.9739,$
      32.9263,$
      32.9746,$
      33.7139,$
; ACS/WFC
      33.8093,$
      34.2620,$
      34.7933,$ ; F555W
      34.7106,$
      34.1231,$
      33.8965,$
      35.4940,$
      33.8511,$
; WFC3/IR
      34.8709,$
      35.2621,$
      34.7311,$
      34.8603,$
      35.1988]
    zpt = zpt - 0.3 ; aperture correction

; add IRAC
    zpt = [zpt,21.56,21.56]
    fact = 10D^(-0.4D*zpt)

; construct maggies and ivarmaggies in each band       
    maggies = dblarr(nbands,ngal)
    ivar = dblarr(nbands,ngal)
    for ib = 0, nbands-1 do begin
       ftag = tag_indx(clash[0],tags[ib])
       utag = tag_indx(clash[0],errtags[ib])

       good = where(clash.(utag) gt 0.0,ngood)
       if (ngood ne 0L) then begin
          maggies[ib,good] = clash[good].(ftag)*fact[ib]
          ivar[ib,good] = 1D/(clash[good].(utag)*fact[ib])^2.0
       endif
    endfor
    
; apply a minimum photometric error
;   if (keyword_set(nominerror) eq 0) then begin
;      minerr = replicate(0.6,nbands)
;      minerr[16:17] = 0.0
;      minerr = replicate(0.02,nbands)
;      isirac = where(strmatch(filterlist,'*irac*',/fold))
;      if isirac[0] ne -1 then minerr[isirac] = 0.1 ; note!
;      k_minerror, maggies, ivar, minerr
;   endif

return   
end
