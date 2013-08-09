;+
; NAME:
;   MASKPOPS_TO_MAGGIES
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

pro maskpops_to_maggies, clash, maggies, ivar, filterlist=filterlist, zpt=zpt

    ngal = n_elements(clash)    
    if (ngal le 0L) then begin
       doc_library, 'maskpops_to_maggies'
       return
    endif

    filterlist = maskpops_filterlist(short_filter=filt)
    nbands = n_elements(filterlist)

    tags = filt+'_flux'
    errtags = filt+'_fluxerr'

; conversion factor to AB mags
    zpt = [$
; WFC3/UVIS
      24.097,$
      24.174,$
      24.645,$
      25.371,$
; ACS/WFC
      25.658,$
      26.059,$
      26.505,$
      25.907,$
      25.665,$
      25.943,$
      24.842,$
; WFC3/IR
      26.271,$
      26.825,$
      26.247,$
      26.465,$
      25.956]
;; add IRAC
;   zpt = [zpt,21.56,21.56]
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
    
return   
end
