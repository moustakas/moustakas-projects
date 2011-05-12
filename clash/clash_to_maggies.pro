;+
; NAME:
;   CLASH_TO_MAGGIES
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

pro clash_to_maggies, clash, maggies, ivar, filterlist=filterlist, $
  nominerror=nominerror

    ngal = n_elements(clash)    
    if (ngal le 0L) then begin
       doc_library, 'clash_to_maggies'
       return
    endif

    filterlist = clash_filterlist(short_filter=filt)
    nbands = n_elements(filterlist)

;; correct for Galactic extinction    
;    weff = k_lambda_eff(filterlist=filterlist)
;    kl = k_lambda(weff,/odonnell,/silent)
;    glactc, clash.i_alpha_j2000, clash.i_delta_j2000, $
;      2000.0, gl, gb, 1, /deg
;    ebv = dust_getval(gl,gb,/interp,/noloop)

    tags = filt+'_mag'
    errtags = filt+'_magerr'

; construct maggies and ivarmaggies in each band       
    maggies = fltarr(nbands,ngal)
    ivar = fltarr(nbands,ngal)
    for ii = 0, nbands-1 do begin
       ftag = tag_indx(clash[0],tags[ii])
       utag = tag_indx(clash[0],errtags[ii])

; ## Both fluxes and magnitudes have been corrected for:
; ##  - galactic extinction: E(B-V) = 0.03118
; ##  - finite apertures (from encircled energy tables)
; ## mag, magerr =  99, 1-sigma limit: non-detection (flux < 0)
; ## mag, magerr = -99, 0: unobserved (outside FOV, in chip gap, etc.)

       limit = where((clash.(ftag) gt 90.0) and $
         (clash.(utag) gt 0.0) and (clash.(utag) lt 90.0),nlimit)
       if (nlimit ne 0L) then begin
          message, 'Deal with me'
;         maggies[ii,limit] = 0.0
;         ivar[ii,limit]= 1.0/(10.^(-0.4*(clash[limit].(ftag[ib])-extinction[ib,ilimit]+v2ab[ib])))^2.0
       endif
         
       good = where($
         (clash.(ftag) gt 0.0) and (clash.(ftag) lt 90.0) and $
         (clash.(utag) gt 0.0) and (clash.(utag) lt 90.0),ngood)
       if (ngood ne 0L) then begin
          magerr = clash[good].(utag)
          mag = clash[good].(ftag); + vega2ab[ii] - kl[ii]*ebv[good] + zpoffset[ii] 
          maggies[ii,good] = 10.0^(-0.4*mag)
          ivar[ii,good] = 1.0/(0.4*alog(10.0)*(maggies[ii,good]*magerr))^2
       endif
    endfor

; apply a minimum photometric error
    if (keyword_set(nominerror) eq 0) then begin
       minerr = replicate(0.02,nbands)
       k_minerror, maggies, ivar, minerr
    endif

return   
end
