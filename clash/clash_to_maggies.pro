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
  nominerror=nominerror, useirac=useirac

    ngal = n_elements(clash)    
    if (ngal le 0L) then begin
       doc_library, 'clash_to_maggies'
       return
    endif

    filterlist = clash_filterlist(short_filter=filt)
    if keyword_set(useirac) then begin
       filterlist = [filterlist,(irac_filterlist())[0:1]]
       filt = [filt,['ch1','ch2']]
    endif
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
    for ib = 0, nbands-1 do begin
       ftag = tag_indx(clash[0],tags[ib])
       utag = tag_indx(clash[0],errtags[ib])

; ## Both fluxes and magnitudes have been corrected for:
; ##  - galactic extinction: E(B-V) = 0.03118
; ##  - finite apertures (from encircled energy tables)
; ## mag, magerr =  99, 1-sigma limit: non-detection (flux < 0)
; ## mag, magerr = -99, 0: unobserved (outside FOV, in chip gap, etc.)

       limit = where((clash.(ftag) gt 90.0) and $
         (clash.(utag) gt 0.0) and (clash.(utag) lt 90.0),nlimit)
       if (nlimit ne 0L) then begin
          maggies[ib,limit] = 0.0
          ivar[ib,limit]= 1.0/(10.0^(-0.4*clash[limit].(utag)))^2.0
       endif
         
       good = where($
         (clash.(ftag) gt 0.0) and (clash.(ftag) lt 90.0) and $
         (clash.(utag) gt 0.0) and (clash.(utag) lt 90.0),ngood)
       if (ngood ne 0L) then begin
          magerr = clash[good].(utag)
          mag = clash[good].(ftag); + vega2ab[ib] - kl[ib]*ebv[good] + zpoffset[ib] 
          maggies[ib,good] = 10.0^(-0.4*mag)
          ivar[ib,good] = 1.0/(0.4*alog(10.0)*(maggies[ib,good]*magerr))^2
       endif
       check = where(finite(maggies[ib,*]) eq 0 or finite(ivar[ib,*]) eq 0)
       if check[0] ne -1 then stop

; jm11nov13ucsd - replace <5-sigma photometry with limits 
       lim = where((maggies[ib,*] gt 0.0) and (maggies[ib,*]*sqrt(ivar[ib,*]) lt 5.0) and $
         strmatch(filterlist[ib],'*irac*') eq 0,nlim)
;      lim = where((maggies gt 0.0) and (maggies*sqrt(ivar) lt 5.0),nlim)
       if (nlim ne 0L) then begin
          ivar[ib,lim] = 1.0/(5.0*maggies[ib,lim])^2.0
          maggies[ib,lim] = 0.0
       endif 
    endfor
    
; apply a minimum photometric error
    if (keyword_set(nominerror) eq 0) then begin
       minerr = replicate(0.02,nbands)
       k_minerror, maggies, ivar, minerr
    endif

return   
end
