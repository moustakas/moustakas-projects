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
  nominerror=nominerror, useirac=useirac, usemag=usemag

    ngal = n_elements(clash)    
    if (ngal le 0L) then begin
       doc_library, 'clash_to_maggies'
       return
    endif

    filterlist = clash_filterlist(short_filter=filt,useirac=useirac,zpt=zpt)
    nbands = n_elements(filterlist)

    if keyword_set(usemag) then begin
       tags = filt+'_mag'
       errtags = filt+'_magerr'
    endif else begin
       tags = filt+'_flux'
       errtags = filt+'_fluxerr'
    endelse

; construct maggies and ivarmaggies in each band       
    maggies = dblarr(nbands,ngal)
    ivar = dblarr(nbands,ngal)
    for ib = 0, nbands-1 do begin
       ftag = tag_indx(clash[0],tags[ib])
       utag = tag_indx(clash[0],errtags[ib])

; ## Both fluxes and magnitudes have been corrected for:
; ##  - galactic extinction: E(B-V) = 0.03118
; ##  - finite apertures (from encircled energy tables)
; ## mag, magerr =  99, 1-sigma limit: non-detection (flux < 0)
; ## mag, magerr = -99, 0: unobserved (outside FOV, in chip gap, etc.)

       if keyword_set(usemag) then begin
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
             mag = clash[good].(ftag) ; + vega2ab[ib] - kl[ib]*ebv[good] + zpoffset[ib] 
             maggies[ib,good] = 10.0^(-0.4*mag)
             ivar[ib,good] = 1.0/(0.4*alog(10.0)*(maggies[ib,good]*magerr))^2
          endif
          check = where(finite(maggies[ib,*]) eq 0 or finite(ivar[ib,*]) eq 0)
          if check[0] ne -1 then stop

; jm11nov13ucsd - replace <5-sigma photometry with 2-sigma limits 
          lim = where((clash.(ftag) lt 90.0) and (maggies[ib,*] gt 0.0) and (maggies[ib,*]*sqrt(ivar[ib,*]) lt 5.0),nlim)
;         lim = where((maggies[ib,*] gt 0.0) and (maggies[ib,*]*sqrt(ivar[ib,*]) lt 2.0) and $
;           strmatch(filterlist[ib],'*irac*') eq 0,nlim)
          if (nlim ne 0L) then begin
             ivar[ib,lim] = 1.0/(2.0*maggies[ib,lim])^2.0
             maggies[ib,lim] = 0.0
          endif
;if ib eq 3 then stop
       endif else begin
          if strmatch(filterlist[ib],'*irac*',/fold) then fact = 1D else fact = 10D^(-0.4D*zpt[ib])

          good = where(clash.(utag) gt 0.0,ngood)
          if (ngood ne 0L) then begin
             maggies[ib,good] = clash[good].(ftag)*fact
             ivar[ib,good] = 1D/(clash[good].(utag)*fact)^2.0
          endif

          lim = where((maggies[ib,*] gt 0.0) and (maggies[ib,*]*sqrt(ivar[ib,*]) lt 5.0),nlim)
;         lim = where((maggies[ib,*] gt 0.0) and (maggies[ib,*]*sqrt(ivar[ib,*]) lt 2.0) and $
;           strmatch(filterlist[ib],'*irac*') eq 0,nlim)
          if (nlim ne 0L) then begin
             ivar[ib,lim] = 1.0/(2.0*maggies[ib,lim])^2.0
             maggies[ib,lim] = 0.0
          endif 
       endelse
    endfor 
    
; apply a minimum photometric error
    if (keyword_set(nominerror) eq 0) then begin
       minerr = replicate(0.02,nbands)
       isirac = where(strmatch(filterlist,'*irac*',/fold))
;      if isirac[0] ne -1 then minerr[isirac] = 0.1 ; note!
       k_minerror, maggies, ivar, minerr
    endif

return   
end
