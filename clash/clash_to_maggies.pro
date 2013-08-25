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
       if keyword_set(useirac) then message, 'Code me'
    endif else begin
       tags = filt+'_flux'
       errtags = filt+'_fluxerr'
       if keyword_set(useirac) then begin
          isirac = where(strmatch(filterlist,'*irac*',/fold))
          tags[isirac] = 

       endif
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
;         lim = where((clash.(ftag) lt 90.0) and (maggies[ib,*] gt 0.0) and $
;           (maggies[ib,*]*sqrt(ivar[ib,*]) lt 5.0),nlim)
;         lim = where((maggies[ib,*] gt 0.0) and (maggies[ib,*]*sqrt(ivar[ib,*]) lt 2.0) and $
;           strmatch(filterlist[ib],'*irac*') eq 0,nlim)
          nlim = 0
          if (nlim ne 0L) then begin
             ivar[ib,lim] = 1.0/maggies[ib,lim]^2.0
             maggies[ib,lim] = 0.0
          endif
       endif else begin
          if strmatch(filterlist[ib],'*irac*',/fold) then begin
; aperture fluxes are at the following *diameters*:
;   print, [2.0,5.0,10.0,16.67,20,40]*0.6 ; [pixels]
;     =[1.2,3.0,6.0,10.0,12.0,24.0]        ; [arcsec]

; we want to use the 3.0" diameter (=1.5 arcsec radius)
; aperture fluxes; for comparison SWIRE and SCOSMOS recommended the
; 1.9 arcsec radius fluxes, but our CLASH fields are typically more
; crowded; interpolate the aperture corrections published here:
; http://irsa.ipac.caltech.edu/data/COSMOS/tables/scosmos/scosmos_irac_200706_colDescriptions.html
; to the correct aperture
             print, interpol([0.610,0.765,0.900,0.950],[1.4,1.9,2.9,4.1],1.5) ; ch1
             print, interpol([0.590,0.740,0.900,0.940],[1.4,1.9,2.9,4.1],1.5) ; ch2
             
             apercor =  [0.736, 0.716, 0.606, 0.543]  ;; from SWIRE DOCS
             apercor *=  [1.021, 1.012, 1.022, 1.014] ;; from SPITZER IRAC handbook
 
             fact = 1D
             2.,5.,10.0,16.67,20,40        # diameter in pixels

; apply an aperture correction              
             
             
             
             
             
stop             
          endif else begin
             fact = 10D^(-0.4D*zpt[ib])
             good = where(clash.(utag) gt 0.0,ngood)
             if (ngood ne 0L) then begin
                maggies[ib,good] = clash[good].(ftag)*fact
                ivar[ib,good] = 1D/(clash[good].(utag)*fact)^2.0
             endif

;             lim = where((maggies[ib,*] gt 0.0) and (maggies[ib,*]*sqrt(ivar[ib,*]) lt 5.0),nlim)
;;            lim = where((maggies[ib,*] gt 0.0) and (maggies[ib,*]*sqrt(ivar[ib,*]) lt 2.0) and $
;;              strmatch(filterlist[ib],'*irac*') eq 0,nlim)
;             if (nlim ne 0L) then begin
;                ivar[ib,lim] = 1.0/(2.0*maggies[ib,lim])^2.0
;                maggies[ib,lim] = 0.0
;             endif 
          endelse
       endelse
    endfor 
    
;; apply a minimum photometric error
;    if (keyword_set(nominerror) eq 0) then begin
;       minerr = replicate(0.02,nbands)
;       isirac = where(strmatch(filterlist,'*irac*',/fold))
;;      if isirac[0] ne -1 then minerr[isirac] = 0.1 ; note!
;       k_minerror, maggies, ivar, minerr
;    endif

return   
end
