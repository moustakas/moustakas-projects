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

function cluster_ebv, cat
; get the cluster E(B-V)
    cluster = ['a209','a383','macs0329','macs0429','macs0744','a611',$
      'macs1115','a1423','macs1206','clj1226','macs1311','rxj1347',$
      'macs1423','rxj1532','macs1720','a2261','macs1931','rxj2129',$
      'ms2137','rxj2248','macs0416','macs0647','macs0717','macs1149',$
      'macs2129']
    sra = ['01:31:52.57','02:48:03.36','03:29:41.68','04:29:36.10',$
      '07:44:52.80','08:00:56.83','11:15:52.05','11:57:17.26','12:06:12.28',$
      '12:26:58.37','13:11:01.67','13:47:30.59','14:23:47.76','15:32:53.78',$
      '17:20:16.95','17:22:27.25','19:31:49.66','21:29:39.94','21:40:15.18',$
      '22:48:44.29','04:16:09.39','06:47:50.03','07:17:31.65','11:49:35.86',$
      '21:29:26.06']
    sdec = ['-13:36:38.8','-03:31:44.7','-02:11:47.7','-02:53:08.0','+39:27:24.4',$
      '+36:03:24.1','+01:29:56.6','+33:36:37.4','-08:48:02.4','+33:32:47.4',$
      '-03:10:39.5','-11:45:10.1','+24:04:40.5','+30:20:58.7','+35:36:23.6',$
      '+32:07:58.6','-26:34:34.0','+00:05:18.8','-23:39:40.7','-44:31:48.4',$
      '-24:04:03.9','+70:14:49.7','+37:45:18.5','+22:23:55.0','-07:41:28.8']
    ra = 15D*hms2dec(sra)
    dec = hms2dec(sdec)
    diff = min(djs_diff_angle(ra,dec,cat[0].ra,cat[0].dec),this)

    glactc, ra[this], dec[this], 2000.0, gl, gb, 1, /deg
    ebv = dust_getval(gl,gb,/interp,/noloop)
return, ebv
end

pro clash_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist, $
  nominerror=nominerror, useirac=useirac, usemag=usemag

    ngal = n_elements(cat)    
    if (ngal le 0L) then begin
       doc_library, 'clash_to_maggies'
       return
    endif

    filterlist = clash_filterlist(short_filter=filt,$
      useirac=useirac,zpt=zpt,weff=weff)
    isirac = where(strmatch(filterlist,'*irac*',/fold),comp=notirac)
    nbands = n_elements(filterlist)
    apcor = fltarr(nbands)+1.0 ; aperture correction

; we need to correct the Spitzer photometry for Galactic extinction,
; but this is what Coe has to say about the HST catalogs:
;
; ## Both fluxes and magnitudes have been corrected for:
; ##  - galactic extinction: E(B-V) = 0.03118
; ##  - finite apertures (from encircled energy tables)
; ## mag, magerr =  99, 1-sigma limit: non-detection (flux < 0)
; ## mag, magerr = -99, 0: unobserved (outside FOV, in chip gap, etc.)
;
; use one extinction value over the whole cluster field 

    ebv = cluster_ebv(cat)
    kl = k_lambda(weff,/odon,/silent)
    kl[notirac] = 0.0 ; don't correct the HST photometry

; the IRAC/Spitzer photometry is on the Vega magnitude system!  the
; SCOSMOS guys
; (http://irsa.ipac.caltech.edu/data/COSMOS/gator_docs/scosmos_irac_colDescriptions.html)
; give the following conversions:
;         ch 1:   K = -2.788
;         ch 2:   K = -3.255
; but I get using Blanton's code:
; print, k_vega2ab(filterlist=irac_filterlist(/warm),/kurucz)
;     2.79981      3.26876
    vega2ab = fltarr(nbands)    ; zero for everything but IRAC
    if keyword_set(useirac) then vega2ab[isirac] = [2.79981,3.26876]

    if keyword_set(usemag) then begin
       tags = filt+'_mag'
       errtags = filt+'_magerr'
       if keyword_set(useirac) then message, 'Code me'
    endif else begin
       tags = filt+'_flux'
       errtags = filt+'_fluxerr'
       if keyword_set(useirac) then begin
          isirac = where(strmatch(filterlist,'*irac*',/fold))
; aperture fluxes are at the following *diameters*:
;   print, [2.0,5.0,10.0,16.67,20,40]*0.6 ; [pixels]
;     =[1.2,3.0,6.0,10.0,12.0,24.0]        ; [arcsec]
          tags[isirac] = ['ch1','ch2']+'_flux_aper1' ; =5 pixel diameter
          errtags[isirac] = ['ch1','ch2']+'_fluxerr_aper1'
; we want to use the 3.0" diameter (=1.5 arcsec radius)
; aperture fluxes; for comparison SWIRE and SCOSMOS recommended the
; 1.9 arcsec radius fluxes, but our CLASH fields are typically more
; crowded; interpolate the aperture corrections published here:
; http://irsa.ipac.caltech.edu/data/COSMOS/tables/scosmos/scosmos_irac_200706_colDescriptions.html
; to the correct aperture
;         print, interpol([0.610,0.765,0.900,0.950],[1.4,1.9,2.9,4.1],1.5), $ ; [ch1]
;           interpol([0.590,0.740,0.900,0.940],[1.4,1.9,2.9,4.1],1.5)         ; [ch2]
;         0.641000     0.620000          
; for comparison the aperure corrections used by SWIRE/SCOSMOS for
; [ch1-4] and a 1.9 arcsec radius aperture were:
;         apercor =  [0.736, 0.716, 0.606, 0.543]  ;; from SWIRE DOCS
;         apercor *=  [1.021, 1.012, 1.022, 1.014] ;; from SPITZER IRAC handbook
          apcor[isirac] = [0.641,0.620]
       endif 
    endelse

; construct maggies and ivarmaggies in each band       
    maggies = dblarr(nbands,ngal)
    ivarmaggies = dblarr(nbands,ngal)
    for ib = 0, nbands-1 do begin
       ftag = tag_indx(cat[0],tags[ib])
       utag = tag_indx(cat[0],errtags[ib])

       if keyword_set(usemag) then begin
          limit = where((cat.(ftag) gt 90.0) and $
            (cat.(utag) gt 0.0) and (cat.(utag) lt 90.0),nlimit)
          if (nlimit ne 0L) then begin
             maggies[ib,limit] = 0.0
             ivarmaggies[ib,limit]= 1.0/(10.0^(-0.4*cat[limit].(utag)))^2.0
          endif
          
          good = where($
            (cat.(ftag) gt 0.0) and (cat.(ftag) lt 90.0) and $
            (cat.(utag) gt 0.0) and (cat.(utag) lt 90.0) and $
            finite(cat.(utag)),ngood)
          if (ngood ne 0L) then begin
             magerr = cat[good].(utag)
             mag = cat[good].(ftag) - kl[ib]*ebv + vega2ab[ib] ; + zpoffset[ib] 
             maggies[ib,good] = 10.0^(-0.4*mag)
             ivarmaggies[ib,good] = 1.0/(0.4*alog(10.0)*(maggies[ib,good]*magerr))^2
          endif

          nlim = 0
          if (nlim ne 0L) then begin
             ivarmaggies[ib,lim] = 1.0/maggies[ib,lim]^2.0
             maggies[ib,lim] = 0.0
          endif
       endif else begin
          fact = 10D^(-0.4D*(zpt[ib]-kl[ib]*ebv+vega2ab[ib]))/apcor[ib]
          good = where(cat.(utag) gt 0.0 and finite(cat.(utag)),ngood)
          if (ngood ne 0L) then begin
             maggies[ib,good] = cat[good].(ftag)*fact
             ivarmaggies[ib,good] = 1D/(cat[good].(utag)*fact)^2.0
          endif
       endelse 

; clean up NAN's; these correspond to masked areas (e.g., chip
; gap, outside the detection image boundaries, etc.) 
       check = where(finite(maggies[ib,*]) eq 0 or finite(ivarmaggies[ib,*]) eq 0)
       if check[0] ne -1 then begin
          maggies[ib,check] = 0.0
          ivarmaggies[ib,check] = 0.0
       endif
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
