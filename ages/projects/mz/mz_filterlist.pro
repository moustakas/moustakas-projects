function mz_filterlist, bands=bands, vega2ab=vega2ab, zpoffset=zpoffset, $
  minerr=minerr, nozpoffset=nozpoffset, sdss=sdss
; jm10oct18ucsd

    if keyword_set(sdss) then begin ; needed by mz_isedfit
       return, [sdss_filterlist(),twomass_filterlist()]
    endif

; now do ages below    
    filterlist = ['lbc_blue_ufilter','ndwfs_'+['Bw','R','I'],$
      'bok_90prime_z','newfirm_'+['J','H','Ks']]+'.par'
    nfilt = n_elements(filterlist)

; -------------------------    
; bandpass (root tag) names
    if arg_present(bands) then begin
       bands = [$
         'U',$
         'Bw','R','I',$
         'z',$
         'J','H','Ks'];,$
;        'ch1','ch2','ch3','ch4']
    endif
    
; -------------------------    
; Vega to AB conversions; note: LBC/U-band and zBootes magnitudes are
; already AB; for IRAC, use Reach et al. 2005 for the IRAC channels
; (see also http://ssc.spitzer.caltech.edu/irac/calib)
    if arg_present(vega2ab) then begin
       vega2ab = k_vega2ab(filterlist=filterlist,/kurucz,/silent)

       other = where(strmatch(filterlist,'*lbc*',/fold) or $
         strmatch(filterlist,'*90prime_z*',/fold),nother)
       if (nother ne 0) then vega2ab[other] = 0.0

;      irac = where(strmatch(filterlist,'*irac*',/fold),nirac)
;      if (nirac ne 0) then vega2ab[irac] = -2.5*alog10([280.9,179.7,115.0,64.13]*1D-23)-48.6
    endif
       
; -------------------------    
; zeropoint offsets; use the results of AGES_DERIVE_ZPTOFFSETS 
    if arg_present(zpoffset) then begin
       zpoffset = fltarr(nfilt)
       if (keyword_set(nozpoffset) eq 0) then begin
          zfile = ages_path(/mycatalogs)+'zptoffsets/'+$
            'ages_zptoffsets.fits.gz'
          if file_test(zfile) then begin
             splog, 'Reading '+zfile
             zpt = mrdfits(zfile,1)

             match, strtrim(filterlist,2), strtrim(zpt[0].filterlist,2), m1, m2

             if keyword_set(absolute) then $
               zpoffset[m1] = total(zpt.zptoffset[m2],2) else $
                 zpoffset[m1] = total(zpt.relative_zptoffset[m2],2)
;            niceprint, zptoffset, filterlist
          endif else begin
             splog, 'Zeropoint corrections file '+zfile+' not found!'
          endelse
       endif
    endif

; -------------------------    
; minimum error
    if arg_present(minerr) then begin
       minerr = [$
         0.05,$                 ; U
         0.02,0.02,0.02,$       ; BwRI
         0.05,$                 ; z
         0.05,0.02,0.02]        ; JHKs
;        0.1,0.1,0.1,0.1]       ; ch1-4
    endif
    
return, filterlist
end
