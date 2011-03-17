function bootes_zpoffset, nozpoffset=nozpoffset
; jm10jan05ucsd - zeropoint offsets for BwRI were derived in
; BOOTES_ZEROPOINTS; assume the JHKs zeropoint corrections are
; negligible

; for IRAC the code below deals with *aperture corrections* (not
; zeropoint corrections, but they're effectively the same);
; because of the wide IRAC PSF, it is generally recommended to
; use a 1.9" radius (3.8" diameter) aperture magnitude and then to
; apply the relevant *point-source* correction to get the total
; magnitude; see, for example:
;   http://ssc.spitzer.caltech.edu/irac/calib
;   http://data.spitzer.caltech.edu/popular/scosmos/20070525_enhanced/docs/irac_info.txt
;   http://swire.ipac.caltech.edu/swire/astronomers/publications/SWIRE2_doc_083105.pdf
;   Surace et al. 2004 and Ilbert et al. 2009
;
; the table of aperture corrections from the SCOSMOS README is:
;
;           1.4''   1.9''   2.9''   4.1''  
;  
;   ch1    0.610   0.765   0.900   0.950
;   ch2    0.590   0.740   0.900   0.940
;   ch3    0.490   0.625   0.840   0.940
;   ch4    0.450   0.580   0.730   0.910
; 
; of course, M. Brown measured aperture magnitudes (in unconvolved
; IRAC images) in a 4" diameter aperture, but we're going to ignore
; the potential difference between the 4" and 3.8" aperture magnitudes 

    irac_apercor = +2.5*alog10([0.765,0.740,0.625,0.580])
;   irac_apercor = +2.5*alog10([0.736,0.716,0.606,0.543]) ; SWIRE/DR2 document
;   irac_apercor = -alog10([1.205,1.221,1.363,1.571])     ; IRAC website

; now deal with BwRIJHKs    
    splog, 'Ignoring zeropoint offsets!'
    bwri_zpoffset = [-0.129,-0.066,-0.029]*0.0
    jhks_zpoffset = [0.0,0.0,0.0]
    if keyword_set(nozpoffset) then begin
       bwri_zpoffset = bwri_zpoffset*0.0
       jhks_zpoffset = jhks_zpoffset*0.0
    endif
    zpoffset = [bwri_zpoffset,jhks_zpoffset,irac_apercor]
;   splog, zpoffset
    
return, zpoffset
end