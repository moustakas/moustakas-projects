;+
; NAME:
;       UNPACK_AGES_CATALOGS
;
; PURPOSE:
;       Parse the AGES (2.0) catalogs.
;
; TODO:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Jul 30, U of A - written, based on the
;         release1.1 code 
;       jm06feb05uofa - use dedicated code to unpack the ZMERGE
;         catalog 
;-

pro unpack_ages_catalogs
    
    catalogs_path = ages_path(/catalogs)
    outpath = catalogs_path
    
; ---------------------------------------------------------------------------
; now unpack the rest of the catalogs    
; ---------------------------------------------------------------------------
    
    catalogfiles = 'catalog.'+['usno.v1.1','ndwfs.short','usno','wsrt',$
      'mips','ebv','zmerge','spectweight','cat.noguidestars',$
      'codes','sdss_photometry','ndwfsb','ndwfsr',$
      'ndwfsi','ndwfsk','kcorr.v3','Vmax.v3']

    sexfiles = catalogfiles+'.sex'
    fitsfiles = catalogfiles+'.fits'
    sexheaderfiles = 'header.'+catalogfiles

    ncat = n_elements(catalogfiles)

    for i = 0, 0 do begin
;   for i = 0, ncat-1 do begin

; read an ASCII version of the full catalog and the sextractor header 

       splog, 'Reading '+catalogs_path+catalogfiles[i]
       catalog = djs_readlines(catalogs_path+catalogfiles[i])
       sexheader = djs_readlines(catalogs_path+sexheaderfiles[i])
       
; remove the existing catalog header and prepend the sextractor header

       headindx = where(strmatch(catalog,'*#*') eq 1B,nhead,comp=keep)
       basecatalog = catalog[keep]

       sexcatalog = [sexheader,basecatalog]

; temporarily write out the sextractor-compatible catalog, read it
; back in and write out a binary FITS table; finally delete the ASCII
; catalog 

       openw, lun, outpath+sexfiles[i], /get_lun
       for j = 0L, n_elements(sexcatalog)-1L do printf, lun, sexcatalog[j]
       free_lun, lun

       splog, 'Reading the temporary sextractor catalog'
       cat = rsex(outpath+sexfiles[i])

       splog, 'Writing '+catalogs_path+fitsfiles[i]
       mwrfits, cat, outpath+fitsfiles[i], /create
       spawn, ['gzip -f '+outpath+fitsfiles[i]], /sh

;      if i eq 4L then stop

       spawn, ['/bin/rm -f '+outpath+sexfiles[i]], /sh
       print
stop
    endfor
stop
;;; ---------------------------------------------------------------------------
;;; OBSOLETE!
;;; unpack the IRAC catalogs (see CATALOGS_PATH+IRAC_README); assume the
;;; same number of elements for all four channels and take the 6"
;;; diameter aperture magnitude
;;; ---------------------------------------------------------------------------
;;
;;    if (n_elements(ch1) eq 0L) then ch1 = rsex(catalogs_path+'spectra_IRAC_ch1_Idet.phot')
;;    if (n_elements(ch2) eq 0L) then ch2 = rsex(catalogs_path+'spectra_IRAC_ch2_Idet.phot')
;;    if (n_elements(ch3) eq 0L) then ch3 = rsex(catalogs_path+'spectra_IRAC_ch3_Idet.phot')
;;    if (n_elements(ch4) eq 0L) then ch4 = rsex(catalogs_path+'spectra_IRAC_ch4_Idet.phot')
;;    nirac = n_elements(ch1)
;;
;;    irac = replicate({ra: 0.0D, dec: 0.0D, $
;;      ch1: -999.0, ch1_err: -999.0, $
;;      ch2: -999.0, ch2_err: -999.0, $
;;      ch3: -999.0, ch3_err: -999.0, $
;;      ch4: -999.0, ch4_err: -999.0},nirac)
;;    irac.ra  = ch1.alpha_j2000
;;    irac.dec = ch1.delta_j2000
;;
;;    good = where((ch1.mag_aper5 gt 0.0) and (ch1.mag_aper5 lt 90.0) and $
;;      (ch1.magerr_aper5 gt 0.0) and (ch1.magerr_aper5 lt 0.3),ngood)
;;    if (ngood ne 0L) then begin
;;       irac[good].ch1     = ch1[good].mag_aper5
;;       irac[good].ch1_err = ch1[good].magerr_aper5
;;    endif
;;
;;    good = where((ch2.mag_aper5 gt 0.0) and (ch2.mag_aper5 lt 90.0) and $
;;      (ch2.magerr_aper5 gt 0.0) and (ch2.magerr_aper5 lt 0.3),ngood)
;;    if (ngood ne 0L) then begin
;;       irac[good].ch2     = ch2[good].mag_aper5
;;       irac[good].ch2_err = ch2[good].magerr_aper5
;;    endif
;;
;;    good = where((ch3.mag_aper5 gt 0.0) and (ch3.mag_aper5 lt 90.0) and $
;;      (ch3.magerr_aper5 gt 0.0) and (ch3.magerr_aper5 lt 0.3),ngood)
;;    if (ngood ne 0L) then begin
;;       irac[good].ch3     = ch3[good].mag_aper5
;;       irac[good].ch3_err = ch3[good].magerr_aper5
;;    endif
;;
;;    good = where((ch4.mag_aper5 gt 0.0) and (ch4.mag_aper5 lt 90.0) and $
;;      (ch4.magerr_aper5 gt 0.0) and (ch4.magerr_aper5 lt 0.3),ngood)
;;    if (ngood ne 0L) then begin
;;       irac[good].ch4     = ch4[good].mag_aper5
;;       irac[good].ch4_err = ch4[good].magerr_aper5
;;    endif
;;
;;    iracoutfile = outpath+'catalog.irac.fits'
;;    splog, 'Writing '+iracoutfile
;;    mwrfits, irac, iracoutfile, /create
;;    spawn, 'gzip -f '+iracoutfile, /sh
;;
;;; ---------------------------------------------------------------------------
;;; unpack the X-ray catalog, which is special; see
;;; CATALOGS_PATH+xbootes_README 
;;; ---------------------------------------------------------------------------
;;
;;    xrayinfile = catalogs_path+'xbootes_cat_xray_opt_IR_21jun_v1.0.txt'
;;    xrayoutfile = outpath+'catalog.xbootes.fits'
;;    if (n_elements(xray) eq 0L) then begin
;;       splog, 'Reading '+xrayinfile
;;       xray = im_read_fmr(xrayinfile)
;;    endif
;;    
;;    splog, 'Writing '+xrayoutfile
;;    mwrfits, xray[where(xray.optrank eq 1)], xrayoutfile, /create
;;    spawn, 'gzip -f '+xrayoutfile, /sh

return
end
    
