pro unpack_bootes, parse_refband=parse_refband
; jm10jan04ucsd - unpack the field-by-field BOOTES catalogs to make a
;   merged, trimmed-down set of catalogs that are easier to manage; do
;   not bother parsing the NDWFS/K-band data

    year = '2010b' ; '2009b'
    iopath = getenv('DATA_DIR')+'/data/bootes/'+year+'/'

    refband = 'I'
    bands = ['J','H','Ks']
;   bands = ['u','Bw','R','z','J','H','Ks','ch1','ch2','ch3','ch4']
;   bands = ['Bw','R','K','J','H','Ks','ch1','ch2','ch3','ch4']
    field1 = ['32','33','34','35']
    field2 = ['33','34','35','36']

    iband_cut = 24.0 ; note!
    reffile = iopath+'bootes_'+refband+'.fits'

    if keyword_set(parse_refband) then begin
; first parse the I-band catalog
       for ii = 0, n_elements(field1)-1 do begin
          file = iopath+'Bootes_'+refband+'_'+year+'_Idet_'+$
            field1[ii]+'_'+field2[ii]+'_phot.fits.gz'
          splog, 'Reading '+file
          refcat1 = mrdfits(file,1)
          ngood = n_elements(refcat1)
          good = lindgen(ngood)
;         good = where((refcat1.mag_auto gt 0.0) and $
;           (refcat1.mag_auto lt 90) and $
;           (refcat1.mag_auto lt iband_cut),ngood)
          splog, 'Good = ', ngood
          refcat1 = struct_trimtags(refcat1[good],$
            except=['*flux_aper*','*fluxerr_aper*']) ; ,'*flags*'])
          moretags = replicate({field: ii, field_object_position: 0L, $
            object_position: 0L},ngood)
          refcat1 = struct_addtags(moretags,temporary(refcat1))
          refcat1.field_object_position = good
          if (ii eq 0) then refcat = temporary(refcat1) else $
            refcat = [temporary(refcat),temporary(refcat1)]
       endfor
       refcat.object_position = lindgen(n_elements(refcat))
; write out       
       im_mwrfits, refcat, reffile, /clobber
       refcat = 0
       return
    endif else begin
       splog, 'Reading '+reffile+'.gz'
       refcat = mrdfits(reffile+'.gz',1)
    endelse

; now loop back through and parse the other bandpasses    
;   for jj = 0, 0 do begin
    for jj = 0, n_elements(bands)-1 do begin
       for ii = 0, n_elements(field1)-1 do begin
          file = iopath+'Bootes_'+bands[jj]+'_'+year+'_Idet_'+$
            field1[ii]+'_'+field2[ii]+'_phot.fits.gz'
          these = refcat[where(refcat.field eq ii)].field_object_position ;
          splog, 'Reading '+file
          cat1 = mrdfits(file,1,rows=these)
          cat1 = struct_trimtags(temporary(cat1),$
            except=['*flux_aper*','*fluxerr_aper*']) ; ,'*flags*'])
          moretags = replicate({field: ii, field_object_position: 0L, $
            object_position: 0L},n_elements(these))
          cat1 = struct_addtags(moretags,temporary(cat1))
          cat1.field_object_position = these
          if (ii eq 0) then cat = temporary(cat1) else $
            cat = [temporary(cat),temporary(cat1)]
       endfor
; write out
       outfile = iopath+'bootes_'+bands[jj]+'.fits'
       cat.object_position = lindgen(n_elements(cat))
       im_mwrfits, cat, outfile, /clobber
       cat = 0
    endfor 

stop    

return    
end
