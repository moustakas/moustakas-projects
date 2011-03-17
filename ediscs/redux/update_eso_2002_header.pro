function update_eso_2002_header, h, hinfo, bias=bias, arclamp=arclamp, $
  domeflat=domeflat, object=object, left=left, center=center, right=right
; jm04jun21uofa
; keywords that need to be added or modified to the ESO header to be
; compatible with ispec2d: OBJECT, IMAGETYP, RA, DEC, EPOCH, OBSERVAT,
; UT, APERTURE, RDNOISE, GAIN

    good = where(strcompress(h,/remove) ne '',ngood)
    if (ngood ne 0L) then newh = h[good] else newh = h
    
    if keyword_set(bias) then begin
       hobj = 'bias'
       imagetyp = 'zero'
       ra = '00:00:00'
       dec = '00:00:00'
    endif

    if keyword_set(arclamp) then begin
       hobj = 'HeArHgCd'
       if keyword_set(left) then hobj = hobj+' left'
       if keyword_set(center) then hobj = hobj+' center'
       if keyword_set(right) then hobj = hobj+' right'

       imagetyp = 'comp'
       ra = '00:00:00'
       dec = '00:00:00'
    endif

    if keyword_set(domeflat) then begin
       hobj = 'dome flat'
       if keyword_set(left) then hobj = hobj+' left'
       if keyword_set(center) then hobj = hobj+' center'
       if keyword_set(right) then hobj = hobj+' right'

       imagetyp = 'flat'
       ra = '00:00:00'
       dec = '00:00:00'
    endif

    if keyword_set(object) then begin
       hobj = strtrim(hinfo.obs_targ_name,2)
       if keyword_set(left) then hobj = hobj+' left'
       if keyword_set(center) then hobj = hobj+' center'
       if keyword_set(right) then hobj = hobj+' right'

       imagetyp = 'object'
       ra = strjoin(strsplit(dec2hms(hinfo.ra/15.0D),' ',/extract),':')
       dec = strjoin(strsplit(dec2hms(hinfo.dec),' ',/extract),':')
    endif

    sxaddpar, newh, 'OBJECT', hobj, after='DATE-OBS'
    sxaddpar, newh, 'IMAGETYP', imagetyp, after='OBJECT'
    sxaddpar, newh, 'RA', ra
    sxaddpar, newh, 'DEC', dec
    sxaddpar, newh, 'EPOCH', '2000.0', after='DEC'
    sxaddpar, newh, 'OBSERVAT', 'ESO', after='EPOCH'
    sxaddpar, newh, 'UT', strmid(hinfo.date_obs,strpos(hinfo.date_obs,'T')+1,$
      strlen(hinfo.date_obs)), after='EPOCH'
    sxaddpar, newh, 'RDNOISE', 5.41, after='UT'
    sxaddpar, newh, 'GAIN', 0.52, after='RDNOISE'
    sxaddpar, newh, 'APERTURE', 1.0, after='GAIN'
    
return, newh
end
