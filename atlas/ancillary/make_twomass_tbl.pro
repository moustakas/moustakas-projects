pro make_twomass_tbl
; jm02mar26uofa

    nedpath = atlas_path(/ned)
    outpath = nedpath

    datainfo = mrdfits(nedpath+'atlas2d_ned.fits.gz',1,/silent)
    ngalaxy = n_elements(datainfo)
    
    ra = 15.0*im_hms2dec(datainfo.ra)
    de = im_hms2dec(datainfo.dec)

    openw, lun, outpath+'twomass.tbl', /get_lun
;   printf, lun, '\ twomass.tbl'
;   printf, lun, "\EQUINOX ='J2000'"
    printf, lun, '| ra          | dec         |'
    printf, lun, '| double      | double      |'
    printf, lun, '| deg         | deg         |'
    printf, lun, '| -999        | -999        |'
;   printf, lun, '| galaxy        | ra          | dec         |'
;   printf, lun, '| char          | double      | double      |'

; parse

    ranew = strcompress(ra,/remove)
    w = where(ra lt 10.0,nw)
    if nw ne 0L then ranew[w] = '0'+strcompress(ranew[w],/remove)

    denew = strcompress(abs(de),/remove)
    w = where(abs(de) lt 10.0,nw)
    if nw ne 0L then denew[w] = '0'+strcompress(denew[w],/remove)

;   pos = where(de gt 0.0,npos,comp=neg,ncomp=nneg)
;   if npos ne 0L then denew[pos] = '+'+strcompress(denew[pos],/remove)
;   if nneg ne 0L then denew[neg] = '-'+strcompress(denew[neg],/remove)

    for j = 0L, ngalaxy-1L do printf, lun, ra[j], de[j], format='(1x,F11.6,5x,F10.6)'
;   for j = 0L, ngalaxy-1L do printf, lun, ranew[j], denew[j], format='(1x,A9,5x,A10)'
;   for j = 0L, ngalaxy-1L do printf, lun, datainfo[j].galaxy, ranew[j], denew[j], $
;     format='(1x,A14,2x,A9,5x,A10)'
    printf, lun, ' '
    free_lun, lun

stop
    
return
end
