pro writeit, cat, outfile
    ngal = n_elements(cat)
    
    splog, 'Writing '+outfile
    openw, lun, outfile, /get_lun
    printf, lun, "\EQUINOX = 'J2000.0'"
    printf, lun, '| ra            | dec           | cntr             |'
    printf, lun, '|double         |double         |long              |'
    for kk = 0, ngal-1L do printf, lun, cat[kk].ra, cat[kk].dec, $
      kk+1, format='(E14.6,3x,E14.6,3x,I14)'
    free_lun, lun
    
return
end

pro write_vagc_wise_input
; jm13jan01siena - match the entire VAGC dr7.2 catalog to IRSA/WISE
; upload here:
; http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query

; verify the formating of the table:
; http://irsa.ipac.caltech.edu/applications/TblCheck/

; search radius: 3"    

    outpath = getenv('VAGC_REDUX')+'/'
    wisepath = outpath+'wise/'
    cat = mrdfits(vagc_name('object_sdss_imaging')+'.gz',$
      1,columns=['ra','dec']);,range=[0,10000])
    ngal = n_elements(cat)

; split into NCHUNK files to avoid IRSA's file size limit
    nchunk = 4
    chunksize = ceil(ngal/nchunk)
    for ii = 0, nchunk-1 do writeit, cat[ii*chunksize:((ii+1)*chunksize-1)<(ngal-1)], $
      wisepath+'vagc_dr72_wise_input'+strtrim(ii+1,2)+'.tbl'

stop    
    
; code to parse the output and write the VAGC object_wise.fits file
    tbl = [im_read_tbl(wisepath+'vagc_dr72_wise_output1.tbl'),$
      im_read_tbl(wisepath+'vagc_dr72_wise_output2.tbl'),$
      im_read_tbl(wisepath+'vagc_dr72_wise_output3.tbl'),$
      im_read_tbl(wisepath+'vagc_dr72_wise_output4.tbl')]

    indx = tbl.cntr_01-1        ; index number in CAT
    
    out = im_empty_structure(tbl[0],ncopies=ngal)
;   out.ra = cat.ra & out.dec = cat.dec ; copy for blank rows
;   out[indx] = im_struct_assign(tbl,out[indx],/nozero)
    out[indx] = tbl

    im_mwrfits, out, outpath+'object_wise.fits', /gzip, /clobber
    
return
end    
