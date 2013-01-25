pro wise_writeit, cat, outfile
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

pro write_nsa_wise_input, write=write, parse=parse
; jm13jan21siena - match the NSA to IRSA/WISE:
; http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-dd

; verify the formating of the table:
; http://irsa.ipac.caltech.edu/applications/TblCheck/

    vv = 'v0_1_2'
    
; search radius: 3"
    outpath = getenv('IM_ARCHIVE_DIR')+'/nsa/'
    wisepath = outpath+'wise/'
    photofile = outpath+'nsa_'+vv+'.fits.gz'
    
    if keyword_set(write) then begin
       cat = mrdfits(photofile,1,columns=['ra','dec'],1)
       wise_writeit, cat, wisepath+'nsa_'+vv+'_wise_input.tbl'
    endif

; code to parse the output and write the SDSS object_wise.fits file
    if keyword_set(parse) then begin
       tbl = im_read_tbl(wisepath+'nsa_'+vv+'_wise_output.tbl')
       indx = tbl.cntr_01-1     ; index number in CAT

       fits_open, photofile, fcb
       ngal = fcb.axis[1,1]
       close, /all
       out = im_empty_structure(tbl[0],ncopies=ngal)
;      out.ra = cat.ra & out.dec = cat.dec ; copy for blank rows
;      out[indx] = im_struct_assign(tbl,out[indx],/nozero)
       out[indx] = tbl

       im_mwrfits, out, outpath+'nsa_'+vv+'_wise.fits', /gzip, /clobber
    endif

stop    
    
return
end    
