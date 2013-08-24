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

pro write_sdss_wise_input, write=write, parse=parse
; jm13jan01siena - match the SDSS/DR9 catalog to IRSA/WISE:
; http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query

; verify the formating of the table:
; http://irsa.ipac.caltech.edu/applications/TblCheck/

; search radius: 3"
    dr = 'dr9'
    outpath = sdss_path()+dr+'/'
    wisepath = outpath+'wise/'
o    photofile = outpath+'photoPosPlate-dr9.fits'
    
    if keyword_set(write) then begin
       cat = hogg_mrdfits(photofile,1,columns=['ra','dec'],nrow=50000L)  ;,range=[0,10000])
       ngal = n_elements(cat)

; split into NCHUNK files to avoid IRSA's file size limit
       nchunk = 4
       chunksize = ceil(ngal/nchunk)
       for ii = 0, nchunk-1 do wise_writeit, cat[ii*chunksize:((ii+1)*chunksize-1)<(ngal-1)], $
         wisepath+'sdss_'+dr+'_wise_input'+strtrim(ii+1,2)+'.tbl'
    endif

; code to parse the output and write the SDSS object_wise.fits file
    if keyword_set(parse) then begin
       tbl = [im_read_tbl(wisepath+'sdss_'+dr+'_wise_output1.tbl'),$
         im_read_tbl(wisepath+'sdss_'+dr+'_wise_output2.tbl'),$
         im_read_tbl(wisepath+'sdss_'+dr+'_wise_output3.tbl'),$
         im_read_tbl(wisepath+'sdss_'+dr+'_wise_output4.tbl')]
       indx = tbl.cntr_01-1     ; index number in CAT

THIS CODE IS WRONG! SEE WRITE_REDMAPPER_WISE_INPUT

       fits_open, photofile, fcb
       ngal = fcb.axis[1,1]
       close, /all
       out = im_empty_structure(tbl[0],ncopies=ngal)
;      out.ra = cat.ra & out.dec = cat.dec ; copy for blank rows
;      out[indx] = im_struct_assign(tbl,out[indx],/nozero)
       out[indx] = tbl

       im_mwrfits, out, outpath+'sdss_'+dr+'_wise.fits', /gzip, /clobber
    endif

stop    
    
return
end    
