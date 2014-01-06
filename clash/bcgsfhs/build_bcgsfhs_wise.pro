pro wise_writeit, cat, outfile
    ngal = n_elements(cat)
    splog, 'Writing '+outfile
    openw, lun, outfile, /get_lun
    printf, lun, "\EQUINOX = 'J2000.0'"
    printf, lun, '| ra            | dec           | cntr             |'
    printf, lun, '|double         |double         |long              |'
    for kk = 0, ngal-1L do printf, lun, 15D*hms2dec(cat[kk].ra), $
      hms2dec(cat[kk].dec), kk+1, format='(E14.6,3x,E14.6,3x,I14)'
    free_lun, lun
return
end

pro build_bcgsfhs_wise, query=query, parse=parse
; jm14jan01siena - match the BCGSFHS sample to IRSA/WISE:
; http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query

; verify the formating of the table:
; http://irsa.ipac.caltech.edu/applications/TblCheck/

; search radius: 2"
    path = bcgsfhs_path(/ancillary)
    cat = read_bcgsfhs_sample()
    ngal = n_elements(cat)

    if keyword_set(query) then begin
       wise_writeit, cat, path+'bcgsfhs_wise_input.tbl'
       return
    endif

;; split into NCHUNK files to avoid IRSA's file size limit
;    nchunk = 3
;    chunksize = ceil(ngal/nchunk)
;    
;    if keyword_set(write) then for ii = 0, nchunk-1 do wise_writeit, $
;      cat[ii*chunksize:((ii+1)*chunksize-1)<(ngal-1)], $
;         wisepath+'bcgsfhs_'+dr+'_wise_input'+strtrim(ii+1,2)+'.tbl'

; code to parse the output and write the BCGSFHS object_wise.fits file
    if keyword_set(parse) then begin
       tbl = im_read_tbl(path+'bcgsfhs_wise_output.tbl')

; remove duplicates the dumbest way possible
       uu = uniq(tbl.cntr_01,sort(tbl.cntr_01))
       tbl = tbl[uu]

       out = im_empty_structure(tbl[0],ncopies=ngal)
       out[tbl.cntr_01-1] = tbl
       
;      for ii = 0, nchunk-1 do begin
;         tbl = im_read_tbl(wisepath+'bcgsfhs_'+ver+'_wise_output'+strtrim(ii+1,2)+'.tbl')
;         if ii eq 0 then out = im_empty_structure(tbl[0],ncopies=ngal,empty_value=-999.0)
;
;         thisout = out[ii*chunksize:((ii+1)*chunksize-1)<(ngal-1)]
;         thisout[tbl.cntr_01-1] = tbl
;
;         out[ii*chunksize:((ii+1)*chunksize-1)<(ngal-1)] = thisout
;      endfor

       out = struct_addtags(struct_trimtags(cat,select=['cluster','shortname']),out)
       im_mwrfits, out, path+'bcgsfhs_wise.fits', /clobber
    endif

return
end    
