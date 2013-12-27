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

pro build_deep2_wise, query=query, parse=parse
; jm13dec23siena - match the DEEP2/DR4 catalog to IRSA/WISE:
; http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query

; verify the formating of the table:
; http://irsa.ipac.caltech.edu/applications/TblCheck/

; search radius: 2"
    dr = 'dr4'
    path = deep2_path(/catalogs)
    cat = read_deep2_zcat(/all)
    ngal = n_elements(cat)

    if keyword_set(query) then begin
       wise_writeit, cat, '~/tmp/deep2_'+dr+'_wise_input.tbl'
       return
    endif

;; split into NCHUNK files to avoid IRSA's file size limit
;    nchunk = 3
;    chunksize = ceil(ngal/nchunk)
;    
;    if keyword_set(write) then for ii = 0, nchunk-1 do wise_writeit, $
;      cat[ii*chunksize:((ii+1)*chunksize-1)<(ngal-1)], $
;         wisepath+'deep2_'+dr+'_wise_input'+strtrim(ii+1,2)+'.tbl'

; code to parse the output and write the DEEP2 object_wise.fits file
    if keyword_set(parse) then begin
       tbl = im_read_tbl('~/tmp/deep2_'+dr+'_wise_output.tbl')

; remove duplicates the dumbest way possible
       uu = uniq(tbl.cntr_01,sort(tbl.cntr_01))
       tbl = tbl[uu]

       out = im_empty_structure(tbl[0],ncopies=ngal)
       out[tbl.cntr_01-1] = tbl
       
;      for ii = 0, nchunk-1 do begin
;         tbl = im_read_tbl(wisepath+'deep2_'+ver+'_wise_output'+strtrim(ii+1,2)+'.tbl')
;         if ii eq 0 then out = im_empty_structure(tbl[0],ncopies=ngal,empty_value=-999.0)
;
;         thisout = out[ii*chunksize:((ii+1)*chunksize-1)<(ngal-1)]
;         thisout[tbl.cntr_01-1] = tbl
;
;         out[ii*chunksize:((ii+1)*chunksize-1)<(ngal-1)] = thisout
;      endfor
       im_mwrfits, out, path+'deep2_'+dr+'_wise.fits', /clobber
    endif

return
end    
