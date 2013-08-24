pro wise_writeit, cat, outfile
    ngal = n_elements(cat)
    splog, 'Writing '+outfile
    openw, lun, outfile, /get_lun
    printf, lun, "\EQUINOX = 'J2000.0'"
    printf, lun, '| ra            | dec           | cntr             |'
    printf, lun, '|double         |double         |long              |'
    for kk = 0L, ngal-1L do printf, lun, cat[kk].ra, cat[kk].dec, $
      kk+1, format='(E14.6,3x,E14.6,3x,I14)'
    free_lun, lun
return
end

pro build_sdss_wise, query=query, parse=parse, dr72=dr72
; jm13aug15siena - match the SDSS catalog to IRSA/WISE:
; http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-scan?mission=irsa&submit=Select&projshort=WISE

; verify the formating of the table:
; http://irsa.ipac.caltech.edu/applications/TblCheck/

; search radius: 2"
    if keyword_set(dr72) then begin
       dr = 'dr72'
       path = getenv('VAGC_REDUX')+'/'
       outfile = path+'object_wise.fits'
       catfile = path+'object_sdss_imaging.fits.gz'
    endif else begin
       dr = 'dr9'
       path = getenv('IM_DATA_DIR')+'/sdss/'+dr+'/'
       outfile = path+'sdss_'+dr+'_wise.fits'
       catfile = path+'photoPosPlate-dr9.fits.gz'
    endelse
    ngal = sxpar(headfits(catfile,ext=1),'NAXIS2')

    if keyword_set(query) then begin
       cat = mrdfits(catfile,1,columns=['ra','dec'])             
       wise_writeit, cat, '~/tmp/sdss_'+dr+'_wise_input.tbl'
       return
    endif

;; split into NCHUNK files to avoid IRSA's file size limit
;    nchunk = 3
;    chunksize = ceil(ngal/nchunk)
;    
;    if keyword_set(write) then for ii = 0, nchunk-1 do wise_writeit, $
;      cat[ii*chunksize:((ii+1)*chunksize-1)<(ngal-1)], $
;         wisepath+'sdss_'+ver+'_wise_input'+strtrim(ii+1,2)+'.tbl'

; code to parse the output and write the SDSS object_wise.fits file
    if keyword_set(parse) then begin
       tbl = im_read_tbl('~/tmp/sdss_'+dr+'_wise_output.tbl')

; remove duplicates the dumbest way possible
       uu = uniq(tbl.cntr_01,sort(tbl.cntr_01))
       tbl = tbl[uu]

       out = im_empty_structure(tbl[0],ncopies=ngal) ; ivar must be zero for no-match
       out[tbl.cntr_01-1] = tbl
       
;      for ii = 0, nchunk-1 do begin
;         tbl = im_read_tbl(wisepath+'sdss_'+ver+'_wise_output'+strtrim(ii+1,2)+'.tbl')
;         if ii eq 0 then out = im_empty_structure(tbl[0],ncopies=ngal,empty_value=-999.0)
;
;         thisout = out[ii*chunksize:((ii+1)*chunksize-1)<(ngal-1)]
;         thisout[tbl.cntr_01-1] = tbl
;
;         out[ii*chunksize:((ii+1)*chunksize-1)<(ngal-1)] = thisout
;      endfor
       im_mwrfits, out, outfile, /clobber
    endif

return
end    
