pro write_sdss_galex_input, gr=gr
; jm13jan19siena - based on WRITE_VAGC_GALEX_CASJOBS_INPUT; modified
; to run on all of SDSS/DR9
;
    if (n_elements(gr) eq 0) then gr = 'gr6'
    
    dr = 'dr9'
    outpath = sdss_path()+dr+'/'
    galexpath = outpath+'galex/'

    cat = hogg_mrdfits(outpath+'photoPosPlate-dr9.fits',1,$
      columns=['ra','dec'],nrow=50000L) ;,range=[0,10000])
    ngal = n_elements(cat)

    out = struct_addtags(replicate({sdss_id: 0L},ngal),cat)
    out.sdss_id = lindgen(ngal)

    outfile = galexpath+'sdss_'+dr+'_galex_'+gr+'_casjobs.dat'
    openw, lun, outfile, /get_lun
    printf, lun, '# sdss_id ra dec'
    struct_print, cat, lun=lun, ddigit=12, /no_head
    free_lun, lun

;; split into NCHUNK files to avoid Casjob's file size limits
;    nchunk = 4
;    chunksize = ceil(ngal/nchunk)
;    for ii = 0, nchunk-1 do begin
;       outfile = galexpath+'sdss_'+dr+'_galex_'+gr+'_input'+strtrim(ii+1,2)+'.dat'
;       openw, lun, outfile, /get_lun
;       printf, lun, '# sdss_id ra dec'
;       struct_print, cat[ii*chunksize:((ii+1)*chunksize-1)<(ngal-1)], $
;         lun=lun, ddigit=12, /no_head
;       free_lun, lun
;    endfor

stop
    
return
end
    
