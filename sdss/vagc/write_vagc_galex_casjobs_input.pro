pro write_vagc_galex_casjobs_input, gr=gr
; jm10apr30ucsd - write the coordinate file that needs to be uploaded
;   to the MAST/casjobs database; see
;   /Users/ioannis/home/research/projects/vagc/catalogs/galex/README
;   for more details
; jm10jul23ucsd - updated to GALEX/GR6
;
; Note: we don't query the full VAGC, just objects with spectroscopic
; redshifts; also note that this routine likely needs to be run on an
; NYU machine because of memory issues

    if (n_elements(gr) eq 0) then gr = 'gr6'
    
    vagcpath = getenv('VAGC_REDUX')+'/'
    imaging1 = hogg_mrdfits(vagcpath+'object_sdss_imaging.fits.gz',1,nrowchunk=20000L)
    spectro1 = hogg_mrdfits(vagcpath+'object_sdss_spectro.fits.gz',1,nrowchunk=20000L)
    ngal = n_elements(imaging1)
    vagc_id = lindgen(ngal)

; just the spectroscopy    
    specz = where(spectro1.sdss_spectro_tag ne -1,nspecz)
    out = struct_addtags(replicate({vagc_id: 0L},nspecz),$
      struct_trimtags(imaging1[specz],select=['ra','dec']))
    out.vagc_id = vagc_id[specz]
    
    outfile = vagcpath+'galex/casjobs_'+gr+'.dat'
    openw, lun, outfile, /get_lun
    printf, lun, '# vagc_id ra dec'
    struct_print, out, lun=lun, ddigit=12, /no_head
    free_lun, lun

; everything!
    out = struct_addtags(replicate({vagc_id: 0L},ngal),$
      struct_trimtags(imaging1,select=['ra','dec']))
    out.vagc_id = vagc_id

    outfile = vagcpath+'galex/casjobs_'+gr+'_big.dat'
    openw, lun, outfile, /get_lun
    printf, lun, '# vagc_id ra dec'
    struct_print, out, lun=lun, ddigit=12, /no_head
    free_lun, lun
    
stop
    
return
end
    
