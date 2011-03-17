pro write_ages_galex_casjobs_input
; jm10apr30ucsd - write the coordinate file that needs to be uploaded
;   to the MAST/casjobs database; see
;   /Users/ioannis/home/research/projects/ages/catalogs/galex/README
;   for more details
; jm10jul22ucsd - code rerun for GR6

    ages = mrdfits(ages_path(/catalogs)+'catalog.cat.noguidestars.fits.gz',1)
    ages.ra = ages.ra*15.0D
    ngal = n_elements(ages)

    out = struct_addtags(replicate({ages_id: 0L},ngal),$
      struct_trimtags(ages,select=['ra','dec']))
    out.ages_id = lindgen(ngal)
    
    outfile = ages_path(/mycatalogs)+'galex/ages_galex_casjobs_input_gr6.dat'
    openw, lun, outfile, /get_lun
    printf, lun, '# ages_id ra dec'
    struct_print, out, lun=lun, ddigit=12, /no_head
    free_lun, lun

return
end
    
