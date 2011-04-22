pro build_wisesfrs_parent, clobber=clobber
; jm11apr19ucsd - build the parent sample for the WISE SFR calibration
; paper 
    
    catpath = wise_path(/catalogs)
    outpath = wise_path(/wisesfrs)
    dimpath = wise_path(/dimage)

    atlas = read_atlas(measure=measure)
    wise = mrdfits(catpath+'atlas_wise_cat.fits.gz',1,row=atlas.nsaid)
    natlas = n_elements(atlas)

; these catalogs are not line-matched to ATLAS!    
    sdss = mrdfits(dimpath+'sdss_atlas.fits',1)
    sdssline = mrdfits(dimpath+'sdssline_atlas.fits',1)

; require at least one solid measurement in either W3 or W4
    wise_to_maggies, wise, maggies, ivarmaggies
    keep = where((maggies[3,*] gt 0),ngal)
;   keep = where((total(maggies[2:3,*] gt 0,1) ge 1),ngal)

    out_sdss = im_empty_structure(sdss[0],ncopies=ngal,empty_value=-999)
    out_sdssline = im_empty_structure(sdssline[0],ncopies=ngal,empty_value=-999)

    good = where(atlas[keep].isdss gt -1)
    out_sdss[good] = sdss[atlas[keep[good]].isdss]
    out_sdssline[good] = sdssline[atlas[keep[good]].isdss]
    
; write out
    im_mwrfits, atlas[keep], outpath+'wisesfrs_atlas.fits', clobber=clobber
    im_mwrfits, measure[keep], outpath+'wisesfrs_atlas_measure.fits', clobber=clobber
    im_mwrfits, wise[keep], outpath+'wisesfrs_wisecat.fits', clobber=clobber
    im_mwrfits, out_sdss, outpath+'wisesfrs_sdss.fits', clobber=clobber
    im_mwrfits, out_sdssline, outpath+'wisesfrs_sdssline.fits', clobber=clobber

return
end
    
