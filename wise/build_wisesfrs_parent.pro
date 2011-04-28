pro build_wisesfrs_parent, clobber=clobber
; jm11apr19ucsd - build the parent sample for the WISE SFR calibration
; paper 
    
    catpath = wise_path(/catalogs)
    outpath = wise_path(/wisesfrs)
    dimpath = wise_path(/dimage)

    atlas = read_atlas(measure=measure)
    wise = mrdfits(catpath+'atlas_wise_cat.fits.gz',1,row=atlas.nsaid)
    natlas = n_elements(atlas)

; get the stellar masses
    splog, 'temporary hack!'
    lowz = mrdfits(sdss_path(/lowz)+'lowz_catalog.dr6.fits.gz',1)
    
    spherematch, atlas.ra, atlas.dec, lowz.ra, lowz.dec, 1.0/3600.0, m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    atlas = atlas[m1]
    measure = measure[m1]
    wise = wise[m1]
    ised = mrdfits(sdss_path()+'lowz_isedfit/lowz_bc03_chab_charlot_sfhgrid04.fits.gz',1,rows=m2)
    galex = mrdfits(sdss_path(/lowz)+'lowz_galex.dr6.fits.gz',1,rows=m2)

; require at least one solid measurement in either W3 or W4
    wise_to_maggies, wise, maggies, ivarmaggies
    keep = where((maggies[3,*] gt 0),ngal)
;   keep = where((total(maggies[2:3,*] gt 0,1) ge 1),ngal)

    atlas = atlas[keep]
    measure = measure[keep]
    wise = wise[keep]
    ised = ised[keep]
    galex = galex[keep]

; compute K-corrections
    sdss_to_maggies, calib=measure, smaggies, sivarmaggies, flux='petro' ; 'sersic'
    galex_to_maggies, galex, gmaggies, givarmaggies
    filterlist = [galex_filterlist(),sdss_filterlist()]

    splog, 'Temporary hack on negative IVAR!'
    kcorr = wisesfrs_do_kcorrect(atlas.zdist,[gmaggies,smaggies],$
      abs([givarmaggies,sivarmaggies]),filterlist=filterlist)
    
; now read the SDSS photometry and spectroscopy tables; note that
; these catalogs are not line-matched to ATLAS!
    sdss = mrdfits(dimpath+'sdss_atlas.fits',1)
    sdssline = mrdfits(dimpath+'sdssline_atlas.fits',1)
    
    good = where(atlas[keep].isdss gt -1)
    out_sdss = im_empty_structure(sdss[0],ncopies=ngal,empty_value=-999)
    out_sdss[good] = sdss[atlas[keep[good]].isdss]

    out_sdssline = im_empty_structure(sdssline[0],ncopies=ngal,empty_value=-999)
    out_sdssline[good] = sdssline[atlas[keep[good]].isdss]
    
; write out
    im_mwrfits, atlas, outpath+'wisesfrs_atlas.fits', clobber=clobber
    im_mwrfits, measure, outpath+'wisesfrs_atlas_measure.fits', clobber=clobber
    im_mwrfits, wise, outpath+'wisesfrs_wisecat.fits', clobber=clobber
    im_mwrfits, ised, outpath+'wisesfrs_isedfit.fits', clobber=clobber
    im_mwrfits, galex, outpath+'wisesfrs_galex.fits', clobber=clobber
    im_mwrfits, kcorr, outpath+'wisesfrs_kcorr.fits', clobber=clobber

    im_mwrfits, out_sdss, outpath+'wisesfrs_sdss.fits', clobber=clobber
    im_mwrfits, out_sdssline, outpath+'wisesfrs_sdssline.fits', clobber=clobber

return
end
    
