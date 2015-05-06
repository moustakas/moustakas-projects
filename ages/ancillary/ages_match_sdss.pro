pro ages_match_sdss, gather=gather
; jm09apr29nyu - match the AGES and SDSS catalogs

    ages = mrdfits(ages_path(/catalogs)+'catalog.cat.noguidestars.fits.gz',1)
    ages.ra = ages.ra*15.0D
    ngal = n_elements(ages)

    cra = mean(ages.ra)
    cdec = mean(ages.dec)
    dra = (max(ages.ra)-min(ages.ra))
    ddec = max(ages.dec)-min(ages.dec)
;   radius = sqrt(dra^2+ddec^2)
    radius = float(ceil(dra>ddec))*1.1
    print, cra, cdec, radius

; grab stars and galaxies
    star = sdss_sweep_circle(cra,cdec,radius,type='star')
    gal = sdss_sweep_circle(cra,cdec,radius,type='gal')
    phot = [star,gal]

stop    
    
; match    
    spherematch, phot.ra, phot.dec, ages.ra, $
      ages.dec, 1.5/3600.0, m1, m2
    splog, 'N = ', n_elements(m1)
    srt = sort(m2)
    m1 = m1[srt]
    m2 = m2[srt]
    
    sdss = im_empty_structure(phot[0],empty_value=-999.0,ncopies=ngal)
    sdss = struct_addtags(temporary(sdss),replicate({sdss_match: 0},ngal))
    sdss[m2] = phot[m1]
    sdss[m2].sdss_match = 1

; write out    
    outfile = ages_path(/mycatalogs)+'ages.sdss.phot.dr72.fits'
    im_mwrfits, sdss, outfile
    
return
end
    
