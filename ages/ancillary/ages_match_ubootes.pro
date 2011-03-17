pro ages_match_ubootes, ubootes
; jm10may01ucsd - build the line-matched UBootes catalog

    catpath = ages_path(/mycatalogs)

    ages = mrdfits(ages_path(/catalogs)+'catalog.cat.noguidestars.fits.gz',1)
    ages.ra = ages.ra*15.0D
    ngal = n_elements(ages)

    splog, 'Reading '+catpath+'ubootes/catalog_u2.fits.gz'
    zband = mrdfits(catpath+'ubootes/catalog_u2.fits.gz',1)

    m1 = im_spherematch(ages,zband,ratagname2='ra',$
      dectagname2='dec',radius=1.0,match2=m2)
    splog, 'Found '+string(n_elements(m1),format='(I0)')+$
      ' objects with good z-band photometry'

    ubootes = im_empty_structure(zband[0],empty_value=-999.0,ncopies=ngal)
    ubootes[m1] = zband[m2]

    moretags = replicate({match: 0, object_position: 0L},ngal)
    ubootes = struct_addtags(moretags,temporary(ubootes))
    ubootes[m1].match = 1
    ubootes[m1].object_position = m2
    
; write out    
    outfile = catpath+'ages_ubootes.fits'
    im_mwrfits, ubootes, outfile, /clobber

return
end
