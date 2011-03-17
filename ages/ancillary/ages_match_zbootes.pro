pro ages_match_zbootes, zbootes
; jm10feb09ucsd - build the line-matched zbootes catalog

    ages = mrdfits(ages_path(/catalogs)+'catalog.cat.noguidestars.fits.gz',1)
    ages.ra = ages.ra*15.0D
    ngal = n_elements(ages)

    catpath = ages_path(/mycatalogs)
    splog, 'Reading '+catpath+'zbootes/zbootes-cat.fits.gz'
    zband = mrdfits(catpath+'zbootes/zbootes-cat.fits.gz',1)

    m1 = im_spherematch(ages,zband,ratagname2='alpha_j2000',$
      dectagname2='delta_j2000',radius=1.0,match2=m2)
    splog, 'Found '+string(n_elements(m1),format='(I0)')+$
      ' objects with good z-band photometry'

    zbootes = im_empty_structure(zband[0],empty_value=-999.0,ncopies=ngal)
    zbootes[m1] = zband[m2]

    moretags = replicate({match: 0, object_position: 0L},ngal)
    zbootes = struct_addtags(moretags,temporary(zbootes))
    zbootes[m1].match = 1
    zbootes[m1].object_position = m2
    
; write out    
    outfile = catpath+'ages_zbootes.fits'
    im_mwrfits, zbootes, outfile, /clobber

return
end
