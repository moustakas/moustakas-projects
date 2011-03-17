pro ages_match_wsrt, ubootes
; jm10oct26ucsd - build the line-matched WSRT catalog

    catpath = ages_path(/mycatalogs)
    ages = mrdfits(ages_path(/catalogs)+'catalog.cat.noguidestars.fits.gz',1)
    ages.ra = ages.ra*15.0D
    ngal = n_elements(ages)

    radio = mrdfits(catpath+'wsrt/wsrt_02devries.fits.gz',1)

    m1 = im_spherematch(ages,radio,ratagname2='raj2000',$
      dectagname2='dej2000',radius=2.0,match2=m2,raoff=raoff,$
      decoff=decoff)
    splog, raoff*3600, decoff*3600
    splog, 'Found '+string(n_elements(m1),format='(I0)')+$
      ' objects with a 1.4-GHz radio flux'

    outradio = im_empty_structure(radio[0],empty_value=-999.0,ncopies=ngal)
    outradio[m1] = radio[m2]

    moretags = replicate({match: 0, object_position: 0L},ngal)
    outradio = struct_addtags(moretags,temporary(outradio))
    outradio[m1].match = 1
    outradio[m1].object_position = m2

; write out    
    outfile = catpath+'ages_wsrt_02devries.fits'
    im_mwrfits, outradio, outfile, /clobber

return
end
