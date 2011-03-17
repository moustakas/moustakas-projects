pro junk

;   gal = 'NGC3627'
;   gal = 'NGC3938'
    gal = 'NGC4254'
    gal = 'NGC4321'
    gal = 'NGC5055'
    gal = 'NGC5194'
;   gal = 'NGC0628'
;   gal = 'NGC2403'
    
    sall = rsex('extranucs_log12oh.dat')
    s1 = sall[speclinefit_locate(sall,gal)]
    hii = s1.galaxy+s1.region & hiira = s1.ra & hiidec = s1.dec
    
    iall = sings_read_info()
    i1 = iall[speclinefit_locate(iall,gal)]
    ra = i1.ra & dec = i1.dec
    
    cosd = cos(im_hms2dec(dec)*!dtor)
    decoffset = (im_hms2dec(hiidec)-im_hms2dec(dec))*3600.0
    raoffset = (im_hms2dec(hiira)-im_hms2dec(ra))*15.0*cosd*3600.0

    rad = im_hiiregion_deproject(i1.inclination,i1.posangle,$
      raoffset,decoffset)
    factor = i1.distance*1E3*!dtor/3600.0
    rad_kpc = rad*factor
    niceprint, hii, hiira, hiidec, raoffset, decoffset, rad, rad_kpc

    print
    print, rad_kpc
    
;   for i = 0L, n_elements(hii)-1L do print, hii[i], hiira[i], hiidec[i], $
;     raoffset[i], decoffset[i], format='(A8,2A12,2I7)'

    stop
    
Return
end
    
