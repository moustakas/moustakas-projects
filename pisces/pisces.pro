pro pisces

    path = pisces_path()

    pushd, path

; read the mask [1 is good, 0 is bad]
    
    mask = byte(readfits('etc/pisces_mask.fits',/silent))
    goodpix = where(mask eq 1B,ngood,comp=badpix,ncomp=nbad)
    
    im = readfits('qj35a000f.fits',h,/silent)
    im2 = readfits('qj35a001f.fits',h2,/silent)

    info = sstretch(im[goodpix],npix=1000,/silent)
    imin = info.min
    imax = imin + 10*info.sigma
    
; compute the sky mode

    sky, im[goodpix], skymode, skysig
    im[badpix] = skymode

    sky, im2[goodpix], skymode2, skysig
    im2[badpix] = skymode2

    atv, im, min=imin, max=imax

    info2 = sstretch(im2[goodpix],npix=1000,/silent)
    imin2 = info2.min
    imax2 = imin2 + 10*info2.sigma
    
    atv, im2, min=imin2, max=imax2

; parallactic angle

    pa1 = sxpar(h,'PA')
    pa2 = sxpar(h2,'PA')
    angle = pa1-pa2
    
    result = rot(im2,-angle)

    find, im, x, y, flux, sharp, round

    


    popd

stop    
    
return
end
