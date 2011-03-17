pro imsim, mipsband=mipsband
; jm01may13uofa

; simulate MIPS images

    if not keyword_set(mipsband) then mipsband = '70' ; choose a bandpass

    sirtfpath = (sirtf_datapath())[0]
    psf = readfits(sirtfpath+'psfs/'+mipsband+'_model.fits',hpsf)
    cirrus = readfits(sirtfpath+'imsim/cirrus_'+mipsband+'.fits',chead)
    path = sirtfpath+'catalogs/SIMULATIONS/'

    smin = 0.001 & smax = 250 ; minmax fluxes (1 microJy to 0.25 Jy)
    
    scale = 2.03/2.46       ; 2.03"/pixel (cirrus); 2.46"/pixel (70 mu PSF)
    cpix = long(1024*scale) ; rebin the cirrus image to match plate scales
    newc = frebin(cirrus,cpix,cpix,/total)

    image = fltarr(cpix,cpix) ; final image
    nsources = 0L
    for i = 1L, 5L do begin     ; pull out subsections

       cmrestore, path+'imsim_'+mipsband+'_34.6_20.idlsave', imflux, names=['imflux_'+strn(i)]
       
       good = where((imflux gt smin) and (imflux lt smax),ngals)
       if ngals ne 0L then begin

          flux = imflux[good]

          nsources = nsources + ngals
          temp = randomu(seed,ngals,2)
          xygal = temp-temp

          xygal[*,0] = (temp[*,0]/max(temp[*,0]))*cpix
          xygal[*,1] = (temp[*,1]/max(temp[*,1]))*cpix

          im = image-image
          xy2image, xygal[*,0], y=xygal[*,1], im, binsize=1.0, weights=flux, /float ; Jy

          delvarx, temp, xygal, good ; clean up memory
          
          image = image + temporary(im)
          
       endif

    endfor

    print, 'Number of sources = '+strn(nsources)+'.'

    imtemp = image + newc*1E3*0.28 ; add the cirrus (scaled to "low" cirrus - mJy)
    imfinal = convolve(imtemp,psf) ; convolve with the PSF

    window, 0, xs=845, ys=845
    info = sstretch(imfinal,npix=5E4)
;   display, bytscl(imfinal,min=info.min,max=info.max)
    plotimage, bytscl(imfinal,min=info.min,max=info.max)
    
stop

    writefits, sirtfpath+'imsim/imsim_'+mipsband+'.fits', imfinal

stop    

return
end
