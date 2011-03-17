pro quick_dss_aperture, galaxy, imsize=imsize, aperture=aperture, $
  scan=scan, posangle=posangle, _extra=extra
; jm05jul08uofa
    
    if (n_elements(galaxy) eq 0L) then return
    if (n_elements(imsize) eq 0L) then imsize = 5.0

    querydss, galaxy, image, h, imsize=imsize, /ned, _extra=extra

    extast, h, astr

    imsize = size(image,/dimension)
    xsize = imsize[0] & xcen = xsize/2.0 & ysize = imsize[1] & ycen = ysize/2.0

    xpixscale = (astr.pltscl*1E-3*astr.xsz)/60 ; [arcmin/pixel]
    ypixscale = (astr.pltscl*1E-3*astr.ysz)/60 ; [arcmin/pixel]

    xaxis = (findgen(xsize)-xcen)*xpixscale ; centered on the image [arcsec]
    yaxis = (findgen(ysize)-ycen)*ypixscale ; centered on the image [arcsec]

    img = imgscl(image,/log)

    im_window, 0, xratio=0.5, /square
    plotimage, img, /normal, position=imposition, margin=0, $
      imgxrange=minmax(xaxis), imgyrange=minmax(yaxis), $
      charthick=postthick, xthick=postthick, ythick=postthick, $
      /noaxes, noerase=noerase, _extra=extra

    if (n_elements(aperture) eq 0L) then aperture = 0.0
    if (n_elements(scan) eq 0L) then scan = 0.0
    if (n_elements(posangle) eq 0L) then posangle = 0.0

    im_oplot_box, aperture/60.0, scan/60.0, posangle, thick=2.0, color=djs_icolor('red')

return
end
