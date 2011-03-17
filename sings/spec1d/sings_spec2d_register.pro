pro sings_spec2d_register
; jm05apr07uofa
; register the 2D SINGS spectra by cross-correlating with the DSS
; images

    datapath = sings_path(/analysis)
    dsspath = sings_path(/dss)
    spec1dpath = sings_path(/spec1d)
    spec2dpath = sings_path(/spec2d)

    sings = mrdfits(datapath+'sings_data.fits.gz',1,/silent)
    spec1d = mrdfits(datapath+'sings_spec1d_info.fits.gz',1,/silent)
    ngalaxy = n_elements(spec1d)

    galaxy = strtrim(strlowcase(spec1d.galaxy),2)
    dssfits = galaxy+'.fits.gz'

    im_window, 0, /square
    
    for k = 8L, ngalaxy-1L do begin
;   for k = 0L, ngalaxy-1L do begin

       match, strlowcase(galaxy[k]), strlowcase(strtrim(sings.galaxy,2)), indx1, indx2
       singsinfo = sings[indx2]
       
; read the DSS image; rebin to the same pixel number of pixels in each
; direction; also rotate such that N is up (zero rotation angle)
       
       image = readfits(dsspath+dssfits[k],h,/silent)
;      hcongrid, image, h, -1, '', min(size(image,/dim)), $
;        min(size(image,/dim)), /half_half
;      hrot, image, h, rimage, rh, 0.0, -1, -1, 0
       extast, h, astr

; sky-subtract

       mmm, image, skymod
       image = image - skymod

       imsize = size(image,/dimension)
       xsize = imsize[0] & xcen = xsize/2.0
       ysize = imsize[1] & ycen = ysize/2.0

       xpixscale = (astr.pltscl*1E-3*astr.xsz)/60.0 ; [arcmin/pixel]
       ypixscale = (astr.pltscl*1E-3*astr.ysz)/60.0 ; [arcmin/pixel]

       xaxis = (findgen(xsize)-xcen)*xpixscale ; centered on the image [arcsec]
       yaxis = (findgen(ysize)-ycen)*ypixscale ; centered on the image [arcsec]

       img = im_imgscl(image)

;      plotimage, img, /preserve_aspect, position=pos, /normal, imgxrange=minmax(xaxis), $
;        imgyrange=minmax(yaxis), charsize=1.8, charthick=postthick, xthick=postthick, $
;        ythick=postthick, xtitle=textoidl('\Delta\alpha [arcmin]'), /noaxes, $
;        ytitle=textoidl('\Delta\delta [arcmin]') ;, title=nicegalaxy[k]
       
       if spec1d[k].drift56 then begin

; read the spectrum and collapse the dispersion axis

          specdata = rd2dspec(strtrim(spec1d[k].drift56_2dfile,2),datapath=spec2dpath,/silent)

          slitpa = spec1d[k].drift56_pa                  ; [degree]
          slitarcsecaxis = make_arcsecaxis(*specdata.header,pscale=slitpixscale)
          slitarcsecaxis = slitarcsecaxis + abs(min(slitarcsecaxis))
          
          nrows = specdata.naxis2
          rowaxis = lindgen(nrows)

;         slitysize = fix(spec1d[k].drift56_ap/slitpixscale)   ; along the slit [pixel]
;         slitxsize = fix(spec1d[k].drift56_scan/slitpixscale) ; perpendicular to the slit [pixel]

          slitprofile = im_normalize(total(specdata.image,1),/max)
;         slitprofile = im_normalize(total(specdata.image,1),rowaxis,normwave=nrows-15L,binsize=14L)

; match the DSS image to the 2D spectrum; first rotate to the
; appropriate angle

          hrot, image, h, rotimage, roth, slitpa, -1, -1, 0, missing=missing
          extast, roth, rotastr

          getrot, roth, angle, cdelt
          splog, 'Desired rotation angle = '+strtrim(string(slitpa,format='(F12.3)'),2)
          splog, 'Actual  rotation angle = '+strtrim(string(angle,format='(F12.3)'),2)

          newxpixscale = abs(cdelt[0])*3600.0 ; [arcsec/pixel]
          newypixscale = cdelt[1]*3600.0      ; [arcsec/pixel]

          slitxsize = spec1d[k].drift56_scan/newxpixscale ; perpendicular to the slit [pixel]
          slitysize = nrows*slitpixscale/newypixscale     ; along the slit [pixel]
;         slitysize = spec1d[k].drift56_ap/newypixscale   ; along the slit [pixel]

; extract a subimage corresponding to the length of the scan (or the
; slit width, in the case of a nuclear spectrum)
          
          x0 = xcen - slitxsize/2.0 & x1 = xcen + slitxsize/2.0
          y0 = ycen - slitysize/2.0 & y1 = ycen + slitysize/2.0

          hextract, rotimage, roth, subimage, subh, x0, x1, y0, y1, /silent
;         getrot, subh, angle, cdelt

          imsize = size(subimage,/dimension)
          newxsize = imsize[0] & newxcen = newxsize/2.0
          newysize = imsize[1] & newycen = newysize/2.0

          plotimage, im_imgscl(rotimage), /preserve_aspect, position=pos, /normal, $
            imgxrange=xsize*[-0.5,0.5], imgyrange=ysize*[-0.5,0.5], charsize=1.8, $
            charthick=postthick, xthick=postthick, ythick=postthick, /noaxes
          im_oplot_box, slitxsize/2.0, slitysize/2.0, 0.0, color='red'

; collapse the X dimension (perpendicular to the slit); interpolate
; the DSS image onto the pixel size of the 2D spectrum

          dssrowaxis = lindgen(newysize)
          dssarcsecaxis = dssrowaxis*newxpixscale

          dssprofile = interpol(im_normalize(total(subimage,1),/max),dssarcsecaxis,slitarcsecaxis)
;         dssprofile = im_normalize(total(subimage,1),dssrowaxis,normwave=newysize-15L,binsize=14L)

          plot, slitarcsecaxis, slitprofile, ps=10, xsty=3, ysty=3
          djs_oplot, slitarcsecaxis, dssprofile, ps=10, color='red', thick=2
          
stop
          
          
       endif

    endfor

    
return
end
    
