pro atlas_display_image, atlas, imagepath=imagepath, postscript=postscript, $
  lcharsize=lcharsize, pcharsize=pcharsize, pspath=pspath, labeltype=labeltype, $
  circle_radius=circle_radius, nobar=nobar, nolabelbar=nolabelbar, with_spectrum=with_spectrum, $
  imposition=imposition, noscanbox=noscanbox, norc3box=norc3box, astr=astr, $
  imageinfo=imageinfo, _extra=extra, postthick=postthick, encapsulated=encapsulated, $
  noblackback=noblackback, arcminlabel=arcminlabel, speclabeltype=speclabeltype, $
  barlabelcolor=barlabelcolor
; jm05may16uofa

    natlas = n_elements(atlas)
    if (natlas eq 0L) then begin
       print, 'Syntax - atlas_display_image, atlas'
       return
    endif

    if (n_elements(imagepath) eq 0L) then imagepath = atlas_path(/dss)
    if (n_elements(lcharsize) eq 0L) then lcharsize = 1.5
    if (n_elements(pcharsize) eq 0L) then pcharsize = 1.0
    if (n_elements(labeltype) eq 0L) then labeltype = 0L

    if (n_elements(pspath) eq 0L) then pspath = cwd()
    
; call this routine recursively

    if (natlas gt 1L) then begin
       for i = 0L, natlas-1L do begin
          atlas_display_image, atlas[i], imagepath=imagepath, postscript=postscript, $
            lcharsize=lcharsize, pcharsize=pcharsize, pspath=pspath, labeltype=labeltype, $
            with_spectrum=with_spectrum, nobar=nobar, nolabelbar=nolabelbar, imposition=imposition, $
            noscanbox=noscanbox, norc3box=norc3box, astr=astr1, imageinfo=imageinfo1, $
            _extra=extra, postthick=postthick, encapsulated=encapsulated, noblackback=noblackback, $
            arcminlabel=arcminlabel, speclabeltype=speclabeltype, barlabelcolor=barlabelcolor
          if (i eq 0L) then begin
             astr = astr1 
             if arg_present(imageinfo) then imageinfo = imageinfo1 
          endif else begin
             astr = [ [astr], [astr1] ]
             if arg_present(imageinfo) then imageinfo = [ [imageinfo], [imageinfo1] ]
          endelse
          if (not keyword_set(postscript)) then cc = get_kbrd(1)
       endfor
       astr = reform(astr)
       if arg_present(imageinfo) then imageinfo = reform(imageinfo)
       return
    endif

    galaxy = strtrim(atlas.galaxy,2)
    nicegalaxy = strtrim(atlas.nice_galaxy,2)
    if (tag_exist(atlas,'ATLAS_ID')) then id = string(atlas.atlas_id,format='(I3.3)') else id = ''
    fitsname = strlowcase(galaxy)+'.fits.gz'

    irlum = strarr(n_elements(galaxy))
    irgood = where(atlas.ir_lum gt -900,ngood)
    if (ngood ne 0L) then irlum[irgood] = string(atlas[irgood].ir_lum,format='(F4.1)')
;   if (ngood ne 0L) then irlum[irgood] = 'log L(IR) = '+string(atlas[irgood].ir_lum,format='(F4.1)')

    if keyword_set(postscript) then begin
       if (n_elements(postthick) eq 0L) then postthick = 8.0
       if keyword_set(encapsulated) then suffix = '.eps' else suffix = '.ps'
       if keyword_set(with_spectrum) then begin
          psname = strlowcase(galaxy)+'_imspec'+suffix
          dfpsplot, pspath+psname, /color, xsize=8.5, ysize=3.5, encapsulated=encapsulated
       endif else begin
          psname = strlowcase(galaxy)+'_image'+suffix
          dfpsplot, pspath+psname, /square, encapsulated=encapsulated, /color
       endelse
    endif else begin
       if (n_elements(postthick) eq 0L) then postthick = 2.0
    endelse

    if keyword_set(with_spectrum) then begin
       pagemaker, nx=2, ny=1, width=[4.9,2.7], height=2.7, xspace=0.0, yspace=0.0, $
         xmargin=[0.8,0.1], ymargin=[0.1,0.7], xpage=8.5, ypage=3.5, position=position, /normal
       imposition = reform(position[*,1])
       spposition = reform(position[*,0])
    endif else begin
       if (n_elements(imposition) eq 0L) then begin
          pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=8.5, height=8.0, $
            xmargin=0.0, ymargin=0.0, xpage=8.5, ypage=8.0, position=position, $
            /normal 
          imposition = position
       endif
    endelse
       
; read and display the DSS image and the integrated spectrum

    image = readfits(imagepath+fitsname,h,/silent)
    gsssextast, h, astr
    
    imsize = size(image,/dimension)
    xsize = imsize[0] & xcen = xsize/2.0 & ysize = imsize[1] & ycen = ysize/2.0

    xpixscale = (astr.pltscl*1E-3*astr.xsz)/60 ; [arcmin/pixel]
    ypixscale = (astr.pltscl*1E-3*astr.ysz)/60 ; [arcmin/pixel]

;   if (tag_exist(astr,'PLTSCL')) then begin
;      xpixscale = (astr.pltscl*1E-3*astr.xsz)/60 ; [arcmin/pixel]
;      ypixscale = (astr.pltscl*1E-3*astr.ysz)/60 ; [arcmin/pixel]
;   endif else begin
;      xpixscale = astr.cdelt[1]*60.0 ; [arcmin/pixel]
;      ypixscale = astr.cdelt[1]*60.0 ; [arcmin/pixel]
;   endelse

    xaxis = (findgen(xsize)-xcen)*xpixscale ; centered on the image [arcsec]
    yaxis = (findgen(ysize)-ycen)*ypixscale ; centered on the image [arcsec]

    imageinfo = {$
      fitsname:  fitsname, $
      xsize:     xsize,    $
      ysize:     ysize,    $
      xcen:      xcen,     $
      ycen:      ycen,     $
      xpixscale: xpixscale,$
      ypixscale: ypixscale,$
      xaxis:     xaxis,    $
      yaxis:     yaxis}

    image = alog10(image>1)
    zz = zscale_range(image,0.1)
    top = 230L
    img = bytscl(image,min=zz[0],max=zz[1],top=top)
    img = bytscl(top-img,min=-40,top=top)
    
;   img = imgscl(image,/log)
;   img = im_imgscl(image,losig=-3.0,hisig=5.0,boxfrac=0.3,/log,minvalue=-10)

    if keyword_set(with_spectrum) and (not keyword_set(noblackback)) then begin
       noerase = 1L
       polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')
    endif

    plotimage, img, /normal, position=imposition, margin=0, $
      imgxrange=minmax(xaxis), imgyrange=minmax(yaxis), $
      charthick=postthick, xthick=postthick, ythick=postthick, $
      /noaxes, noerase=noerase, _extra=extra

    case labeltype of
       -1L: label = ''
       0L: label = nicegalaxy
       1L: label = [nicegalaxy,strtrim(atlas.lit_type,2)]
       2L: label = [nicegalaxy,strtrim(atlas.lit_type,2),irlum]
       3L: label = [nicegalaxy+' ('+id+')',strtrim(atlas.lit_type,2)]
       else: label = nicegalaxy
    endcase

    if (labeltype ge 0L) then begin
       lpos = [imposition[0,0]*1.1,imposition[3,0]*0.95]
       legend, textoidl(label), box=0, charsize=lcharsize, $
         charthick=postthick, textcolor=djs_icolor('black'), $
         /normal, _extra=extra; /left, /top, 
    endif

; visualize the RC3 box             

    if (not keyword_set(norc3box)) then begin
    
;      adxy, astr, 15.0*im_hms2dec(atlas.ra), im_hms2dec(atlas.dec), xrc3, yrc3
       gsssadxy, astr, 15.0*im_hms2dec(atlas.ra), im_hms2dec(atlas.dec), xrc3, yrc3
       xrc3 = (xrc3 - xcen)*xpixscale & yrc3 = (yrc3 - ycen)*ypixscale

;      if (atlas.posangle gt -900) then pa = atlas.posangle else pa = 90.0
       if (atlas.d25_maj gt -900) and (atlas.posangle gt -900) then begin

          pa = atlas.posangle
          
          dmaj = atlas.d25_maj
          if (atlas.d25_min gt -900) then dmin = atlas.d25_min else dmin = dmaj

; plot the RC3 aperture as a box       
          
;         im_oplot_box, dmin, dmaj, pa, /noplot, corners=corners
;     
;         indx = [0,1,2,3,0]
;         djs_oplot, [corners[0,indx],corners[0,indx+1]]+xrc3, $
;           [corners[1,indx],corners[1,indx+1]]+yrc3, line=2, $
;           color=djs_icolor('blue'), thick=postthick

; plot the RC3 aperture as an ellipse

          tvellipse, dmaj/2.0, dmin/2.0, 0.0, 0.0, 90+pa, color=djs_icolor('purple'), $
            line=2, thick=postthick, /data
          
       endif

    endif
       
; visualize the spectroscopic apertures

    if (not keyword_set(noscanbox)) then begin
       
;      adxy, astr, 15.0*im_hms2dec(atlas.ra), im_hms2dec(atlas.dec), xspec, yspec
       gsssadxy, astr, 15.0*im_hms2dec(atlas.ra), im_hms2dec(atlas.dec), xspec, yspec
       xspec = (xspec - xcen)*xpixscale & yspec = (yspec - ycen)*ypixscale

       im_oplot_box, atlas.drift_scan/60.0, atlas.drift_ap/60.0, $
         atlas.drift_posangle, /noplot, corners=corners

       indx = [0,1,2,3,0]
       djs_oplot, [corners[0,indx],corners[0,indx+1]]+xspec, $
         [corners[1,indx],corners[1,indx+1]]+yspec, line=0, $
         color=djs_icolor('red'), thick=postthick

    endif
       
; overlay a bar representing 30"

    if (not keyword_set(nobar)) then begin
       if (n_elements(barlabelcolor) eq 0L) then barlabelcolor = ''
       if keyword_set(arcminlabel) then begin
          xbar = 60.0
          barlabel = '1 arcmin'
       endif else begin
          xbar = 30.0
          barlabel = '30"'
       endelse
       im_oplot_box, 0.0, xbar/60.0, 90.0, xoffset=-xsize*xpixscale*0.8/2.0, $
         yoffset=+ysize*ypixscale*0.7/2.0, line=0, $
         thick=postthick, corners=corners, /noplot
       oplot, [corners[0,1],corners[0,0]], corners[1,0]*[1,1], $
         line=0, thick=postthick+1, color=djs_icolor(barlabelcolor);, color=djs_icolor('')
    endif
    if (not keyword_set(nolabelbar)) then begin
;   if (not keyword_set(nobar)) or (not keyword_set(nolabelbar)) then begin
       xyouts, (corners[0,1]-corners[0,0])/2.0+corners[0,0], corners[1,0]*1.15, $
         barlabel, charsize=1.1, charthick=2.0, align=0.5, /data, $
         color=djs_icolor(barlabelcolor)
    endif
       
; overlay a circle if requested

    if (n_elements(circle_radius) ne 0L) then begin
       tvcircle, circle_radius, 0.0, 0.0, /data, $
         line=0, color=djs_icolor('dark green'), thick=postthick
    endif    

; plot the spectrum     

    if (n_elements(speclabeltype) eq 0L) then speclabeltype = labeltype
    
    if keyword_set(with_spectrum) then begin
       atlas_display_spectrum, atlas, position=spposition, xstyle=3, ystyle=3, $
         pcharsize=pcharsize, /noerase, postscript=postscript, _extra=extra, $
         labeltype=speclabeltype
    endif
    
    if keyword_set(postscript) then begin
       dfpsclose
;      spawn, ['convert '+pspath+psname+' '+pspath+repstr(psname,'.ps','.png')], /sh
    endif

return
end
    
