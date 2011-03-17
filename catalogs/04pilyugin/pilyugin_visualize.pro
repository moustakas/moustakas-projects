function pilyugin_imgscl, image

    imsize = size(image,/dimension)
    xsize = imsize[0] & xcen = xsize/2.0
    ysize = imsize[1] & ycen = ysize/2.0
    
    lo = -2.0
    hi = 3.0
    
    xbox = fix(xsize*0.20) & ybox = fix(ysize*0.20)
    djs_iterstat, image[xcen-xbox:xcen+xbox,ycen-ybox:ycen+ybox], sigma=rms, mean=mean

    topvalue = !d.table_size-2
    minvalue = -40L
    
    img = imgscl(image,min=(mean+lo*rms)>min(image),max=(mean+hi*rms)<max(image),top=topvalue)
    img = bytscl(topvalue-img,min=minvalue,top=topvalue)

return, img
end    

pro overlay_boxes, boxes, scanlabel, postthick=postthick

    colors = ['orange','dark red']
    
    nbox = n_elements(boxes)
    if (boxes[0].scanlen gt 0.0) then begin

       boxindx = lindgen(nbox)
       
       scanlabel = [$
         'Scan Length = '+strtrim(string(boxes[0].scanlen,format='(I3)'),2)+'"', $
         'Slit PA = '+strtrim(string(boxes[0].scanlen_pa,format='(I3)'),2), $
         'Offset Perp = '+strtrim(string(boxes[0].offset_perp,format='(I5)'),2)+'"', $
         'Offset'+string(boxindx+1L,format='(I0)')+' Along = '+$
         strtrim(string(boxes[boxindx].offset_along,format='(I5)'),2)+'"']
       
       for ibox = 0L, nbox-1L do begin
          im_oplot_box, boxes[ibox].scanlen/60.0, 3.5, $
            boxes[ibox].scanlen_pa, line=0, color=djs_icolor(colors[ibox]), $
            thick=postthick+2, xoffset=boxes[ibox].offset_perp/60.0, $
            yoffset=boxes[ibox].offset_along/60.0
       endfor

    endif else scanlabel = ''

return
end    

pro pilyugin_visualize, postscript=postscript
; jm05mar01uofa
; visualize the galaxy images

;   Y is "along" the slit, X is perpendicular to the slit. 

    root = '04pilyugin'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'

    dsspath = path+'DSS/'

    psname = 'pilyugin_visualize.ps'
    if keyword_set(postscript) then begin
       dfpsplot, path+psname, /square, /color
       postthick = 6.0
    endif else begin
       im_window, 0, xratio=0.6, /square
       postthick = 2.0
    endelse

; read the catalog; read the atlas and the 11HUGS galaxies; exclude
; southern galaxies, really big objects, Sy1's, and any objects in the
; atlas or in 11HUGS

    pilyugin = read_04pilyugin()
    boxes = rsex(path+'pilyugin_driftscan_input.txt')

    keep = where((pilyugin.atlas eq 0L) and (pilyugin.hugs eq 0L) and $
      (strmatch(pilyugin.morph,'*Sy1*') eq 0B) and (im_hms2dec(pilyugin.dec) gt -20.0) $
      and (pilyugin.dmaj lt 15.0) and (im_hms2dec(pilyugin.ra) gt 0.0) and $
      (im_hms2dec(pilyugin.ra) lt 20.0))
    pilyugin = pilyugin[keep]
    boxes = boxes[keep]

;   niceprint, pilyugin.galaxy, pilyugin.nedgalaxy, pilyugin.ra, pilyugin.dec
    
    galaxy = pilyugin.galaxy
    nedgalaxy = pilyugin.nedgalaxy
    ngalaxy = n_elements(galaxy)
    
    dssfits = dsspath+strlowcase(galaxy)+'.fits.gz'
    
    for i = 0L, ngalaxy-1L do begin

       i = i > 0L
       
       if (galaxy[i] eq nedgalaxy[i]) then name = galaxy[i] else $
         name = galaxy[i]+'/'+nedgalaxy[i]

       dssimage = readfits(dssfits[i],hdss,/silent)
       extast, hdss, astr

       imsize = size(dssimage,/dimension)
       xsize = imsize[0] & xcen = xsize/2.0 & ysize = imsize[1] & ycen = ysize/2.0

       xpixscale = (astr.pltscl*1E-3*astr.xsz)/60 ; [arcmin/pixel]
       ypixscale = (astr.pltscl*1E-3*astr.ysz)/60 ; [arcmin/pixel]

       xaxis = (findgen(xsize)-xcen)*xpixscale ; centered on the image [arcsec]
       yaxis = (findgen(ysize)-ycen)*ypixscale ; centered on the image [arcsec]

       img = pilyugin_imgscl(dssimage)
       
       plotimage, img, /preserve_aspect, position=pos, /normal, imgxrange=minmax(xaxis), $
         imgyrange=minmax(yaxis), charsize=1.8, charthick=postthick, xthick=postthick, $
         ythick=postthick, xtitle=textoidl('\Delta\alpha [arcmin]'), $
         ytitle=textoidl('\Delta\delta [arcmin]'), title=name
       im_oplot_box, pilyugin[i].dmin, pilyugin[i].dmaj, pilyugin[i].pa, line=2, $
         color=djs_icolor('red'), thick=postthick

       overlay_boxes, boxes[i], scanlabel, postthick=postthick

       label1 = [$
         '\alpha = '+pilyugin[i].ra, $
         '\delta = '+pilyugin[i].dec, $
         'Dmaj = '+strtrim(string(pilyugin[i].dmaj,format='(F5.1)'),2)+"'", $
         'Dmin = '+strtrim(string(pilyugin[i].dmin,format='(F5.1)'),2)+"'", $
         '\Delta\theta = '+strtrim(string(pilyugin[i].pa,format='(I3)'),2), $
         'NHII = '+strtrim(string(pilyugin[i].n_hii,format='(I0)'),2)+"'"]

       legend, textoidl(label1), /left, /top, box=0, charsize=1.5, $
         charthick=postthick, textcolor=djs_icolor('black'), /clear

       if (not keyword_set(postscript)) then begin
          print, galaxy[i]+' - '+string(i+1L,format='(I0)')+' [Options: b,g,q]'
          cc = get_kbrd(1)
          case strlowcase(strcompress(cc,/remove)) of
             'b': begin
                i = (i-2L);>0L   ; back
             end
             'g': begin         ; goto 
                number = ''
                read, number, prompt='Goto galaxy (1-'+string(ngalaxy,format='(I0)')+'): '
                number = 0 > long(number) < (ngalaxy-1L)
                i = (number-2L);>0L
             end
             'q': return        ; quit
                else: 
          endcase
             
       endif

       if keyword_set(postscript) then begin ; generate the airmass plot

          airmass_plots, '2005-03-11', pilyugin[i].ra, pilyugin[i].dec, $
            object=pilyugin[i].galaxy, obsname='KPNO'

       endif

    endfor

    if keyword_set(postscript) then begin
       dfpsclose
       spawn, ['gzip -f '+path+psname], /sh
    endif
    
stop
    
return
end
    
