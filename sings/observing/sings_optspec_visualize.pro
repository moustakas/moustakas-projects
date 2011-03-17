pro draw_legends, name, data, boxes, irsinfo, scanlabel, postscript=postscript

    if keyword_set(postscript) then begin
       postthick = 6.0
       charsize = 1.0
    endif else begin
       postthick = 2.0
       charsize = 1.3
    endelse

    label1 = [$
      '\alpha = '+data.ra, $
      '\delta = '+data.dec, $
      'Dmaj = '+strtrim(string(data.d25_maj,format='(F5.1)'),2)+"'",$
      'Dmin = '+strtrim(string(data.d25_min,format='(F5.1)'),2)+"'",$
      'PA = '+strtrim(string(data.posangle,format='(I3)'),2)]

;   label3 = [$
;     'B = '+strtrim(string(hugs.B,format='(F5.1)'),2),$
;     'M_{B} = '+strtrim(string(hugs.M_B,format='(F5.1)'),2),$
;     'T total = '+strtrim(string(hugs.expt_tot,format='(I10)'),2)+' min']
;     'T eff = '+strtrim(string(hugs.expt_eff,format='(F4.2)'),2)+' min',$
;     'S/N = '+strtrim(string(hugs.sn_hb,format='(I10)'),2)]

    legend, textoidl(label1), /left, /top, box=0, charsize=charsize, $
      charthick=postthick, textcolor=djs_icolor('black'), /clear
;   legend, textoidl(label3), /left, /bottom, box=0, charsize=charsize, $
;     charthick=postthick, textcolor=djs_icolor('black'), /clear

    comments = ''
    if (strmatch(boxes[0].comments1,'*...*') eq 0B) then comments = [comments,repstr(boxes[0].comments1,'_',' ')]

    if (n_elements(comments) gt 1L) then $
      legend, comments[1L:n_elements(comments)-1L], /left, /bottom, box=0, charsize=charsize, $
        charthick=postthick, textcolor=djs_icolor('black'), /clear

    if (scanlabel[0] ne '') then $
      legend, textoidl(scanlabel), /right, /top, box=0, charsize=charsize, $
        charthick=postthick, textcolor=djs_icolor('black'), /clear

    if (size(irsinfo,/type) eq 8L) then begin
       irslabel = [$
         '\alpha = '+irsinfo.ra, $
         '\delta = '+irsinfo.dec, $
;        ' ', $
         'IRS Pa = '+strtrim(string(im_angle_format(round(irsinfo.pa)),format='(I3)'),2), $
         'Width = '+strtrim(string(irsinfo.width*60.0,format='(I0)'),2)+'"', $
         'Length (2x) = '+strtrim(string(irsinfo.length_2x,format='(F5.1)'),2)+"'"]
;        '\Delta\theta = '+strtrim(string(im_angle_format(round(irsinfo.pa)),format='(I3)'),2)]
;        '1x = '+strtrim(string(irsinfo.length_1x,format='(F5.1)'),2)+"'", $
    endif else irslabel = ['NOT SCHEDULED']
    
    legend, textoidl(irslabel), /right, /bottom, box=0, charsize=charsize, $
      charthick=postthick, textcolor=djs_icolor('black'), /clear

return
end    
    
pro overlay_boxes, boxes, irsinfo, scanlabel, postthick=postthick

    color1 = 'orange'
    color2 = 'dark red'
    color3 = 'light blue'
;   colors = ['orange','dark red']
    
    nbox = n_elements(boxes)
    if (boxes[0].scanlen gt 0.0) then begin

       boxindx = lindgen(nbox)

       scanlabel = [$
         'Slit PA = '+strtrim(string(boxes[0].scanlen_pa,format='(I3)'),2), $
         'Rotator PA = '+strtrim(string(181.2-boxes[0].scanlen_pa,format='(I3)'),2), $
         'Scan Length = '+strtrim(string(boxes[0].scanlen,format='(I3)'),2)+'"', $
         'Offset Perp = '+strtrim(string(boxes[0].offset_perp,format='(I5)'),2)+'"', $
         'Number of Slits = '+string(nbox,format='(I0)')]

       if (nbox gt 1L) then begin
          slit_length = 3.5*60.0   ; [arcsec]
          noffsets = nbox / 2.0
          if odd(nbox) then $
            overlap = slit_length - abs(boxes[boxindx[0]].offset_along) else $
            overlap = slit_length - 2.0*min(abs(boxes[boxindx].offset_along))
          total_distance = noffsets*slit_length - overlap/2.0 - (noffsets-1)*overlap
          scanlabel = [scanlabel,$
            'Slit Overlap = '+strtrim(string(overlap,format='(I5)'),2)+'"',$
            'Diameter = '+strtrim(string(2*total_distance/60.0,format='(F5.1)'),2)+"'"]
       endif else scanlabel = [scanlabel,'Offset Along = '+strtrim(string(boxes[0].offset_along,format='(I5)'),2)+'"']
       
;      scanlabel = [$
;        'Scan Length = '+strtrim(string(boxes[0].scanlen,format='(I3)'),2)+'"', $
;        'Slit PA = '+strtrim(string(boxes[0].scanlen_pa,format='(I3)'),2), $
;        'Offset Perp = '+strtrim(string(boxes[0].offset_perp,format='(I5)'),2)+'"', $
;        'Offset'+string(boxindx+1L,format='(I0)')+' Along = '+$
;        strtrim(string(boxes[boxindx].offset_along,format='(I5)'),2)+'"']
       
       for ibox = 0L, nbox-1L do begin

          if (boxes[ibox].offset_along ge 0.0) then color = color1 else color = color2
          if (odd(nbox)) and (nbox gt 1L) then if (ibox eq nbox/2L) then color = color3

          im_oplot_box, boxes[ibox].scanlen/60.0, 3.5, $
            boxes[ibox].scanlen_pa, line=0, color=djs_icolor(color), $
            thick=postthick+2, xoffset=boxes[ibox].offset_perp/60.0, $
            yoffset=boxes[ibox].offset_along/60.0, corners=corners
          if (boxes[ibox].offset_along eq 0.0) or (nbox eq 1L) then begin
             arrow, corners[0,1], corners[1,1], corners[0,2], corners[1,2], /data, $
               hsize=-0.5, thick=postthick, hthick=postthick, color=djs_icolor(color)
             arrow, corners[0,0], corners[1,0], corners[0,3], corners[1,3], /data, $
               hsize=-0.5, thick=postthick, hthick=postthick, color=djs_icolor(color)
          endif else begin
             if (boxes[ibox].offset_along gt 0.0) then $
               arrow, corners[0,1], corners[1,1], corners[0,2], corners[1,2], /data, $
                 hsize=-0.5, thick=postthick, hthick=postthick, color=djs_icolor(color) else $
               arrow, corners[0,0], corners[1,0], corners[0,3], corners[1,3], /data, $
                 hsize=-0.5, thick=postthick, hthick=postthick, color=djs_icolor(color)
          endelse
       endfor

    endif else scanlabel = ''

; visualize the IRS 2x coverage

    if (size(irsinfo,/type) eq 8L) then $
      im_oplot_box, irsinfo.width, irsinfo.length_2x, im_angle_format(irsinfo.pa), $
        line=2, color=djs_icolor('yellow'), thick=postthick
;   im_oplot_box, boxes[0].scanlen/60.0, boxes[0].irs_length, boxes[0].scanlen_pa, $
;     line=2, color=djs_icolor('yellow'), thick=postthick

return
end    

pro sings_optspec_visualize, sings, postscript=postscript
; jm05mar22uofa

; visualize the 2005 April 4-8 SINGS galaxy images, some relevant
; information, and the scan parameters

; grep "LL/SL nuc" mips_irac_irs_schedule_05mar22.txt > irs_schedule_05mar22.txt
; mpage -2 -l -W140 irs_schedule_05mar22.txt > irs_schedule_05mar22.ps

; Y is "along" the slit, X is perpendicular to the slit. 

; slit_length = 3.5*60; slit_overlap = 40.0; the first slit offset is
; first_offset = (slit_length-slit_overlap)/2.0 = 85"; subsequently,
; slit_offsets = first_offset + (lindgen(noffsets)+1)*slit_length - slit_overlap/2.0

    if (n_elements(sings) eq 0L) then sings = sings_read_info()
    
    obspath = sings_path(/observing)
    dsspath = sings_path()+'DSS/'
    analysis_path = sings_path(/analysis)

; read the April spectroscopy sample, driftscan parameters, and the
; NED data; SAMPLE and DATA both contain the NED galaxy centers; the
; IRS centers, when available (see IRS_GETDSS), are in the DSS
; headers; also read the IRS parameters for other quantities of
; interest 

    irs = rsex(obspath+'irs_ll_params.txt')
    boxes = rsex(obspath+'sings_driftscan_params.txt')        

    galaxy = boxes.galaxy
    unique = uniq(galaxy)
    ugalaxy = galaxy[unique]
    altgalaxy = boxes[unique].alt_galaxy
    nicegalaxy = boxes[unique].nice_galaxy
    ngalaxy = n_elements(ugalaxy)

;   nopa = where(sings.posangle lt -900.0,nnopa) & if (nnopa ne 0L) then sings[nopa].posangle = 90.0
    
    dssfits = dsspath+strlowcase(nicegalaxy)+'.fits.gz'
    
    if keyword_set(postscript) then begin
       dfpsplot, obspath+'sings_visualize.ps', /square, /color
       postthick = 6.0
    endif else begin
       im_window, 0, xratio=0.6, /square
       postthick = 2.0
    endelse

    for i = 0L, ngalaxy-1L do begin

       redraw: 
       
       boxes = rsex(obspath+'sings_driftscan_params.txt')
       match = where(galaxy eq ugalaxy[i],nmatch)

       if (ugalaxy[i] eq altgalaxy[i]) then name = ugalaxy[i] else $
         name = ugalaxy[i]+'/'+altgalaxy[i]

       if (file_test(dssfits[i],/regular) eq 0L) then begin
          splog, 'DSS image for '+ugalaxy[i]+' not found!'
          return
       endif

       dssimage = readfits(dssfits[i],hdss,/silent)
       gsssextast, hdss, astr

;         irs_ra = repstr(strtrim(im_dec2hms(astr.crval[0]/15.0),2),' ',':')
;         irs_dec = repstr(strtrim(im_dec2hms(astr.crval[1]),2),' ',':')

       irs_ra = repstr(strtrim(sxpar(hdss,'OBJCTRA'),2),' ',':')
       irs_dec = repstr(strtrim(sxpar(hdss,'OBJCTDEC'),2),' ',':')
       
       imsize = size(dssimage,/dimension)
       xsize = imsize[0] & xcen = xsize/2.0 & ysize = imsize[1] & ycen = ysize/2.0

       xpixscale = (astr.pltscl*1E-3*astr.xsz)/60 ; [arcmin/pixel]
       ypixscale = (astr.pltscl*1E-3*astr.ysz)/60 ; [arcmin/pixel]

       xaxis = (findgen(xsize)-xcen)*xpixscale ; centered on the image [arcsec]
       yaxis = (findgen(ysize)-ycen)*ypixscale ; centered on the image [arcsec]

       img = logscl(dssimage,exponent=1.5,negative=keyword_set(postscript),omin=35,omax=255)
       
       plotimage, img, /preserve_aspect, position=pos, /normal, imgxrange=minmax(xaxis), $
         imgyrange=minmax(yaxis), charsize=1.8, charthick=postthick, xthick=postthick, $
         ythick=postthick, xtitle=textoidl('\Delta\alpha [arcmin]'), $
         ytitle=textoidl('\Delta\delta [arcmin]'), title=name

; overplot the RC3 box

       gsssadxy, astr, 15.0*im_hms2dec(sings[i].ra), im_hms2dec(sings[i].dec), xrc3, yrc3
       xrc3 = (xrc3 - xcen)*xpixscale & yrc3 = (yrc3 - ycen)*ypixscale
       
;      im_oplot_box, sings[i].d25_min, sings[i].d25_maj, sings[i].pa, /noplot, corners=corners
;      indx = [0,1,2,3,0]
;      djs_oplot, [corners[0,indx],corners[0,indx+1]]+xrc3, $
;        [corners[1,indx],corners[1,indx+1]]+yrc3, line=2, $
;        color=djs_icolor('red'), thick=postthick

       if (sings[i].posangle gt -900.0) then tvellipse, sings[i].d25_maj/2.0, $
         sings[i].d25_min/2.0, xrc3, yrc3, 90+sings[i].posangle, color=djs_icolor('red'), $
         thick=postthick, /data, line=2

       boxinfo = boxes[match]

; is there IRS data on this object?

       match, ugalaxy[i], irs.nice_galaxy, indx1, indx2
       if (indx1[0] ne -1L) then irsinfo = irs[indx2] else begin
          splog, 'No IRS information on '+ugalaxy[i]+'.'
          irsinfo = -1L
       endelse
       
; finally  visualize everything
       
       overlay_boxes, boxes[match], irsinfo, scanlabel, postthick=postthick
       draw_legends, name, sings[i], boxinfo, irsinfo, scanlabel, postscript=postscript

       if (not keyword_set(postscript)) then begin
          print, ugalaxy[i]+' - '+string(i+1L,format='(I0)')+' [Options: b,g,q]'
;            print, 'Galaxy '+string(i+1L,format='(I0)')+' [Options: b,g,q]'
          cc = get_kbrd(1)
          case strlowcase(strcompress(cc,/remove)) of
             'b': begin
                i = (i-1L)>0L   ; back
                goto, redraw
             end
             'g': begin         ; goto 
                number = ''
                read, number, prompt='Goto galaxy (1-'+string(ngalaxy,format='(I0)')+'): '
                number = 0 > long(number) < (ngalaxy-1L)
                i = (number-1L)>0L
                goto, redraw
             end
             'q': return        ; quit
             else: 
          endcase
          
       endif

;      if keyword_set(postscript) then begin ; generate the airmass plot
;         airmass_plots, '2005-04-5', sample[i].ra, sample[i].dec, $
;           object=sample[i].galaxy, obsname='KPNO'
;      endif

    endfor 

    if keyword_set(postscript) then begin
       dfpsclose
;      spawn, ['gzip -f '+path+'sings_05apr_visualize.ps'], /sh
    endif
    
stop
    
return
end
    
