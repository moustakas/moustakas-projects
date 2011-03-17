pro sings_visualize_boxes, boxes, irsinfo, scanlabel, postthick=postthick
; jm05apr05uofa - see SINGS_HTML
    
    color1 = 'orange'
    color2 = 'dark red'
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

          im_oplot_box, boxes[ibox].scanlen/60.0, 3.5, $
            boxes[ibox].scanlen_pa, line=0, color=djs_icolor(color), $
            thick=postthick+2, xoffset=boxes[ibox].offset_perp/60.0, $
            yoffset=boxes[ibox].offset_along/60.0, corners=corners
          if (boxes[ibox].offset_along eq 0.0) then begin
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

