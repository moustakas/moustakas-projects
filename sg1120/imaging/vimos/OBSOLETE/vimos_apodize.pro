;+
; NAME:
;       VIMOS_APODIZE
;
; PURPOSE:
;       Apodize the VIMOS images, immediately after CCD processing.
;
; INPUTS: 
;       imagelist - full path name to the images to process
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2007 Jan 25, NYU - written
;
; Copyright (C) 2007, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

pro vimos_apodize, imagefile, weight_suffix=weight_suffix, wfits=wfits

    nimage = n_elements(imagefile) 
    if (nimage eq 0L) then begin
       doc_library, 'vimos_apodize'
       return
    endif

    for ii = 0L, nimage-1L do begin

       if (file_test(imagefile[ii],/regular) eq 0L) then begin
          splog, 'Image '+imagefile[ii]+' not found.'
          continue
       endif
       
       splog, 'Apodizing '+imagefile[ii]

       fits_info, imagefile[ii], n_ext=next, /silent
       image = mrdfits(imagefile[ii],0,hdr,/silent)
       imsize = size(image,/dim)
       nx = imsize[0] & ny = imsize[1]
       
; this code could break easily       
       
       if (n_elements(weight_suffix) eq 0L) then begin
          weight = mrdfits(imagefile[ii],1,weighthdr,/silent)
       endif else begin
          weightfile = repstr(imagefile[ii],'.fits',weight_suffix+'.fits')
          weight = mrdfits(weightfile,0,weighthdr,/silent)
       endelse

       ximage = findgen(nx) # replicate(1.0,ny)
       yimage = replicate(1.0,nx) # findgen(ny)
       
       if strmatch(imagefile[ii],'*Q1*') then begin
          if strmatch(imagefile[ii],'*2003*') then begin
             xcen = 1115.0-1.0 & ycen = 1220.0-1.0
             xsize = 1930.0    & ysize = 2280.0
          endif else begin
             xcen = 1115.0-1.0 & ycen = 1250.0-1.0
             xsize = 1930.0    & ysize = 2200.0
;            xcen = 1115.0-1.0 & ycen = 1220.0-1.0
;            xsize = 1930.0    & ysize = 2280.0
          endelse
       endif
       
       if strmatch(imagefile[ii],'*Q2*') then begin
          if strmatch(imagefile[ii],'*2003*') then begin
             xcen = 1020.0-1.0 & ycen = 1210.0-1.0
             xsize = 1900.0    & ysize = 2240.0
          endif else begin
             xcen = 1020.0-1.0 & ycen = 1240.0-1.0
             xsize = 1900.0    & ysize = 2160.0
;            xcen = 1020.0-1.0 & ycen = 1210.0-1.0
;            xsize = 1900.0    & ysize = 2240.0
          endelse
       endif
       
       if strmatch(imagefile[ii],'*Q3*') then begin
          if strmatch(imagefile[ii],'*2003*') then begin
             xcen = 1080.0-1.0 & ycen = 1400.0-1.0
             xsize = 1990.0    & ysize = 1890.0
          endif else begin
             xcen = 1080.0-1.0 & ycen = 1400.0-1.0
             xsize = 1990.0    & ysize = 1760.0
;            xcen = 1080.0-1.0 & ycen = 1420.0-1.0
;            xsize = 1990.0    & ysize = 1840.0
          endelse
       endif
       
       if strmatch(imagefile[ii],'*Q4*') then begin
          if strmatch(imagefile[ii],'*2003*') then begin
             xcen = 1110.0-1.0 & ycen = 1225.0-1.0
             xsize = 1920.0    & ysize = 2220.0
          endif else begin
             xcen = 1110.0-1.0 & ycen = 1200.0-1.0
             xsize = 1920.0    & ysize = 2110.0
;            xcen = 1110.0-1.0 & ycen = 1225.0-1.0
;            xsize = 1920.0    & ysize = 2220.0
          endelse
       endif
       
       apimage = ((ximage gt (xcen-xsize/2.0)) and (ximage lt (xcen+xsize/2.0)) and $
         (yimage gt (ycen-ysize/2.0)) and (yimage lt (ycen+ysize/2.0))) gt 0

       if keyword_set(wfits) then begin
          if (n_elements(weight_suffix) eq 0L) then begin
             mwrfits, float(apimage*image), imagefile[ii], hdr, /create
             mwrfits, float(apimage*weight), imagefile[ii], weighthdr
          endif else begin
             mwrfits, float(apimage*image), imagefile[ii], hdr, /create
             mwrfits, byte(float(apimage*weight)), weightfile, weighthdr, /create
          endelse

       endif 
    endfor
       
return    
end


   
