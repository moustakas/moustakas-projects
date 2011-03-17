; ###########################################################################
; OBSOLETE!!  SUPERSEDED BY UNPACK_LDSS3 AND LDSS3_MAKE_BADPIX!
; ###########################################################################
;
;+
; NAME:
;       LDSS3_INIT_WCS
;
; PURPOSE:
;       Add the preliminary WCS to the LDSS3 headers, and reflect 
;       chip 1 to account for the way the chip was read out.  Also
;       optionally apodize the images.
;
; INPUTS: 
;       imagelist - full path name to the images to process
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;       apodize - apodize
;
; OUTPUTS: 
;       Only the headers for chip 2 images are modified; chip 1 images
;       are reflected and written back out as FITS files.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       This routine should be run immediately after LDSS3_REDUX, /CCDPROC
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2007 Jan 17, NYU - written from a previous code 
;       jm07jan24nyu - also apodize the images on output
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

pro ldss3_init_wcs, imagefile, weight_suffix=weight_suffix, apodize=apodize, wfits=wfits

    nimage = n_elements(imagefile) 
    if (nimage eq 0L) then begin
       doc_library, 'ldss3_init_wcs'
       return
    endif

    orientation = 270.0         ; NOTE!

    for ii = 0L, nimage-1L do begin

       if (file_test(imagefile[ii],/regular) eq 0L) then begin
          splog, 'Image '+imagefile[ii]+' not found.'
          return
       endif
       
       splog, 'Initializing astrometry for image '+imagefile[ii]+'.'

       fits_info, imagefile[ii], n_ext=next, /silent
       image = mrdfits(imagefile[ii],0,hdr,/silent)
       imsize = size(image,/dim)
       
       ra = sxpar(hdr,'RA',count=racount)
       dec = sxpar(hdr,'DEC',count=deccount)
       pixscale = sxpar(hdr,'SCALE',count=scalecount)

       if ((racount+deccount+scalecount) ne 3L) then begin
          splog, 'Unexpected LDSS3 header!'
          return
       endif
       
       racen = 15.0*hms2dec(ra)
       deccen = hms2dec(dec)
       dra1 = imsize[0]*pixscale/3600.0
       ddec1 = imsize[1]*pixscale/3600.0

       astr = hogg_make_astr(racen,deccen,dra1,ddec1,$
         pixscale=pixscale/3600.0,orientation=orientation)

       putast, hdr, astr
       hdr = hdr[where((strmatch(strcompress(hdr,/remove),'*HISTORY*PUTAST*',/fold) eq 0B))]
       sxaddhist, "'Initial WCS astrometry added "+hogg_iso_date()+"'", hdr

; this code could break easily       
       
       if (n_elements(weight_suffix) eq 0L) then begin
          weight = mrdfits(imagefile[ii],1,weighthdr,/silent)
       endif else begin
          weightfile = repstr(imagefile[ii],'.fits',weight_suffix+'.fits')
          weight = mrdfits(weightfile,0,weighthdr,/silent)
       endelse

       putast, weighthdr, astr
       weighthdr = weighthdr[where((strmatch(strcompress(weighthdr,/remove),'*HISTORY*PUTAST*',/fold) eq 0B))]
       sxaddhist, "'Initial WCS astrometry added "+hogg_iso_date()+"'", weighthdr

; write out       
       
;      if keyword_set(wfits) then begin
;         modfits, imagefile[ii], 0, hdr, exten_no=0
;         if (n_elements(weight) ne 0L) then modfits, imagefile[ii], 0, weighthdr, exten_no=1
;      endif

; reflect chip 1 to account for the way it was read out, and apodize
       
       if keyword_set(wfits) then begin
          if strmatch(imagefile[ii],'*c1*') then begin
             image = reverse(temporary(image))
             weight = reverse(temporary(weight))
          endif
          if keyword_set(apodize) then begin
             apimage = ldss3_apodize(size(image,/dim),file=imagefile[ii],$
               chip1=strmatch(imagefile[ii],'*c1*'))
          endif else begin
             apimage = image*0.0+1.0 ; no apodization
          endelse
          if (n_elements(weight_suffix) eq 0L) then begin
             mwrfits, float(apimage*image), imagefile[ii], hdr, /create
             mwrfits, float(apimage*weight), imagefile[ii], weighthdr
          endif else begin
             mwrfits, float(apimage*image), imagefile[ii], hdr, /create
             mwrfits, byte(float(apimage*weight)), weightfile, weighthdr, /create
          endelse

;         modfits, imagefile[ii], 0, hdr, exten_no=0
;         if (n_elements(weight) ne 0L) then modfits, imagefile[ii], 0, weighthdr, exten_no=1

       endif

    endfor
       
return    
end


   
