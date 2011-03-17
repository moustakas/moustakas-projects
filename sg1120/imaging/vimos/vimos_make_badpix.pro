;+
; NAME:
;       VIMOS_MAKE_BADPIX
;
; PURPOSE:
;       Generate generic bad pixel masks for each VIMOS quadrant.
;       Differentiate between bad columns and vignetted regions.
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2007 Jun 20, NYU - based on VIMOS_APODIZE 
;       jm08jul30nyu - use IM_FLAGVAL() to distinguish between bad
;         columns and vignetted ("no data") regions
;
; Copyright (C) 2007-2008, John Moustakas
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

pro vimos_make_badpix, suffix=suffix, nx=nx, ny=ny, $
  trimsec=trimsec, outpath=outpath, wfits=wfits

    if (strmatch(suffix,'2003') eq 0B) and (strmatch(suffix,'2006') eq 0B) then begin
       splog, 'Unrecognized SUFFIX keyword!'
       return
    endif

    if (n_elements(nx) eq 0L) or (n_elements(ny) eq 0L) then begin
       doc_library, 'vimos_make_badpix'
       return
    endif

    if (n_elements(trimsec) eq 0L) then trimsec = [0L,nx-1L,0L,ny-1L] ; no trim

    ximage = findgen(nx) # replicate(1.0,ny)
    yimage = replicate(1.0,nx) # findgen(ny)
    
    quadrant = ['Q1','Q2','Q3','Q4']

    for iq = 1L, n_elements(quadrant) do begin

       badpix = lonarr(nx,ny)
;      badpix = bytarr(nx,ny) + 1B ; 1 = good, 0 = bad

       case iq of
          1L: begin
             if (suffix eq '2003') then begin
                xcen = 1115.0-1.0 & ycen = 1220.0-1.0
                xsize = 1930.0    & ysize = 2280.0
             endif
             if (suffix eq '2006') then begin
                xcen = 1115.0-1.0 & ycen = 1250.0-1.0
                xsize = 1930.0    & ysize = 2200.0
             endif
; bad columns
             badpix[743:743,0:ny-1] = badpix[743:743,0:ny-1] or $
               im_flagval('PIXMASK','BAD_COLUMN')
             badpix[1404:1404,1696:ny-1] = badpix[1404:1404,1696:ny-1] or $
               im_flagval('PIXMASK','BAD_COLUMN')
; vignetted corners
             badpix[0:380,2177:ny-1] = badpix[0:380,2177:ny-1] or $
               im_flagval('PIXMASK','NODATA') ; 0B
             badpix[1928:nx-1,2250:ny-1] = badpix[1928:nx-1,2250:ny-1] or $
               im_flagval('PIXMASK','NODATA') ; 0B
          end
          2L: begin
             if (suffix eq '2003') then begin
                xcen = 1020.0-1.0 & ycen = 1210.0-1.0
                xsize = 1900.0    & ysize = 2240.0
             endif
             if (suffix eq '2006') then begin
                xcen = 1020.0-1.0 & ycen = 1240.0-1.0
                xsize = 1900.0    & ysize = 2160.0
             endif
             badpix[150:150,0:ny-1] = badpix[150:150,0:ny-1] or im_flagval('PIXMASK','BAD_COLUMN')
             badpix[840:840,0:ny-1] = badpix[840:840,0:ny-1] or im_flagval('PIXMASK','BAD_COLUMN')
          end
          3L: begin
             if (suffix eq '2003') then begin
                xcen = 1080.0-1.0 & ycen = 1400.0-1.0
                xsize = 1990.0    & ysize = 1890.0
             endif
             if (suffix eq '2006') then begin
                xcen = 1080.0-1.0 & ycen = 1400.0-1.0
                xsize = 1990.0    & ysize = 1760.0
             endif
;            badpix = badpix + ((ximage gt (xcen-xsize/2.0)) and (ximage lt (xcen+xsize/2.0)) and $
;              (yimage gt (ycen-ysize/2.0)) and (yimage lt (ycen+ysize/2.0))) eq 2B
             badpix[1328:1329,2282:ny-1] = badpix[1328:1329,2282:ny-1] or im_flagval('PIXMASK','BAD_COLUMN')
          end
          4L: begin
             if (suffix eq '2003') then begin
                xcen = 1110.0-1.0 & ycen = 1325.0-1.0
                xsize = 1920.0    & ysize = 2020.0
             endif
             if (suffix eq '2006') then begin
                xcen = 1110.0-1.0 & ycen = 1200.0-1.0 ; NOTE DIFFERENT REGION RELATIVE TO 2003!
                xsize = 1920.0    & ysize = 2110.0
             endif
;            badpix = badpix + ((ximage gt (xcen-xsize/2.0)) and (ximage lt (xcen+xsize/2.0)) and $
;              (yimage gt (ycen-ysize/2.0)) and (yimage lt (ycen+ysize/2.0))) eq 2B
          end
       endcase 

; mask vignetted/no data regions       
       
       edgemask = where((ximage lt (xcen-xsize/2.0)) or (ximage gt (xcen+xsize/2.0)) or $
         (yimage lt (ycen-ysize/2.0)) or (yimage gt (ycen+ysize/2.0)),nedgemask)
       if (nedgemask ne 0L) then badpix[edgemask] = badpix[edgemask] or $
         im_flagval('PIXMASK','NODATA')

;      flatfiles = file_search(outpath+'flat_*_'+quadrant[iq-1L]+'.fits',count=nflat)
;      flats = im_fits_cube(flatfiles) ; now read the flat-fields
;      flatbadpix = total((flats.image gt 0.9),nflat) eq nflat ; must be a good pixel in all flats

;      window, 0, xs=450, ys=450 & plotimage, logscl(badpix), /preserve
;      window, 1, xs=450, ys=450 & plotimage, logscl(flatbadpix), /preserve
;      window, 2, xs=450, ys=450 & plotimage, logscl(flatbadpix*badpix), /preserve

;      badpix = badpix*flatbadpix ; this looks crappy so don't do it

; now trim the bad pixel mask

       mkhdr, hdr, badpix
;      mkhdr, hdr, byte(badpix)
       sxdelpar, hdr, 'COMMENT'
       sxdelpar, hdr, 'DATE'
       im_trim, badpix, trimsec=trimsec, hdr=hdr, /silent

; finally reset the "no data" bad pixel values to just "no data"

       nodata = where((badpix and im_flagval('PIXMASK','NODATA')),nnodata)
       if (nnodata ne 0L) then $
         badpix[nodata] = im_flagval('PIXMASK','NODATA')
       
       if keyword_set(wfits) then begin
          outfile = outpath+'badpix_'+quadrant[iq-1L]+'.fits'
          splog, 'Writing '+outfile
          mwrfits, badpix, outfile, hdr, /create
;         mwrfits, byte(badpix), outfile, hdr, /create
       endif

    endfor  

return
end
