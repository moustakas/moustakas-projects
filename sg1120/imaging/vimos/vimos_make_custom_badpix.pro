;+
; NAME:
;       VIMOS_MAKE_BADPIX
;
; PURPOSE:
;       Generate customized bad pixel masks for each image, adding it
;       to the generic bad pixel mask generated for each QUADRANT.
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
;       jm07aug28nyu - updated to generate customized bad pixel masks
;                      for each object
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

pro vimos_make_custom_badpix, flist=flist, datapath=datapath, trimsec=trimsec, $
  quadrant=quadrant, runsex=runsex, wfits=wfits

    if (n_elements(datapath) eq 0L) then datapath = './'
    if (n_elements(quadrant) eq 0L) then quadrant = 'Q1'

    pixscale = 0.205
    
    if (n_elements(flist) eq 0L) then begin
       flist = file_search(datapath+'a.sg1120*_[B,V,R]_'+quadrant+'.fits',count=nimage)
       if (nimage eq 0L) then begin
          splog, 'No files found in '+datapath
          return
       endif
    endif
    catlist = repstr(flist,'.fits','.cat')
    seglist = repstr(flist,'.fits','.seg.fits')
    nimage = n_elements(flist)
    
    badpixlist = repstr(flist,'.fits','.badpix.fits')

    masterbadpixfile = datapath+'badpix_'+quadrant+'.fits'
    if (file_test(masterbadpixfile,/regular) eq 0L) then begin
       splog, 'Master bad pixel mask '+masterbadpixfile+' not found.'
       return
    endif else masterbadpix = mrdfits(masterbadpixfile,0,/silent)
    
    sexpath = sg1120_path(/sex)
    sexconfig = sexpath+'default.sex'
    wwconfig = sexpath+'default.ww'
    sexparam = sexpath+'sg1120.sex.param.badpix'
    sexconv = sexpath+'default.conv'
    sexnnw = sexpath+'default.nnw'

; now loop on each image    
    
    for ii = 0L, nimage-1L do begin

       if (file_test(flist[ii],/regular) eq 0L) then begin
          splog, 'Image '+flist[ii]+' not found.'
          return
       endif
       
       hdr = headfits(flist[ii])
       imsize = [sxpar(hdr,'NAXIS1'),sxpar(hdr,'NAXIS2')]
       badpix = masterbadpix
       
       if (n_elements(trimsec) eq 0L) then trimsec = [0L,imsize[0]-1L,0L,imsize[1]-1L] ; no trim

; ---------------------------------------------------------------------------
       if strmatch(flist[ii],'*a.sg1120_2006-02-04T08:45:50.646_R_Q2.fits*') then begin ; satellite trail

          xcen = 1030.0 & ycen = 958.0 & xsize = 2300.0 & ysize = 20.0 & angle = 152.3 ; ds9 coordinates
          xx = [-xsize/2.0,+xsize/2.0,+xsize/2.0,-xsize/2.0]
          yy = [-ysize/2.0,-ysize/2.0,+ysize/2.0,+ysize/2.0]

          corners = im_offset_and_rotate(transpose([[xx],[yy]]),angle-180.0,$
            xoffset=0.0,yoffset=0.0)+rebin([xcen,ycen],2,4)
          corners[0,*] = (corners[0,*]>0L)<(imsize[0]-1L)
          corners[1,*] = (corners[1,*]>0L)<(imsize[1]-1L)

          morebad = polyfillv(corners[0,*],corners[1,*],imsize[0],imsize[1])
          badpix[morebad] = badpix[morebad] or im_flagval('PIXMASK','SATTRAIL')

;         morebadpix = byte(badpix*0.0)+1B
;         morebadpix[morebad] = 0B
;         badpix = (badpix + morebadpix) eq 2B

       endif
       
; ---------------------------------------------------------------------------
       if strmatch(flist[ii],'*a.sg1120_2006-02-04T08:45:49.664_R_Q4.fits*') then begin ; satellite trail

          xcen = 1750.0 & ycen = 2225.0 & xsize = 900.0 & ysize = 20.0 & angle = 152.3 ; ds9 coordinates
          xx = [-xsize/2.0,+xsize/2.0,+xsize/2.0,-xsize/2.0]
          yy = [-ysize/2.0,-ysize/2.0,+ysize/2.0,+ysize/2.0]

          corners = im_offset_and_rotate(transpose([[xx],[yy]]),angle-180.0,$
            xoffset=0.0,yoffset=0.0)+rebin([xcen,ycen],2,4)
          corners[0,*] = (corners[0,*]>0L)<(imsize[0]-1L)
          corners[1,*] = (corners[1,*]>0L)<(imsize[1]-1L)

          morebad = polyfillv(corners[0,*],corners[1,*],imsize[0],imsize[1])
          badpix[morebad] = badpix[morebad] or im_flagval('PIXMASK','SATTRAIL')

;         badpix[morebad] = 0B
;         morebadpix = byte(badpix*0.0)+1B
;         morebadpix[morebad] = 0B
;         badpix = (badpix + morebadpix) eq 2B

       endif
       
; ###########################################################################
; TEST CODE!!
; ###########################################################################
       
       if keyword_set(runsex) then begin

;         spawn, 'sex '+flist[ii]+' -c '+sexconfig+' -CATALOG_NAME '+catlist[ii]+' -CATALOG_TYPE ASCII_HEAD'+$
;           ' -DETECT_MINAREA 200 -DETECT_THRESH 30.0 -ANALYSIS_THRESH 0.5 -DEBLEND_MINCONT 0.0001 -CLEAN Y -FILTER Y'+$
;           ' -FILTER_NAME '+sexconv+' -STARNNW_NAME '+sexnnw+' -PARAMETERS_NAME '+sexparam+$
;           ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+masterbadpixfile+' -WEIGHT_THRESH 0'+$
;           ' -VERBOSE_TYPE NORMAL -NTHREADS 4 -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME '+seglist[ii], /sh
;         cat = rsex(catlist[ii])
;         struct_print, cat
;
;         seg = mrdfits(seglist[ii],0,/silent)
;         badpix = seg eq 0
;
;         bad = where(cat.flux_radius*pixscale gt 10.0,nbad)
;         badpix = byte(seg*0.0)+1B
;
;         for ibad = 0L, nbad-1L do badpix[where(seg eq cat[bad[ibad]].number)] = 0B
;         badpix = (badpix + masterbadpix) eq 2B

; try to use DIMAGE - not so great at detecting artifacts!
          
          image = mrdfits(flist[ii],0,hdr,/silent)
          sky, image[where(masterbadpix ne 0)], skymode, skysig, /silent
          invvar = image*0.0+1.0/skysig^2
;         invvar = (image*0.0+1.0/skysig^2)*masterbadpix

          splog, 'Median smoothing and detecting objects.'
          msimage = dmedsmooth(image,invvar,box=41L)
          simage = image - msimage
;         writefits, 'junk.fits', simage
;         writefits, 'junk.fits', invvar

          dobjects, simage, objects=obj, plim=10.0

          splog, 'Sorting and measuring object list.'
          id = obj[uniq(obj,sort(obj))]
          id = id[where(id ne -1L)]
          nobj = n_elements(id)
          xsize = fltarr(nobj)
          ysize = fltarr(nobj)
          for iobj = 0L, nobj-1L do begin
             pix = where(obj eq id[iobj])
             xpix = pix mod imsize[0]
             ypix = pix / imsize[0]
             xsize[iobj] = max(xpix)-min(xpix)
             ysize[iobj] = max(ypix)-min(ypix)
          endfor

          xysize = sqrt(xsize*ysize)

          badpix = byte(obj*0.0)+1B
          bad = where(xysize gt 60.0,nbad)
          for ibad = 0L, nbad-1L do badpix[where(obj eq id[bad[ibad]])] = 0B
          
;         writefits, seglist[ii], fix(obj)
;         writefits, seglist[ii], fix(obj*masterbadpix);, hdr

          badpix = (badpix + masterbadpix) eq 2B
          
       endif
; ###########################################################################

; now trim and write out

       mkhdr, hdr, badpix
;      mkhdr, hdr, byte(badpix)
       sxdelpar, hdr, 'COMMENT'
       sxdelpar, hdr, 'DATE'
       im_trim, badpix, trimsec=trimsec, hdr=hdr, /silent
       
       if keyword_set(wfits) then begin
          splog, 'Writing '+badpixlist[ii]
          mwrfits, badpix, badpixlist[ii], hdr, /create
;         mwrfits, byte(badpix), badpixlist[ii], hdr, /create
       endif
       
    endfor

return
end
