;+
; NAME:
;       LDSS3_MAKE_BADPIX
;
; PURPOSE:
;       Generate bad pixel masks for each LDSS3 image.
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; COMMENTS:
;       Generate the bad pixel mask using an elliptical aperture that
;       I defined "by hand" using DS9; this apodization image should
;       be applied to the inverse variance map *and* the data because
;       it also crops nicely.
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2007 Aug 21, NYU - based on LDSS3_APODIZE 
;       jm08sep03nyu - use IM_FLAGVAL() to distinguish between bad
;         columns and vignetted ("no data") regions
;       jm09mar23nyu - add some bad columns from c1 to the bad pixel
;         mask 
;
; Copyright (C) 2007-2009, John Moustakas
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

pro ldss3_make_badpix, flist=flist, datapath=datapath, trimsec=trimsec, wfits=wfits

    if (n_elements(datapath) eq 0L) then datapath = './'

    if (n_elements(flist) eq 0L) then begin
;      flist = file_search(datapath+'a.????_sg1120_?_[g,r]_c?.fits',count=nimage) ; SDSS fields ignored!!
       flist = file_search(count=nimage,datapath+$
         ['a.????_sg1120_?_[g,r]_c?.fits','a.????_sdss_1_[g,r]_c?.fits'])
       if (nimage eq 0L) then begin
          splog, 'No files found in '+datapath
          return
       endif
    endif
    nimage = n_elements(flist)

    badpixlist = repstr(flist,'.fits','.badpix.fits')

    for ii = 0L, nimage-1L do begin

       if (file_test(flist[ii],/regular) eq 0L) then begin
          splog, 'Image '+flist[ii]+' not found.'
          return
       endif
       
       hdr = headfits(flist[ii])
       imsize = [sxpar(hdr,'NAXIS1'),sxpar(hdr,'NAXIS2')]
       badpix = lonarr(imsize)
       
       if (n_elements(trimsec) eq 0L) then trimsec = [0L,imsize[0]-1L,0L,imsize[1]-1L] ; no trim
       
       if strmatch(flist[ii],'*c1*') then rmajor = 1250.0 else rmajor = 1280.0 ; determined using DS9
       e_xcen = trimsec[1]
       e_ycen = fix((trimsec[3]-trimsec[2])/2.0)+trimsec[2]-10L

       dist_circle, distimage, imsize, e_xcen, e_ycen
       nodata = where(distimage gt rmajor)
       badpix[nodata] = badpix[nodata] or im_flagval('PIXMASK','NODATA')
       
;      badpix = (distimage lt rmajor) ; 1 = good, 0 = bad

; bad columns
       if strmatch(flist[ii],'*c1*') then begin
          badpix[1733:1734,1301:*] = badpix[1733:1734,1301:*] or $
            im_flagval('PIXMASK','BAD_COLUMN')
       endif

; customized bad pixel masks; depends on the trim!!

; ---------------------------------------------------------------------------    
       if strmatch(flist[ii],'*2016_sg1120_2_g_c1*') then begin
          rmajor = 430.0 & rminor = 200.0 & e_pa = 0.0 & e_xcen = 815.0 & e_ycen = 2105.0
          dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
          morebad = where(mask lt rmajor)
          badpix[morebad] = badpix[morebad] or im_flagval('PIXMASK','NODATA')
;         badpix = temporary(badpix) * (mask gt rmajor)
       endif
       if strmatch(flist[ii],'*2017_sg1120_2_r_c1*') then begin
          rmajor = 430.0 & rminor = 200.0 & e_pa = 0.0 & e_xcen = 815.0 & e_ycen = 2105.0
          dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
          morebad = where(mask lt rmajor)
          badpix[morebad] = badpix[morebad] or im_flagval('PIXMASK','NODATA')
;         badpix = temporary(badpix) * (mask gt rmajor)
       endif
; ---------------------------------------------------------------------------    
       if strmatch(flist[ii],'*2018_sg1120_2_r_c1*') then begin
          rmajor = 492.0 & rminor = 270.0 & e_pa = 155.0 & e_xcen = 935.0 & e_ycen = 2510.0
          dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
          morebad = where(mask lt rmajor)
          badpix[morebad] = badpix[morebad] or im_flagval('PIXMASK','NODATA')
;         badpix = temporary(badpix) * (mask gt rmajor)
       endif

       if strmatch(flist[ii],'*2019_sg1120_2_g_c1*') then begin
          rmajor = 492.0 & rminor = 270.0 & e_pa = 155.0 & e_xcen = 935.0 & e_ycen = 2510.0
          dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
          morebad = where(mask lt rmajor)
          badpix[morebad] = badpix[morebad] or im_flagval('PIXMASK','NODATA')
;         badpix = temporary(badpix) * (mask gt rmajor)
       endif
       
; ---------------------------------------------------------------------------    
       if strmatch(flist[ii],'*2018_sg1120_2_r_c2*') then begin
          rmajor = 545.0 & rminor = 190.0 & e_pa = 330.0 & e_xcen = 1432.0 & e_ycen = 1135.0
          dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa-90.0 ; note angle offset
          morebad = where(mask lt rmajor)
          badpix[morebad] = badpix[morebad] or im_flagval('PIXMASK','NODATA')
;         badpix = temporary(badpix) * (mask gt rmajor)
       endif

       if strmatch(flist[ii],'*2019_sg1120_2_g_c2*') then begin
          rmajor = 545.0 & rminor = 190.0 & e_pa = 330.0 & e_xcen = 1432.0 & e_ycen = 1135.0
          dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa-90.0 ; note angle offset
          morebad = where(mask lt rmajor)
          badpix[morebad] = badpix[morebad] or im_flagval('PIXMASK','NODATA')
;         badpix = temporary(badpix) * (mask gt rmajor)
       endif

; ---------------------------------------------------------------------------    
       if strmatch(flist[ii],'*2020_sg1120_2_r_c1*') then begin
          rmajor = 820.0 & rminor = 525.0 & e_pa = 0.0 & e_xcen = 1145.0 & e_ycen = 1995.0
          dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
          morebad = where(mask lt rmajor)
          badpix[morebad] = badpix[morebad] or im_flagval('PIXMASK','NODATA')
;         badpix = temporary(badpix) * (mask gt rmajor)
       endif

       if strmatch(flist[ii],'*2021_sg1120_2_g_c1*') then begin
          rmajor = 820.0 & rminor = 525.0 & e_pa = 0.0 & e_xcen = 1145.0 & e_ycen = 1995.0
          dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
          morebad = where(mask lt rmajor)
          badpix[morebad] = badpix[morebad] or im_flagval('PIXMASK','NODATA')
;         badpix = temporary(badpix) * (mask gt rmajor)
       endif

; ---------------------------------------------------------------------------    
       if strmatch(flist[ii],'*2022_sg1120_3_g_c2*') then begin
          rmajor = 260.0 & rminor = 52.0 & e_pa = 35.0 & e_xcen = 967.0 & e_ycen = 1400.0
          dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
          morebad = where(mask lt rmajor)
          badpix[morebad] = badpix[morebad] or im_flagval('PIXMASK','NODATA')
;         badpix = temporary(badpix) * (mask gt rmajor)
       endif

       if strmatch(flist[ii],'*2023_sg1120_3_r_c2*') then begin
          rmajor = 260.0 & rminor = 52.0 & e_pa = 35.0 & e_xcen = 967.0 & e_ycen = 1400.0
          dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
          morebad = where(mask lt rmajor)
          badpix[morebad] = badpix[morebad] or im_flagval('PIXMASK','NODATA')
;         badpix = temporary(badpix) * (mask gt rmajor)
       endif

       if strmatch(flist[ii],'*2024_sg1120_3_r_c2*') then begin
          rmajor = 260.0 & rminor = 52.0 & e_pa = 35.0 & e_xcen = 967.0 & e_ycen = 1400.0
          dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
          morebad = where(mask lt rmajor)
          badpix[morebad] = badpix[morebad] or im_flagval('PIXMASK','NODATA')
;         badpix = temporary(badpix) * (mask gt rmajor)
       endif

       if strmatch(flist[ii],'*2025_sg1120_3_g_c2*') then begin
          rmajor = 260.0 & rminor = 52.0 & e_pa = 35.0 & e_xcen = 967.0 & e_ycen = 1400.0
          dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
          morebad = where(mask lt rmajor)
          badpix[morebad] = badpix[morebad] or im_flagval('PIXMASK','NODATA')
;         badpix = temporary(badpix) * (mask gt rmajor)
       endif

; ---------------------------------------------------------------------------    
       if strmatch(flist[ii],'*2032_sg1120_4_g_c2*') then begin
          rmajor = 550.0 & rminor = 320.0 & e_pa = 0.0 & e_xcen = 962.0 & e_ycen = 2430.0
          dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
          morebad = where(mask lt rmajor)
          badpix[morebad] = badpix[morebad] or im_flagval('PIXMASK','NODATA')
;         badpix = temporary(badpix) * (mask gt rmajor)
       endif

       if strmatch(flist[ii],'*2033_sg1120_4_r_c2*') then begin
          rmajor = 550.0 & rminor = 320.0 & e_pa = 0.0 & e_xcen = 962.0 & e_ycen = 2430.0
          dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
          morebad = where(mask lt rmajor)
          badpix[morebad] = badpix[morebad] or im_flagval('PIXMASK','NODATA')
;         badpix = temporary(badpix) * (mask gt rmajor)
       endif

;; ---------------------------------------------------------------------------    
;       if strmatch(flist[ii],'*2034_sdss_1_r_c1*') then begin
;          rmajor = 500.0 & rminor = 240.0 & e_pa = 40.0 & e_xcen = 1155.0 & e_ycen = 2835.0
;          dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
;          badpix = temporary(badpix) * (mask gt rmajor)
;       endif
;
;       if strmatch(flist[ii],'*2035_sdss_1_g_c1*') then begin
;          rmajor = 500.0 & rminor = 240.0 & e_pa = 40.0 & e_xcen = 1155.0 & e_ycen = 2835.0
;          dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
;          badpix = temporary(badpix) * (mask gt rmajor)
;       endif
;
;; ---------------------------------------------------------------------------    
;       if strmatch(flist[ii],'*2038_sdss_1_r_c1*') then begin
;          rmajor = 730.0 & rminor = 350.0 & e_pa = 35.0 & e_xcen = 1270.0 & e_ycen = 2740.0
;          dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
;          badpix = temporary(badpix) * (mask gt rmajor)
;       endif
;
;       if strmatch(flist[ii],'*2039_sdss_1_g_c1*') then begin
;          rmajor = 730.0 & rminor = 350.0 & e_pa = 35.0 & e_xcen = 1270.0 & e_ycen = 2740.0
;          dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
;          badpix = temporary(badpix) * (mask gt rmajor)
;       endif

; now trim the bad pixel mask

       mkhdr, hdr, fix(badpix)
;      mkhdr, hdr, byte(badpix)
       sxdelpar, hdr, 'COMMENT'
       sxdelpar, hdr, 'DATE'
       im_trim, badpix, trimsec=trimsec, hdr=hdr, /silent
       
       if keyword_set(wfits) then begin
          splog, 'Writing '+badpixlist[ii]
          mwrfits, fix(badpix), badpixlist[ii], hdr, /create
;         mwrfits, byte(badpix), badpixlist[ii], hdr, /create
       endif

    endfor
       
return
end
    
