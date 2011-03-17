;+
; NAME:
;       SINGS2D_STITCH
;
; PURPOSE:
;       Stitch together, e.g., East-West observations of the SINGS 56"
;       spectral scans.
;
; CALLING SEQUENCE:
;       sings2d_stitch, speclist, user_pixshift=, wmapname=, $
;         minwave=, dwave=, stitchrange=, scalefactor=, refindx=, $
;         outname=, object=, datapath=, /debug, /wfits, /gzip
;
; INPUTS:
;       speclist - list of 2D spectra to stitch together; the
;                  reference spectrum is determined by the minimum
;                  value of USER_PIXSHIFT [NSPEC]
;       user_pixshift - shift each 2D spectrum other than the zeroth
;                       spectrum by this amount [NSPEC-1, PIXEL]; must
;                       be positive
;       wmapname - use this wavelength map to rectify
;                  (wavelength-calibrate) the spectra before combining 
;
; OPTIONAL INPUTS:
;       minwave     - see ICALIBRATE
;       dwave       - see ICALIBRATE
;       stitchrange - set the YRANGE manually when DEBUG=1
;       scalefactor - manually set the scale factor to match the 2D
;                     spectra in the overlap region [NSPEC]
;       refindx     - reference index number (default 0)
;       outname     - output file name
;       object      - output object name
;       datapath    - I/O path name
;
; KEYWORD PARAMETERS:
;       debug - show a plot useful for debuggin
;       wfits - write out the stitched FITS file
;       gzip  - see WRT2DSPEC
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; COMMENTS:
;       Generalization required if the total exposure time is
;       different for each spectrum.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Jul 13, U of A - written, based on
;          SINGS_STITCH 
;       jm06jan24uofa - generalized to more than two 2D spectra; many
;                       of the input variables were renamed;
;                       simplified in the sense that user inputs are
;                       *required* 
;
; Copyright (C) 2005-2006, John Moustakas
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

pro sings2d_stitch, speclist, outname=outname, object=object, $
  user_pixshift=user_pixshift, wmapname=wmapname, minwave=minwave, $
  dwave=dwave, stitchrange=stitchrange, scalefactor=scalefactor, $
  refindx=refindx, datapath=datapath, debug=debug, wfits=wfits, gzip=gzip

    nspec = n_elements(speclist)
    if (nspec eq 0L) then begin
       print, 'Syntax - sings2d_stitch, speclist, outname=, object=, $'
       print, '   user_pixshift=, wmapname=, minwave=, dwave=, $'
       print, '   stitchrange=, scalefactor=, refindx=, datapath=, /debug, $'
       print, '   /wfits, /gzip'
       return
    endif

    if (n_elements(outname) eq 0L) then outname = 'sings2d_stitch.fits'
    if (n_elements(datapath) eq 0L) then datapath = cwd()

    if (n_elements(scalefactor) eq 0L) then begin
       scalefactor = fltarr(nspec)+1.0
       findscale = 1L
    endif else begin
       if (n_elements(scalefactor) ne nspec) then begin
          splog, 'Dimensions of SCALEFACTOR and SPEC1LIST do not agree.'
          return
       endif
       findscale = 0L
    endelse

    if (n_elements(wmapname) eq 0L) then begin
       splog, 'WMAPNAME must be provided.'
       return
    endif

    colorlist = ['red','blue','green','purple','cyan','yellow']
    if (nspec gt n_elements(colorlist)) then message, 'Number of colors exceeded!'

    if keyword_set(wfits) then debug = 0L
    
; forage header information and read the 2D spectra; assume that the
; spectra are the same dimensions
    
    forage = iforage(datapath+speclist)
    data = rd2dspec(speclist,datapath=datapath)
    header = *data[0].header

    if (n_elements(object) eq 0L) then object = forage[0].object
    
    ncols = data[0].naxis1
    nrows = data[0].naxis2

    if (n_elements(user_pixshift) eq 0L) then user_pixshift = fltarr(nspec-1)+nrows else begin
       if (n_elements(user_pixshift) ne nspec-1L) then begin
          splog, 'Dimensions of USER_PIXSHIFT must be NSPEC-1.'
          return
       endif
    endelse

    neg = where(user_pixshift lt 0.0,nneg)
    if (nneg ne 0L) then begin
       splog, 'Negative USER_PIXSHIFT not supported!'
       return
    endif
    
    delta_pixshift = nrows-user_pixshift ; overlap spacing between adjacent images
    bignrows = nrows*nspec-total(delta_pixshift)    
    bigrowaxis = lindgen(bignrows)
    rowaxis = lindgen(nrows)

    xstart = lindgen(nspec)*nrows
    for i = 1L, nspec-1L do xstart[i] = xstart[i] - total(delta_pixshift[0L:i-1L])
    xend = xstart+nrows-1L

;   niceprint, xstart, xend
    
; wavelength-calibrate each object before stitching them together 

    if keyword_set(wfits) then begin
       icalibrate, speclist, newdata, datapath=datapath, wmapname=wmapname, $
         sensname='', tracename='', minwave=minwave, dwave=dwave
    endif else newdata = data

    dims = size(newdata.image,/dimension)
    if (dims[0] ne ncols) or (dims[1] ne nrows) then begin
       splog, 'Dimensions changed!'
       return
    endif
    
    profile = total(newdata.image,1)

    if keyword_set(debug) then begin
       im_window, 0, xratio=0.9, yratio=0.6
       pagemaker, nx=3, ny=1, position=pos, /normal, xspace=0.0, $
         xmargin=[1.3,0.2], ymargin=[0.3,1.3]
       if (n_elements(stitchrange) eq 0L) then profrange = minmax(profile) else $
         profrange = stitchrange
       djs_plot, rowaxis+xstart[0], profile[*,0], xsty=3, ysty=3, xthick=2, $
         ythick=2, charsize=1.7, charthick=2.0, xtitle='Row Number', $
         ytitle='Flux', xrange=minmax(bigrowaxis), ps=10, position=pos[*,0], $
         yrange=profrange, color=colorlist[0]
       for i = 1L, nspec-1L do djs_oplot, rowaxis+xstart[i], $
         profile[*,i], color=colorlist[i], ps=10
    endif

; form the output arrays and co-add the spectra; construct the average
; image and sky spectrum; only divide by two in the overlapping rows

    imagecube = fltarr(ncols,bignrows,nspec)
    sigmacube = imagecube*0.0
    maskcube = imagecube*0.0
    wavecube = imagecube*0.0

    for i = 0L, nspec-1L do begin
       if (i ge 1L) then begin
          im1 = newdata[i-1L].image[*,nrows-delta_pixshift[i-1L]:nrows-1L]
          im2 = newdata[i].image[*,0L:delta_pixshift[i-1L]-1L]
          npix = (size(im1,/dimension))[1]
          imratio = im1/reverse(im2)
          if findscale then scalefactor[i] = median(imratio)
;         print, nrows-delta_pixshift[i-1L], nrows-1L, 0, delta_pixshift[i-1L]-1L, scalefactor[i]
          if keyword_set(debug) then begin
             im_window, 1, xratio=0.5, yratio=0.4
             p1 = total(im1,1)/(npix-1.0) & p2 = total(im2,1)/(npix-1.0)
             djs_plot, lindgen(npix), p1*scalefactor[i-1], xsty=3, ysty=3, xthick=2, $
               ythick=2, charsize=1.7, charthick=2.0, xtitle='Row Number', line=0, $
               ytitle='Flux', ps=10, yrange=[min(p1*scalefactor[i-1])<min(p2*scalefactor[i]),$
               max(p1*scalefactor[i-1])>max(p2*scalefactor[i])]
             djs_oplot, lindgen(npix), p2*scalefactor[i], color='cyan', ps=10, line=0
;            djs_oplot, lindgen(npix), p1, ps=10, line=0
;            djs_oplot, lindgen(npix), p2, ps=10, line=2, color='cyan'
             if (nspec gt 2L) then cc = get_kbrd(1)
          endif
       endif
       imagecube[*,xstart[i]:xend[i],i] = newdata[i].image*scalefactor[i]
    endfor
    for i = 0L, nspec-1L do sigmacube[*,xstart[i]:xend[i],i] = newdata[i].sigmamap*scalefactor[i]
    for i = 0L, nspec-1L do maskcube[*,xstart[i]:xend[i],i] = newdata[i].mask
    for i = 0L, nspec-1L do wavecube[*,xstart[i]:xend[i],i] = newdata[i].wavemap

    denom = total(imagecube ne 0.0,3) ; denominator term

    outimage = total(imagecube,3)/denom
    outsigmap = sqrt(total(sigmacube^2,3)/denom)
    outmask = total(maskcube,3)
    outskyimage = outimage*0.0

    outwavemap = total(wavecube,3)/denom
    if keyword_set(wfits) then outwset = {func: newdata[0].func, $
      xmin: newdata[0].xmin, xmax: newdata[0].xmax, $
      coeff: rebin(newdata[0].coeff[*,0],2,bignrows)}
    
    if keyword_set(debug) then begin
       wset, 0
; scaled, shifted profiles
       djs_plot, rowaxis+xstart[0], total(newdata[0].image,1)*scalefactor[0], $
         xsty=3, ysty=3, xthick=2, ythick=2, charsize=1.7, charthick=2.0, $
         xtitle='Row Number', xrange=minmax(bigrowaxis), ps=10, position=pos[*,1], $
         yrange=profrange, ytickname=replicate(' ',10), /noerase, color=colorlist[0]
       for i = 1L, nspec-1L do djs_oplot, rowaxis+xstart[i], $
         total(newdata[i].image,1)*scalefactor[i], color=colorlist[i]
; combined image       
;      djs_oplot, bigrowaxis, total(outimage,1), color='green'
       djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=2, ythick=2, $
         charsize=1.7, charthick=2.0, xtitle='Row Number', $
         xrange=minmax(bigrowaxis), yrange=profrange, $
         position=pos[*,2], /noerase, ytickname=replicate(' ',10)
       oplot, bigrowaxis, total(outimage,1), ps=10;, color=djs_icolor('green')
;      oploterror, bigrowaxis, total(outimage,1), total(outsigmap,1), $
;        ps=10, color=djs_icolor('green'), errcolor=djs_icolor('green')
; wavemap
;      djs_plot, rowaxis, djs_median(wavemap,1), xsty=3, ysty=3, xthick=2, $
;        ythick=2, charsize=1.7, charthick=2.0, xtitle='Row Number', $
;        ytitle='Flux', xrange=minmax(bigrowaxis), ps=10, position=pos[*,0]
;      for i = 1L, nspec-1L do djs_oplot, rowaxis+user_pixshift[i], $
;        djs_median(wavemap,1), color=colorlist[i]
    endif

; update the header; assume that EXPTIME, DATE-OBS, SCANLEN, POSANGLE,
; CCD noise parameters and the spectrophotometric errors are the same;
; as the default, adopt the header with the minimum USER_PIXSHIFT;
; retain the time parameters of this spectrum (UT, UTMIDDLE, JD, and
; ST) and all the ISPEC2D keywords (IRAW, etc.); compute the mean
; AIRMASS, ZD, PARANGLE, RA, and DEC; compute the new reference pixel
; CRPIX2; add a history note and a STITCHED boolean keyword

    if (n_elements(refindx) eq 0L) then minshift = min(user_pixshift,refindx)
    outheader = *data[refindx].header

    mean_ra = repstr(strtrim(im_dec2hms(djs_mean(im_hms2dec(forage.ra))),2),' ',':')
    mean_dec = repstr(strtrim(im_dec2hms(djs_mean(im_hms2dec(forage.dec))),2),' ',':')
    mean_airmass = djs_mean(forage.airmass)   ; mean airmass
    mean_zd = djs_mean(forage.zd)             ; zenith angle
    mean_parangle = djs_mean(forage.parangle) ; mean parallactic angle
    mean_jd = djs_mean(forage.jd)             ; mean Julian date
    mean_mjd = djs_mean(forage.mjd)           ; mean modified Julian date
    
    sxaddpar, outheader, 'OBJECT', object
    sxaddpar, outheader, 'NAXIS2', long(bignrows)
    sxaddpar, outheader, 'CRPIX2', long(bignrows/2.0+1.0), ' reference pixel number'
    sxaddpar, outheader, 'RA', mean_ra, ' right ascension [HMS]'
    sxaddpar, outheader, 'DEC', mean_dec, ' declination [DMS]', after='RA'
    sxaddpar, outheader, 'AIRMASS', float(mean_airmass), ' airmass at mid-exposure', format='(F12.4)'
    sxaddpar, outheader, 'ZD', float(mean_zd), after='AIRMASS', $
      ' zenith distance at mid-exposure [degrees]', format='(F12.2)'
    sxaddpar, outheader, 'PARANGLE', float(mean_parangle), after='AIRMASS', $
      ' parallactic angle at mid-exposure [degrees]', format='(F12.2)'

    for i = 0L, nspec-1L do sxaddpar, outheader, 'SSCALE'+string(i+1L,format='(I2.2)'), $
      float(scalefactor[i]), ' SINGS2D_STITCH scale factor', before='HISTORY'
    sxaddpar, outheader, 'STITCHED', 1L, before='HISTORY', ' stitched spectrum'

; wavelength-calibration parameters    

    if keyword_set(wfits) then begin
       
       sxaddpar, outheader, 'CRVAL1', float(minwave), ' wavelength at CRPIX1', before='HISTORY'
       sxaddpar, outheader, 'CRPIX1', float(1.0), ' reference pixel number', before='HISTORY'
       sxaddpar, outheader, 'CD1_1', float(dwave), ' dispersion [Angstrom/pixel]', before='HISTORY'
       sxaddpar, outheader, 'CDELT1', float(dwave), ' dispersion [Angstrom/pixel]', before='HISTORY'
       sxaddpar, outheader, 'CTYPE1', 'LINEAR', ' projection type', before='HISTORY'
       sxaddpar, outheader, 'WMAPNAME', wmapname, ' wavelength map name', before='HISTORY'
       sxaddhist, "'Linearized wavelength scale applied "+im_today()+"'", header

    endif
       
    sxaddhist, "'Images "+strjoin(repstr(repstr(repstr(forage.file,'.fits.gz',''),$
      '.fits',''),'cra.',''),',')+" stitched together "+im_today()+"'", outheader

; clean up memory
    
    icleanup, data
    icleanup, newdata

; write out

    if keyword_set(wfits) then begin

       splog, 'Writing '+datapath+outname+'.'
       wrt2dspec, outname, outimage, outsigmap, outmask, outheader, $
         skyimage=outskyimage, wset=outwset, telluric_spec=telluric_spec, $
         telluric_head=telluric_head, datapath=datapath, gzip=gzip

    endif

return
end    
