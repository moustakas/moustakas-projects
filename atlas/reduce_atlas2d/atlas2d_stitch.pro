;+
; NAME:
;       ATLAS2D_STITCH
;
; PURPOSE:
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; COMMENTS:
;       Only works with two images so far.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Jul 13, U of A - written, based on
;          SINGS_STITCH 
;
;-

pro atlas2d_stitch, spec1list, spec2list, stitchlist, objectlist, $
  nlags=nlags, user_pixshift=user_pixshift, datapath=datapath, $
  debug=debug, wfits=wfits

    ngalaxy = n_elements(spec1list)
    if (n_elements(spec2list) ne ngalaxy) or (n_elements(stitchlist) ne ngalaxy) then begin
       print, 'Syntax - sings_stitch, spec1list, spec2list, stitchlist, $'
       print, '   nlags=, user_pixshift=, datapath=, /debug, /wfits'
       return
    endif

    if (n_elements(datapath) eq 0L) then datapath = cwd()

    if (n_elements(user_pixshift) eq 0L) then user_pixshift = fltarr(ngalaxy)-999.0 else begin
       if (n_elements(user_pixshift) ne ngalaxy) then begin
          splog, 'Dimensions of USER_PIXSHIFT and SPEC1LIST do not agree.'
          return
       endif
    endelse
    
    if (n_elements(nlags) eq 0L) then nlags = 101L

; call this routine recursively    
    
    if (ngalaxy gt 1L) then begin
       for igalaxy = 0L, ngalaxy-1L do $
         atlas2d_stitch, spec1list[igalaxy], spec2list[igalaxy], stitchlist[igalaxy], $
           objectlist[igalaxy], nlags=nlags, user_pixshift=user_pixshift[igalaxy], $
           datapath=datapath, debug=debug, wfits=wfits
       return
    endif

; forage header information and read the 2D spectra; assume the
; spectra have the same wavelength solution and telluric spectrum 
    
    forage = iforage([spec1list,spec2list],datapath=datapath)
    spec1 = rd2dspec(spec1list,datapath=datapath,wset=wset1)
    spec2 = rd2dspec(spec2list,datapath=datapath,wset=wset2)

    telluric_spec = spec1.telluric_spec
    telluric_head = spec1.telluric_head

    ncols = spec1.naxis1
    nrows = spec1.naxis2
    bignrows = 2*nrows
    bigrowaxis = lindgen(bignrows)
    rowaxis = lindgen(nrows)

; cross-correlate the normalized spatial profiles to determine the
; pixel shift, unless USER_PIXSHIFT has been specified; shift SPEC2 to
; match SPEC1 
    
    bigrefspec = bigrowaxis*0.0       ; reference spectrum
    refspec = total(spec2.image,1)
    bigrefspec[0L:nrows-1L] = im_normalize(refspec,rowaxis,const=specnorm,/median)

    bigspec = bigrowaxis*0.0
    spec = total(spec1.image,1)
    bigspec[nrows:bignrows-1L] = spec/specnorm

; cross-correlate    

    if (user_pixshift le -999.0) then begin
    
       lags = -reverse(lindgen(nlags))
       result = djs_correlate(bigspec,bigrefspec,lags)
       maxcorr = max(result,bestlag)
       pixshift = lags[bestlag]
       splog, 'Best pixel shift = '+string(pixshift,format='(I0)')

       if keyword_set(debug) then begin

          im_window, 0, xratio=0.4, /square
          djs_plot, lags, result, xsty=3, ysty=3, xthick=2, ythick=2, charsize=2.0, $
            charthick=2.0, xtitle='Pixel Shift', ytitle='Cross-Correlation'
          oplot, pixshift*[1,1], !y.crange, line=2, thick=2

       endif
       
    endif else pixshift = user_pixshift
       
; form the output arrays and co-add the spectra; construct the average
; image and sky spectrum; only divide by two in the overlapping rows

    outnrows = bignrows + pixshift

    outimage1 = fltarr(ncols,outnrows)
    outimage2 = outimage1*0.0
    outsigmap1 = outimage1*0.0
    outsigmap2 = outimage1*0.0
    outskyimage1 = outimage1*0.0
    outskyimage2 = outimage1*0.0
    outwavemap1 = outimage1*0D0 ; must be double!
    outwavemap2 = outimage1*0D0
    outmask1 = outimage1*0B
    outmask2 = outimage1*0B

    outimage1[*,0L:nrows-1L] = spec2.image
    outimage2[*,nrows+pixshift:outnrows-1L] = spec1.image

    outskyimage1[*,0L:nrows-1L] = spec2.sky
    outskyimage2[*,nrows+pixshift:outnrows-1L] = spec1.sky

    outwavemap1[*,0L:nrows-1L] = spec2.wavemap
    outwavemap2[*,nrows+pixshift:outnrows-1L] = spec1.wavemap

    denom = outimage1*0.0+1.0
    denom[where((abs(outimage1) gt 0.0) and (abs(outimage2) gt 0.0))] = 2.0
    
    outimage = float(total([ [ [outimage1] ], [ [outimage2] ] ],3)/denom)
    outskyimage = float(total([ [ [outskyimage1] ], [ [outskyimage2] ] ],3)/denom)

; error map
    
    outsigmap1[*,0L:nrows-1L] = spec2.sigmamap
    outsigmap2[*,nrows+pixshift:outnrows-1L] = spec1.sigmamap
    outsigmap = float(sqrt(total([ [ [outsigmap1^2.0] ], [ [outsigmap2^2.0] ] ],3)/denom))

; bad pixel mask
    
    outmask1[*,0L:nrows-1L] = spec2.mask
    outmask2[*,nrows+pixshift:outnrows-1L] = spec1.mask
    outmask = fix(total([ [ [outmask1] ], [ [outmask2] ] ],3))

; wavelength map

    outwavemap = total([ [ [outwavemap1] ], [ [outwavemap2] ] ],3) / denom
    xxmap = rebin(findgen(ncols),ncols,outnrows)
    ncoeff = (size(wset1.coeff,/dimension))[0]
    xy2traceset, xxmap, outwavemap, outwset, ncoeff=ncoeff, /silent
    
; update the header; assume that EXPTIME, DATE-OBS, SCANLEN, POSANGLE,
; CCD noise parameters and the spectrophotometric errors are the same;
; retain the time parameters of the zeroth exposure (UT, UTMIDDLE, JD,
; and ST) and all the ISPEC2D keywords (IRAW, etc.); compute the mean
; AIRMASS, ZD, PARANGLE, RA, and DEC; compute the new reference pixel
; CRPIX2; add a history note
    
    outheader = *spec1.header

    mean_ra = repstr(strtrim(im_dec2hms(djs_mean(im_hms2dec(forage.ra))),2),' ',':')
    mean_dec = repstr(strtrim(im_dec2hms(djs_mean(im_hms2dec(forage.dec))),2),' ',':')
    mean_airmass = djs_mean(forage.airmass)   ; mean airmass
    mean_zd = djs_mean(forage.zd)             ; zenith angle
    mean_parangle = djs_mean(forage.parangle) ; mean parallactic angle
    mean_jd = djs_mean(forage.jd)             ; mean Julian date
    mean_mjd = djs_mean(forage.mjd)           ; mean modified Julian date
    
    sxaddpar, outheader, 'OBJECT', objectlist
    sxaddpar, outheader, 'CRPIX2', outnrows/2L+1L, ' reference pixel number'
    sxaddpar, header, 'RA', mean_ra, ' right ascension [HMS]'
    sxaddpar, header, 'DEC', mean_dec, ' declination [DMS]', after='RA'
    sxaddpar, header, 'AIRMASS', float(mean_airmass), ' airmass at mid-exposure', format='(F12.4)'
    sxaddpar, header, 'ZD', float(mean_zd), after='AIRMASS', $
      ' zenith distance at mid-exposure [degrees]', format='(F12.2)'
    sxaddpar, header, 'PARANGLE', float(mean_parangle), after='AIRMASS', $
      ' parallactic angle at mid-exposure [degrees]', format='(F12.2)'

    sxaddhist, "'Images "+strjoin(repstr(forage.file,'.fits.gz',''),' and ')+$
      " stitched "+im_today()+"'", outheader

; write out

    if keyword_set(wfits) then begin

       splog, 'Writing '+datapath+stitchlist+'.'

       wrt2dspec, stitchlist, outimage, outsigmap, outmask, outheader, $
         skyimage=outskyimage, wset=outwset, telluric_spec=telluric_spec, $
         telluric_head=telluric_head, datapath=datapath, gzip=gzip

    endif

; debugging plots    
    
    if keyword_set(debug) then begin

       im_window, 2, /square

       profile = djs_median(outimage[(ncols/2L-10L)>0L:(ncols/2L+10L)<(ncols-1L),*],1)
       profile1 = djs_median(outimage1[(ncols/2L-10L)>0L:(ncols/2L+10L)<(ncols-1L),*],1)
       profile2 = djs_median(outimage2[(ncols/2L-10L)>0L:(ncols/2L+10L)<(ncols-1L),*],1)

       good1 = where(profile1 gt 0.0)
       good2 = where(profile2 gt 0.0)

       yrange = minmax(profile)
;      yrange = median(profile)+djsig(profile)*[-1,5]
       yrange[0] = yrange[0] > 0
       
       scale = 1E15
       djs_plot, lindgen(outnrows), scale*profile, ps=10, xsty=3, ysty=3, thick=2.0, $
         xthick=2, ythick=2, charsize=1.8, charthick=2.0, xtitle='Row', $
         ytitle='Median Spatial profile [10^{-15} '+flam_units()+']', $
         yrange=scale*yrange
       djs_oplot, (lindgen(outnrows))[good1], scale*profile1[good1], ps=10, thick=1.0, color='red'
       djs_oplot, (lindgen(outnrows))[good2], scale*profile2[good2], ps=10, thick=1.0, color='yellow'

       legend, ['Combined','Spectrum 1','Spectrum 2'], /right, /top, box=0, $
         charthick=2.0, charsize=1.5, textcolor=djs_icolor(['white','red','yellow'])

       splog, 'Press any key to continue.'
       cc = get_kbrd(1)
       
    endif

; clean up memory
    
    icleanup, spec1 
    icleanup, spec2

return
end    
