;+
; NAME:
;       ATLAS_LIGHTFRACTION
;
; PURPOSE:
;       Compute the fraction of galaxy light enclosed by our
;       integrated spectroscopic apertures by using the DSS-2 FITS
;       images.
;
; CALLING SEQUENCE:
;       atlas_lightfraction, object, /find, /noplot
;
; INPUTS:
;
; OPTIONAL INPUTS:
;       object - galaxy name
;
; KEYWORD PARAMETERS:
;       find   - identify stars in the image using FIND and
;                cross-match them with the GSC results 
;       noplot - do not generate a visualization of the apertures 
;       write  - write out a binary FITS file of the results 
;
; OUTPUTS:
;       A compressed binary fits file ATLAS_LIGHTFRACTION.FITS is
;       written if WRITE is set.
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       ATLAS_PATH(), MRDFITS(), MATCH_STRING(), VISUALIZE_ATLAS,
;       VISUALIZE_APERTURES, MMM, QUERYGSC(), FIND, GSSSXYAD,
;       IM_DJS_ANGLE_MATCH(), GSSSADXY, TVCIRCLE, DJS_ICOLOR(),
;       IM_OPLOT_BOX,  SPLOG, MWRFITS
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 April 2, U of A
;-

pro atlas_lightfraction, object, find=find, noplot=noplot, write=write

    atlaspath = atlas_path(/analysis1d) ; ATLAS path
    dsspath = atlas_path(/dss)          ; DSS FITS path

; read the atlas structure and only keep integrated spectra 
    
    atlas = mrdfits(atlaspath+'atlas_info.fits.gz',1,/silent)
    keep = where(atlas.drift eq 'Y',natlas)
    if natlas ne 0L then atlas = atlas[keep]
    
    if n_elements(object) ne 0L then begin

       doit = match_string(object,atlas.galaxy,index=match) ; match by galaxy name
       if match[0] eq -1L then begin
          splog, 'Object '+object+' not found!'
          return
       endif
       atlas = atlas[match]

    endif
    natlas = n_elements(atlas)

; initialize the output structure

    result= {$
      galaxy:     '', $
      lfraction: 0.0}
    result = replicate(result,natlas)

; set some FIND parameters

    if keyword_set(find) then begin
    
       sharplim = [0.2, 1.0]
       roundlim = [-1.0, 1.0]   ;*1.5
       findsig   = 50.0         ; DAO find sigma cut

    endif

; conduct aperture photometry on each image

    for k = 0L, natlas-1L do begin

; visualize the galaxy and the apertures

       galaxy = strcompress(atlas[k].galaxy,/remove)
       visualize_atlas, galaxy, dsspath=dsspath, visinfo=visinfo, $
         gsssastr=gsssastr, noplot=noplot, silent=write
       visualize_apertures, atlas[k], visinfo, apinfo=apinfo, /silent, $
         noplot=noplot

; compute the sky value in the image

       image = visinfo.image
       mmm, image, skymode, skysigma, skyskew

; query the HST-GSC and also FIND stars in the image; replace the
; stellar pixels with sky

       if keyword_set(find) then begin

          searchrad = visinfo.xcen*visinfo.pixscale/60.0 ; [arcmin]

          stars = querygsc(visinfo.crval,searchrad,magrange=[0,25])
          keep = where(stars.c eq 0L,nkeep) ; star flag
       
          if nkeep ne 0L then begin

             hmin = 3*djsig(image)*findsig
             find, image, imx, imy, flux, sharp, round, hmin, 3.0, $
               roundlim, sharplim, /silent
             gsssxyad, imx, imy, gsssastr, ima, imd

             stars = stars[keep]

; cross-match GSC stars and FIND stars             
             
             ntot = im_djs_angle_match(stars.ra,stars.dec,ima,imd,$
               dtheta=10.0/3600.0,units='degrees',mcount=mcount,$
               mindx=mindx,mdist=mdist)
             good = where(mindx ne -1L,ngood)

; convert from (RA,DEC) to (x,y)             

             gsssadxy, gsssastr, stars[good].ra, stars[good].dec, xstar, ystar 
;            xstar = (xstar-visinfo.xcen)*visinfo.pixscale ; convert to arcsec
;            ystar = (ystar-visinfo.ycen)*visinfo.pixscale 

             starrad = 10.0     ; star "radius" []

;            if not keyword_set(noplot) then tvcircle, starrad, $
;              xstar, ystar, /data, color=djs_icolor('red')

; replace "star" pixels with sky                

             im = image
             for j = 0L, ngood-1L do begin

                im_oplot_box, starrad, starrad, 0.0, xoffset=xstar[j], yoffset=ystar[j], $
                  color='red', corners=scorners, noplot=noplot
                indx = polyfillv(scorners[0,*],scorners[1,*],visinfo.xsize,visinfo.ysize)
                im[indx] = skymode
             
             endfor
             
          endif else begin

             splog, 'No GSC stars found in the field around '+galaxy+'.'
             im = image

          endelse
                   
       endif

; do rectangular photometry on the image

       wtrue = polyfillv(apinfo.rc3box[0,*],apinfo.rc3box[1,*],visinfo.xsize,visinfo.ysize)
       wspec = polyfillv(apinfo.specbox[0,*],apinfo.specbox[1,*],visinfo.xsize,visinfo.ysize)
       
       truecounts = total(image[wtrue]-skymode)
       speccounts = total(image[wspec]-skymode)

       result[k].galaxy = galaxy
       result[k].lfraction = (speccounts/truecounts) < 1.0
       
       splog, galaxy+': fraction = ', result[k].lfraction
       if (not keyword_set(noplot)) and (not keyword_set(write)) then $
         if (natlas gt 1L) then cc = get_kbrd(1)

    endfor 
    
; write a binary FITS file of the results

    if keyword_set(write) then begin

       outfile = 'atlas_lightfraction.fits'
       splog, 'Writing '+atlaspath+outfile+'.'
       mwrfits, result, atlaspath+outfile, /create
       spawn, ['gzip -f '+atlaspath+outfile], /sh
    
    endif
    
return
end 
