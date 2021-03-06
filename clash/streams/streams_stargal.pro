;+
; NAME:
;   dstargal_atlas
; PURPOSE:
;   identify stars and galaxies in NASA-Sloan Atlas image
; CALLING SEQUENCE:
;   dstargal_atlas 
; COMMENTS:
;   Requires dparents_atlas, dpsf_atlas to have been run.
;   Agressive about calling things stars (ie. tuned to large galaxies)
;   Outputs files in:
;     atlases/#/[base]-#-nimage.fits - star-subtracted images
;     atlases/#/[base]-#-sgset.fits - locations of stars, galaxies
; REVISION HISTORY:
;   11-Jan-2006  Written by Blanton, NYU
;-
;------------------------------------------------------------------------------
pro streams_stargal, base, plot=plot, gsmooth=gsmooth, glim=glim, $
  gsaddle=gsaddle, nsigma=nsigma, noclobber=noclobber, $
  psffiles=psffiles, maxnstar=maxnstar
    
    if(NOT keyword_set(gsmooth)) then gsmooth=10.
    if(NOT keyword_set(glim)) then glim=3.0 ; 25.
    if(NOT keyword_set(gsaddle)) then gsaddle=10.0 ; 50.
    if(NOT keyword_set(nsigma)) then nsigma=3.0 ; 20.
    if(NOT keyword_set(maxnstar)) then maxnstar=3000L
    if(NOT keyword_set(maxmem)) then maxmem=2.e+9
    
; default to use base name same as directory name
    if n_elements(base) eq 0 then begin
       spawn, 'pwd', cwd
       base=(file_basename(cwd))[0]
    endif

; read in pset
    pset= gz_mrdfits(base+'-pset.fits',1)
    imfiles=strtrim(pset.imfiles,2)
    puse=pset.puse
    nim= n_elements(imfiles)
    nx=lonarr(nim)
    ny=lonarr(nim)
    images=ptrarr(nim)
    nimages=ptrarr(nim)
    ivars=ptrarr(nim)
    hdrs=ptrarr(nim)
    psfs=ptrarr(nim)

    ref = pset.ref

; set up star and galaxy locations
    sgset1={base: base, $
      ref: pset.ref, $
      iparent: -1L, $
      nstars: 0L, $
      ra_stars: dblarr(maxnstar), $
      dec_stars: dblarr(maxnstar), $
      ngals: 0L, $
      ra_gals: dblarr(maxnstar), $
      dec_gals: dblarr(maxnstar) }

; read in the psf information
    if n_elements(psffiles) eq 0 then begin
       for k=0L, nim-1L do begin
          psfs[k]=ptr_new(gz_mrdfits(psffiles[k]))
       endfor
    endif else begin
       for k=0L, nim-1L do begin
          bimfile=(stregex(imfiles[k], '(.*)\.fits.*', /sub, /extr))[1]
          psfs[k]=ptr_new(gz_mrdfits(bimfile+'-bpsf.fits'))
       endfor
    endelse
    
; loop on each parent
    pim=gz_mrdfits(base+'-pimage.fits')
    pcat=gz_mrdfits(base+'-pcat.fits',1)

    subdir= 'atlases'
    spawn, 'mkdir -p '+subdir, /sh

    dhdr= dimage_hdr()
    
;   splog, 'HACK!'
;   for iparent = 50, 50 do begin
    for iparent = 0L, n_elements(pcat)-1 do begin
       splog, 'Parent ', iparent
       sgsetfile=subdir+'/'+strtrim(string(iparent),2)+'/'+base+'-'+ $
         strtrim(string(iparent),2)+'-sgset.fits'

       sgset = sgset1
       sgset.iparent=iparent

; setup output directory
       spawn, 'mkdir -p '+subdir+'/'+strtrim(string(iparent),2), /sh

       if(gz_file_test(sgsetfile) ne 0 AND $
         keyword_set(noclobber) ne 0) then begin
          nimfile=subdir+'/'+strtrim(string(iparent),2)+'/'+base+ $
            '-nimage-'+strtrim(string(iparent),2)+'.fits'
          if(gz_file_test(nimfile) eq 0) then begin
             message, 'sgset file exists, but not nimfile'
          endif
          return
       endif
    
; check if this parent is junk by computing the fraction of masked
; pixels; if greater than zero then skip it (usually because it's on
; the edge) and move on!
       ivarref = gz_mrdfits('parents/'+base+'-parent-'+ $
         strtrim(string(iparent),2)+'.fits',1+ref*2L,hdr,/silent)
       fracmask = total(ivarref eq 0)/cmproduct(size(ivarref,/dim))
       if fracmask gt 0.01 then begin
          im_mwrfits, sgset, sgsetfile, dhdr, /clobber
          continue
       endif

; read in the images
       for k=0L, nim-1L do begin
          images[k]=ptr_new(gz_mrdfits('parents/'+base+'-parent-'+ $
            strtrim(string(iparent),2)+'.fits',0+k*2L,hdr,/silent))
          hdrs[k]=ptr_new(hdr)
          nx[k]=(size(*images[k],/dim))[0]
          ny[k]=(size(*images[k],/dim))[1]
          ivars[k]=ptr_new(gz_mrdfits('parents/'+base+'-parent-'+ $
            strtrim(string(iparent),2)+'.fits',1+k*2L,/silent))
       endfor
    
; find and subtract stars in the reference image
       splog, 'Skipping PSF-subtraction of stars for now!'
       nstars = 0
       nimages[ref] = ptr_new(*images[ref])
;      adxy, *hdrs[ref], raex, decex, xex, yex
;      nimages[ref]= ptr_new(streams_psfsub(image=*images[ref], ivar=*ivars[ref], $
;        psf=*psfs[ref], x=xstar_r, y=ystar_r, $
;        flux=fluxstar_r, nstars=nstars, exx=xex, $
;        exy=yex, exr=1.5, nsigma=nsigma))

       if(nstars gt 0) then begin
          xyad, *hdrs[ref], xstar_r, ystar_r, ra_stars, dec_stars
          cirrange, ra_stars
          fluxstar= fltarr(n_elements(fluxstar_r), nim)
          fluxstar[*,ref]= fluxstar_r
       endif
    
; subtract same stars in all other images
       if (nstars gt 0) then begin
          for k=0L, nim-1L do begin
             if(k ne ref) then begin
                adxy, *hdrs[k], ra_stars, dec_stars, xstar, ystar
                nimages[k]= ptr_new(dpsfsub_atlas(image=*images[k], ivar=*ivars[k], $
                  psf=*psfs[k], x=xstar, y=ystar, $
                  flux=tmp_fluxstar))
                fluxstar[*,k]=tmp_fluxstar
             endif
          endfor
       endif else begin
          for k=0L, nim-1L do $
            if(k ne ref) then $
              nimages[k]= ptr_new(*images[k])
       endelse
       
; output subtracted images
       nimfile=subdir+'/'+strtrim(string(iparent),2)+'/'+base+ $
         '-nimage-'+strtrim(string(iparent),2)+'.fits'
       for k=0L, nim-1L do mwrfits, *nimages[k], nimfile, $
         *hdrs[k], create=(k eq 0), /silent
       spawn, 'gzip -f '+nimfile, /sh
    
; find galaxies in reference image
       simage=dsmooth(*nimages[ref], gsmooth)
       ssig=dsigma(simage, sp=long(gsmooth*5.))
       saddle=gsaddle*ssig
       dpeaks, simage, xc=xc, yc=yc, sigma=ssig, minpeak=glim*ssig, $
         npeaks=ngals, saddle=gsaddle, /refine, /check

;; if no galaxies found, let's just assume there is one at center
;       if(ngals eq 0) then begin
;          ngals=1L
;          ra_gals = pcat[iparent].cra
;          dec_gals = pcat[iparent].cdec
;          adxy, *hdrs[ref], ra_gals, dec_gals, xc, yc
;       endif
       
       if(ngals gt 0) then begin
          xgals=float(xc)
          ygals=float(yc)
          
; refine the centers
          drefine, *nimages[ref], xgals, ygals, smooth=2., $
            xr=xrgals, yr=yrgals, box=long(5)
          xyad, *hdrs[ref], xrgals, yrgals, ra_gals, dec_gals
          cirrange, ra_gals
       endif
       
; only store up to maxnstars stars
       if(nstars gt maxnstar) then begin
          nstars=maxnstar
          ra_stars= ra_stars[0L:nstars-1L]
          dec_stars= dec_stars[0L:nstars-1L]
       endif
    
; store locations in sgset
       sgset.nstars= nstars
       if(nstars gt 0) then begin
          sgset.ra_stars[0:nstars-1]= ra_stars
          sgset.dec_stars[0:nstars-1]= dec_stars
       endif
       sgset.ngals= ngals
       if(ngals gt 0) then begin
          sgset.ra_gals[0:ngals-1]= ra_gals
          sgset.dec_gals[0:ngals-1]= dec_gals
       endif 
;if iparent eq 54 then stop       
; get memory use
       maxpix= float(max(nx*ny))
       memuse= 4.*maxpix*float(ngals)
       if(memuse gt maxmem) then begin
          splog, 'LIMITING MEMORY USE'
          sgset.ngals= long(maxmem/(4.*maxpix))
          limitedmem=1
       endif
       
; output star and galaxy info
       if(keyword_set(limitedmem)) then $
         sxaddpar, dhdr, 'LIMMEM', 1, '1 if limited NGALS due to memory'
       im_mwrfits, sgset, sgsetfile, dhdr, /clobber
       
       if(keyword_set(plot)) then begin
          atv, *nimages[ref]
          if(ngals gt 0) then begin
             adxy, *hdrs[ref], sgset.ra_gals[0:ngals-1], sgset.dec_gals[0:ngals-1], $
               xg, yg
             atvplot, xg, yg, psym=4, color='blue', th=4
          endif
          if(nstars gt 0) then begin
             adxy, *hdrs[ref], sgset.ra_stars[0:nstars-1], sgset.dec_stars[0:nstars-1], $
               xs, ys
             atvplot, xs, ys, psym=4, color='red'
          endif
       endif
    endfor 
    
    heap_free, psfs
    heap_free, images
    heap_free, nimages
    heap_free, ivars
    heap_free, hdrs
    
return 
end
;------------------------------------------------------------------------------
