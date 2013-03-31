function get_zpt, hdr
; get the image zeropoint
    photzpt = sxpar(hdr,'PHOTZPT')
    photflam = sxpar(hdr,'PHOTFLAM')
    photplam = sxpar(hdr,'PHOTPLAM')
    zpt = photzpt - 2.5*alog10(photflam) - 5.0*alog10(photplam/5475.4)
return, zpt
end

pro bcgimf_phot, bcgmodel=bcgmodel, clobber=clobber, debug=debug
; jm13mar15siena
    
    bcgimfpath = getenv('BCGIMF_DATA')+'/'
    mosaicpath = clash_path('abell_2261',/mosaics,mas30=keyword_set(bcgmodel))
    modelpath = clash_path('abell_2261',/bcgmodels,mas30=keyword_set(bcgmodel))

    filt = bcgimf_filterlist(short=short,instr=instr,weff=weff)
    nfilt = n_elements(filt)
    
    ra_bcg = 15D*hms2dec('17:22:27.18')
    dec_bcg = hms2dec('32:07:57.30')
    glactc, ra_bcg, dec_bcg, 2000.0, gl, gb, 1, /deg
    ebv = dust_getval(gl,gb,/interp,/noloop)
    kl = k_lambda(weff,/odon)

    if keyword_set(bcgmodel) then begin
       scale = '030mas'
       pixscale = 0.030D        ; [arcsec/pixel]
    endif else begin
       scale = '065mas'
       pixscale = 0.065D        ; [arcsec/pixel]
    endelse

; read the F606W SE catalog to find and mask the foreground stars
; using the segmentation map
    if keyword_set(bcgmodel) eq 0 then begin
       cat = read_clash_catalog('abell_2261',path=sexpath)
       segm = mrdfits(sexpath+'SExtractor/detectionImage_SEGM.fits.gz',$
         0,seghdr,/silent)

       istar = where(cat.stel gt 0.9 and cat.f606w_mag lt 28.5,nstar)
       splog, 'Masking '+strtrim(nstar,2)+' stars...'
       if nstar ne 0 then begin
          starmask = abs(segm*0.0)
          for is = 0, nstar-1 do starmask += (segm eq cat[istar[is]].id)
          starmask = starmask eq 0 ; invert the mask [1=good,0=bad]
       endif
    endif 
    
; innermost radius: 3 pixels = 0.195 arcsec
; outermost radius: 35 arcsec = 540 pixels
    nring = 100
    rmin = 0.1D ; [arcsec]
;   rmin = 3D*pixscale ; [arcsec]
    rmax = 25D  ; [arcsec]
    
    rad = range(rmin,rmax,nring+1,/asinh)/pixscale ; [pixels]
    rin = rad[0:n_elements(rad)-2]  ; [pixels]
    rout = rad[1:n_elements(rad)-1] ; [pixels]

; initialize the output structure
    phot = {rin: 0.0, rout: 0.0, radius: 0.0, sb: fltarr(nfilt), $
      maggies: fltarr(nfilt), ivarmaggies: fltarr(nfilt)}
    phot = replicate(phot,nring)
    phot.rin = rin
    phot.rout = rout
    phot.radius = rout*pixscale ; [arcsec]
    
; loop on each band and then on each ring and do the photometry
;   for ib = 4, 4 do begin
    for ib = 0, nfilt-1 do begin
;      splog, 'Working on band '+short[ib]

; use Marc's model of the BCG       
       if keyword_set(bcgmodel) then begin
          image = mrdfits(modelpath+'a2261_mosaic_'+scale+'_'+strtrim(instr[ib],2)+$
            '_'+strtrim(short[ib],2)+'_cln_drz_model.fits.gz',0,hdr,/silent)
;         image = mrdfits(modelpath+'a2261_mosaic_'+scale+'_'+strtrim(instr[ib],2)+$
;           '_'+strtrim(short[ib],2)+'_cln_drz_????????_model.fits.gz',0,hdr,/silent)
          starmask = image*0+1
       endif else begin
          image = mrdfits(mosaicpath+'a2261_mosaic_'+scale+'_'+strtrim(instr[ib],2)+$
            '_'+strtrim(short[ib],2)+'_drz_????????.fits.gz',0,hdr,/silent)
       endelse
       inf = where(finite(image) eq 0)
       if inf[0] ne -1L then image[inf] = min(image)

       invvar = mrdfits(mosaicpath+'a2261_mosaic_'+scale+'_'+strtrim(instr[ib],2)+$
         '_'+strtrim(short[ib],2)+'_wht_????????.fits.gz',/silent)
       invvar = invvar*starmask
       inf = where(finite(invvar) eq 0)
       if inf[0] ne -1L then invvar[inf] = 0

       sz = size(image,/dim)
       zpt = get_zpt(hdr)-kl[ib]*ebv ; extinction-corrected (matches Coe's numbers)

       extast, hdr, astr
       ad2xy, ra_bcg, dec_bcg, astr, xcen, ycen

       dispimage = im_imgscl((image*(invvar gt 0))[(xcen-rmax/pixscale):(xcen+rmax/pixscale),$
         (ycen-rmax/pixscale):(ycen+rmax/pixscale)],/neg,/sqrroot)

       if keyword_set(debug) then begin
          if !d.window ne 0 then window, 0
          plotimage, dispimage, /noaxes, /preserve_aspect, /norm
       endif
       sz2 = size(dispimage,/dim)
;      tvcircle, 6.0, sz2[0]/2, sz2[1]/2, color=im_color('firebrick') ; innermost lens

       for ir = 0L, nring-1 do begin
          if keyword_set(debug) then begin
             if (ir mod 10) eq 0 then tvcircle, rout[ir], $
               sz2[0]/2, sz2[1]/2, color=im_color('orange')
          endif

;         djs_photfrac, xcen, ycen, [rin[ir],rout[ir]], fracs=fracs, $
          djs_photfrac, xcen, ycen, [0D,rout[ir]], fracs=fracs, $
            pixnum=pix, xdimen=sz[0], ydimen=sz[1], xpix=xpix, ypix=ypix

          area = (total(fracs*(invvar[pix] gt 0))*pixscale)^2 ; [arcsec^2]
          counts = total(image[pix]*(invvar[pix] gt 0)*fracs,/nan)
          var = 1.0/(invvar[pix]+(invvar[pix] eq 0))*(invvar[pix] ne 0)
          sigma2 = total(fracs*var,/nan)

;         counts2 = djs_phot(xcen,ycen,rout[ir],skyrad,image,$
;           invvar,calg='none',salg='none',flerr=flerr)
;         print, counts, counts2

          maggies = counts*10.0^(-0.4*zpt) ; AB maggies
          errmaggies = sqrt(sigma2)*10.0^(-0.4*zpt)
          ivarmaggies = 1.0/errmaggies^2

          phot[ir].maggies[ib] = maggies
          phot[ir].ivarmaggies[ib] = ivarmaggies
          phot[ir].sb[ib] = -2.5*alog10(counts)+zpt-2.5*alog10(area)      ; SB [mag/arcsec^2]
       endfor
;      window, 2
;      plot, phot.radius, phot.sb[ib], psym=8, ysty=3, /xlog
    endfor 

;   mag = maggies2mag(phot.maggies[0],ivarmaggies=phot.ivarmaggies[0],magerr=magerr)
;   plot, phot.radius, -2.5*alog10(phot.maggies[3]/phot.maggies[4]), psym=8 ; V-R

; write out
    if keyword_set(bcgmodel) then $
      outfile = bcgimfpath+'bcgimf_phot_bcgmodel.fits' else $
        outfile = bcgimfpath+'bcgimf_phot.fits'
    im_mwrfits, phot, outfile, clobber=clobber

stop    
    
return
end
    
