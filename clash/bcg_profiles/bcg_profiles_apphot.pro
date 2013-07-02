pro bcg_profiles_apphot
; jm13jul01siena - perform circular aperture photometry on the *model*
; BCGs 
    
    clash = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    ncl = n_elements(clash)
    
    filt = bcgimf_filterlist(short=short,weff=weff,zpt=zpt,instr=instr)
;   srt = reverse(sort(weff))
;   filt = filt[srt]
;   short = short[srt]
;   weff = weff[srt]
;   zpt = zpt[srt]
;   instr = instr[srt]
    
    kl = k_lambda(weff,/odon)
    nfilt = n_elements(filt)

; establish the radii and initialize the output structure 
    nring = 100
    rmin = 0.1D ; [kpc]
    rmax = 200D ; [kpc]
    radius = range(rmin,rmax,nring,/asinh) ; [kpc]

    inrad = radius-(radius-shift(radius,1))/2.0 ; inner radius [kpc]
    outrad = radius+(shift(radius,-1)-radius)/2 ; outer radius [kpc]
    inrad[0] = 0D
    outrad[nring-1] = radius[nring-1]+(radius[nring-1]-radius[nring-2])/2.0

    phot = replicate({radius: 0.0, sb: fltarr(nfilt), $
      maggies: fltarr(nfilt), ivarmaggies: fltarr(nfilt)},nring)
    phot.radius = radius
    
 ; loop on each cluster    
    for ic = 0, ncl-1 do begin
       cluster = strtrim(clash[ic].shortname,2)
       this =  where(cluster[ic] eq strtrim(clash.shortname,2))
       arcsec2kpc = dangular(clash[this].z,/kpc)/206265D ; [kpc/arcsec]
       ebv = clash[this].ebv
       
       if cluster eq 'a2261' then begin
          scale = '030mas'
          pixscale = 0.030D     ; [arcsec/pixel]
          mas30 = 1
       endif else begin
          scale = '065mas'
          pixscale = 0.065D     ; [arcsec/pixel]
          mas30 = 0
       endelse
       mosaicpath = clash_path(clash[ic].dirname,/mosaics,mas30=mas30)
       modelpath = clash_path(clash[ic].dirname,/bcgmodels,mas30=mas30)

; for GAIN see notes at
; http://www.ifa.hawaii.edu/~rgal/science/sextractor_notes.html
       if file_test(modelpath,/dir) then begin
          splog, 'Working on cluster '+cluster
          delvarx, xcen_lum, ycen_lum
          for ib = nfilt-1, 0, -1 do begin
;         for ib = 0, nfilt-1 do begin
             imfile = modelpath+cluster+'_mosaic_'+scale+'_'+strtrim(instr[ib],2)+$
               '_'+strtrim(short[ib],2)+'_drz_????????_BCG.fits.gz'
             weightfile = mosaicpath+cluster+'_mosaic_'+scale+'_'+strtrim(instr[ib],2)+$
               '_'+strtrim(short[ib],2)+'_wht_????????.fits.gz'
             if file_test(imfile) and file_test(weightfile) then begin
                image = mrdfits(imfile,0,hdr)
                invvar = mrdfits(weightfile,0,ivarhdr)
                sz = size(image,/dim)
                exptime = sxpar(hdr,'EXPTIME') ; [sec]

; find the BCG but only in the reference filter, F160W
                if short[ib] eq 'f160w' then begin
                   splog, 'Searching for the BCG...'
                   find_galaxy, image, size, ellipticity, posangle, xcen, ycen, $
                     xcen_lum, ycen_lum, fraction=fraction, index=index, level=level, $
                     nblob=1, plot=plot ;, /quiet
                   extast, hdr, astr
                   xy2ad, xcen_lum, ycen_lum, astr, ra_bcg, dec_bcg
                endif
                
; do the photometry                
                rin = inrad/pixscale/arcsec2kpc ; [pixels]
                rout = outrad/pixscale/arcsec2kpc ; [pixels]

                for ir = 0L, nring-1 do begin
                   djs_photfrac, xcen_lum, ycen_lum, [rin[ir],rout[ir]], fracs=fracs, $
                     pixnum=pix, xdimen=sz[0], ydimen=sz[1], xpix=xpix, ypix=ypix

                   area = (total(fracs*(invvar[pix] gt 0))*pixscale)^2 ; [arcsec^2]
                   counts = total(image[pix]*(invvar[pix] gt 0)*fracs,/nan)

; add the shot noise of the object in quadrature to the inverse
; variance map; poisson error = sqrt(electrons), converted back to
; native electrons/pixel/sec units
                   shotvar = image[pix]/exptime ; [electron/pixel/sec]
                   var = (shotvar + 1.0/(invvar[pix]+(invvar[pix] eq 0)))*(invvar[pix] ne 0)
;                  var = 1.0/(invvar[pix]+(invvar[pix] eq 0))*(invvar[pix] ne 0)
                   sigma2 = total(fracs*var,/nan)

                   fact = 10.0^(-0.4*(zpt[ib]-kl[ib]*ebv))
                   maggies = counts*fact          ; [AB maggies]
                   errmaggies = sqrt(sigma2)*fact ; [AB maggies]
                   ivarmaggies = 1.0/errmaggies^2
                   
                   phot[ir].maggies[ib] = maggies
                   phot[ir].ivarmaggies[ib] = ivarmaggies
                   phot[ir].sb[ib] = -2.5*alog10(counts)-2.5*alog10(area)+zpt[ib]-kl[ib]*ebv ; SB [mag/arcsec^2]
                endfor
stop                
             endif else splog, '  No photometry in band '+short[ib]
          endfor
       endif else splog, 'No BCG data for cluster '+cluster
    endfor 
    

return
end
    
