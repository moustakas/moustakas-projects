pro bcg_profiles_apphot, debug=debug
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
    rmin = 0.1D ; [kpc]
    rmax = 150D ; [kpc]
    nring = round(12.0*alog10(rmax-rmin)) ; 10 isophotes per decade with 1.1 spacing
    radius = range(rmin,rmax,nring,/log) ; [kpc]
;   radius = range(rmin,rmax,nring,/asinh) ; [kpc]

    inrad = radius-(radius-shift(radius,1))/2.0 ; inner radius [kpc]
    outrad = radius+(shift(radius,-1)-radius)/2 ; outer radius [kpc]
    inrad[0] = 0D
    outrad[nring-1] = radius[nring-1]+(radius[nring-1]-radius[nring-2])/2.0

    phot = replicate({radius: 0.0, sb: fltarr(nfilt), $
      maggies: fltarr(nfilt), ivarmaggies: fltarr(nfilt)},nring)
    phot.radius = radius

    info = replicate({cluster: '', z: 0.0, ra_bcg: 0D, dec_bcg: 0D, $
      ellipticity: 0.0, posangle: 0.0, size: 0.0, size_kpc: 0.0},ncl)
    info.cluster = clash.shortname
    info.z = clash.z
    
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
             modelfile = modelpath+cluster+'_mosaic_'+scale+'_'+strtrim(instr[ib],2)+$
               '_'+strtrim(short[ib],2)+'_drz_????????_BCG.fits.gz'
             imfile = mosaicpath+cluster+'_mosaic_'+scale+'_'+strtrim(instr[ib],2)+$
               '_'+strtrim(short[ib],2)+'_drz_????????.fits.gz'
             weightfile = mosaicpath+cluster+'_mosaic_'+scale+'_'+strtrim(instr[ib],2)+$
               '_'+strtrim(short[ib],2)+'_wht_????????.fits.gz'
             if file_test(imfile) and file_test(weightfile) then begin
                model = mrdfits(modelfile,0)
                image = mrdfits(imfile,0,hdr)
                invvar = mrdfits(weightfile,0,ivarhdr)
                extast, hdr, astr
                exptime = sxpar(hdr,'EXPTIME') ; [sec]
                sz = size(image,/dim)

; add the shot noise of the object in quadrature to the inverse
; variance map; poisson error = sqrt(electrons), converted back to
; native electrons/pixel/sec units
                shotvar = image/exptime ; [electron/pixel/sec]
                totalvar = (shotvar + 1.0/(invvar+(invvar eq 0)))*(invvar ne 0)
                totalivar = 1.0/(totalvar+(totalvar eq 0))*(totalvar ne 0)

; conversion from electrons/second/pixel to AB mag                
                zpt1 = zpt[ib]-kl[ib]*ebv
                fact = 10.0^(-0.4*zpt1)
                
; find the BCG but only in the reference filter, F160W
                if short[ib] eq 'f160w' then begin
                   splog, 'Searching for the BCG...'
                   wid = floor(max(radius)*0.5/pixscale/arcsec2kpc)
                   x0 = sz[0]/2-wid/2 & x1 = sz[0]/2+wid/2
                   y0 = sz[1]/2-wid/2 & y1 = sz[1]/2+wid/2
                   model1 = model[x0:x1,y0:y1]
                   find_galaxy, model1, size, ellipticity, posangle, xcen, ycen, $
                     xcen_lum, ycen_lum, fraction=fraction, index=index, level=level
                   xcen += x0 & xcen_lum += x0
                   ycen += y0 & ycen_lum += y0
                   xy2ad, xcen_lum, ycen_lum, astr, ra_bcg, dec_bcg
                   info[ic].ra_bcg = ra_bcg
                   info[ic].dec_bcg = dec_bcg
                   info[ic].ellipticity = ellipticity
                   info[ic].posangle = posangle
                   info[ic].size = size*pixscale ; [arcsec]
                   info[ic].size_kpc = info[ic].size*arcsec2kpc
                endif else begin
                   ad2xy, info[ic].ra_bcg, info[ic].dec_bcg, $
                     astr, xcen_lum, ycen_lum
                endelse

; elliptical photometry; ellipticity = 1-b/a
;               xcen_lum1 = xcen_lum
;               ycen_lum1 = ycen_lum
;               sectors_photometry, image, info[ic].ellipticity, info[ic].posangle, $
;                 xcen_lum1, ycen_lum1, ellradius, phi, ellcounts, ellsigmacounts, n_sectors=1, $
;                 sector_width=90, badpixels=badpixels, minlevel=minlevel

;; initialize the output photometry structure                
;                if short[ib] eq 'f160w' then begin
;                   ellphot = replicate({radius: 0.0, sb: fltarr(nfilt), maggies: fltarr(nfilt), $
;                     ivarmaggies: fltarr(nfilt)},n_elements(ellradius))
;                   ellphot.radius = ellradius*pixscale*arcsec2kpc ; [kpc]
;                endif

;               ellphot.maggies[ib] = ellcounts*fact ; [AB maggies]
;               phot[ir].ivarmaggies[ib] = ivarmaggies
;               ellphot.sb[ib] = -2.5*alog10(ellcounts)+zpt[ib]-kl[ib]*ebv ; [mag/arcsec^2]
;               ellphot.sb[ib] = -2.5*alog10(ellcounts)+5*alog10(pixscale)+zpt[ib]-kl[ib]*ebv ; [mag/arcsec^2]

; circular and elliptical aperture photometry
                wid = max(radius)/pixscale/arcsec2kpc
                dispimage = alog10(model[xcen-wid:xcen+wid,ycen-wid:ycen+wid])
;               dispimage = image[xcen-wid:xcen+wid,ycen-wid:ycen+wid]
                sz2 = size(dispimage,/dim)
                if keyword_set(debug) then begin
                   if !d.window ne 0 then window, 0
;                  plotimage, im_imgscl(dispimage,/neg,/sqr), /noaxes, $
;                    /preserve_aspect, /norm
                   cgimage, dispimage, stretch=2, /keep_aspect, clip=2.0, /save
                endif
                
                rin = inrad/pixscale/arcsec2kpc   ; [pixels]
                rout = outrad/pixscale/arcsec2kpc ; [pixels]
;               rout = rin*1.10
                
                dist_ellipse, ellrad, sz, xcen_lum, ycen_lum, $
                  1.0-info[ic].ellipticity, info[ic].posangle
                
                for ir = 0, nring-1 do begin
                   if keyword_set(debug) then begin
                      if (ir mod 5) eq 0 then begin
                         tvellipse, rout[ir], rout[ir]*(1-info[ic].ellipticity), sz2[0]/2, $
                           sz2[1]/2, info[ic].posangle+90, color=im_color('orange'), /data
                         tvellipse, rin[ir], rin[ir]*(1-info[ic].ellipticity), sz2[0]/2, $
                           sz2[1]/2, info[ic].posangle+90, color=im_color('orange'), /data
                      endif
                   endif

                   pix = where(ellrad ge rin[ir] and ellrad lt rout[ir],npix)
                   totnpix = (npix-total(totalivar[pix] eq 0))
;                  area = totnpix*pixscale^2 ; [arcsec^2]
                   area = (rout[ir]-rin[ir])*!pi*(1-info[ic].ellipticity)
                   
                   counts = total(model[pix]*(totalivar[pix] gt 0))
                   errcounts = sqrt(total(totalvar[pix]*(totalivar[pix] gt 0)))
                   print, rin[ir], rout[ir], counts, area
                   
;                  djs_photfrac, xcen_lum, ycen_lum, [rin[ir],rout[ir]], fracs=fracs, $
;                    pixnum=pix, xdimen=sz[0], ydimen=sz[1], xpix=xpix, ypix=ypix
;                  area = (total(fracs*(totalivar[pix] gt 0))*pixscale)^2 ; [arcsec^2]
;                  counts = total(image[pix]*(totalivar[pix] gt 0)*fracs,/nan)
;                  errcounts = sqrt(total(fracs*totalvar[pix]*(totalivar[pix] gt 0)))

                   phot[ir].maggies[ib] = counts*fact ; [AB maggies]
                   phot[ir].ivarmaggies[ib] = 1.0/(errcounts*fact)^2
                   phot[ir].sb[ib] = -2.5*alog10(counts/totnpix)+5.0*$ ; [mag/arcsec^2]
                     alog10(pixscale)+zpt1 
;if ir eq 10 then stop
                endfor
;               jj = read_bcg_profiles(cluster,these_filters=filt[ib])
;               djs_plot, phot.radius, phot.sb[ib], psym=8, yr=[25,15], ysty=3, /xlog
;               djs_oplot, jj.sma, jj.mu, psym=8, color='red'

                mag = maggies2mag(phot.maggies[ib],ivarmaggies=phot.ivarmaggies[ib],$
                  magerr=magerr)
                ploterror, phot.radius, mag, magerr, xsty=3, ysty=3, /xlog, psym=8, /trad
                
                

                
stop                
             endif else splog, '  No photometry in band '+short[ib]
          endfor
       endif else splog, 'No BCG data for cluster '+cluster
    endfor 
    

return
end
    
