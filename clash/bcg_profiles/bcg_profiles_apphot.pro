pro bcg_profiles_apphot, getinfo=getinfo, photometry=photometry, $
  debug=debug, clobber=clobber
; jm13jul01siena - perform circular aperture photometry on the *model*
; BCGs 

    outpath = getenv('CLASH_DATA')+'/bcg_apphot/'
    infofile = outpath+'bcg_apphot_info.fits'

    filt = bcgimf_filterlist(short=short,weff=weff,zpt=zpt,instr=instr)
    kl = k_lambda(weff,/odon)
    nfilt = n_elements(filt)

; ---------------------------------------------------------------------------
; first loop through and get some basic info on each BCG in the F160W
; band 
    if keyword_set(getinfo) then begin
       clash = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
       ncl = n_elements(clash)

       info = struct_addtags(clash,replicate({arcsec2kpc: 0.0, $
         ellipticity: -99.0, posangle: -99.0, size: -99.0, $
         size_kpc: -99.0},ncl))
       info.arcsec2kpc = dangular(clash.z,/kpc)/206265D ; [kpc/arcsec]
       reffilt = (where(short eq 'f160w'))[0]

;      for ic = 20, 20 do begin
       t0 = systime(1)
       for ic = 0, ncl-1 do begin
          cluster = strtrim(info[ic].shortname,2)
          if cluster eq 'a2261' then begin
             scale = '030mas'
             pixscale = 0.030D  ; [arcsec/pixel]
             mas30 = 1
          endif else begin
             scale = '065mas'
             pixscale = 0.065D  ; [arcsec/pixel]
             mas30 = 0
          endelse
          mosaicpath = clash_path(info[ic].dirname,/mosaics,mas30=mas30)
          modelpath = clash_path(info[ic].dirname,/bcgmodels,mas30=mas30)

          if file_test(modelpath,/dir) then begin
             splog, 'Getting info on cluster '+cluster
             delvarx, xcen_lum, ycen_lum
             modelfile = file_search(modelpath+cluster+'_mosaic_'+scale+'_'+strtrim(instr[reffilt],2)+$
               '_'+strtrim(short[reffilt],2)+'_drz_????????_BCG.fits.gz')
             imfile = file_search(mosaicpath+cluster+'_mosaic_'+scale+'_'+strtrim(instr[reffilt],2)+$
               '_'+strtrim(short[reffilt],2)+'_drz_????????.fits.gz')
             weightfile = file_search(mosaicpath+cluster+'_mosaic_'+scale+'_'+strtrim(instr[reffilt],2)+$
               '_'+strtrim(short[reffilt],2)+'_wht_????????.fits.gz')
             splog, modelfile
             if file_test(imfile) and file_test(weightfile) and file_test(modelfile) then begin
                model = mrdfits(modelfile,0,hdr)
                extast, hdr, astr
                sz = size(model,/dim)
                
; find the BCG and measure its basic properties
                splog, 'Searching for the BCG...'
;               wid = floor(150.0*0.5/pixscale/info[ic].arcsec2kpc) ; [+/-150 kpc in pixels]
;               x0 = sz[0]/2-wid/2 & x1 = sz[0]/2+wid/2
;               y0 = sz[1]/2-wid/2 & y1 = sz[1]/2+wid/2
;               model1 = model[x0:x1,y0:y1]
                find_galaxy, model, size, ellipticity, posangle, xcen, ycen, $
                  xcen_lum, ycen_lum, fraction=fraction, index=index, level=level
;               xcen += x0 & xcen_lum += x0
;               ycen += y0 & ycen_lum += y0
                xy2ad, xcen_lum, ycen_lum, astr, ra_bcg, dec_bcg
                info[ic].ra_bcg = ra_bcg ; update the coordinates here!
                info[ic].dec_bcg = dec_bcg
                info[ic].ellipticity = ellipticity ; ellipticity = 1-b/a
                info[ic].posangle = posangle
                info[ic].size = size*pixscale ; [arcsec]
                info[ic].size_kpc = info[ic].size*info[ic].arcsec2kpc
             endif else splog, '  No BCG model found!'
          endif else splog, 'No BCG data for cluster '+cluster
       endfor
       keep = where(info.posangle gt -99)
       im_mwrfits, info[keep], infofile, clobber=clobber
       splog, 'Total time to get info = ', (systime(1)-t0)/60.
    endif 

; ---------------------------------------------------------------------------
; now get the photometry for each cluster in INFO
    if keyword_set(photometry) then begin
       info = gz_mrdfits(infofile,1)
       ncl = n_elements(info)
       
; establish the radii and initialize the output structure 
       rmin = 0.3D                      ; [kpc]
       rmax = 120D                      ; [kpc]
       nring = round(10*alog10(rmax-rmin)) ; total number of rings
;      nring = round(24.2*alog10(rmax-rmin)) ; 10 isophotes per decade with 1.1 spacing
;      radius = range(rmin,rmax,nring,/log) ; [kpc]
;      radius = range(rmin,rmax,nring,/asinh) ; [kpc]
       radius = 10^range(alog10(rmin),alog10(rmax),nring)

       inrad = radius-(radius-shift(radius,1))/2.0 ; inner radius [kpc]
       inrad[0] = 0D
       outrad = radius+(shift(radius,-1)-radius)/2 ; outer radius [kpc]
       outrad[nring-1] = radius[nring-1]+(radius[nring-1]-radius[nring-2])/2.0

; loop on each cluster    
       t1 = systime(1)
;      for ic = 15, 15 do begin
       for ic = 0, ncl-1 do begin
          t0 = systime(1)
          cluster = strtrim(info[ic].shortname,2)
          splog, 'Working on cluster '+cluster
          arcsec2kpc = info[ic].arcsec2kpc ; [kpc/arcsec]

          if cluster eq 'a2261' then begin
             scale = '030mas'
             pixscale = 0.030D  ; [arcsec/pixel]
             mas30 = 1
          endif else begin
             scale = '065mas'
             pixscale = 0.065D  ; [arcsec/pixel]
             mas30 = 0
          endelse
          mosaicpath = clash_path(strtrim(info[ic].dirname,2),/mosaics,mas30=mas30)
          modelpath = clash_path(strtrim(info[ic].dirname,2),/bcgmodels,mas30=mas30)

; initialize the output structure          
          phot = struct_addtags(info[ic],{radius: radius, inrad: inrad, outrad: outrad, $
            mu: fltarr(nfilt,nring), mu_err: fltarr(nfilt,nring), $
            abmag: fltarr(nfilt,nring)-99.0, abmag_err: fltarr(nfilt,nring)-99.0, $
            maggies: fltarr(nfilt,nring), ivarmaggies: fltarr(nfilt,nring)})
          
; for GAIN see notes at
; http://www.ifa.hawaii.edu/~rgal/science/sextractor_notes.html
          delvarx, xcen_lum, ycen_lum
;         for ib = 12, nfilt-1 do begin
          for ib = 0, nfilt-1 do begin
             splog, 'Working on band '+short[ib]
             modelfile = file_search(modelpath+cluster+'_mosaic_'+scale+'_'+strtrim(instr[ib],2)+$
               '_'+strtrim(short[ib],2)+'_drz_????????_BCG.fits.gz')
             imfile = file_search(mosaicpath+cluster+'_mosaic_'+scale+'_'+strtrim(instr[ib],2)+$
               '_'+strtrim(short[ib],2)+'_drz_????????.fits.gz')
             weightfile = file_search(mosaicpath+cluster+'_mosaic_'+scale+'_'+strtrim(instr[ib],2)+$
               '_'+strtrim(short[ib],2)+'_wht_????????.fits.gz')
             if file_test(imfile) and file_test(weightfile) and file_test(modelfile) then begin
                model = mrdfits(modelfile,0)
                image = mrdfits(imfile,0,hdr)
                invvar = mrdfits(weightfile,0,ivarhdr)
                var = 1.0/(invvar+(invvar eq 0))*(invvar ne 0)
                extast, hdr, astr
                exptime = sxpar(hdr,'EXPTIME') ; [sec]
                sz = size(image,/dim)
             
; add the shot noise of the object in quadrature to the inverse
; variance map; poisson error = sqrt(electrons), converted back to
; native electrons/pixel/sec units
                shotvar = image/exptime ; [electron/pixel/sec]
                totalvar = shotvar + var
                totalivar = 1.0/(totalvar+(invvar eq 0))*(invvar ne 0)
             
; conversion from electrons/second/pixel to AB mag                
                zpt1 = zpt[ib]-kl[ib]*info[ic].ebv
                fact = 10.0^(-0.4*zpt1)

; get the coordinates of the BCG and then its surface brightness
; profile 
                ad2xy, info[ic].ra_bcg, info[ic].dec_bcg, $
                  astr, xcen_lum, ycen_lum
                xcen_lum1 = xcen_lum
                ycen_lum1 = ycen_lum
                model1 = model
                sectors_photometry, model1, info[ic].ellipticity, info[ic].posangle, $
                  xcen_lum1, ycen_lum1, ellradius, phi, ellcounts, ellsigmacounts, $
                  n_sectors=1, sector_width=90.0, badpixels=where(totalivar eq 0)
                
                mu = -2.5*alog10(ellcounts)+5.0*alog10(pixscale)+zpt1
                mu_err = 2.5*ellsigmacounts/ellcounts/alog(10)
                
                fnd = findex(ellradius,radius/pixscale/arcsec2kpc)
                phot.mu[ib,*] = interpolate(mu,fnd)
                phot.mu_err[ib,*] = interpolate(mu_err,fnd)
                
;               jj = read_bcg_profiles(cluster,these_filters=filt[ib])
;               djs_plot, phot.radius, phot.mu[ib,*], psym=8, yr=[25,15], ysty=3, /xlog
;               djs_oplot, jj.sma, jj.mu, psym=8, color='red'

; elliptical aperture photometry (with code to also perform circular
; aperture photometry)
                wid = max(radius)/pixscale/arcsec2kpc
                dispimage = alog10(model[xcen_lum-wid:xcen_lum+wid,ycen_lum-wid:ycen_lum+wid])
;               dispimage = image[xcen-wid:xcen+wid,ycen-wid:ycen+wid]
                sz2 = size(dispimage,/dim)
                if keyword_set(debug) then begin
                   if !d.window ne 0 then window, 0
;                  plotimage, im_imgscl(dispimage,/neg,/sqr), /noaxes, $
;                    /preserve_aspect, /norm
                   cgimage, dispimage, stretch=2, /keep_aspect, clip=2.0, /save
                endif
                
                rin = inrad/pixscale/arcsec2kpc ; [pixels]
                rout = outrad/pixscale/arcsec2kpc ; [pixels]
;               rout = rin*1.10
                
                dist_ellipse, ellrad, sz, xcen_lum, ycen_lum, $
                  1.0-info[ic].ellipticity, info[ic].posangle
                for ir = 0, nring-1 do begin
                   if keyword_set(debug) then begin
                      if (ir mod 4) eq 0 then begin
                         tvellipse, rout[ir], rout[ir]*(1-info[ic].ellipticity), sz2[0]/2, $
                           sz2[1]/2, info[ic].posangle+90, color=im_color('orange'), /data
                         tvellipse, rin[ir], rin[ir]*(1-info[ic].ellipticity), sz2[0]/2, $
                           sz2[1]/2, info[ic].posangle+90, color=im_color('orange'), /data
                      endif
                      if ir eq nring-1 then cc = get_kbrd(1)
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
                   if counts gt 0.0 then begin
                      phot.abmag[ib,ir] = -2.5*alog10(counts*fact) ; [AB maggies]
                      phot.abmag_err[ib,ir] = 2.5*errcounts/counts/alog(10)
                   endif
                   
                   phot.maggies[ib,ir] = counts*fact ; [AB maggies]
                   phot.ivarmaggies[ib,ir] = 1.0/(errcounts*fact)^2
;                  phot.sb[ib,ir] = -2.5*alog10(counts/totnpix)+5.0*$ ; [mag/arcsec^2]
;                    alog10(pixscale)+zpt1 
                endfor 

;               if keyword_set(debug) then begin
;                  if !d.window ne 0 then window, 2
;                  mag = maggies2mag(phot.maggies[ib,*],magerr=magerr,$
;                    ivarmaggies=phot.ivarmaggies[ib,*])
;                  ploterror, phot.radius, mag, magerr, xsty=3, ysty=3, $
;                    /xlog, psym=8, /trad
;               endif
             endif else splog, '  No photometry in band '+short[ib]
          endfor
; write out
          outfile = outpath+cluster+'_bcg_apphot.fits'
          im_mwrfits, phot, outfile, clobber=clobber
          splog, 'Total time for this cluster = ', (systime(1)-t0)/60.
       endfor 
       splog, 'Total time for all clusters = ', (systime(1)-t1)/60.
    endif 

; write out txt files for Jack
    
    
    
stop    

return
end
    


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

;               bgt_ellipse, model, imivar=invvar, ini_xcen=xcen_lum, $
;                 ini_ycen=ycen_lum, ini_e0=info[ic].ellipticity, $
;                 ini_pa0=info[ic].posangle, ellipse=ell
;               bgt_ellipse_sersic, ell, outellipse=final
                
