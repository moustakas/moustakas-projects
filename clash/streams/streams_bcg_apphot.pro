pro convert_sb, ellipse, pixscale=pixscale
; convert the surface brightnesses in GZ's code
    ff = 0 & ll = ellipse.na-1
    pixarea = 5.0*alog10(pixscale)

    ellipse.sb0_ivar[ff:ll] *= (ellipse.sb0[ff:ll]*alog(10)/2.5)^2.0
    ellipse.sb0fit_ivar[ff:ll] *= (ellipse.sb0fit[ff:ll]*alog(10)/2.5)^2.0
    
    ellipse.sb0[ff:ll] = -2.5*alog10(ellipse.sb0[ff:ll])+pixarea
    ellipse.sb0fit[ff:ll] = -2.5*alog10(ellipse.sb0fit[ff:ll])+pixarea
return
end    

function get_radius, rmin=rmin, rmax=rmax, inrad=inrad, outrad=outrad
    if n_elements(rmin) eq 0 then rmin = 0.3D ; [kpc]
    if n_elements(rmax) eq 0 then rmax = 120D ; [kpc]
    nring = round(10*alog10(rmax-rmin))    ; total number of rings
;   nring = round(24.2*alog10(rmax-rmin)) ; 10 isophotes per decade with 1.1 spacing
;   radius = range(rmin,rmax,nring,/log) ; [kpc]
;   radius = range(rmin,rmax,nring,/asinh) ; [kpc]
    radius = 10^range(alog10(rmin),alog10(rmax),nring)

    inrad = radius-(radius-shift(radius,1))/2.0 ; inner radius [kpc]
    inrad[0] = 0D
    outrad = radius+(shift(radius,-1)-radius)/2 ; outer radius [kpc]
    outrad[nring-1] = radius[nring-1]+(radius[nring-1]-radius[nring-2])/2.0
return, radius
end

pro streams_bcg_apphot, debug=debug, clobber=clobber
; jm13sep16siena - get aperture photometry from the outputs of
; STREAMS_BCG_POSTMAN

; note! images in units of [10^-12 erg/s/cm^2/Hz] (pico-maggies)

; read the sample and then match against the info structure written by
; streams_find_bcg 
    sample = rsex(streams_path(/propath)+'streams_sample.sex')
    splog, 'IGNORING A2261!!!'
    keep = where(strtrim(sample.shortname,2) ne 'a2261')
    sample = sample[keep]
;   struct_print, sample
    ncl = n_elements(sample)

    fact = 1D-12                           ; conversion from picomaggies to maggies
    pixscale = 0.065D                      ; [arcsec/pixel]
    pixarea = 5.0*alog10(pixscale)         ; 2.5*log10(pixscale^2)
    rmax = 150.0                           ; [kpc]

; wrap on each cluster    
    for ic = 0, 0 do begin
;   for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       splog, 'Working on cluster '+cluster
       outpath = streams_path(/postman_bcg)+cluster+'/'

; read the info structure to get the filters
       info = mrdfits(outpath+cluster+'-mgeskyinfo.fits.gz',1,/silent)
       short = strtrim(info.band,2)
       reffilt = where(short eq 'f160w') ; reference filter
       nfilt = n_elements(info)

       arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
       namax = long(rmax/pixscale/arcsec2kpc)           ; [pixels]
       
; initialize the output structure          
       phot = struct_addtags(struct_trimtags(sample[ic],select=[$
         'cluster','shortname','z']),struct_trimtags(info[reffilt],$
         select=['ra','dec','ellipticity','posangle',$
         'size','size_kpc']))
       phot = struct_addtags(phot,{sigma: fltarr(nfilt), sblimit: fltarr(nfilt)})
       phot.sigma = info.sigma
       phot.sblimit = info.sblimit

       n_sectors = 2
       sector_width = 90/(n_sectors-1.0)
       
       for ib = nfilt-1, 0, -1 do begin
;      for ib = nfilt-1, nfilt-1 do begin
          imfile = outpath+cluster+'-'+short[ib]+'.fits.gz'
          splog, 'Reading '+file_basename(imfile)

          image = mrdfits(imfile,0,hdr,/silent)*1D-12
          model = mrdfits(imfile,1,/silent)*1D-12
          invvar = mrdfits(imfile,2,/silent)*1D24
          mask = mrdfits(imfile,3,/silent)
          var = 1.0/(invvar+(invvar eq 0))*(invvar ne 0)
          sz = size(model,/dim)

          adxy, hdr, phot.ra, phot.dec, xcen, ycen
;         print, xcen, ycen

; do the ellipse-fitting; what units should MODEL be in?!? ask GZ!
          splog, 'Fitting the model ...'
;         t0 = systime(1)
          bgt_ellipse, model, imivar=invvar, ini_xcen=xcen, $
            ini_ycen=ycen, ini_e0=phot.ellipticity, /silent, $
            ini_pa0=phot.posangle, ellipse=ellipse1, namax=namax
          muradius = ellipse1.majora[0:ellipse1.na-1]*$
            sqrt(1.0-ellipse1.ellipticityfit[0:ellipse1.na-1])
          muradius_kpc = muradius*pixscale*arcsec2kpc
          
; fit various Sersic models: (1) single Sersic; (2) Sersic + disk; (3)
; double-Sersic 
          bgt_ellipse_sersic, ellipse1, outellipse=ellipse ; (1)
          
          streams_sersic_fitbd, muradius, ellipse.sb0fit[0:ellipse.na-1], $ ; (2)
            bdcoeff, sb_ivar=ellipse.sb0fit_ivar[0:ellipse.na-1], $
            sersicfit=sersicfit
          sersicfit = -2.5*alog10(sersicfit)+pixarea

          convert_sb, ellipse, pixscale=pixscale
          
          mu = ellipse.sb0fit[0:ellipse.na-1]
          mu_err = 1.0/sqrt(ellipse.sb0fit_ivar[0:ellipse.na-1])

          djs_plot, muradius_kpc, mu, /xlog, psym=8, xsty=1, ysty=1, $
            yr=[28,15], xrange=[pixscale*arcsec2kpc,300]
          djs_oplot, 10^!x.crange, phot.sblimit[ib]*[1,1], line=5

          ellparams = [ellipse.sersic_lnsb0,ellipse.sersic_k,ellipse.sersic_n]
          djs_oplot, muradius_kpc, -2.5*alog10(exp(sersic_func(muradius,ellparams)))+pixarea, $
            color='blue'
          djs_oplot, muradius_kpc, sersicfit, color='red'
          cc = get_kbrd(1)
          
;         bgt_ellipse_sersicradius, ellipse, outradius=outrad
;         bgt_ellipse_radius() ; get the half-light radius
          
;         splog, 'Time to fit ', t0-systime(1)
;         bgt_ellipse_show, model, ellipse
          
;; --------------------------------------------------
;; all this code is good but I don't think I'm going to use it!
;          
;; photometer the model          
;          model1 = model
;          xcen1 = xcen & ycen1 = ycen
;          sectors_photometry, model1, phot.ellipticity, phot.posangle, $
;            xcen1, ycen1, muradius, phi, counts, errcounts, $
;            n_sectors=n_sectors, sector_width=sector_width, $
;            badpixels=where(invvar eq 0)
;          srt = sort(muradius)
;          muradius = muradius[srt] & phi = phi[srt]
;          counts = counts[srt] & errcounts = errcounts[srt]
;          muradius_kpc = muradius*pixscale*arcsec2kpc ; [kpc]
;
;          mu = -2.5*alog10(counts)+5.0*alog10(pixscale)
;          mu_err = 2.5*errcounts/counts/alog(10)
;
;; photometer the data (both masked & unmasked); my conclusion is that
;; something funky is going on with the masked data at large radii 
;          image1 = image
;          xcen1 = xcen & ycen1 = ycen
;          sectors_photometry, image1, phot.ellipticity, phot.posangle, $
;            xcen1, ycen1, muobsradius, obsphi, obscounts, obserrcounts, $
;            n_sectors=n_sectors, sector_width=sector_width, $
;            badpixels=where(invvar eq 0); or mask eq 0)
;          srt = sort(muobsradius)
;          muobsradius = muobsradius[srt] & obsphi = phi[srt]
;          obscounts = obscounts[srt] & obserrcounts = obserrcounts[srt]
;
;          muobsradius_kpc = muobsradius*pixscale*arcsec2kpc ; [kpc]
;          mu_obs = -2.5*alog10(obscounts)+5.0*alog10(pixscale)
;          mu_obs_err = 2.5*obserrcounts/obscounts/alog(10)
;
;          image1 = image
;          xcen1 = xcen & ycen1 = ycen
;          sectors_photometry, image1, phot.ellipticity, phot.posangle, $
;            xcen1, ycen1, muobs2radius, obs2phi, obs2counts, obs2errcounts, $
;            n_sectors=n_sectors, sector_width=sector_width, $
;            badpixels=where(invvar eq 0 or mask eq 1)
;          srt = sort(muobs2radius)
;          muobs2radius = muobs2radius[srt] & obs2phi = phi[srt]
;          obs2counts = obs2counts[srt] & obs2errcounts = obs2errcounts[srt]
;
;          muobs2radius_kpc = muobs2radius*pixscale*arcsec2kpc ; [kpc]
;          mu_obs2 = -2.5*alog10(obs2counts)+5.0*alog10(pixscale)
;          mu_obs2_err = 2.5*obs2errcounts/obs2counts/alog(10)
;
;;         if n_elements(muradius) ne n_elements(muobsradius) then message, 'Problem here!'
;          if keyword_set(debug) then begin
;             djs_plot, muobsradius_kpc, mu_obs, psym=8, /xlog, xrange=[0.06,500], $
;               yrange=[32,15], xsty=3, ysty=1, xtitle='Semi-Major Axis (kpc)', $
;               ytitle='Surface Brightness ('+short[ib]+')'
;             djs_oplot, muobs2radius_kpc, mu_obs2, psym=6, color=cgcolor('orange')
;             djs_oplot, muradius_kpc, mu, color=cgcolor('dodger blue'), psym=8
;             djs_oplot, 10^!x.crange, phot.sblimit[ib]*[1,1], line=5
;             cc = get_kbrd(1)
;          endif
;          
;; build the radial annuli and the output structure
;          if ib eq reffilt then begin ; not general!
;             radius = get_radius(rmin=min(muradius_kpc),$
;               rmax=max(muradius_kpc),inrad=inrad,outrad=outrad)
;             nrad = n_elements(radius)
;             phot = struct_addtags(phot,{$
;               radius:     radius, $
;               inrad:      inrad, $
;               outrad:     outrad, $
;               mu:         fltarr(nfilt,nrad), $
;               mu_err:     fltarr(nfilt,nrad), $
;               mu_obs:     fltarr(nfilt,nrad), $
;               mu_obs_err: fltarr(nfilt,nrad), $
;               abmag:      fltarr(nfilt,nrad)-99.0, $
;               abmag_err:  fltarr(nfilt,nrad)-99.0, $
;               maggies:    fltarr(nfilt,nrad), $
;               ivarmaggies: fltarr(nfilt,nrad)})
;          endif else nrad = n_elements(phot.radius)
;
;          fnd = findex(muradius,phot.radius/pixscale/arcsec2kpc)
;          phot.mu[ib,*] = interpolate(mu,fnd)
;          phot.mu_err[ib,*] = sqrt(interpolate(mu_err^2,fnd))
;; --------------------------------------------------

          continue
          
; now do elliptical and circular aperture photometry (with code to
; also perform circular aperture photometry)
          rin = phot.inrad/pixscale/arcsec2kpc   ; [pixels]
          rout = phot.outrad/pixscale/arcsec2kpc ; [pixels]

          dist_ellipse, ellrad, sz, xcen, ycen, $
            1.0-phot.ellipticity, phot.posangle

          if keyword_set(debug) then begin
             window, 0
             cgimage, image, stretch=3, /keep_aspect, /save, $
               clip=3, /negative, minvalue=0.0, maxvalue=3.0
          endif
          
          for ir = 0, nrad-1 do begin
             if keyword_set(debug) then begin
                if (ir mod 4) eq 0 then begin
                   tvellipse, rout[ir], rout[ir]*(1-phot.ellipticity), sz[0]/2, $
                     sz[1]/2, phot.posangle+90, color=im_color('orange'), /data
                   tvellipse, rin[ir], rin[ir]*(1-phot.ellipticity), sz[0]/2, $
                     sz[1]/2, phot.posangle+90, color=im_color('orange'), /data
                endif
                if ir eq nrad-1 then cc = get_kbrd(1)
             endif
             
             pix = where(ellrad ge rin[ir] and ellrad lt rout[ir],npix)
             totnpix = (npix-total(invvar[pix] eq 0))
;            area = totnpix*pixscale^2 ; [arcsec^2]
             area = (rout[ir]-rin[ir])*!pi*(1-phot.ellipticity)
             
             counts = total(model[pix]*(invvar[pix] gt 0))
             errcounts = sqrt(total(var[pix]*(invvar[pix] gt 0)))
;            print, rin[ir], rout[ir], counts, area
                
             if counts gt 0.0 then begin
                phot.abmag[ib,ir] = -2.5*alog10(counts*fact) ; [AB maggies]
                phot.abmag_err[ib,ir] = 2.5*errcounts/counts/alog(10)
             endif
             
             phot.maggies[ib,ir] = counts*fact ; [AB maggies]
             phot.ivarmaggies[ib,ir] = 1.0/(errcounts*fact)^2
          endfor 
       endfor                   ; close filter loop
; write out
       outfile = outpath+cluster+'-apphot.fits'
       im_mwrfits, phot, outfile, clobber=clobber
    endfor                      ; close cluster loop

return
end
    
