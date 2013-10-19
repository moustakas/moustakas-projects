function get_radius, rmin=rmin, rmax=rmax, inrad=inrad, outrad=outrad
    if n_elements(rmin) eq 0 then rmin = 5D   ; [pixels]
    if n_elements(rmax) eq 0 then rmax = 500D ; [pixels]
;   nring = round(10*alog10(rmax-rmin))    ; total number of rings
    nring = round(10.0*alog10(rmax-rmin)) ; 10 isophotes per decade with 1.1 spacing
;   radius = range(rmin,rmax,nring,/log) ; [kpc]
    radius = range(rmin,rmax,nring,/asinh) ; [kpc]
;   radius = 10.0^range(alog10(rmin),alog10(rmax),nring)

    inrad = radius-(radius-shift(radius,1))/2.0 ; inner radius [kpc]
    inrad[0] = 0D
    outrad = radius+(shift(radius,-1)-radius)/2 ; outer radius [kpc]
    outrad[nring-1] = radius[nring-1]+(radius[nring-1]-radius[nring-2])/2.0
return, radius
end

function get_median_ellipse, ellipse, sblimit=sblimit
; concatenate the two structures plus add the median ellipse
; parameters; the position angle is measured counter-clockwise with
; respect to the image Y-axis
    good = where(ellipse.sb0fit lt sblimit,ngood)
    out = {ellipticity_median: 0.0, posangle_median: 0.0}
    out.ellipticity_median = djs_median(ellipse.ellipticityfit)
    out.posangle_median = djs_median(ellipse.pafit)*!radeg+90.0 ; [deg]
return, out
end

function streams_sersic_func, fit
; return the magnitude version of the best-fitting Sersic model
    params = [fit.sersic_lnsb0,fit.sersic_k,fit.sersic_n]
    model = sersic_func(fit.radius,params)
    model = -2.5*alog10(exp(model)/alog(10))
return, model    
end

function do_ellipse, image, invvar=invvar, ini_xcen=ini_xcen, $
  ini_ycen=ini_ycen, ini_e0=ini_e0, ini_pa0=ini_pa0, namax=namax, $
  pixscale=pixscale, arcsec2kpc=arcsec2kpc  
; do the ellipse-fitting    

    bgt_ellipse, image, imivar=invvar, ini_xcen=ini_xcen, ini_ycen=ini_ycen, $
      ini_e0=ini_e0, ini_pa0=ini_pa0, ellipse=ellipse1, namax=namax, /silent
    if n_elements(ellipse1.xcen) ne namax then message, 'Generalize below...'

    radius = ellipse1.majora*sqrt(1.0-ellipse1.ellipticityfit)
    ellipse1 = struct_addtags(ellipse1, {radius: radius, $
      radius_kpc: radius*pixscale*arcsec2kpc})

; single Sersic fitting    
;   bgt_ellipse_sersic, ellipse1, outellipse=ellipse
    
return, ellipse1
end

pro convert_sb, ellipse, pixscale=pixscale
; convert the surface brightnesses in GZ's code from intensity
; to mag/arcsec^2 
    ff = 0 & ll = ellipse.na-1
    pixarea = 5.0*alog10(pixscale)

    ellipse.sb0_ivar[ff:ll] *= (ellipse.sb0[ff:ll]*alog(10)/2.5)^2.0
    ellipse.sb0fit_ivar[ff:ll] *= (ellipse.sb0fit[ff:ll]*alog(10)/2.5)^2.0
    
    ellipse.sb0[ff:ll] = -2.5*alog10(ellipse.sb0[ff:ll])+pixarea
    ellipse.sb0fit[ff:ll] = -2.5*alog10(ellipse.sb0fit[ff:ll])+pixarea
return
end    

pro streams_bcg_ellipse, debug=debug, clobber=clobber
; jm13sep16siena - perform ellipse fitting on the BCG cutouts written
; out by STREAMS_BCG_POSTMAN

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
;   for ic = 0, 0 do begin
    for ic = 0, ncl-1 do begin
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

       delvarx, imellipse, modellipse, phot
       for ib = nfilt-1, 0, -1 do begin
;      for ib = nfilt-1, nfilt-1 do begin
          imfile = outpath+cluster+'-'+short[ib]+'.fits.gz'
          splog, 'Reading '+file_basename(imfile)

; read the data and convert to intensity (maggies surface brightness)
          image = mrdfits(imfile,0,hdr,/silent)*1D-12 ; [maggies]
          model = mrdfits(imfile,1,/silent)*1D-12     ; [maggies]
          invvar = mrdfits(imfile,2,/silent)*1D24     ; [1/maggies^2]
;         image = mrdfits(imfile,0,hdr,/silent)*1D-12/pixscale^2 ; [maggies/arcsec^2]
;         model = mrdfits(imfile,1,/silent)*1D-12/pixscale^2     ; [maggies/arcsec^2]
;         invvar = mrdfits(imfile,2,/silent)*1D24*pixscale^2
          mask = mrdfits(imfile,3,/silent)
          var = 1.0/(invvar+(invvar eq 0))*(invvar ne 0)
          sz = size(image,/dim)

          adxy, hdr, info[reffilt].ra, info[reffilt].dec, xcen, ycen

; ellipse-fit the model and then use the identical ellipse parameters
; to generate the surface-brightness profile from the (unmasked)
; image; note that the ellipse PA (in radians) is measured positive
; from the X-axis, not the *Y-axis* as we assume below
          modellipse1 = do_ellipse(model,invvar=invvar,ini_xcen=xcen, $
            ini_ycen=ycen,ini_e0=info[reffilt].mge_ellipticity,$
            ini_pa0=info[reffilt].mge_posangle,namax=namax,$
            pixscale=pixscale,arcsec2kpc=arcsec2kpc)

          imellipse1 = modellipse1
          for ia = 0, modellipse1.na-1 do begin
             sbcoord = bgt_ellipse_sb(image,imivar=invvar,xcen=modellipse1.xcenfit[ia],$
               ycen=modellipse1.xcenfit[ia],pa=modellipse1.pafit[ia],$
               ee=modellipse1.ellipticityfit[ia],majora=modellipse1.majora[ia],$
               ea=ea,sbcoord_ivar=sbcoord_ivar)
             bgt_ellipse_fit, ea, sbcoord, ecoeff, ini_value=ini_value, $
               sb_ivar=sbcoord_ivar, status=status, fixed=[0,1,1,1,1]
             imellipse1.sb0[ia] = ecoeff.sb0 ; only free parameter is the SB
             imellipse1.sb0_ivar[ia] = median(sbcoord_ivar)
          endfor
          bgt_ellipse_allfit, image, imellipse1, imivar=invvar

; convert from intensity to mag/arcsec^2          
          convert_sb, modellipse1, pixscale=pixscale
          convert_sb, imellipse1, pixscale=pixscale

; get the median ellipticity and position angle          
          median_ellipse = get_median_ellipse(modellipse1,sblimit=info[ib].sblimit)
          modellipse1 = struct_addtags(struct_addtags(info[ib],modellipse1),median_ellipse)
          imellipse1 = struct_addtags(struct_addtags(info[ib],imellipse1),median_ellipse)

          if n_elements(modellipse) eq 0 then modellipse = modellipse1 else $
            modellipse = [modellipse,modellipse1]
          if n_elements(imellipse) eq 0 then imellipse = imellipse1 else $
            imellipse = [imellipse,imellipse1]
          
;;         splog, 'Ellipse-fitting the data...'
;          imellipse1 = do_ellipse(image,invvar=invvar,$;*mask,$
;            ini_xcen=xcen,ini_ycen=ycen,ini_e0=info[reffilt].mge_ellipticity,$
;            ini_pa0=info[reffilt].mge_posangle,namax=namax,$
;            pixscale=pixscale,arcsec2kpc=arcsec2kpc)
;          convert_sb, imellipse1, pixscale=pixscale
;          imellipse1 = struct_addtags(info[ib],imellipse1)
;          if n_elements(imellipse) eq 0 then imellipse = imellipse1 else $
;            imellipse = [imellipse,imellipse1]

; now do elliptical-aperture photometry on the model images; build the
; radial annuli and the output structure
          if ib eq reffilt then begin ; not general!
             radius = get_radius(rmin=min(imellipse1.radius_kpc/pixscale/arcsec2kpc)>5,$
               rmax=max(imellipse1.radius_kpc/pixscale/arcsec2kpc),inrad=rin,$
               outrad=rout)
          endif
          nrad = n_elements(radius)

          phot1 = {$
            radius_kpc:    radius*pixscale*arcsec2kpc, $
            inradius_kpc:  rin*pixscale*arcsec2kpc, $
            outradius_kpc: rout*pixscale*arcsec2kpc, $
            npix:          fltarr(nrad), $
            area:          fltarr(nrad), $
            abmag:         fltarr(nrad)-99.0, $
            abmag_err:     fltarr(nrad)-99.0, $
            maggies:       fltarr(nrad), $
            ivarmaggies:   fltarr(nrad)}
          
          dist_ellipse, ellrad, sz, xcen, ycen, $
            1.0/(1.0-modellipse1.ellipticity_median), $
            modellipse1.posangle_median
;         atv, ellrad, /bl

          for ir = 0, nrad-1 do begin
             pix = where(ellrad ge rin[ir] and ellrad lt rout[ir],npix)
             totnpix = (npix-total(invvar[pix] eq 0))
;            area = totnpix*pixscale^2 ; [arcsec^2]
             area = !pi*(rout[ir]^2*(1-modellipse1.ellipticity_median)-$
               rin[ir]^2*(1-modellipse1.ellipticity_median)) ; [pixel^2]
             phot1.npix[ir] = totnpix
             phot1.area[ir] = area*pixscale^2 ; [arcsec^2]
             
             maggies = total(model[pix]*(invvar[pix] gt 0))
             errmaggies = sqrt(total(var[pix]*(invvar[pix] gt 0)))
             print, rin[ir], rout[ir], maggies, area
                
             if maggies gt 0.0 then begin
                phot1.abmag[ir] = -2.5*alog10(maggies) ; [AB maggies]
                phot1.abmag_err[ir] = 2.5*errmaggies/maggies/alog(10)
             endif

             if errmaggies le 0 then message, 'Problem here!'
             phot1.maggies[ir] = maggies ; [AB maggies]
             phot1.ivarmaggies[ir] = 1.0/errmaggies^2
          endfor
          if n_elements(phot) eq 0 then phot = phot1 else phot = [phot,phot1]
;         djs_plot, phot1.radius_kpc, phot1.abmag, ysty=3, psym=8, /xlog, xsty=3
          
;          djs_plot, [0], [0], /nodata, /xlog, xsty=1, ysty=1, $
;            yr=[28,15], xrange=[pixscale*arcsec2kpc,300]
;          oploterror, imellipse.radius_kpc, imellipse.sb0fit, $
;            1.0/sqrt(imellipse.sb0fit_ivar), psym=symcat(9), $
;            symsize=0.5, color=cgcolor('orange')
;          oploterror, modellipse.radius_kpc, modellipse.sb0fit, $
;            1.0/sqrt(modellipse.sb0fit_ivar), psym=symcat(9), $
;            symsize=0.5, color=cgcolor('dodger blue')
;          djs_oplot, imellipse.radius_kpc, streams_sersic_func(imellipse), $
;            color=cgcolor('firebrick')
;          djs_oplot, modellipse.radius_kpc, -2.5*alog10(exp(sersic_func(modellipse.radius,$
;            [modellipse.sersic_lnsb0,modellipse.sersic_k,modellipse.sersic_n]))/alog(10)), $
;            color=cgcolor('cyan')
;          djs_oplot, 10^!x.crange, info[ib].sblimit*[1,1], line=5
;          
;stop          
;          
;;         t0 = systime(1)
;          bgt_ellipse, model/pixscale^2, imivar=invvar*pixscale^2, $
;            ini_xcen=xcen, ini_ycen=ycen, ini_e0=out.ellipticity, $            
;            ini_pa0=out.posangle, ellipse=ellipse1, namax=namax, $
;            /silent
;          if n_elements(ellipse1.xcen) ne ellipse1.namax then message, 'Generalize below...'
;
;          radius = ellipse1.majora*sqrt(1.0-ellipse1.ellipticityfit)
;          ellipse1 = struct_addtags(ellipse1, {radius: radius, $
;            radius_kpc: radius*pixscale*arcsec2kpc})
;          
;; fit various Sersic models: (1) single Sersic; (2) Sersic + disk; (3)
;; double-Sersic 
;          
;          streams_sersic_fitbd, muradius, ellipse.sb0fit[0:ellipse.na-1], $ ; (2)
;            bdcoeff, sb_ivar=ellipse.sb0fit_ivar[0:ellipse.na-1], $
;            sersicfit=sersicfit
;          sersicfit = -2.5*alog10(sersicfit)
;
;          convert_sb, ellipse, pixscale=pixscale
;          
;          mu_err = 1.0/sqrt(ellipse.sb0fit_ivar[0:ellipse.na-1])
;
;          djs_plot, muradius_kpc, mu, /xlog, psym=8, xsty=1, ysty=1, $
;            yr=[28,15], xrange=[pixscale*arcsec2kpc,300]
;          djs_oplot, 10^!x.crange, out.sblimit[ib]*[1,1], line=5
;
;          djs_oplot, muradius_kpc, sersicfit, color='red'
;          cc = get_kbrd(1)
;          
;;         bgt_ellipse_sersicradius, ellipse, outradius=outrad
;;         bgt_ellipse_radius() ; get the half-light radius
;          
;;         splog, 'Time to fit ', t0-systime(1)
;;         bgt_ellipse_show, model, ellipse
;          
       endfor                   ; close filter loop
; write everything out
       im_mwrfits, imellipse, outpath+cluster+'-ellipse-image.fits', clobber=clobber
       im_mwrfits, modellipse, outpath+cluster+'-ellipse-model.fits', clobber=clobber
       im_mwrfits, phot, outpath+cluster+'-ellipse-ellphot.fits', clobber=clobber
stop
    endfor                      ; close cluster loop

stop    
    
return
end
    
