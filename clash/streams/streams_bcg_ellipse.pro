function streams_sersic_func, fit.
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
      ini_e0=ini_e0, ini_pa0=ini_pa0, ellipse=ellipse1, namax=namax, $
      /silent
    if n_elements(ellipse1.xcen) ne namax then message, 'Generalize below...'

    radius = ellipse1.majora*sqrt(1.0-ellipse1.ellipticityfit)
    ellipse1 = struct_addtags(ellipse1, {radius: radius, $
      radius_kpc: radius*pixscale*arcsec2kpc})

; single Sersic fitting    
    bgt_ellipse_sersic, ellipse1, outellipse=ellipse
    
return, ellipse
end

pro convert_sb, ellipse, pixscale=pixscale
; convert the surface brightnesses in GZ's code
    ff = 0 & ll = ellipse.na-1
    pixarea = 1.0 ; 5.0*alog10(pixscale)

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

       for ib = nfilt-1, 0, -1 do begin
;      for ib = nfilt-1, nfilt-1 do begin
          imfile = outpath+cluster+'-'+short[ib]+'.fits.gz'
          splog, 'Reading '+file_basename(imfile)

; read the data and convert to intensity (maggies surface brightness)
          image = mrdfits(imfile,0,hdr,/silent)*1D-12/pixscale^2 ; [maggies/arcsec^2]
          model = mrdfits(imfile,1,/silent)*1D-12/pixscale^2     ; [maggies/arcsec^2]
          invvar = mrdfits(imfile,2,/silent)*1D24*pixscale^2
          mask = mrdfits(imfile,3,/silent)

          adxy, hdr, info[reffilt].ra, info[reffilt].dec, xcen, ycen

; ellipse-fit the data and then the model
          splog, 'Ellipse-fitting the data...'
          imellipse = do_ellipse(image,invvar=invvar*mask,$
            ini_xcen=xcen,ini_ycen=ycen,ini_e0=info[reffilt].ellipticity,$
            ini_pa0=info[reffilt].posangle,namax=namax,$
            pixscale=pixscale,arcsec2kpc=arcsec2kpc)
          convert_sb, imellipse, pixscale=pixscale

          splog, 'Ellipse-fitting the Postman model...'
          modellipse = do_ellipse(model,invvar=invvar,ini_xcen=xcen, $
            ini_ycen=ycen,ini_e0=info[reffilt].ellipticity,$
            ini_pa0=info[reffilt].posangle,namax=namax,$
            pixscale=pixscale,arcsec2kpc=arcsec2kpc)
          convert_sb, modellipse, pixscale=pixscale

          djs_plot, [0], [0], /nodata, /xlog, xsty=1, ysty=1, $
            yr=[28,15], xrange=[pixscale*arcsec2kpc,300]
          oploterror, imellipse.radius_kpc, imellipse.sb0fit, $
            1.0/sqrt(imellipse.sb0fit_ivar), psym=symcat(9), $
            symsize=0.5, color=cgcolor('orange')
          oploterror, modellipse.radius_kpc, modellipse.sb0fit, $
            1.0/sqrt(modellipse.sb0fit_ivar), psym=symcat(9), $
            symsize=0.5, color=cgcolor('dodger blue')
          djs_oplot, imellipse.radius_kpc, streams_sersic_func(imellipse), $
            color=cgcolor('firebrick')
          djs_oplot, modellipse.radius_kpc, -2.5*alog10(exp(sersic_func(modellipse.radius,$
            [modellipse.sersic_lnsb0,modellipse.sersic_k,modellipse.sersic_n]))/alog(10)), $
            color=cgcolor('cyan')
          djs_oplot, 10^!x.crange, info[ib].sblimit*[1,1], line=5
          
stop          
          
;         t0 = systime(1)
          bgt_ellipse, model/pixscale^2, imivar=invvar*pixscale^2, $
            ini_xcen=xcen, ini_ycen=ycen, ini_e0=out.ellipticity, $            
            ini_pa0=out.posangle, ellipse=ellipse1, namax=namax, $
            /silent
          if n_elements(ellipse1.xcen) ne ellipse1.namax then message, 'Generalize below...'

          radius = ellipse1.majora*sqrt(1.0-ellipse1.ellipticityfit)
          ellipse1 = struct_addtags(ellipse1, {radius: radius, $
            radius_kpc: radius*pixscale*arcsec2kpc})
          
; fit various Sersic models: (1) single Sersic; (2) Sersic + disk; (3)
; double-Sersic 
          
          streams_sersic_fitbd, muradius, ellipse.sb0fit[0:ellipse.na-1], $ ; (2)
            bdcoeff, sb_ivar=ellipse.sb0fit_ivar[0:ellipse.na-1], $
            sersicfit=sersicfit
          sersicfit = -2.5*alog10(sersicfit)

          convert_sb, ellipse, pixscale=pixscale
          
          mu_err = 1.0/sqrt(ellipse.sb0fit_ivar[0:ellipse.na-1])

          djs_plot, muradius_kpc, mu, /xlog, psym=8, xsty=1, ysty=1, $
            yr=[28,15], xrange=[pixscale*arcsec2kpc,300]
          djs_oplot, 10^!x.crange, out.sblimit[ib]*[1,1], line=5

          djs_oplot, muradius_kpc, sersicfit, color='red'
          cc = get_kbrd(1)
          
;         bgt_ellipse_sersicradius, ellipse, outradius=outrad
;         bgt_ellipse_radius() ; get the half-light radius
          
;         splog, 'Time to fit ', t0-systime(1)
;         bgt_ellipse_show, model, ellipse
          
       endfor                   ; close filter loop
; write out
       outfile = outpath+cluster+'-apphot.fits'
       im_mwrfits, out, outfile, clobber=clobber
    endfor                      ; close cluster loop

return
end
    
