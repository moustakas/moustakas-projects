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

function do_ellipse, image, invvar=invvar, ini_xcen=ini_xcen, $
  ini_ycen=ini_ycen, ini_e0=ini_e0, ini_pa0=ini_pa0, namax=namax, $
  pixscale=pixscale, arcsec2kpc=arcsec2kpc, fixcenter=fixcenter
; do the ellipse-fitting    

    bgt_ellipse, image, imivar=invvar, ini_xcen=ini_xcen, ini_ycen=ini_ycen, $
      ini_e0=ini_e0, ini_pa0=ini_pa0, ellipse=ellipse1, namax=namax, $
      fixcenter=fixcenter, /silent
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

pro bcgmstar_ellipse, debug=debug, clobber=clobber
; jm13sep16siena - perform ellipse fitting on the BCG cutouts written
; out by BCGMSTAR_GET_BCGS

; note! images in units of [10^-12 erg/s/cm^2/Hz] (pico-maggies)

    clobber = 1
    
; read the sample
    sample = read_bcgmstar_sample()
    ncl = n_elements(sample)

    fact = 1D-12                           ; conversion from picomaggies to maggies
    pixscale = 0.065D                      ; [arcsec/pixel]
    pixarea = 5.0*alog10(pixscale)         ; 2.5*log10(pixscale^2)
;   rmax = 10.0                            ; [kpc]
;   rmax = 25.0                            ; [kpc]
    rmax = 200.0                           ; [kpc]

    ellpath = bcgmstar_path(/ellipse)
    skyinfopath = bcgmstar_path()+'skyinfo/'

; wrap on each cluster    
    for ic = 4, 4 do begin
;   for ic = 5, ncl-1 do begin
;   for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       splog, 'Working on cluster '+cluster
       datapath = bcgmstar_path(/bcg)+cluster+'/'

; read the info structure to get the filters; exclude certain
; bandpasses here
       info = mrdfits(skyinfopath+cluster+'-mgeskyinfo.fits.gz',1,/silent)
       case cluster of
          'macs1311': info = info[where(strtrim(info.band,2) ne 'f435w')]
          'macs0744': info = info[where(strtrim(info.band,2) ne 'f475w')]
          'macs1149': info = info[where(strtrim(info.band,2) ne 'f475w')]
          'macs2129': info = info[where(strtrim(info.band,2) ne 'f435w' and strtrim(info.band,2) ne 'f475w')]
          else:
       endcase

       short = strtrim(info.band,2)
       reffilt = where(short eq 'f160w') ; reference filter
       nfilt = n_elements(info)

       arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
       namax = long(rmax/pixscale/arcsec2kpc)           ; [pixels]

; ellipse-fit the data and the model with all parameters free **in the
; reference band (NFILT-1)**; then perform a *constrained* fit to all
; the other bands; note that the ellipse PA (in radians) in GZ's code
; is measured positive from the X-axis, not the *Y-axis* as we assume
; below
       imfile = datapath+cluster+'-'+short[reffilt]+'.fits.gz'
       splog, 'Reading '+file_basename(imfile)

; read the data and convert to intensity (maggies surface brightness)
;      image = mrdfits(imfile,0,hdr,/silent)*1D-12 ; [maggies]
;      model = mrdfits(imfile,1,/silent)*1D-12     ; [maggies]
;      invvar = mrdfits(imfile,2,/silent)*1D24     ; [1/maggies^2]
       image = mrdfits(imfile,0,hdr,/silent)*1D-12/pixscale^2 ; [maggies/arcsec^2]
       invvar = mrdfits(imfile,1,/silent)*(1D12*pixscale^2)^2
       mask = mrdfits(imfile,2,/silent)
       model = mrdfits(imfile,3,/silent)*1D-12/pixscale^2 ; [maggies/arcsec^2]
       var = 1.0/(invvar+(invvar eq 0))*(invvar ne 0)
       sz = size(image,/dim)
       
       adxy, hdr, info[reffilt].ra, info[reffilt].dec, xcen, ycen

; fit in the reference band          
       refmodellipse = do_ellipse(model,invvar=invvar,ini_xcen=xcen, $
         ini_ycen=ycen,ini_e0=info[reffilt].mge_ellipticity,$
         ini_pa0=info[reffilt].mge_posangle,namax=namax,$
         pixscale=pixscale,arcsec2kpc=arcsec2kpc)
       mod_median_ellipse = get_median_ellipse(refmodellipse,sblimit=info[reffilt].sblimit)
       refmodellipse = struct_addtags(refmodellipse,mod_median_ellipse)

; also fit the actual data in the reference band
       if cluster eq 'a611' then begin
          refimellipse = do_ellipse(image,invvar=invvar*mask,ini_xcen=xcen, $
            ini_ycen=ycen,ini_e0=info[reffilt].mge_ellipticity,$
            ini_pa0=info[reffilt].mge_posangle,namax=namax,$
            pixscale=pixscale,arcsec2kpc=arcsec2kpc)
          im_median_ellipse = get_median_ellipse(refimellipse,sblimit=info[reffilt].sblimit)
          refimellipse = struct_addtags(refimellipse,im_median_ellipse)
          im_mwrfits, refimellipse, ellpath+cluster+'-ellipse-refimage.fits', clobber=clobber
       endif

; now loop through and do a constrained fit of the *model* images in
; every other band 
       delvarx, imellipse, modellipse, phot
       for ib = nfilt-1, 0, -1 do begin
;      for ib = nfilt-2, nfilt-1 do begin
;      for ib = 0, nfilt-1 do begin
          imfile = datapath+cluster+'-'+short[ib]+'.fits.gz'
          splog, 'Reading '+file_basename(imfile)

; read the data and convert to intensity (maggies surface brightness)
;         image = mrdfits(imfile,0,hdr,/silent)*1D-12 ; [maggies]
;         model = mrdfits(imfile,1,/silent)*1D-12     ; [maggies]
;         invvar = mrdfits(imfile,2,/silent)*1D24     ; [1/maggies^2]
          image = mrdfits(imfile,0,hdr,/silent)*1D-12/pixscale^2 ; [maggies/arcsec^2]
          invvar = mrdfits(imfile,1,/silent)*(1D12*pixscale^2)^2
          mask = mrdfits(imfile,2,/silent)
          model = mrdfits(imfile,3,/silent)*1D-12/pixscale^2     ; [maggies/arcsec^2]
          var = 1.0/(invvar+(invvar eq 0))*(invvar ne 0)
          sz = size(image,/dim)

          adxy, hdr, info[reffilt].ra, info[reffilt].dec, xcen, ycen

; ellipse-fit the data and the model with all parameters free **in the
; reference band (NFILT-1)**; then perform a *constrained* fit to all
; the other bands; note that the ellipse PA (in radians) in GZ's code
; is measured positive from the X-axis, not the *Y-axis* as we assume
; below

;;         imellipse1 = do_ellipse(image,invvar=invvar*mask,ini_xcen=xcen,$
;          ii = do_ellipse(image,invvar=invvar*mask,ini_xcen=xcen,$
;            ini_ycen=ycen,ini_e0=info[reffilt].mge_ellipticity,$
;            ini_pa0=info[reffilt].mge_posangle,namax=namax,$
;            pixscale=pixscale,arcsec2kpc=arcsec2kpc);,/fixcenter)
;
;;         modellipse1 = do_ellipse(model,invvar=invvar,ini_xcen=xcen, $
;          mm = do_ellipse(model,invvar=invvar,ini_xcen=xcen, $
;            ini_ycen=ycen,ini_e0=info[reffilt].mge_ellipticity,$
;            ini_pa0=info[reffilt].mge_posangle,namax=namax,$
;            pixscale=pixscale,arcsec2kpc=arcsec2kpc)
;
;          djs_plot, mm.majora, mm.pafit
;          djs_oplot, ii.majora, ii.pafit, color='red'
;          djs_plot, mm.majora, mm.ellipticityfit
;          djs_oplot, ii.majora, ii.ellipticityfit, color='red'
;          djs_plot, mm.majora, mm.sb0fit, xsty=3, ysty=3, /xlog, /ylog
;          djs_oplot, ii.majora, ii.sb0fit, color='red'

          modellipse1 = refmodellipse
          for ia = 0, refmodellipse.na-1 do begin
             print, format='("BGT_ELLIPSE: Major-axis pixel ",I0,"/",I0, A10,$)', $
               ia+1, refmodellipse.na, string(13b)
             sbcoord = bgt_ellipse_sb(model,imivar=invvar,xcen=refmodellipse.xcenfit[ia],$
               ycen=refmodellipse.xcenfit[ia],pa=refmodellipse.pafit[ia],$
               ee=refmodellipse.ellipticityfit[ia],majora=refmodellipse.majora[ia],$
               ea=ea,sbcoord_ivar=sbcoord_ivar)
             bgt_ellipse_fit, ea, sbcoord, ecoeff, ini_value=ini_value, $
               sb_ivar=sbcoord_ivar, status=status, fixed=[0,1,1,1,1]
             modellipse1.sb0[ia] = ecoeff.sb0 ; only free parameter is the SB
             modellipse1.sb0_ivar[ia] = median(sbcoord_ivar)
          endfor
          bgt_ellipse_allfit, model, modellipse1, imivar=invvar
          modellipse1 = struct_addtags(info[ib],modellipse1)
          if n_elements(modellipse) eq 0 then modellipse = modellipse1 else $
            modellipse = [modellipse,modellipse1]

;; debugging plot
;          djs_plot, refmodellipse.majora, refmodellipse.sb0fit, $
;            xsty=3, ysty=3, /xlog, /ylog, yrange=[min(refmodellipse.sb0fit)<$
;            min(modellipse.sb0fit),max(refmodellipse.sb0fit)>$
;            max(modellipse.sb0fit)]
;          djs_oplot, modellipse1.majora, modellipse1.sb0fit, color='green'
;          im_legend, short[ib], /right, /top, box=0
;;         cc = get_kbrd(1)

       endfor                   ; close filter loop
; write everything out
;      im_mwrfits, imellipse, ellpath+cluster+'-ellipse-image.fits', clobber=clobber
       im_mwrfits, modellipse, ellpath+cluster+'-ellipse-model.fits', clobber=clobber
;      im_mwrfits, phot, ellpath+cluster+'-ellipse-ellphot.fits', clobber=clobber
    endfor                      ; close cluster loop

return
end
    
;;; optionally use the identical ellipse parameters to generate the
;;; surface-brightness profile from the data themselves
;;          constrained_fit = 0
;;          if constrained_fit then begin
;;             imellipse1 = modellipse1
;;             for ia = 0, modellipse1.na-1 do begin
;;                sbcoord = bgt_ellipse_sb(image,imivar=invvar,xcen=modellipse1.xcenfit[ia],$
;;                  ycen=modellipse1.xcenfit[ia],pa=modellipse1.pafit[ia],$
;;                  ee=modellipse1.ellipticityfit[ia],majora=modellipse1.majora[ia],$
;;                  ea=ea,sbcoord_ivar=sbcoord_ivar)
;;                bgt_ellipse_fit, ea, sbcoord, ecoeff, ini_value=ini_value, $
;;                  sb_ivar=sbcoord_ivar, status=status, fixed=[0,1,1,1,1]
;;                imellipse1.sb0[ia] = ecoeff.sb0 ; only free parameter is the SB
;;                imellipse1.sb0_ivar[ia] = median(sbcoord_ivar)
;;             endfor
;;             bgt_ellipse_allfit, image, imellipse1, imivar=invvar
;;          endif else begin
;;             imellipse1 = do_ellipse(image,invvar=invvar*mask,ini_xcen=xcen,$
;;               ini_ycen=ycen,ini_e0=info[reffilt].mge_ellipticity,$
;;               ini_pa0=info[reffilt].mge_posangle,namax=namax,$
;;               pixscale=pixscale,arcsec2kpc=arcsec2kpc)
;;          endelse
;
;; get the median ellipticity and position angle from the model 
;          median_ellipse = get_median_ellipse(modellipse1,sblimit=info[ib].sblimit)
;          modellipse1 = struct_addtags(struct_addtags(info[ib],modellipse1),median_ellipse)
;          imellipse1 = struct_addtags(struct_addtags(info[ib],imellipse1),median_ellipse)
;
;          if n_elements(modellipse) eq 0 then modellipse = modellipse1 else $
;            modellipse = [modellipse,modellipse1]
;          if n_elements(imellipse) eq 0 then imellipse = imellipse1 else $
;            imellipse = [imellipse,imellipse1]
;          
;; now do elliptical-aperture photometry on the model images; build the
;; radial annuli and the output structure
;          if ib eq reffilt then begin ; not general!
;             radius = get_radius(rmin=min(imellipse1.radius_kpc/pixscale/arcsec2kpc)>5,$
;               rmax=max(imellipse1.radius_kpc/pixscale/arcsec2kpc),inrad=rin,$
;               outrad=rout)
;          endif
;          nrad = n_elements(radius)
;
;          phot1 = {$
;            radius_kpc:    radius*pixscale*arcsec2kpc, $
;            inradius_kpc:  rin*pixscale*arcsec2kpc, $
;            outradius_kpc: rout*pixscale*arcsec2kpc, $
;            npix:          fltarr(nrad), $
;            area:          fltarr(nrad), $
;;           abmag:         fltarr(nrad)-99.0, $
;;           abmag_err:     fltarr(nrad)-99.0, $
;            maggies:       fltarr(nrad), $
;            ivarmaggies:   fltarr(nrad)}
;          
;          dist_ellipse, ellrad, sz, xcen, ycen, $
;            1.0/(1.0-modellipse1.ellipticity_median), $
;            modellipse1.posangle_median
;;         atv, ellrad, /bl
;
;          for ir = 0, nrad-1 do begin
;;            atv, model*pixscale^2*(ellrad ge rin[ir] and ellrad lt rout[ir]), /bl
;             pix = where(ellrad ge rin[ir] and ellrad lt rout[ir],npix)
;             totnpix = (npix-total(invvar[pix] eq 0))
;;            area = totnpix*pixscale^2 ; [arcsec^2]
;             area = !pi*(rout[ir]^2*(1-modellipse1.ellipticity_median)-$
;               rin[ir]^2*(1-modellipse1.ellipticity_median)) ; [pixel^2]
;             phot1.npix[ir] = totnpix
;             phot1.area[ir] = area*pixscale^2 ; [arcsec^2]
;             
;             maggies = total(model[pix]*pixscale^2*(invvar[pix] gt 0))
;             errmaggies = sqrt(total(var[pix]*pixscale^4*(invvar[pix] gt 0)))
;             if errmaggies le 0 then message, 'Problem here!'
;
;             print, rin[ir], rout[ir], maggies, errmaggies, totnpix, area
;                
;;            if maggies gt 0.0 then begin
;;               phot1.abmag[ir] = -2.5*alog10(maggies) ; [AB maggies]
;;               phot1.abmag_err[ir] = 2.5*errmaggies/maggies/alog(10)
;;            endif
;
;             phot1.maggies[ir] = maggies ; [AB maggies]
;             phot1.ivarmaggies[ir] = 1.0/errmaggies^2
;          endfor
;          if n_elements(phot) eq 0 then phot = phot1 else phot = [phot,phot1]
;;         djs_plot, phot1.radius_kpc, -2.5*alog10(phot1.maggies), ysty=3, psym=8, /xlog, xsty=3
;       endfor                   ; close filter loop
;; write everything out
;;      im_mwrfits, imellipse, ellpath+cluster+'-ellipse-image.fits', clobber=clobber
;       im_mwrfits, modellipse, ellpath+cluster+'-ellipse-model.fits', clobber=clobber
;       im_mwrfits, phot, ellpath+cluster+'-ellipse-ellphot.fits', clobber=clobber
;    endfor                      ; close cluster loop

;          djs_plot, [0], [0], /nodata, /xlog, xsty=1, ysty=1, $
;            yr=[28,15], xrange=[pixscale*arcsec2kpc,300]
;          oploterror, imellipse.radius_kpc, imellipse.sb0fit, $
;            1.0/sqrt(imellipse.sb0fit_ivar), psym=symcat(9), $
;            symsize=0.5, color=cgcolor('orange')
;          oploterror, modellipse.radius_kpc, modellipse.sb0fit, $
;            1.0/sqrt(modellipse.sb0fit_ivar), psym=symcat(9), $
;            symsize=0.5, color=cgcolor('dodger blue')
;          djs_oplot, imellipse.radius_kpc, bcgmstar_sersic_func(imellipse), $
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
;          
;;         splog, 'Time to fit ', t0-systime(1)
;;         bgt_ellipse_show, model, ellipse
;          





;;; ellipse-fit the model; note that the ellipse PA (in radians) in GZ's code is
;;; measured positive from the X-axis, not the *Y-axis* as we assume
;;; below
;;          modellipse1 = do_ellipse(model,invvar=invvar,ini_xcen=xcen, $
;;            ini_ycen=ycen,ini_e0=info[reffilt].mge_ellipticity,$
;;            ini_pa0=info[reffilt].mge_posangle,namax=namax,$
;;            pixscale=pixscale,arcsec2kpc=arcsec2kpc)
;;
;;; optionally use the identical ellipse parameters to generate the
;;; surface-brightness profile from the data themselves
;;          constrained_fit = 0
;;
;;          if constrained_fit then begin
;;             imellipse1 = modellipse1
;;             for ia = 0, modellipse1.na-1 do begin
;;                sbcoord = bgt_ellipse_sb(image,imivar=invvar,xcen=modellipse1.xcenfit[ia],$
;;                  ycen=modellipse1.xcenfit[ia],pa=modellipse1.pafit[ia],$
;;                  ee=modellipse1.ellipticityfit[ia],majora=modellipse1.majora[ia],$
;;                  ea=ea,sbcoord_ivar=sbcoord_ivar)
;;                bgt_ellipse_fit, ea, sbcoord, ecoeff, ini_value=ini_value, $
;;                  sb_ivar=sbcoord_ivar, status=status, fixed=[0,1,1,1,1]
;;                imellipse1.sb0[ia] = ecoeff.sb0 ; only free parameter is the SB
;;                imellipse1.sb0_ivar[ia] = median(sbcoord_ivar)
;;             endfor
;;             bgt_ellipse_allfit, image, imellipse1, imivar=invvar
;;          endif else begin
;;             imellipse1 = do_ellipse(image,invvar=invvar*mask,ini_xcen=xcen,$
;;               ini_ycen=ycen,ini_e0=info[reffilt].mge_ellipticity,$
;;               ini_pa0=info[reffilt].mge_posangle,namax=namax,$
;;               pixscale=pixscale,arcsec2kpc=arcsec2kpc)
;;          endelse
;;
