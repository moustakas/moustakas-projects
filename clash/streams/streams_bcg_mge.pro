function mge_bcg_fit, image, invvar=invvar, minlevel=minlevel, $
  badpixels=badpixels, pixscale=pixscale, nprint=nprint, twist=twist, $
  nofit=nofit, mgeinfo=mgeinfo, ngauss=ngauss

    forward_function multi_gauss, multi_gauss_twist
    resolve_routine, 'mge_print_contours', /compile_full_file, /no_recompile
    resolve_routine, 'mge_print_contours_twist', /compile_full_file, /no_recompile

    n_sectors = 25
;   sector_width = 10.0
    sigmapsf = 0.0
    
    if n_elements(nprint) eq 0 then npring = 100
    
    if keyword_set(nofit) then begin
; scale by the inverse variance to get the new model
       stop       
    endif else begin
; (re)find the galaxy
       find_galaxy, image, size, ellipticity, posangle, xcen, ycen, $
         xcen_lum, ycen_lum, fraction=0.1, index=index, level=level, $
         nblob=nblob, plot=plot;, /quiet
; get the photometry
       xcen_lum1 = xcen_lum
       ycen_lum1 = ycen_lum
       if keyword_set(twist) then begin
          sectors_photometry_twist, image, posangle, xcen_lum1, ycen_lum1, $
            radius, phi, counts, n_sectors=n_sectors, sector_width=sector_width, $
            badpixels=badpixels, minlevel=minlevel
       endif else begin
          sectors_photometry, image, ellipticity, posangle, xcen_lum1, ycen_lum1, $
            radius, phi, counts, n_sectors=n_sectors, sector_width=sector_width, $
            badpixels=badpixels, minlevel=minlevel
       endelse

; do the fit
       if keyword_set(twist) then begin
          mge_fit_sectors_twist, radius, phi, counts, ellipticity, ngauss=ngauss, $
            scale=pixscale, negative=negative, normpsf=normpsf, print=print, $
            qbounds=qbounds, rbounds=rbounds, sigmapsf=sigmapsf, $
            outer_slope=outer_slope, sol=sol, absdev=absdev
          model = multi_gauss_twist(sol,image,0.0,1.0,xcen_lum1,ycen_lum1,posangle)
       endif else begin
          mge_fit_sectors, radius, phi, counts, ellipticity, ngauss=ngauss, $
            scale=pixscale, bulge_disk=bulge_disk, fastnorm=fastnorm, $
            linear=linear, negative=negative, normpsf=normpsf, print=print, $
            qbounds=qbounds, rbounds=rbounds, /quiet, sigmapsf=sigmapsf, $
            outer_slope=outer_slope, sol=sol, absdev=absdev, nprint=nprint
          model = float(multi_gauss(sol,image,0.0,1.0,xcen_lum1,ycen_lum1,posangle))
       endelse

; pack everything into a structure
       mgeinfo = {$
         xcen:         float(xcen),$ ; luminosity-weighted x,y centroid 
         ycen:         float(ycen),$
         xcen_lum:     float(xcen_lum),$ ; luminosity-weighted x,y centroid 
         ycen_lum:     float(ycen_lum),$
         size:         float(size*pixscale),$ ; approximate galaxy "size" [arcsec]
;        size_kpc:     float(size*pixscale*arcsec2kpc),$ ; approximate galaxy "size" [kpc]
         ellipticity:  float(ellipticity),$ ; average ellipticity = 1-b/a
         posangle:     float(posangle),$    ; position angle measured from the image Y-axis
         minlevel:     minlevel,$ 
         absdev:       absdev,$ ; mean absolute deviation of the data from the model
         sol:          sol}
    endelse
return, model
end

pro streams_bcg_mge, debug=debug
; jm13sep04siena - do the preliminary BCG subtraction using MGE 

; this code doesn't do what it should; the MGE model needs to
; be fitten in each band independently; i should also remove the
; dependence on BCG_INFO.FITS.

; note! images in units of [10^-12 erg/s/cm^2/Hz] (pico-maggies)

    qapath = streams_path(/bcg)+'qaplots/'
    
; read the sample and then match against the info structure written by
; streams_find_bcg 
    sample = rsex(streams_path(/propath)+'streams_sample.sex')
    splog, 'IGNORING A2261!!!'
    keep = where(strtrim(sample.shortname,2) ne 'a2261')
    sample = sample[keep]
;   sample = sample[5] ; =MACS1206
    struct_print, sample

    allinfo = mrdfits(streams_path()+'bcg_info.fits.gz',1)
    match, strtrim(sample.shortname,2), strtrim(allinfo.cluster,2), m1, m2
    srt = sort(m1)
    sample = sample[m1[srt]]
    info = allinfo[m2[srt]]
    struct_print, info

    ncl = n_elements(sample)

    pixscale = 0.065D ; [arcsec/pixel]
    rmaxkpc = 30D     ; [kpc]

; wrap on each cluster    
    for ic = 0, ncl-1 do begin
;   for ic = 5, 5 do begin
       cluster = strtrim(sample[ic].shortname,2)
       splog, 'Working on cluster '+cluster
       skypath = streams_path(/skysub)+cluster+'/'
       outpath = streams_path(/bcg)+cluster+'/'
       if file_test(outpath,/dir) eq 0 then file_mkdir, outpath

       bcgqafile = qapath+cluster+'_bcg.ps'
       qafile = qapath+cluster+'_all.ps'

; get a cutout that scales with the size of the galaxy       
       arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
;      rmax = ceil(rmaxkpc/arcsec2kpc/pixscale)         ; [pixels]
;      rmax = rmax+odd(rmax)                            ; force even
       rmax = long(15.0*info[ic].size/pixscale) ; [pixels]
       rmax = rmax+odd(rmax)                    ; force even

; read the skyinfo structure to get the filters
       skyinfo = mrdfits(streams_path(/skysub)+'skyinfo-'+$
         cluster+'.fits.gz',1,/silent)
       short = strtrim(skyinfo.band,2)
       reffilt = where(short eq 'f160w') ; reference filter
       nfilt = n_elements(skyinfo)
       scales = fltarr(nfilt)+1 ; this could be generalized
       niceprint, skyinfo.band, skyinfo.sigma, skyinfo.sblimit

       refsig = skyinfo[reffilt].sigma
       
; read the PSFs and the reference image (where we will do the Sersic
; modeling) 
       psffiles = file_search(streams_path(/psfs)+short+'-bpsf.fits*')
       psf = gz_mrdfits(psffiles[reffilt],/silent)

       imfile = skypath+cluster+'-'+short[reffilt]+'.fits.gz'
       splog, 'Reading '+imfile
       rimage = mrdfits(imfile,0,rhdr,/silent)
       rinvvar = mrdfits(imfile,1,ivarrhdr,/silent)

; get a RMAXKPC by RMAXKPC cutout centered on the BCG
       xcen = info[ic].xcen
       ycen = info[ic].ycen
       
       hextract, rimage, rhdr, cutrimage, cutrhdr, xcen-rmax, $
         xcen+rmax, ycen-rmax, ycen+rmax, /silent
       hextract, rinvvar, ivarrhdr, cutrinvvar, cutivarrhdr, $
         xcen-rmax, xcen+rmax, ycen-rmax, ycen+rmax, /silent
       cutsz = size(cutrimage,/dim)

; find other galaxies in the field (except the BCG!) and mask them out 
;      simage = dsmooth(cutrimage,5L)
;      dpeaks, simage, xc=xc, yc=yc, sigma=refsig, saddle=10.0, $
;        minpeak=10.0*refsig, /refine, /check
;      dist = sqrt((xc-rmax)^2+(yc-rmax)^2)
;      wmask = where(dist gt min(dist),nmask)
;      for im = 0, nmask-1 do begin ; more than one?
;         dmeasure, cutrimage, cutrinvvar, xcen=xc[wmask[im]], $
;           ycen=yc[wmask[im]], measure=measure1, /fixcen
          
;         bgt_ellipse, model, imivar=invvar, ini_xcen=xcen_lum, $
;           ini_ycen=ycen_lum, ini_e0=info[ic].ellipticity, $
;           ini_pa0=info[ic].posangle, ellipse=ell
;         bgt_ellipse_sersic, ell, outellipse=final
;      endfor
       
;      cgimage, cutrimage, /keep, stretch=2, /save, /neg
;      djs_oplot, xc, yc, psym=symcat(9), color=cgcolor('red'), symsize=1.5

; do the free-parameter modeling in the reference band; iterate the
; fit twice, each time improving the masking of nearby objects
       dobjects, cutrimage, objects=obj, dpsf=dpsf, plim=25.0, nlevel=1L
       ignore = where((obj ne -1 and obj ne obj[cutsz[0]/2,cutsz[1]/2]) or $
         (cutrinvvar eq 0),nignore)
       if nignore eq 0L then delvarx, ignore
       
; initial fit       
       cutrimage_mge = cutrimage ; copy because MGE messes with the image
       cutrimage_mge[ignore] = 0.0
       refmodel = mge_bcg_fit(cutrimage_mge,invvar=cutrinvvar,$
         badpixels=ignore,minlevel=refsig,pixscale=pixscale,$
         nprint=100,nofit=nofit,mgeinfo=mgeinfo)

; now detect objects in the (smoothed) residual image but "put back"
; masked pixels located in regions where the BCG *model* is bright
;      dobjects, dsmooth(cutrimage-refmodel,5L), objects=obj, plim=25, nlevel=1L
;      restore = where(refmodel gt 50.0*refsig,nrestore)
;      obj[restore] = -1

; now refit       
       ignore = where((obj ne -1 and obj ne obj[cutsz[0]/2,cutsz[1]/2]) or $
         (cutrinvvar eq 0),nignore)
       if nignore eq 0L then delvarx, ignore
       
       cutrimage_mge = cutrimage ; copy because MGE messes with the image
       cutrimage_mge[ignore] = 0.0
       refmodel = mge_bcg_fit(cutrimage_mge,invvar=cutrinvvar,$
         badpixels=ignore,minlevel=refsig,pixscale=pixscale,$
         nprint=100,nofit=nofit,mgeinfo=mgeinfo) ;,/twist)
       
       mgeinfo = struct_addtags({cluster: cluster, z: info[ic].z, $
         ra: info[ic].ra, dec: info[ic].dec},mgeinfo)
       
       outfile = outpath+cluster+'-refmodel-mgeinfo.fits'
       im_mwrfits, mgeinfo, outfile, /clobber
       
       cutrimage_mge[ignore] = randomn(seed,nignore)*refsig

; build a 4-panel QAplot
       im_plotconfig, 5, pos, psfile=bcgqafile, xspace=0.02, yspace=0.02, xmargin=[0.4,0.4]
       cgloadct, 0, /silent
       cgimage, cutrimage, clip=3, /negative, stretch=2, $
         margin=0, /keep_aspect, position=pos[*,0], /save
;      cgplots, mgeinfo.xcen_lum, mgeinfo.ycen_lum, psym=symcat(7,thick=2), color='red'
       xyouts, pos[0,0]+0.02, pos[3,0]-0.05, strupcase(cluster), /norm
       cgimage, cutrinvvar, clip=3, /negative, /noerase, $
         stretch=2, margin=0, /keep_aspect, position=pos[*,1]
;      cgimage, cutrimage_mge, clip=3, /negative, /noerase, $
;        stretch=2, margin=0, /keep_aspect, position=pos[*,1]
       cgimage, refmodel, clip=3, /negative, /noerase, $
         stretch=2, margin=0, /keep_aspect, position=pos[*,2]
       cgimage, cutrimage-refmodel, clip=3, /negative, /noerase, $
         stretch=2, margin=0, /keep_aspect, position=pos[*,3]
       im_plotconfig, /psclose, psfile=bcgqafile, /pdf

; write out
       outfile = outpath+cluster+'-refmodel.fits'
       splog, 'Writing '+file_basename(outfile)

       mwrfits, cutrimage, outfile, cutrhdr, /create, /silent
       mwrfits, refmodel, outfile, cutrhdr, /silent
       mwrfits, cutrinvvar, outfile, cutivarrhdr, /silent
       spawn, 'gzip -f '+outfile
       
;      im_plotconfig, 0, pos, psfile=qafile
;      pos = reform(im_getposition(nx=2,ny=nfilt,xspace=0.0,$
;        yspace=0.0,xmargin=[0.2,0.2],ymargin=[0.2,0.2],$
;        height=replicate(2.0,nfilt),width=2*[1,1]),4,2,nfilt)
       
; now loop through each band and scale the fit

; this is WRONG!  it doesn't give us the color gradients; the
; bands need to be fitted independently...       
       
       for ib = 0, nfilt-1 do begin
          imfile = skypath+cluster+'-'+short[ib]+'.fits.gz'
          splog, 'Reading '+imfile
          image = gz_mrdfits(imfile,0,hdr,/silent)
          invvar = gz_mrdfits(imfile,1,ivarhdr,/silent)

          hextract, image, hdr, cutimage, cuthdr, xcen-rmax, $
            xcen+rmax, ycen-rmax, ycen+rmax, /silent
          hextract, invvar, ivarhdr, cutinvvar, cutivarhdr, $
            xcen-rmax, xcen+rmax, ycen-rmax, ycen+rmax, /silent

          scale = total(cutinvvar*cutimage*refmodel)/total(refmodel^2*cutinvvar)
          model = scale*refmodel

;          cgimage, cutimage, clip=3, /negative, noerase=ib gt 0, $
;            stretch=2, margin=0, keep_aspect=1, position=pos[*,0,ib]
;;         cgimage, model, clip=3, /negative, /noerase, $
;;           stretch=2, margin=0, /keep_aspect, position=pos[*,ib,1]
;          cgimage, cutimage-model, clip=3, /negative, /noerase, $
;            stretch=2, margin=0, keep_aspect=1, position=pos[*,1,ib]
;          xyouts, pos[0,0,ib]-0.02, (pos[3,0,ib]-pos[1,0,ib])/2.0+pos[1,0,ib], $
;            strupcase(short[ib]), align=0.5, orientation=90, /norm, charsize=1.0

;         cgimage, cutimage-model, clip=3, /negative, stretch=2, $
;           margin=0, keep_aspect=1
;         cc = get_kbrd(1)
          
; write out
          outfile = outpath+cluster+'-'+short[ib]+'.fits'
          splog, 'Writing '+file_basename(outfile)

          mwrfits, cutimage, outfile, cuthdr, /create, /silent
          mwrfits, model, outfile, cuthdr, /silent
          mwrfits, cutinvvar, outfile, cutivahdr, /silent
          spawn, 'gzip -f '+outfile
       
;          dsersic, image, invvar, xcen=xcen, ycen=ycen, sersic=curr_sersic, $
;            /onlyflux, /fixcen, /fixsky, psf=psf, model=model
;;         dsersic2, image, invvar, xcen=xcen, ycen=ycen, sersic=curr_sersic, $
;;           /nofit, /fixcen, /fixsky, psf=psf, model=model, bulge=bulge, disk=disk
;             
;          mge1 = streams_mge(cutimage,badpixels=badpixels,$
;            pixscale=pixscale,model=cutmodel,twist=twist)
;;          help, mge1, /str
;
;; rebuild the original image with the BCG subtracted and write out 
;          model = image*0.0
;          embed_stamp, model, cutmodel, xcen-rmax, ycen-rmax
;;         atv, image-model, /bl
;
;          outfile = outpath+cluster+'-'+short[ib]+'.fits'
;          splog, 'Writing '+file_basename(outfile)
;
;          mwrfits, image-model, outfile, hdr, /create
;          mwrfits, invvar, outfile, ivarhdr, /silent
;          spawn, 'gzip -f '+outfile
       endfor                   ; close filter loop
;      im_plotconfig, psfile=qafile, /psclose, /pdf
    endfor                      ; close cluster loop

return
end
    
