pro flamingos_distort, npix=npix, order=order, xpix=xpix, ypix=ypix, xdist=xdist, ydist=ydist
; distortion model for FLAMINGOS
    
; if "align=center" then offset = NPIX/2L+1L; if "align=corner" then
; offset = NPIX/2L+0.5 

    if (n_elements(npix) eq 0L) then npix = 2048L
    if (n_elements(order) eq 0L) then order = 6L 

    crpix1 = npix/2L + 0.5      ; reference pixels [1-indexed]
    crpix2 = npix/2L + 0.5

    if (n_elements(xpix) eq 0L) then xpix = findgen(npix)+1.0 ; native pixel coordinates [1-indexed]
    if (n_elements(ypix) eq 0L) then ypix = findgen(npix)+1.0
    
    xdist_coeff = [$
      -3.13385180345,$
      1.004072924,$
      -0.00183808434959,$
      1.90899554759D-06,$
      -7.02600909987D-07,$
      1.31772038336D-06,$
      -1.46253633834D-08,$
      8.74646975868D-10,$
      -1.27477572409D-08,$
      3.84044359946D-10,$
      -9.23163130406D-13,$
      -3.02699688831D-13,$
      -2.17804944819D-12,$
      7.13140243591D-13,$
      -4.22992633579D-13,$
      7.50779676661D-15,$
      -1.05067728167D-16,$
      1.32865593246D-14,$
      -1.39547603335D-15,$
      4.79680525157D-15,$
      5.27901629365D-16,$
      -5.28209885681D-19,$
      -2.15829018477D-19,$
      -1.69082032126D-19,$
      2.98401155856D-19,$
      1.42989656823D-19,$
      -2.34978528457D-19,$
      -1.4139870516D-19]

    ydist_coeff = [$
      -1.59317700345,$
      -7.58792414804D-05,$
      1.00403420094,$
      -2.87917603493D-08,$
      1.4431725501D-08,$
      -4.00634996454D-06,$
      -2.46081186836D-10,$
      -1.53044126738D-08,$
      8.91239390236D-10,$
      -1.53282151607D-08,$
      -3.25515890748D-14,$
      -3.40198566384D-13,$
      2.21451114579D-12,$
      -8.13094365249D-13,$
      3.39460425228D-12,$
      2.9734493965D-16,$
      7.34769704289D-15,$
      2.74279596123D-16,$
      1.33306619088D-14,$
      -1.57274385218D-15,$
      7.21586951996D-15,$
      5.0088021963D-19,$
      3.53436803775D-19,$
      -4.51355862217D-19,$
      -1.471392533D-19,$
      2.80336330843D-19,$
      -2.1021589111D-19,$
      -3.01518328277D-19]

    xdist = xpix*0.0            ; distorted pixel coordinates, relative to CRPIX? [1-indexed]
    ydist = ypix*0.0

    n = 0L
    for j = 0L, order do for l = 0L, j do begin
       xdist = xdist + xdist_coeff[n]*(xpix-crpix1)^(j-l)*(ypix-crpix2)^l
       ydist = ydist + ydist_coeff[n]*(xpix-crpix1)^(j-l)*(ypix-crpix2)^l
;      print, n, j, l, j-l, xdist_coeff[n], ydist_coeff[n]
       n = n+1L
    endfor

;   plot, xpix, xpix-(xdist+crpix1), xsty=3
;   oplot, ypix, ypix-(ydist+crpix2), line=2

    xdist = xdist + crpix1      ; NOTE!
    ydist = ydist + crpix2
       
return
end

pro flamingos_ahead_global
; jm07may09nyu - generate the global .ahead file for FLAMINGOS for
;                input into SCAMP

    sexpath = sg1120_path(/sex)
    rootpath = flamingos_path()
    datapath = flamingos_path(/mar06)

; use SCAMP to convert the distortion model into an appropriate
; "global" external (.ahead) FITS header

    order = 6L
    npix = 2048L

    factor = 128.0              ; 64.0 ; 32.0 ; 16.0
    nstar = npix/factor
    xstar = findgen(nstar)/nstar*npix+factor/2.0 ; undistorted coordinates
    ystar = xstar

    flamingos_distort, npix=npix, order=order, xpix=xstar, ypix=ystar, xdist=xstar_dist, ydist=ystar_dist

; write out distorted and undistorted images
    
    dim = 100L & sigma = 5L
    star = mkgaussian(dim,sigma)
    star = star/total(star)
    
    image = fltarr(npix,npix)   ; undistorted image
    for ixstar = 0L, nstar-1L do for iystar = 0L, nstar-1L do embed_stamp, image, star, xstar[ixstar]-dim/2.0, ystar[iystar]-dim/2.0
    mkhdr, hdr, image, image=0
    astr = hogg_make_astr(0.0,30.0,npix/3600.0,npix/3600.0,pixscale=1.0/3600.0)
    putast, hdr, astr       
;      writefits, rootpath+'pinhole.undistorted.fits', image, hdr

    image_dist = fltarr(npix,npix) ; undistorted image
    for ixstar = 0L, nstar-1L do for iystar = 0L, nstar-1L do embed_stamp, image_dist, star, xstar_dist[ixstar]-dim/2.0, ystar_dist[iystar]-dim/2.0
;      writefits, rootpath+'pinhole.distorted.fits', image_dist, hdr

; write out the undistorted coordinates as a FITS-LDAC type catalog
; for use as the astrometric reference catalog in SCAMP

    cat = replicate({x_world: 0.0D, y_world: 0.0D, erra_world: 0.0D, errb_world: 0.0D, $
      mag: 10.0D, magerr: 0.0D},nstar*nstar)
    xy2ad, reform(xstar # (fltarr(nstar)+1.0),nstar*nstar), reform(ystar # (fltarr(nstar)+1.0),nstar*nstar), astr, aa, dd
    cat.x_world = aa
    cat.y_world = dd
;      struct_print, cat
    mkhdr, cathdr, image, /extend
    sxaddpar, cathdr, 'BITPIX', 0
    sxdelpar, cathdr, 'COMMENT'
    sxdelpar, cathdr, 'DATE'
    putast, cathdr, astr
    wfitsldac, cat, fits_header=cathdr, outfile=rootpath+'pinhole.undistorted.cat.fits'
    
    mkhdr, cathdr, image, /extend
    sxaddpar, cathdr, 'BITPIX', 0
    putast, cathdr, astr
    mwrfits, 0, rootpath+'pinhole.undistorted.cat.fits', /create
    mwrfits, {field_header_card: strjoin(cathdr)}, rootpath+'pinhole.undistorted.cat.fits'
    mwrfits, cat, rootpath+'pinhole.undistorted.cat.fits'

;   cat_dist = replicate({xwin_image: 0.0D, ywin_image: 0.0D, errawin_image: 0.0D, errbwin_image: 0.0D, $
;     errthetawin_image: 0.0D, flux_auto: 1.0D, fluxerr_auto: 1D-6},nstar*nstar)
;   cat_dist.xwin_image = reform(xstar_dist # (fltarr(nstar)+1.0),nstar*nstar)
;   cat_dist.ywin_image = reform(ystar_dist # (fltarr(nstar)+1.0),nstar*nstar)
;   wsex, cat_dist, outfile=rootpath+'pinhole.distorted.cat'
;   mwrfits, 0, rootpath+'pinhole.distorted.cat.fits', /create
;   mwrfits, {field_header_card: 'This is a FITS-LDAC file'}, rootpath+'pinhole.distorted.cat.fits'
;   mwrfits, cat_dist, rootpath+'pinhole.distorted.cat.fits'
       
       spawn, 'sex '+rootpath+'pinhole.distorted.fits -c '+sexpath+'sg1120.sex -CATALOG_NAME '+rootpath+'pinhole.distorted.cat.fits'+$
         ' -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME '+sexpath+'sg1120.sex.param -STARNNW_NAME '+sexpath+'default.nnw'+$
         ' -FILTER N -CLEAN N -WEIGHT_TYPE NONE -BACK_TYPE MANUAL -BACK_VALUE 0.0 -CHECKIMAGE_TYPE NONE'

;      spawn, 'sex '+rootpath+'pinhole.undistorted.fits -c '+sexpath+'sg1120.sex -CATALOG_NAME '+rootpath+'pinhole.undistorted.cat.fits'+$
;        ' -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME '+sexpath+'sg1120.sex.param -STARNNW_NAME '+sexpath+'default.nnw'+$
;        ' -FILTER N -CLEAN N -WEIGHT_TYPE NONE -BACK_TYPE MANUAL -BACK_VALUE 0.0 -CHECKIMAGE_TYPE NONE'

       spawn, 'scamp '+rootpath+'pinhole.distorted.cat.fits -ASTREF_CATALOG FILE -ASTREFCAT_NAME '+rootpath+'pinhole.undistorted.cat.fits'+$
         ' -ASTREFCENT_KEYS XWIN_IMAGE,YWIN_IMAGE -CHECKPLOT_DEV NULL -DISTORT_DEGREES '+strtrim(order,2)+$
         ' -HEADER_TYPE NORMAL -MERGEDOUTCAT_TYPE NONE -SAVE_REFCATALOG N -SOLVE_ASTROM Y -SOLVE_PHOTOM N -WRITE_XML N', /sh

; test       
       spawn, 'scamp '+rootpath+'pinhole.distorted.cat.fits -ASTREF_CATALOG USNO-B1 -ASTREFCAT_NAME '+rootpath+'pinhole.undistorted.cat.fits'+$
         ' -ASTREFCENT_KEYS XWIN_IMAGE,YWIN_IMAGE -CHECKPLOT_DEV NULL -DISTORT_DEGREES '+strtrim(order,2)+$
         ' -HEADER_TYPE NORMAL -MERGEDOUTCAT_TYPE NONE -SAVE_REFCATALOG Y -SOLVE_ASTROM Y -SOLVE_PHOTOM N -WRITE_XML N', /sh
       
stop       
       
    endif
    
    catlist = file_search(datapath+'ss_sg1120[a,b]_ks.????.cat',count=nimage)
;   catlist = file_search(datapath+'sra.sg1120*.cat')
    
    root = 'flamingos'
    scampconfig = sexpath+'sg1120.scamp'

    if keyword_set(usno) then astref_catalog = 'USNO-B1' else astref_catalog = '2MASS'
    refout_catpath = datapath
    mergedoutcat_name = root+'.scamp.cat'
    astrinstru_key = 'FILTER'
    photinstru_key = 'FILTER'
    magzero_key = 'PHOT_C'
    extinct_key = 'PHOT_K'
    photomflag_key = 'PHOTFLAG'
    degree = '1'

    for iter = 0L, 1L do begin 

       if (iter eq 0L) then aheader_suffix = '.ahead' else aheader_suffix = '.head'
       case iter of
          0L: begin
             position_maxerr = '5.0'
             mosaic_type = 'UNCHANGED'
          end
          1L: begin
             position_maxerr = '0.5'
             mosaic_type = 'FIX_FOCALPLANE'
          end
          2L: begin
             position_maxerr = '0.1'
          end
          else: position_maxerr = '0.1'
       endcase

       spawn, 'scamp '+strjoin(catlist,',')+' -c '+scampconfig+$
         ' -ASTREF_CATALOG '+astref_catalog+' -SAVE_REFCATALOG N -REFOUT_CATPATH '+refout_catpath+$
         ' -MERGEDOUTCAT_NAME '+mergedoutcat_name+' -MERGEDOUTCAT_TYPE NONE'+$
         ' -PIXSCALE_MAXERR 1.2 -POSANGLE_MAXERR 2.0 -POSITION_MAXERR '+position_maxerr+' -MOSAIC_TYPE '+mosaic_type+$
         ' -ASTRINSTRU_KEY '+astrinstru_key+' -DISTORT_DEGREES '+degree+$
         ' -SOLVE_PHOTOM Y -PHOTINSTRU_KEY '+photinstru_key+' -MAGZERO_KEY '+magzero_key+$
         ' -EXTINCT_KEY '+extinct_key+' -PHOTOMFLAG_KEY '+photomflag_key+$
         ' -WRITE_XML N -VERBOSE_TYPE NORMAL -AHEADER_SUFFIX '+aheader_suffix, /sh ; NOTE!

    endfor

; 
; ----- Astrometric stats (internal) :
; 
;                   All detections         |           High S/N           
;            dAXIS1  dAXIS2   chi2   ndets | dAXIS1  dAXIS2   chi2   ndets
; Group  1:  0.138"   0.12"     58   23023   0.095" 0.0527" 2.8e+02    4123
;   
; ----- Astrometric stats (external):
; 
;                   All detections         |           High S/N           
;            dAXIS1  dAXIS2   chi2  nstars | dAXIS1  dAXIS2   chi2  nstars
; Group  1:  0.269"  0.186"    1.2     142   0.226"  0.124"    1.8      38
; 
       
return
end
    
