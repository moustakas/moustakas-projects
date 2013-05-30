pro streams_bcg_mge
; jm13may23siena - fit the BCGs

; rsync -auv coma2:"/clash-archive/clash_archive/macs1206/HST/images/mosaicdrizzle_image_pipeline/scale_65mas/macs1206_mosaic_065mas_wfc3ir_f160w_drz_20110815.fits.gz" ./
; rsync -auv coma2:"/clash-archive/clash_archive/macs1206/HST/images/mosaicdrizzle_image_pipeline/scale_65mas/macs1206_mosaic_065mas_wfc3ir_f160w_wht_20110815.fits.gz" ./
; rsync -auv coma2:"/clash-archive/clash_archive/macs1206/HST/catalogs/mosaicdrizzle_image_pipeline/IR_detection/SExtractor/detectionImage_SEGM.fits.gz" ./
; rsync -auv coma2:"/clash-archive/clash_archive/macs1206/HST/catalogs/mosaicdrizzle_image_pipeline/IR_detection/SExtractor/macs1206_f160w.cat" ./

    cluster = 'macs1206'
    mosaicpath = getenv('IM_ARCHIVE_DIR')+'/'+cluster+'/images/'
    catpath = getenv('IM_ARCHIVE_DIR')+'/'+cluster+'/catalogs/'

; cluster and filter properties    
    clash = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    this =  where('macs1206' eq strtrim(clash.cluster_short,2))
    arcsec2kpc = dangular(clash[this].z,/kpc)/206265D ; [kpc/arcsec]
    ebv = clash[this].ebv

    filt = bcgimf_filterlist(short=short,instr=instr,weff=weff,zpt=zpt)
    kl = k_lambda(weff,/odon)
    nfilt = n_elements(filt)

    pixscale = 0.065D           ; [arcsec/pixel]
    rmaxkpc = 50D ; [kpc]
    rmax1 = rmaxkpc/arcsec2kpc/pixscale ; [pixels]
    
    image = mrdfits(mosaicpath+'macs1206_mosaic_065mas_wfc3ir_f160w_drz_20110815.fits.gz',0,hdr)
    invvar = mrdfits(mosaicpath+'macs1206_mosaic_065mas_wfc3ir_f160w_wht_20110815.fits.gz')
    ib = 12 ; [=F160W]
    
    extast, hdr, astr
    ad2xy, clash[this].ra_bcg, clash[this].dec_bcg, astr, xcen, ycen
    
    dispimage = im_imgscl((image*(invvar gt 0))[(xcen-rmax1):(xcen+rmax1),$
      (ycen-rmax1):(ycen+rmax1)],/neg,/sqrroot)
    
    plotimage, dispimage, /noaxes, /preserve_aspect, /norm

; mask everything but the BCG

; get initial estimates of the galaxy shape
    cutimage = image[(xcen-rmax1):(xcen+rmax1),(ycen-rmax1):(ycen+rmax1)]
    sz = size(cutimage,/dim)

    find_galaxy, cutimage, majoraxis, eps, ang, xpeak, ypeak, xmed, ymed, $
      fraction=fraction, index=index, level=level

; perform photometry    
    sectors_photometry, cutimage, eps, ang, xpeak, ypeak, radius, angle, counts, $
      n_sectors=n_sectors, sector_width=sector_width, badpixels=badpixels, $
      minlevel=minlevel

; do the MGE fitting
    mge_fit_sectors, radius, angle, counts, eps, ngauss=ngauss, scale=pixscale, $
      /print, sol=sol, absdev=absdev, outer_slope=outer_slope

; convert the solution to surface brightness and build a QAplot
    mge_print_contours, cutimage, ang, xpeak, ypeak, sol, model=model

    sbprofile = sol[0,*]/(2.0*!dpi*sol[1,*]^2*sol[2,*])
    sbmodel = -2.5*alog10(model)+5*alog10(pixscale)+zpt[ib]-kl[ib]*ebv  
    sbimage = -2.5*alog10(cutimage)+5*alog10(pixscale)+zpt[ib]-kl[ib]*ebv  

    axis = (rebin(reform(findgen(sz[0]),sz[0],1),sz[0],sz[1])-xpeak)*pixscale*arcsec2kpc ; [kpc]
    djs_plot, axis[*,135], sbimage[*,135], psym=8, ysty=3, yr=[24,18]
    djs_oplot, axis[*,135], sbmodel[*,135], line=0, color='orange'

    jj = read_bcg_profiles('macs1206')
    djs_oplot, jj[15].sma, jj[15].mu, psym=8, color='blue'
    
    
    
stop    

return
end
