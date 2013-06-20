function recenter, cutout, invvar, rinvvar=rinvvar, hdr=hdr
; recenter the star using sinc interpolation
    dpeaks, cutout, xcen=xx, ycen=yy, maxnpeaks=1, /refine
    xx = xx[0]+x0-nbigpix/2D
    yy = yy[0]+y0-nbigpix/2D
    xy2ad, xx, yy, astr, ra, dec
; remap the cutout so that it's centered in the stamp
    newastr = astr
    newastr.naxis = npix*[1,1]
    newastr.crval = [ra,dec]
    newastr.crpix = npix*[1,1]/2D ; center of the mosaic
    smosaic_remap, image1, astr, newastr, cutimage, offset, $
      weight=invvar, outweight=cutivar
return, rimage
end

function get_cutouts, cat, imfile=imfile, npix=npix
; get postage stamp cutouts without recentering
    nbigpix = npix*5
    pixel_scale = 0.065D
    wid = npix*pixel_scale/3600D ; cutout width [deg]
    nstar = n_elements(cat)
    ivarfile = repstr(imfile,'drz','wht')
    image = mrdfits(imfile,0,hdr,/silent)
    invvar = mrdfits(ivarfile,0,ivarhdr,/silent)
    extast, hdr, astr
    out = replicate({xcen: 0D, ycen: 0D, ra: 0D, dec: 0D, fracmask: 0.0, $
      image: fltarr(npix,npix), invvar: fltarr(npix,npix)},nstar)
    for ii = 0, nstar-1 do begin
       print, ii+1, nstar, string(13b), format='("Cutting out star ",I0,"/",I0,A10,$)'
       xx = cat[ii].x_image
       yy = cat[ii].y_image
       xy2ad, xx, yy, astr, ra, dec
; do a local sky subtraction
       x0 = long(cat[ii].x_image)
       y0 = long(cat[ii].y_image)
       hextract, image, hdr, bigcutout, bignewhdr, x0-nbigpix/2, $
         x0+nbigpix/2, y0-nbigpix/2, y0+nbigpix/2, /silent
       hextract, image, hdr, bigcutoutivar, bignewivarhdr, x0-nbigpix/2, $
         x0+nbigpix/2, y0-nbigpix/2, y0+nbigpix/2, /silent
       dobjects, bigcutout, obj=obj
       good = where(obj eq -1 and bigcutoutivar gt 0,ngood)
       if ngood eq 0 then stop
       mmm, bigcutout[good], skymode, /silent
; now cutout a smaller stamp
       hextract, image, hdr, cutout, newhdr, x0-npix/2, $
         x0+npix/2, y0-npix/2, y0+npix/2, /silent
       hextract, invvar, hdr, cutoutivar, newivarhdr, x0-npix/2, $
         x0+npix/2, y0-npix/2, y0+npix/2, /silent
       cutoutsub = cutout-skymode
       
       out[ii].xcen = xx
       out[ii].ycen = yy
       out[ii].ra = ra
       out[ii].dec = dec
       out[ii].image = cutoutsub
       out[ii].invvar = cutoutivar
       out[ii].fracmask = total(cutoutivar le 0.0)/(1.0*npix*npix)
;      plotimage, im_imgscl(cutimage>1E-5,/log,min=10,top=220), /preserve
    endfor
return, out
end

pro clash_psf_sex, imfile, sexpath=sexpath, catalog_name=catalog_name, $
  detect_imfile=detect_imfile, zpt=zpt, gain=gain
; build a basic SE catalog using the F160W band in double-image mode 
;   rmsfile = repstr(imfile,'drz','rms')
    weightfile = repstr(imfile,'drz','wht')
    detect_weightfile = repstr(detect_imfile,'drz','wht')

    config = init_sex_config()
    config.catalog_name = catalog_name
    config.parameters_name = sexpath+'psf_sex.param'
    config.filter_name = sexpath+'psfex.sex.conv'
    config.starnnw_name = sexpath+'default.nnw'
    
    config.detect_minarea = 4
    config.detect_thresh = 1.5
    config.analysis_thresh = 4
    config.phot_apertures = 1
    
    config.catalog_type = 'FITS_LDAC'
;   config.weight_type = 'MAP_RMS'
;   config.weight_image = rmsfile
    config.weight_type = 'MAP_WEIGHT'
    config.weight_image = detect_weightfile+','+weightfile
    
    config.mag_zeropoint = zpt
    config.satur_level = 249500
    config.pixel_scale = 0.065
    config.gain = gain

    im_sex, imfile, config, detect_imagelist=detect_imfile

return
end
    
pro build_clash_psfs, sexcatalogs=sexcatalogs, choose_stars=choose_stars, $
  build_psf=build_psf
; jm13may30siena - build the PSFs for CLASH in all 16 bands using all
; the clusters 

    splog, 'Do not forget to re-gzip the fullfield directory!!'
    
    clash = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    clash = clash[8] ; use just MACS1206 for now
    ncl = n_elements(clash)

    filt = bcgimf_filterlist(short=short,weff=weff,instr=instr,zpt=zpt)
    srt = reverse(sort(weff))
    filt = filt[srt]
    short = short[srt]
    instr = instr[srt]
    zpt = zpt[srt]
    nfilt = n_elements(filt)

    reffilt = (where(short eq 'f160w'))[0] ; reference filter
    pixel_scale = 0.065                    ; [arcsec/pixel]
    npix = 41                              ; cutout size
    
    path = getenv('IM_ARCHIVE_DIR')+'/projects/clash/psfs/' ; main I/O path
    
;   archive = getenv('CLASH_ARCHIVE')
    archive = getenv('IM_ARCHIVE_DIR')+'/'

; --------------------------------------------------
; build the SExtractor catalogs we need
       if keyword_set(sexcatalogs) then begin
          for ic = 0, ncl-1 do begin
             cluster = strtrim(clash[ic].shortname,2)
             if file_test(archive+strtrim(clash[ic].dirname,2),/dir) then begin
                mosaicpath = archive+strtrim(clash[ic].dirname,2)+'/HST/images/'+$
                  'mosaicdrizzle_image_pipeline/fullfield/'

; for GAIN see notes at
; http://www.ifa.hawaii.edu/~rgal/science/sextractor_notes.html
                detect_imfile = file_search(mosaicpath+cluster+'-fullfield_mosaic_065mas_'+$
                  instr[reffilt]+'_'+short[reffilt]+'_drz_????????.fits')
;               detect_imfile = file_search(mosaicpath+cluster+'_mosaic_065mas_'+$
;                 instr[reffilt]+'_'+short[reffilt]+'_drz_????????.fits')
                
                for ib = 0, nfilt-1 do begin
                   catalog_name = path+cluster+'-'+short[ib]+'.cat'
                   imfile = file_search(mosaicpath+cluster+'-fullfield_mosaic_065mas_'+$
                     instr[ib]+'_'+short[ib]+'_drz_????????.fits')
;                  imfile = file_search(mosaicpath+cluster+'_mosaic_065mas_'+$
;                    instr[ib]+'_'+short[ib]+'_drz_????????.fits')
                   hdr = headfits(imfile)
                   gain = sxpar(hdr,'ccdgain')*sxpar(hdr,'exptime')        
                   clash_psf_sex, imfile, sexpath=path, catalog_name=catalog_name, $
                     detect_imfile=detect_imfile, zpt=zpt[ib], gain=gain
                endfor
             endif
          endfor 
       endif
       
; --------------------------------------------------
; select candidate PSF stars
       if keyword_set(choose_stars) then begin
          magpivot = 23.0
          int = 26.5
          slope = 0.95
          magfaint = 25.0
          magaxis = range(15,magfaint,500)
; select stars from the F160W filter in each cluster
          for ic = 0, ncl-1 do begin
             cluster = strtrim(clash[ic].shortname,2)
             splog, 'Working on cluster '+cluster
             mosaicpath = archive+strtrim(clash[ic].dirname,2)+'/HST/images/'+$
               'mosaicdrizzle_image_pipeline/fullfield/'

             catfile = path+cluster+'-f160w.cat'
             cat = mrdfits(catfile,2,/silent)
             
             istar = where(cat.mag_auto gt 0.0 and cat.mag_auto lt magfaint and $
               cat.mag_aper lt poly(cat.mag_auto-magpivot,[int,slope]),nstar)
;            istar = where(cat.mag_auto gt 0.0 and cat.mag_auto lt magfaint and $
;              cat.mag_aper lt poly(cat.mag_auto-magpivot,[int,slope]) and $
;              cat.flags eq 0,nstar)
;            medpsf = djs_median(cat[istar].vignet,3)
             write_ds9_regionfile, cat.alpha_j2000, cat.delta_j2000, $
               file=path+cluster+'-'+short[reffilt]+'-allstars.reg', $
               color='green', symbol='box'
             write_ds9_regionfile, cat[istar].alpha_j2000, cat[istar].delta_j2000, $
               file=path+cluster+'-'+short[reffilt]+'-stars.reg', color='red'

;; stars that should have been flagged             
;             xshould
;             6164.0, 7057.0

; build a QAplot
             xrange = [15,28]
             psfile = path+'qaplot_'+cluster+'_stars.ps'
             im_plotconfig, 6, pos, psfile=psfile

             for ib = 0, nfilt-1 do begin
                ff = strupcase(short[ib])
                splog, ff
                catfile = path+cluster+'-'+short[ib]+'.cat'
                cat = mrdfits(catfile,2,/silent)
             
;               istar = where(cat.mag_auto gt 0.0 and cat.mag_auto lt magfaint and $
;                 cat.mag_aper lt poly(cat.mag_auto-magpivot,[int,slope]) and $
;                 cat.flags eq 0,nstar)
;               keep = where(cat[istar].flags eq 0,nstar)
                keep = where(cat[istar].mag_aper gt 0 and cat[istar].mag_aper lt 99.0,nstar)
                istar1 = istar[keep]

; get cutouts and write out                
                imfile = file_search(mosaicpath+cluster+'-fullfield_mosaic_065mas_'+$
                  instr[ib]+'_'+short[ib]+'_drz_????????.fits')
                cut = get_cutouts(cat[istar1],imfile=imfile,npix=npix)
                out = struct_addtags(cat[istar1],cut)

                im_mwrfits, out, repstr(catfile,'.cat','-stars.fits'), /clobber
;               medpsf = djs_median(cat[istar1].vignet,3)

; mag_aper vs mag_auto
                djs_plot, [0], [0], /nodata, position=pos[*,0], xrange=xrange, $
                  yrange=[32,15], xtickname=replicate(' ',10), ytitle=ff+' apmag', $
                  xsty=1, ysty=1
                djs_oplot, cat.mag_auto, cat.mag_aper, psym=8, symsize=0.7
                djs_oplot, cat[istar1].mag_auto, cat[istar1].mag_aper, psym=8, $
                  symsize=1.0, color='red'
                djs_oplot, magaxis, poly(magaxis-magpivot,[int,slope]), line=0
                djs_oplot, magfaint*[1,1], [poly(magfaint-magpivot,[int,slope]),!y.crange[1]]
                legend, [ff+' (N='+strtrim(nstar,2)+')'], $
                  /left, /top, box=0, margin=0
;               legend, [cluster,'N='+strtrim(nstar,2)], /left, /top, box=0, margin=0
; flux_radius vs mag_auto
                djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xrange=xrange, $
                  yrange=[0,15], xtitle=ff+' total mag', ytitle='r_{h} (pixels)', $
                  xsty=1, ysty=1
                djs_oplot, cat.mag_auto, cat.flux_radius, psym=8, symsize=0.7
                djs_oplot, cat[istar1].mag_auto, cat[istar1].flux_radius, $
                  psym=8, symsize=1.0, color='red'
             endfor
             im_plotconfig, psfile=psfile, /psclose, /pdf
          endfor 
       endif 

; --------------------------------------------------
; build the PSF; see atlas_fitpsf
       if keyword_set(build_psf) then begin
          rejsigma = [10.0,5.0,3.0]
          splog, 'Cut on FLAGS!'
stop          
          for ic = 0, ncl-1 do begin
             cluster = strtrim(clash[ic].shortname,2)
             splog, 'Working on cluster '+cluster
             for ib = 0, nfilt-1 do begin
                ff = strupcase(short[ib])
                splog, ff
                catfile = path+cluster+'-'+short[ib]+'-stars.fits.gz'
                cat = mrdfits(catfile,1,/silent)
                nstar = n_elements(cat)
                
; normalize to unity                
                norm = total(total(cat.image,1),1)
                norm1 = rebin(reform(norm,1,1,nstar),npix,npix,nstar)
                cat.image = cat.image/norm1
                cat.invvar = cat.invvar*norm1^2

                nkeep = nstar
                for iter = 0, n_elements(rejsigma)-1 do begin
; do initial fit
                   bpsf = djs_median(cat.image,3)
;                  bpsf = total(cat.image,3)
                   bpsf = bpsf/total(bpsf)
;                  dfit_mult_gauss, bpsf, 5, amp, psfsig, model=model, /quiet
                   
; clip the worst outliers
                   chi2 = fltarr(nkeep)
                   scale = fltarr(nkeep)
                   for ii = 0, nkeep-1 do begin  
                      scale[ii] = total(bpsf*cat[ii].image)/total(bpsf^2)
                      chi2[ii] = total((cat[ii].image-scale[ii]*bpsf)^2)
;                     chi2[ii] = total((cat[ii].image-scale[ii]*bpsf)^2*cat[ii].invvar)
;                     atv, cat[ii].image-scale[ii]*bpsf, /blo
                   endfor

stop                   
                   
; reject REJSIGMA outliers
                   dof = float(npix)^2
                   keep = where(chi2 lt dof+rejsigma[iter]*sqrt(2.0*dof),nkeep,comp=rej)
                   splog, 'REJSIGMA = ', rejsigma[iter], '; keeping ', nkeep, ' stars'
                   if (nkeep lt 3L) then message, 'Not enough good stars!'
stop
                   cat = cat[keep]
                endfor
                
stop                
             endfor
          endfor
       endif
       
return
end
    

;
;       if file_test(archive+strtrim(clash[ic].dirname,2),/dir) then begin
;          splog, 'Working on cluster '+clash[ic].cluster
;          catpath = archive+strtrim(clash[ic].dirname,2)+'/HST/catalogs/'+$
;            'mosaicdrizzle_image_pipeline/IR_detection/SExtractor/'
;          mosaicpath = archive+strtrim(clash[ic].dirname,2)+'/HST/images/'+$
;            'mosaicdrizzle_image_pipeline/scale_65mas/'
;          cat = rsex(catpath+'detectionImage.cat')
;          cat.mag_auto += 25
;          cat.mag_aper += 25
;
;          istar = where(cat.mag_auto gt 0.0 and cat.mag_auto lt magfaint and $
;            cat.mag_aper lt poly(cat.mag_auto-magpivot,[int,slope]) and $
;            cat.flags ne 0,nstar)
;          plot, [0], [0], /nodata, xrange=[15,30], yrange=[35,15], $
;            xtitle='mag_auto (AB mag)', ytitle='mag_aper (AB mag)'
;          djs_oplot, cat.mag_auto, cat.mag_aper, psym=8
;          djs_oplot, cat[istar].mag_auto, cat[istar].mag_aper, psym=8, color='green'
;          djs_oplot, magaxis, poly(magaxis-magpivot,[int,slope]), line=0
;          djs_oplot, magfaint*[1,1], [poly(magfaint-magpivot,[int,slope]),!y.crange[1]]
;;         legend, strupcase(short[ii]), /right, /top, box=0
;
;stop          
;          
