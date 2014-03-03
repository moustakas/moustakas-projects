; function isstar, mag_auto, mu_max
; ; choose stars
; return, indx
; end

pro write_thumbnail, outfile, pos=pos
    x0 = fix((convert_coord(pos[0:1],/normal,/to_device))[0])
    nx = fix((convert_coord(pos[2:3],/normal,/to_device))[0])-X0
    y0 = fix((convert_coord(pos[0:1],/normal,/to_device))[1])
    ny = fix((convert_coord(pos[2:3],/normal,/to_device))[1])-y0
    img = tvrd(x0,y0,nx,ny)
    tvlct, rr, gg, bb, /get
    write_png, outfile, img, rr, gg, bb
return
end

pro make_montage, cat, qafile=qafile, title=title
    nstar = n_elements(cat)
    ncols = 8
    nrows = ceil(nstar/float(ncols))
    npix = (size(cat[0].image,/dim))[0]
    
    for ii = 0, nstar-1 do begin
       outfile = '/tmp/star_'+string(cat[ii].number,format='(I4.4)')+'.png'
       loadct, 0, /silent
       image = cat[ii].image
;      crap = where(cat[ii].invvar eq 0,ncrap)
;      if ncrap ne 0L then image[crap] = djsig(image,sigrej=3.0)
       im = asinhscl(image,beta=0.005,omax=254,/negative)
;      im = gmascl(image,gamma=5.0,omax=254,/negative)
;      im = clipscl(image,3.0,omax=254);,percent=99.0,top=250)
;      im = sdevscl(image,mult=5.0,omax=254);,percent=99.0,top=250)
;      im = hist_equal(image,top=254);,percent=99.0,top=250)
       set_plot, 'Z'
       plotimage, im, /noaxes, /preserve_aspect, position=pos, /norm
       if cat[ii].x_image gt 5000 and cat[ii].x_image lt 10000 and $
         cat[ii].y_image gt 5000 and cat[ii].y_image lt 10000 then $
           primary = 1 else primary = 0
       if primary then plots, 3, npix-8-1, psym=symcat(16,thick=10), $
         symsize=2.5, color=im_color('red',255), /data
;      im_legend, string(cat[ii].number,format='(I0)'), /left, /top, $
;        margin=0, charsize=4, charthick=5, box=primary
       xyouts, 6, npix-10-1, string(cat[ii].number,format='(I0)'), $
         color=im_color('red',255), /data, charsize=4, charthick=5
;      if (ii mod 3) eq 0 then xyouts, 5, npix-10-1, $
;        string(cat[ii].number,format='(I0)'), color=im_color('yellow',255), /data, $
;        charsize=4, charthick=6
;      im_legend, 'ID='+string(cat[ii].number,format='(I4.4)'), /left, $
;        /bottom, box=0, margin=0, color=im_color('black',255)
       x0 = fix((convert_coord(pos[0:1],/normal,/to_device))[0])
       nx = fix((convert_coord(pos[2:3],/normal,/to_device))[0])-X0
       y0 = fix((convert_coord(pos[0:1],/normal,/to_device))[1])
       ny = fix((convert_coord(pos[2:3],/normal,/to_device))[1])-y0
       img = tvrd(x0,y0,nx,ny)
       tvlct, rr, gg, bb, /get
       write_png, outfile, img, rr, gg, bb
       set_plot, 'X'
    endfor

; now make a montage
    allpng = file_search('/tmp/star_'+string(cat.number,format='(I4.4)')+'.png')
    cmd = 'montage -bordercolor black -borderwidth 1'+ $
      ' -tile '+string(ncols,format='(I0)')+'x'+string(nrows,format='(I0)')+$
      ' -geometry +0+0 -quality 100 -resize 200x200 -title '+title+' '+$
      strjoin(allpng,' ')+' '+qafile
    spawn, cmd, /sh

; clean up    
    rmfile, file_search('/tmp/star_????.png') 
    
return
end

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

function get_cutouts, cat, image=image, invvar=invvar, hdr=hdr, npix=npix
; get postage stamp cutouts without recentering
    sz = size(image,/dim)
    nbigpix = npix*4
    pixel_scale = 0.065D
    wid = npix*pixel_scale/3600D ; cutout width [deg]
    nstar = n_elements(cat)
;   ivarfile = repstr(imfile,'drz','wht')
;   image = mrdfits(imfile,0,hdr,/silent)
;   invvar = mrdfits(ivarfile,0,ivarhdr,/silent)
    exptime = sxpar(hdr,'EXPTIME')
    extast, hdr, astr
    out = replicate({xcen: 0D, ycen: 0D, ra: 0D, dec: 0D, fracmask: 0.0, $
      nobjpix: 0, sigma: 0.0, image: fltarr(npix,npix), invvar: fltarr(npix,npix)},nstar)
    for ii = 0, nstar-1 do begin
       print, ii+1, nstar, string(13b), format='("Cutting out star ",I0,"/",I0,A10,$)'
       xx = cat[ii].x_image
       yy = cat[ii].y_image
; do a local sky subtraction
       x0 = long(cat[ii].x_image)
       y0 = long(cat[ii].y_image)
       xst = (x0-nbigpix/2L)>0L
       xnd = (x0+nbigpix/2L)<(sz[0]-1)
       yst = (y0-nbigpix/2L)>0L
       ynd = (y0+nbigpix/2L)<(sz[1]-1)
       hextract, image, hdr, bigcutout, bignewhdr, xst, xnd, yst, ynd, /silent
       hextract, invvar, hdr, bigcutoutivar, bignewivarhdr, xst, xnd, yst, ynd, /silent
; need to calculate an offset in case the thumbnail was cropped
; (should equal NBIGPIX by NBIGPIX but may not necessarily)
       sz2 = size(bigcutout,/dim)
       xoff = (nbigpix-sz2[0])/2L
       yoff = (nbigpix-sz2[1])/2L
; identify objects, compute the total number of pixels belonging to
; the star, and improve the sky-subtraction
       dobjects, bigcutout, obj=obj
       out[ii].nobjpix = total(obj eq obj[nbigpix/2-xoff,nbigpix/2-yoff]) 
       good = where(obj eq -1 and bigcutoutivar gt 0,ngood)
       if ngood eq 0 then stop
       mmm, bigcutout[good], skymode, skysig, /silent
; refine the center of the star; need to make sure the big
; cutout wasn't cropped
       drefine, bigcutout, nbigpix/2L-xoff, nbigpix/2L-yoff, xr=xcen, yr=ycen
       xcen += x0 - (nbigpix/2L-xoff)
       ycen += y0 - (nbigpix/2L-yoff)

; now cut out a smaller stamp of the star and center it
       cutout = fltarr(npix,npix)
       cutoutivar = fltarr(npix,npix)
       embed_stamp, cutout, image, npix/2L-xcen, npix/2L-ycen
       embed_stamp, cutoutivar, invvar, npix/2L-xcen, npix/2L-ycen
       cutoutivar = cutoutivar>0
;      hextract, image, hdr, cutout, newhdr, x0-npix/2, $
;        x0+npix/2, y0-npix/2, y0+npix/2, /silent
;      hextract, invvar, hdr, cutoutivar, newivarhdr, x0-npix/2, $
;        x0+npix/2, y0-npix/2, y0+npix/2, /silent
       cutoutsub = cutout-skymode
; add the shot noise of the object in quadrature to the inverse
; variance map; poisson error = sqrt(electrons), converted back to
; native electrons/pixel/sec units; only do this for pixels that are
; detected at >3*sigma
       these = where(cutout gt 3.0*skysig,nthese)
       if (nthese eq 0L) then begin
          splog, 'No significant pixels!'
          totalivar = cutoutivar
       endif else begin
          shotvar = cutout*0.0
          shotvar[these] = cutout[these]/exptime ; [electron/pixel/sec]
          totalvar = (shotvar + 1.0/(cutoutivar+(cutoutivar eq 0)))*(cutoutivar ne 0)
          totalivar = 1.0/(totalvar+(totalvar eq 0))*(totalvar ne 0)
       endelse
; pack it in!
       xy2ad, xcen, ycen, astr, ra, dec
;      xy2ad, xx, yy, astr, ra, dec
       out[ii].xcen = xcen ; xx
       out[ii].ycen = ycen ; yy
       out[ii].ra = ra
       out[ii].dec = dec
       out[ii].sigma = skysig
       out[ii].image = cutoutsub
       out[ii].invvar = totalivar ; cutoutivar
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
    
pro build_clash_psfs, sexcatalogs=sexcatalogs, candidates=candidates, $
  qaplots=qaplots, build_psf=build_psf
; jm13may30siena - build the PSFs for CLASH in all 16+1 bands using
; all the clusters

; unzip and rezip everything
; find /Volumes/Archive_the_First/clash-archive/ -name "*-fullfield_mosaic_065mas_*.fits.gz" -print | xargs gunzip -v
    
    clash = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
; no full-field mosaics yet
    these = where(strtrim(clash.shortname,2) eq 'macs1311' or strtrim(clash.shortname,2) eq 'clj1226')
    clash = clash[these]
;   keep = where(strtrim(clash.shortname,2) ne 'macs1311' and strtrim(clash.shortname,2) ne 'clj1226')
;   clash = clash[keep]

;   ww = where(strtrim(clash.shortname,2) eq 'macs0717' or strtrim(clash.shortname,2) eq 'macs0744' or $
;     strtrim(clash.shortname,2) eq 'rxj1347' or strtrim(clash.shortname,2) eq 'macs2129' or $
;     strtrim(clash.shortname,2) eq 'macs1149' or strtrim(clash.shortname,2) eq 'rxj2129')
;   ww = where(strtrim(clash.shortname,2) ne 'clj1226')
;   ww = lindgen(25)
;   ww = ww[where(ww gt 14)]
;   clash = clash[ww]
    ncl = n_elements(clash)
    struct_print, clash

    filt = clash_filterlist(short=short,weff=weff,instr=instr,zpt=zpt)
    srt = reverse(sort(weff))
    filt = filt[srt]
    short = short[srt]
    instr = instr[srt]
    zpt = zpt[srt]
    nfilt = n_elements(filt)

    reffilt = (where(short eq 'f160w'))[0] ; reference filter
    pixel_scale = 0.065                    ; [arcsec/pixel]
    npix = 51                              ; cutout size

; parameters for selecting candidate PSF stars    
    magpivot = 23.0
    magfaint = 25.0
; mu_max vs mag_auto selection
    slope = 0.95
    int1 = 19.6
    int2 = 20.3
; these coefficients work fine when using the mag_aper vs mag_auto method          
;   int = 26.5
;   slope = 0.95
    magaxis = range(15,magfaint,500)
    
    rootpath = getenv('IM_ARCHIVE_DIR')+'/projects/clash/psfs/' ; main I/O path
    archive = getenv('CLASH_ARCHIVE')
;   archive = getenv('IM_ARCHIVE_DIR')+'/'

; --------------------------------------------------
; build the SExtractor catalogs we need
       if keyword_set(sexcatalogs) then begin
;         for ic = 12, ncl-1 do begin
          for ic = 0, ncl-1 do begin
             cluster = strtrim(clash[ic].shortname,2)
             path = rootpath+cluster+'/'
             if file_test(path,/dir) eq 0 then file_mkdir, path

             if file_test(archive+strtrim(clash[ic].dirname,2),/dir) then begin
                mosaicpath = archive+strtrim(clash[ic].dirname,2)+'/HST/images/'+$
                  'mosaicdrizzle_image_pipeline/fullfield/'
;               splog, 'Unzipping '+mosaicpath
;               spawn, 'gunzip '+mosaicpath+'*.fits.gz', /sh
; for GAIN see notes at
; http://www.ifa.hawaii.edu/~rgal/science/sextractor_notes.html
                detect_imfile = file_search(mosaicpath+cluster+'-fullfield_mosaic_065mas_'+$
                  instr[reffilt]+'_'+short[reffilt]+'_drz_????????.fits')
;               detect_imfile = file_search(mosaicpath+cluster+'_mosaic_065mas_'+$
;                 instr[reffilt]+'_'+short[reffilt]+'_drz_????????.fits*')
                
                for ib = 0, nfilt-1 do begin
                   catalog_name = path+cluster+'-'+short[ib]+'.cat'
                   imfile = file_search(mosaicpath+cluster+'-fullfield_mosaic_065mas_'+$
                     instr[ib]+'_'+short[ib]+'_drz_????????.fits',count=nn)
;                  imfile = file_search(mosaicpath+cluster+'_mosaic_065mas_'+$
;                    instr[ib]+'_'+short[ib]+'_drz_????????.fits')

                   if nn ne 0 then begin
                      hdr = headfits(imfile)
                      gain = sxpar(hdr,'ccdgain')*sxpar(hdr,'exptime')        
                      clash_psf_sex, imfile, sexpath=rootpath, catalog_name=catalog_name, $
                        detect_imfile=detect_imfile, zpt=zpt[ib], gain=gain
                   endif else splog, 'No observations in band '+short[ib]
                endfor
;               splog, 'Zipping '+mosaicpath
;               spawn, 'gzip '+mosaicpath+'*.fits', /sh
             endif
          endfor 
       endif
       
; --------------------------------------------------
; select candidate PSF stars from the F160W filter in each cluster 
       if keyword_set(candidates) then begin
          splog, 'Selecting candidate PSF stars'
;         for ic = 3, 5 do begin
          for ic = 0, ncl-1 do begin
             cluster = strtrim(clash[ic].shortname,2)
             path = rootpath+cluster+'/'

             splog, 'Working on cluster '+cluster
             mosaicpath = archive+strtrim(clash[ic].dirname,2)+'/HST/images/'+$
               'mosaicdrizzle_image_pipeline/fullfield/'

             catfile = path+cluster+'-f160w.cat'
             refcat = mrdfits(catfile,2,/silent)

             istar = where(refcat.mag_auto gt 0.0 and refcat.mag_auto lt magfaint and $
               refcat.mu_max gt poly(refcat.mag_auto-magpivot,[int1,slope]) and $
               refcat.mu_max lt poly(refcat.mag_auto-magpivot,[int2,slope]),nstar)
;            istar = where(refcat.mag_auto gt 0.0 and refcat.mag_auto lt magfaint and $
;              refcat.mag_aper lt poly(refcat.mag_auto-magpivot,[int,slope]),nstar)
;            istar = where(refcat.mag_auto gt 0.0 and refcat.mag_auto lt magfaint and $
;              refcat.mag_aper lt poly(refcat.mag_auto-magpivot,[int,slope]) and $
;              refcat.flags eq 0,nstar)
;            medpsf = djs_median(refcat[istar].vignet,3)
             write_ds9_regionfile, refcat.alpha_j2000, refcat.delta_j2000, $
               file=path+cluster+'-'+short[reffilt]+'-allstars.reg', $
               color='green', symbol='box'
             write_ds9_regionfile, refcat[istar].alpha_j2000, refcat[istar].delta_j2000, $
               file=path+cluster+'-'+short[reffilt]+'-stars.reg', color='red'

; sort by brightness
             brightsort = sort(refcat[istar].mag_auto)
;            refid = lindgen(nstar)
             
; loop through each band  
             for ib = 0, nfilt-1 do begin
                catfile = path+cluster+'-'+short[ib]+'.cat'
                if file_test(catfile) then begin ; not all filters exist
; just read the F160W stars                   
                   cat = mrdfits(catfile,2,/silent,rows=istar[brightsort]) 
;                  cat = struct_addtags(replicate({id: 0L},nstar),cat)
;                  cat.id = refid
             
;                  istar = where(cat.mag_auto gt 0.0 and cat.mag_auto lt magfaint and $
;                    cat.mag_aper lt poly(cat.mag_auto-magpivot,[int,slope]) and $
;                    cat.flags eq 0,nstar)
;                  keep = where(cat[istar].flags eq 0,nstar)
                   keep = where(cat.mag_aper gt 0 and cat.mag_aper lt 99.0,nkeep)
                   cat = cat[keep]

; get cutouts and write out the final catalog of stellar candidates
                   imfile = file_search(mosaicpath+cluster+'-fullfield_mosaic_065mas_'+$
                     instr[ib]+'_'+short[ib]+'_drz_????????.fits*',count=nn)
                   if nn eq 0 then message, 'This should not happen!'
                   ivarfile = repstr(imfile,'drz','wht')
                   image = mrdfits(imfile,0,hdr,/silent)
                   invvar = mrdfits(ivarfile,0,ivarhdr,/silent)
                   exptime = sxpar(hdr,'EXPTIME')
                   extast, hdr, astr

                   cut = get_cutouts(cat,image=image,invvar=invvar,hdr=hdr,npix=npix)
                   out = struct_addtags(cat,temporary(cut))
                   im_mwrfits, out, repstr(catfile,'.cat','-stars.fits'), /clobber
                endif
             endfor 
          endfor
       endif

; --------------------------------------------------
; build QAplots
       if keyword_set(qaplots) then begin
          splog, 'Building QAplots'
;         for ic = 0, 1 do begin
          for ic = 0, ncl-1 do begin
             cluster = strtrim(clash[ic].shortname,2)
             path = rootpath+cluster+'/'
             splog, 'Building montage for cluster '+cluster

; (1) build a montage showing the individual stars in every band
             refcatfile = path+cluster+'-f160w-stars.fits.gz'
             refcat = mrdfits(refcatfile,1,/silent)
             nstar = n_elements(refcat)
             outcat_template = im_empty_structure(refcat,$
               ncopies=nstar,empty_value=-999)
             filters = 'f160w'

             outcat = refcat
             for ib = 1, nfilt-1 do begin
                catfile = path+cluster+'-'+short[ib]+'-stars.fits.gz'
                if file_test(catfile) then begin
                   filters = [filters,short[ib]]
                   cat1 = mrdfits(catfile,1,/silent)
                   match, refcat.number, cat1.number, m1, m2
                   outcat1 = outcat_template
                   outcat1[m1] = cat1[m2]
                   outcat = [[outcat],[outcat1]]
                endif
             endfor

; require that a star be observed in at least 3 bands; for the most
; part, this removes the stars observed in the parallel fields (mostly
; in F160W and F125W)
             these = where(total(outcat.number gt -900,2) gt 3 and $
               outcat[*,0].flags eq 0) ; require flagless F160W
             outcat = outcat[these,*]

             dim = size(outcat,/dim)
             nstar = dim[0]
             nfilters = dim[1]

             loadct, 0, /silent
             count = 0
             for ii = 0, nstar-1 do begin
                for jj = 0, nfilters-1 do begin
                   outfile = '/tmp/star_'+string(count,format='(I4.4)')+'.png'
                   if outcat[ii,jj].number eq -999 then begin
                      im = byte(outcat[ii,jj].image*0) ; no data
                   endif else begin
                      image = outcat[ii,jj].image
                      im = asinhscl(image,beta=0.005,omax=254,/negative)
                   endelse
                   set_plot, 'Z'
                   plotimage, im, /noaxes, /preserve_aspect, position=pos, /norm
                   if jj eq 0 then begin
                      xyouts, 6, 6, string(outcat[ii,jj].number,format='(I0)'), $
                        color=im_color('red',255), /data, charsize=4.0, charthick=5
                   endif
                   if ii eq 0 then begin
                      xyouts, 6, npix-10-1, strupcase(filters[jj]), /data, $
                        color=im_color('red',255), charsize=4.0, charthick=5
                   endif
                   write_thumbnail, outfile, pos=pos                   
                   set_plot, 'X'
                   count++
                endfor
             endfor

; make the montage then clean up             
             qafile = path+'qa_'+cluster+'_allstars_montage.pdf'
             allpng = file_search('/tmp/star_????.png')
             cmd = 'montage -bordercolor black -borderwidth 1'+ $
               ' -tile '+string(nfilters,format='(I0)')+'x'+string(nstar,format='(I0)')+$
               ' -geometry +0+0 -quality 100 -resize 200x200 '+$;-title '+strupcase(cluster)+' '+$
               strjoin(allpng,' ')+' '+qafile
             spawn, cmd, /sh
             rmfile, file_search('/tmp/star_????.png') 
             
;; (1) make a postage stamp montage and then pack all the montages
;; together 
;             allqafile = strarr(nfilt)
;             for ib = 0, nfilt-1 do begin
;                catfile = path+cluster+'-'+short[ib]+'-stars.fits.gz'
;                if file_test(catfile) then begin ; not all filters exist
;                   cat = mrdfits(catfile,1,/silent)
;                   qafile = path+'qa_'+cluster+'_allstars_montage_'+short[ib]+'.pdf'
;                   allqafile[ib] = qafile
;                   splog, 'Building '+qafile
;                   make_montage, cat, qafile=qafile, title=strupcase(cluster)+'/'+$
;                     strupcase(short[ib])
;                endif
;             endfor
;             allqafile = allqafile[where(strtrim(allqafile,2) ne '')]
;
;             qafile = path+'qa_'+cluster+'_allstars_montage.pdf'
;             spawn, 'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+$
;               qafile+' '+strjoin(allqafile,' '), /sh
;             rmfile, allqafile
          endfor 
;         model = (total(outcat.image,3))[*,*,0]
          
; (2) show how the stars were selected
;         for ic = 0, 1 do begin
          for ic = 0, ncl-1 do begin
             cluster = strtrim(clash[ic].shortname,2)
             path = rootpath+cluster+'/'
             splog, 'Building selection plots for cluster '+cluster

             xrange = [15,28]
             psfile = path+'qa_'+cluster+'_allstars.ps'
             im_plotconfig, 4, pos, psfile=psfile

;            for ib = 0, 0 do begin
             for ib = 0, nfilt-1 do begin
                ff = strupcase(short[ib])
                allcatfile = path+cluster+'-'+short[ib]+'.cat'
                catfile = path+cluster+'-'+short[ib]+'-stars.fits.gz'
                if file_test(catfile) then begin
                   cat = mrdfits(catfile,1,/silent)
                   allcat = mrdfits(allcatfile,2,/silent)
                   nstar = n_elements(cat)
; mu_max vs mag_auto
                   djs_plot, [0], [0], /nodata, position=pos[*,0], xrange=xrange, $
                     yrange=[11,26], xtickname=replicate(' ',10), ytitle='\mu_{max} (mag arcsec^{-2})', $
                     xsty=1, ysty=1
                   djs_oplot, allcat.mag_auto, allcat.mu_max, psym=symcat(6,thick=1), symsize=0.5
                   djs_oplot, cat.mag_auto, cat.mu_max, psym=symcat(16), $
                     symsize=0.5, color='red'
                   djs_oplot, magaxis, poly(magaxis-magpivot,[int1,slope]), line=0, color='blue', thick=4
                   djs_oplot, magaxis, poly(magaxis-magpivot,[int2,slope]), line=0, color='blue', thick=4
                   djs_oplot, magfaint*[1,1], [poly(magfaint-magpivot,[int1,slope]),$
                     poly(magfaint-magpivot,[int2,slope])], line=0, color='blue', thick=4
;                  djs_oplot, magaxis, poly(magaxis-magpivot,[int,slope]), line=0
;                  djs_oplot, magfaint*[1,1], [poly(magfaint-magpivot,[int,slope]),!y.crange[1]]
                   im_legend, [ff+' (N='+strtrim(nstar,2)+')'], $
                     /left, /top, box=0, margin=0
; mag_aper vs mag_auto
                   djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xrange=xrange, $
                     yrange=[32,14], xtickname=replicate(' ',10), ytitle='mag_{aper} (1 pixel)', $
                     xsty=1, ysty=1
                   djs_oplot, allcat.mag_auto, allcat.mag_aper, psym=symcat(6,thick=1), symsize=0.5
                   djs_oplot, cat.mag_auto, cat.mag_aper, psym=symcat(16), $
                     symsize=0.5, color='red'
;                  djs_oplot, magaxis, poly(magaxis-magpivot,[int,slope]), line=0
;                  djs_oplot, magfaint*[1,1], [poly(magfaint-magpivot,[int,slope]),!y.crange[1]]
; flux_radius vs mag_auto
                   djs_plot, [0], [0], /nodata, /noerase, position=pos[*,2], xrange=xrange, $
                     yrange=[0,15], xtitle='mag_{auto}', ytitle='r_{h} (pixels)', $
                     xsty=1, ysty=1
                   djs_oplot, allcat.mag_auto, allcat.flux_radius, psym=symcat(6,thick=1), symsize=0.5
                   djs_oplot, cat.mag_auto, cat.flux_radius, $
                     psym=symcat(16), symsize=0.5, color='red'
                endif
             endfor 
             im_plotconfig, psfile=psfile, /psclose, /pdf
          endfor 
       endif 

; --------------------------------------------------
; build the PSF; see also DFITPSF_ATLAS
       if keyword_set(build_psf) then begin
          rejsigma = [10.0,5.0]

          ngrid = 5
          pos = im_getposition(nx=ngrid,ny=ngrid)
          
          for ic = 0, ncl-1 do begin
             cluster = strtrim(clash[ic].shortname,2)
             path = rootpath+cluster+'/'

             splog, 'Working on cluster '+cluster
             mosaicpath = archive+strtrim(clash[ic].dirname,2)+'/HST/images/'+$
               'mosaicdrizzle_image_pipeline/fullfield/'

; loop on each band             
             for ib = 0, 0 do begin
;            for ib = 0, nfilt-1 do begin
                imfile = file_search(mosaicpath+cluster+'-fullfield_mosaic_065mas_'+$
                  instr[ib]+'_'+short[ib]+'_drz_????????.fits*')
                ivarfile = repstr(imfile,'drz','wht')
                image = mrdfits(imfile,0,hdr,/silent)
                invvar = mrdfits(ivarfile,0,ivarhdr,/silent)
                exptime = sxpar(hdr,'EXPTIME')
                extast, hdr, astr
                
; mask objects and then improve the sky subtraction; also compute how
; many pixels each star occupies, so we can identify objects that
; might be blended (i.e., too close) with other objects
;               splog, 'Sky-subtracting band '+short[ib]
;               dobjects, image, obj=obj
;               good = where(obj eq -1 and invvar gt 0,ngood)
;               mmm, image[good], skymode, skysig, /silent

; read the catalog of candidate stars
                catfile = path+cluster+'-'+short[ib]+'-stars.fits.gz'
                cat = mrdfits(catfile,1,/silent)
                nstar = n_elements(cat)

; normalize the peak to unity
;               norm = max(max(cat.image,dim=1),dim=1)
;               norm = total(total(cat.image,1),1)
;               norm1 = rebin(reform(norm,1,1,nstar),npix,npix,nstar)
;               cat.image = cat.image/norm1
;               cat.invvar = cat.invvar*norm1^2

                nkeep = nstar
                for iter = 0, n_elements(rejsigma)-1 do begin
                   ikeep = lindgen(nkeep)
;                  ikeep = where(cat.x_image gt 5000 and cat.x_image lt 10000 and $
;                    cat.y_image gt 5000 and cat.y_image lt 10000,nkeep)

; get cutouts and then build the mean PSF using a simple stack
                   cut = get_cutouts(cat[ikeep],image=image,$
                     invvar=invvar,hdr=hdr,npix=npix)
                   thiscat = struct_addtags(cat[ikeep],cut)

stop                   
                   
                   bpsf = djs_median(cut.image,dim=3)
;                  bpsf = total(cut.image,3)
                   bpsf = bpsf/total(bpsf)
;                  dfit_mult_gauss, bpsf, 1, amp, psfsig, model=model, /quiet

stop                   
                   
; do it!                   
                   ivaruse = invvar < 1.0/(abs(image)*0.1)^2
                   chi2 = dpsfselect_atlas(image,ivaruse,cat.xcen,cat.ycen,$
                     amp=amp,psf=bpsf,flux=flux,dof=dof,subimage=subimage,/noclip)
                   istar = where(chi2 lt dof+rejsigma[iter]*sqrt(2.*dof),comp=notstar,nstar)
                   

                
                   
                   
stop                   
                   
; clip the worst outliers
                   chi2 = fltarr(nkeep)
                   scale = fltarr(nkeep)
                   for ii = 0, nkeep-1 do begin  
                      scale[ii] = total(bpsf*cat[ii].image)/total(bpsf^2)
                      chi2[ii] = total((cat[ii].image-scale[ii]*bpsf)^2)
;                     chi2[ii] = total((cat[ii].image-scale[ii]*bpsf)^2*cat[ii].invvar)
;                     atv, cat[ii].image-scale[ii]*bpsf, /blo
                   endfor

; reject REJSIGMA outliers
                   dof = float(npix)^2
                   keep = where(chi2 lt dof+rejsigma[iter]*sqrt(2.0*dof),nkeep,comp=rej)
                   splog, 'REJSIGMA = ', rejsigma[iter], '; keeping ', nkeep, ' stars'
                   if (nkeep lt 3L) then message, 'Not enough good stars!'

; make a QAplot
                   for pp = 0, npage-1 do begin
;                     plotimage, asin
                   endfor
                   
                   
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
