pro bcgmstar_skysubtract
; jm13sep04siena - sky-subtract the CLASH clusters for the BCGs/Mstar
; project 

; note! images are converted to [10^-12 erg/s/cm^2/Hz] (pico-maggies)
    sample = read_bcgmstar_sample()
    ncl = n_elements(sample) 

    pixscale = 0.065D ; [arcsec/pixel]
    xsize = 5000L
    ysize = 5000L

    maskpath = bcgmstar_path(/objectmask)
    skypath = bcgmstar_path(/skysub)
    skyinfopath = bcgmstar_path()+'skyinfo/'
    qapath = bcgmstar_path()+'qaplots-skyinfo/'
    
; specifiy the filters and some other handy info    
    filt = bcgmstar_filterlist(short=short,instr=instr,$
      weff=weff,zpt=zpt)
    allfiltinfo = replicate({filt: '', short: '', instr: '', $
      weff: 0.0, zpt: 0.0},n_elements(filt))
    allfiltinfo.filt = filt
    allfiltinfo.short = short
    allfiltinfo.instr = instr
    allfiltinfo.weff = weff
    allfiltinfo.zpt = zpt

;   splog, 'HACK!!!'
;   allfiltinfo = allfiltinfo[[2,7,13]]
    nfilt = n_elements(allfiltinfo)
    struct_print, allfiltinfo

; initialize the sky structure    
    skyinfo1 = {file: '', band: '', weff: 0.0, exptime: 0.0, factor: 0.0, $
      sblimit: 0.0, mode: 0.0, sigma: 0.0, skew: 0.0, nsky: 0L}

; initialize the sky-subtraction QAplot
    ncol = 3 ; number of columns
    
; wrap on each cluster    
;   for ic = ncl-1, ncl-1 do begin
    for ic = 0, ncl-1 do begin
;   for ic = 1, 1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       splog, 'Working on cluster '+cluster
       outpath = skypath+cluster+'/'
       if file_test(outpath,/dir) eq 0 then file_mkdir, outpath

; not all clusters have all filters, so check for that case here
       bcgmodelpath = getenv('CLASH_ARCHIVE')+'/'+strtrim(sample[ic].dirname,2)+$
         '/HST/galaxy_subtracted_images/marc/'
       mosaicpath = getenv('CLASH_ARCHIVE')+'/'+strtrim(sample[ic].dirname,2)+$
         '/HST/images/mosaicdrizzle_image_pipeline/scale_65mas/'
       these = where(file_test(mosaicpath+cluster+'_mosaic_065mas_'+$
         strtrim(allfiltinfo.instr,2)+'_'+strtrim(allfiltinfo.short,2)+'_drz_????????.fits*'),nfilt)
       if nfilt gt n_elements(allfiltinfo) or nfilt eq 0 then message, 'Problem here!'
       filtinfo = allfiltinfo[these]
       reffilt = where(filtinfo.short eq 'f160w') ; reference filter

       skyinfo = replicate(skyinfo1,nfilt)
       skyinfo.band = filtinfo.short
       skyinfo.weff = filtinfo.weff

; read the original mosaics and inverse variance maps in a datacube;
; we need to do this so we can crop the mosaics 
       splog, 'Reading the data...'
       drzfiles = file_search(mosaicpath+cluster+'_mosaic_065mas_'+$
         filtinfo.instr+'_'+filtinfo.short+'_drz_????????.fits*')
       whtfiles = file_search(mosaicpath+cluster+'_mosaic_065mas_'+$
         filtinfo.instr+'_'+filtinfo.short+'_wht_????????.fits*')

;      bcgmodelfiles = file_search(bcgmodelpath+cluster+$
;        '_mosaic_065mas_*_'+filtinfo.short+'_drz_*_model.fits.gz')
;      bcgmodelfiles = file_search(bcgmodelpath+cluster+$
;        '_mosaic_065mas_*_'+filtinfo.short+'_drz_*_BCG.fits.gz')

       imagecube = fltarr(xsize,ysize,nfilt)
       ivarcube = fltarr(xsize,ysize,nfilt)
;      bcgcube = fltarr(xsize,ysize,nfilt)

       hdrcube = strarr(4000,nfilt)
       ivarhdrcube = strarr(4000,nfilt)
       exptime = fltarr(nfilt)
       
       for ib = 0, nfilt-1 do begin
          imagecube[*,*,ib] = mrdfits(drzfiles[ib],0,hdr,/silent)
          ivarcube[*,*,ib] = mrdfits(whtfiles[ib],0,ivarhdr,/silent)

          skyinfo[ib].exptime = sxpar(hdr,'EXPTIME')
          hcrop = where(strcompress(hdr,/remove) ne '',nhdr)
          ihcrop = where(strcompress(ivarhdr,/remove) ne '',nihdr)
          hdrcube[0:nhdr-1,ib] = hdr[hcrop]
          ivarhdrcube[0:nihdr-1,ib] = ivarhdr[ihcrop]
       endfor

; crop the mosaics to the reference band
       xcrop = minmax(where(total(imagecube[*,*,reffilt],2) ne 0))
       ycrop = minmax(where(total(imagecube[*,*,reffilt],1) ne 0))
;      xcrop = ceil(minmax(where(total(imagecube[*,*,reffilt],2) gt 0))/10.0)*10
;      ycrop = ceil(minmax(where(total(imagecube[*,*,reffilt],1) gt 0))/10.0)*10
       xcrop = (xcrop>0)<(xsize-1)
       ycrop = (ycrop>0)<(ysize-1)

       xcrop = xcrop+(odd(xcrop) eq 0) ; force odd because hextract, below, adds a pixel
       ycrop = ycrop+(odd(ycrop) eq 0) ; force odd because hextract, below, adds a pixel

; read the object mask (see BCGMSTAR_OBJECTMASK)       
       maskfile = maskpath+cluster+'/'+cluster+'-f160w-objectmask.fits'
       if file_test(maskfile) eq 0 then message, 'Missing '+maskfile
       objmask = mrdfits(maskfile,0,objmaskhdr) eq 0 ; 0 = masked, 1 = unmasked
       
; aggresively build an object mask using all the bands 
       splog, 'Building object mask'
       dobjects, imagecube[*,*,nfilt-1]*objmask, plim=20.0, fobject=fobj, dpsf=dpsf
;      dobjects, imagecube*rebin(reform(objmask,xsize,ysize,1),xsize,ysize,nfilt), $
;        plim=20.0, fobject=fobj, dpsf=dpsf
       
; finally sky-subtract each band independently; convert from
; [electron/s] to [1D-12 erg/s/cm^2/hz] (picomaggies)
       kl = k_lambda(filtinfo.weff,/odon)
       factor = 1D12*10D^(-0.4*(filtinfo.zpt-kl*sample[ic].ebv))
       skyinfo.factor = factor ; keep the conversion factor!

       nrow = ceil(nfilt/float(ncol))
       qafile = qapath+'qa_'+cluster+'_skyinfo.ps'
       im_plotconfig, 0, psfile=qafile, charsize=1.3
       pos = im_getposition(nx=ncol,ny=nrow,yspace=0.0,xspace=0.0,$
         xmargin=[0.9,0.4],width=2.4)
       xx = range(-0.05,0.05,1000)
       binsize = 0.001
       
       skyinfo.file = file_basename(drzfiles)

;      for ib = nfilt-1, nfilt-3, -1 do begin
       for ib = nfilt-1, 0, -1 do begin
;      for ib = 0, nfilt-1 do begin
          band = strupcase(filtinfo[ib].short)
          splog, 'Sky-subtracting band '+band
          good = where(objmask eq 1 and ivarcube[*,*,ib] gt 0 and fobj eq -1,ngood)
;         good = where(fobj eq -1 and ivarcube[xcrop[0]:xcrop[1],$
;           ycrop[0]:ycrop[1],ib] gt 0,ngood)
          mmm, (imagecube[*,*,ib])[good], mode1, sigma1, $
            skew1, nsky=nsky1, /debug;, /silent
          skyinfo[ib].mode = mode1*factor[ib]
          skyinfo[ib].sigma = sigma1*factor[ib]
          skyinfo[ib].skew = skew1*factor[ib]
          skyinfo[ib].nsky = nsky1

; make the QAplot          
;         med = djs_median((imagecube[*,*,ib])[good])
;         mode = 2.5*djs_median((imagecube[*,*,ib])[good])-1.5*djs_mean((imagecube[*,*,ib])[good])
;         sig = djsig((imagecube[*,*,ib])[good])

          if ib eq 1 then title = strupcase(cluster) else delvarx, title
          if ib ge nfilt-3 then begin
             delvarx, xtickname
          endif else begin
             xtickname = replicate(' ',10)
          endelse
          if (ib mod 3) eq 0 then begin
             delvarx, ytickname
          endif else begin
             ytickname = replicate(' ',10)
          endelse
          
          djs_plot, [0], [0], /nodata, noerase=ib lt nfilt-1, xsty=1, ysty=1, $
            xrange=[-0.12,0.04], yrange=[0,1.2], position=pos[*,ib], $
            xtickname=xtickname, ytickname=ytickname, title=title, ytickinterval=0.5
          im_plothist, (imagecube[*,*,ib])[good], bin=binsize, /peak, /overplot, thick=6
          djs_oplot, [0,0], !y.crange, line=0, color='grey'
          djs_oplot, mode1*[1,1], !y.crange, line=2
          djs_oplot, xx, gauss1(xx,[mode1,sigma1,1.0],/peak), color='red'
          im_legend, [band,'NSky='+strtrim(string(100*nsky1/5000.0^2,format='(F12.1)'),2)+'%',$
            'Mode='+strtrim(string(mode1,format='(F12.5)'),2),$
            'Sigma='+strtrim(string(sigma1,format='(F12.5)'),2)], $
            /left, /top, box=0, margin=0, charsize=0.8
       endfor 
       xyouts, min(pos[0,*])-0.06, (max(pos[3,*])-min(pos[1,*]))/2.0+min(pos[1,*]), $
         textoidl('Fraction of Pixels'), orientation=90, align=0.5, charsize=1.3, /norm
       xyouts, (max(pos[2,*])-min(pos[0,*]))/2.0+min(pos[0,*]), min(pos[1,*])-0.06, $
         textoidl('Pixel Values after Masking (counts/s)'), align=0.5, charsize=1.3, /norm
       im_plotconfig, psfile=qafile, /psclose, /pdf
       
; compute the 1-sigma surface brightness limit; to check, do:
; IDL> atv, im*(im gt 3*ss[12].sigma), /log
       skyinfo.sblimit = -2.5*alog10(skyinfo.sigma)+5*alog10(pixscale)-2.5*alog10(1D-12)
       struct_print, skyinfo

       im_mwrfits, skyinfo, skyinfopath+'skyinfo-'+cluster+'.fits', /clobber

; add the Poisson noise of the objects in quadrature to the
; 'intrinsic' inverse variance maps and convert to physical units
       for ib = 0, nfilt-1 do begin
          mask = (imagecube[*,*,ib]*(imagecube[*,*,ib] gt 3.0*skyinfo[ib].sigma/factor[ib])) and $
            (ivarcube[*,*,ib] gt 0)
; Poisson error = sqrt(electrons), converted back to native
; electrons/s then added in quadrature to INVVAR
          intvar = 1.0/(ivarcube[*,*,ib]+(ivarcube[*,*,ib] eq 0))*$ ; intrinsic var [electron/s]^2
            (ivarcube[*,*,ib] gt 0) 
          objvar = (imagecube[*,*,ib]*skyinfo[ib].exptime)*mask/skyinfo[ib].exptime^2 ; extrinsic var [electron/s]^2
          var = intvar + objvar                                       ; quadrature sum [electron/s]^2
          ivarcube[*,*,ib] = 1.0/(var+(var eq 0))*(var gt 0)          ; final map [electron/s]^(-2)

; convert to physical units
          imagecube[*,*,ib] *= factor[ib] ; [electron/s]->[erg/s/cm^2/Hz]
          ivarcube[*,*,ib] = ivarcube[*,*,ib]/factor[ib]^2
       endfor
          
; finally crop the images and write out!
       for ib = 0, nfilt-1 do begin
          outfile = outpath+cluster+'-'+filtinfo[ib].short+'.fits'
          splog, 'Writing '+outfile
          hcrop = where(strcompress(hdrcube[*,ib],/remove) ne '')
          ihcrop = where(strcompress(ivarhdrcube[*,ib],/remove) ne '')
          hextract, imagecube[*,*,ib]-skyinfo[ib].mode, hdrcube[hcrop,ib], $
            newim, newhdr, xcrop[0], xcrop[1], ycrop[0], ycrop[1], /silent
          hextract, ivarcube[*,*,ib], ivarhdrcube[ihcrop,ib], $
            newivar, newivarhdr, xcrop[0], xcrop[1], ycrop[0], ycrop[1], /silent
          hextract, objmask, objmaskhdr, mask, maskhdr, xcrop[0], xcrop[1], $
            ycrop[0], ycrop[1], /silent

          mwrfits, newim, outfile, newhdr, /create
          mwrfits, newivar, outfile, newivarhdr, /silent
          mwrfits, mask, outfile, maskhdr, /silent
          spawn, 'gzip -f '+outfile
       endfor                   ; close filter loop
    endfor                      ; close cluster loop

return
end
    
