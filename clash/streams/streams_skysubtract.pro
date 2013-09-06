pro streams_skysubtract
; jm13sep04siena - sky-subtract the CLASH clusters containing tidal
; streams 

; note! images are converted to [10^-12 erg/s/cm^2/Hz] (pico-maggies)
    
    sample = rsex(streams_path(/propath)+'streams_sample.sex')
    struct_print, sample
    ncl = n_elements(sample) 

    pixscale = 0.065D ; [arcsec/pixel]
    xsize = 5000L
    ysize = 5000L
    
; specifiy the filters and some other handy info    
    filt = clash_filterlist(short=short,instr=instr,$
      weff=weff,zpt=zpt,/dropbluest)
    allfiltinfo = replicate({filt: '', short: '', instr: '', $
      weff: 0.0, zpt: 0.0},n_elements(filt))
    allfiltinfo.filt = filt
    allfiltinfo.short = short
    allfiltinfo.instr = instr
    allfiltinfo.weff = weff
    allfiltinfo.zpt = zpt
    struct_print, allfiltinfo

; initialize the sky structure    
    skyinfo1 = {file: '', band: '', sblimit: 0.0, $
      mode: 0.0, sigma: 0.0, skew: 0.0, nsky: 0L}

; wrap on each cluster    
    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       splog, 'Working on cluster '+cluster
       outpath = streams_path(/skysub)+cluster+'/'
       if file_test(outpath,/dir) eq 0 then file_mkdir, outpath

; not all clusters have all filters, so check for that case here
       mosaicpath = getenv('CLASH_ARCHIVE')+'/'+strtrim(sample[ic].dirname,2)+$
         '/HST/images/mosaicdrizzle_image_pipeline/scale_65mas/'
       these = where(file_test(mosaicpath+cluster+'_mosaic_065mas_'+$
         instr+'_'+short+'_drz_????????.fits*'),nfilt)
       if nfilt gt n_elements(allfiltinfo) or nfilt eq 0 then message, 'Problem here!'
       filtinfo = allfiltinfo[these]
       reffilt = where(filtinfo.short eq 'f160w') ; reference filter

       skyinfo = replicate(skyinfo1,nfilt)
       skyinfo.band = filtinfo.short
       
; read the original mosaics and inverse variance maps
       splog, 'Reading the data...'
       drzfiles = file_search(mosaicpath+cluster+'_mosaic_065mas_'+$
         filtinfo.instr+'_'+filtinfo.short+'_drz_????????.fits*')
       whtfiles = file_search(mosaicpath+cluster+'_mosaic_065mas_'+$
         filtinfo.instr+'_'+filtinfo.short+'_wht_????????.fits*')

       imagecube = fltarr(xsize,ysize,nfilt)
       ivarcube = fltarr(xsize,ysize,nfilt)
       hdrcube = strarr(4000,nfilt)
       ivarhdrcube = strarr(4000,nfilt)
       exptime = fltarr(nfilt)
       
       for ib = 0, nfilt-1 do begin
          imagecube[*,*,ib] = mrdfits(drzfiles[ib],0,hdr,/silent)
          ivarcube[*,*,ib] = mrdfits(whtfiles[ib],0,ivarhdr,/silent)
          exptime[ib] = sxpar(hdr,'EXPTIME')
          hcrop = where(strcompress(hdr,/remove) ne '',nhdr)
          ihcrop = where(strcompress(ivarhdr,/remove) ne '',nihdr)
          hdrcube[0:nhdr-1,ib] = hdr[hcrop]
          ivarhdrcube[0:nihdr-1,ib] = ivarhdr[ihcrop]
       endfor
       
; crop the mosaics to the reference band, rounding to the nearest
; hundred pixels
       xcrop = ceil(minmax(where(total(imagecube[*,*,reffilt],1) gt 0))/100.0)*100
       ycrop = ceil(minmax(where(total(imagecube[*,*,reffilt],2) gt 0))/100.0)*100
       xcrop = (xcrop>0)<(xsize-1)
       ycrop = (ycrop>0)<(ysize-1)
       
; aggresively build an object mask using all the bands
       splog, 'Building object mask'
       dobjects, imagecube[xcrop[0]:xcrop[1],ycrop[0]:ycrop[1],*], $
         object=oimage, plim=5.0, fobject=fobj, dpsf=dpsf
       sz = size(fobj,/dim)

; finally sky-subtract each band independently; convert from
; [electron/s] to [1D-12 erg/s/cm^2/hz] (picomaggies)
       kl = k_lambda(filtinfo.weff,/odon)
       factor = 1D12*10D^(-0.4*(filtinfo.zpt-kl*sample[ic].ebv)) 
       
       skyinfo.file = file_basename(drzfiles)
       for ib = 0, nfilt-1 do begin
          splog, 'Sky-subtracting band '+filtinfo[ib].short
          good = where(fobj eq -1 and ivarcube[xcrop[0]:xcrop[1],$
            ycrop[0]:ycrop[1],ib] gt 0,ngood)
          mmm, (imagecube[xcrop[0]:xcrop[1],ycrop[0]:ycrop[1],ib])[good], $
            mode1, sigma1, skew1, nsky=nsky1, /silent ;, /debug
          skyinfo[ib].mode = mode1*factor[ib]
          skyinfo[ib].sigma = sigma1*factor[ib]
          skyinfo[ib].skew = skew1*factor[ib]
          skyinfo[ib].nsky = nsky1
       endfor

; compute the 1-sigma surface brightness limit; to check, do:
; IDL> atv, im*(im gt 3*ss[12].sigma), /log
       skyinfo.sblimit = -2.5*alog10(skyinfo.sigma)+5*alog10(pixscale)-2.5*alog10(1D-12)
       struct_print, skyinfo

       im_mwrfits, skyinfo, streams_path(/skysub)+'skyinfo-'+cluster+'.fits', /clobber

; add the Poisson noise of the objects in quadrature to the
; 'intrinsic' inverse variance maps and convert to physical units
       for ib = 0, nfilt-1 do begin
          mask = (imagecube[*,*,ib]*(imagecube[*,*,ib] gt 3.0*skyinfo[ib].sigma/factor[ib])) and $
            (ivarcube[*,*,ib] gt 0)
; Poisson error = sqrt(electrons), converted back to native
; electrons/s then added in quadrature to INVVAR
          intvar = 1.0/(ivarcube[*,*,ib]+(ivarcube[*,*,ib] eq 0))*$ ; intrinsic var [electron/s]^2
            (ivarcube[*,*,ib] gt 0) 
          objvar = (imagecube[*,*,ib]*exptime[ib])*mask/exptime[ib]^2 ; extrinsic var [electron/s]^2
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

          mwrfits, newim, outfile, newhdr, /create
          mwrfits, newivar, outfile, newivarhdr, /silent
          spawn, 'gzip -f '+outfile
       endfor                   ; close filter loop
    endfor                      ; close cluster loop

return
end
    
