pro bcgmstar_parallelsky
; jm14aug28siena - analyze the statistics of the sky in the parallel
; fields 

; note! images are converted to [10^-12 erg/s/cm^2/Hz] (pico-maggies)
    sample = read_bcgmstar_sample()
    ncl = n_elements(sample) 

    pixscale = 0.065D ; [arcsec/pixel]
    
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

    nfilt = n_elements(allfiltinfo)
    struct_print, allfiltinfo

    kl = k_lambda(allfiltinfo.weff,/odon)

; initialize the template sky structure    
    skyinfo_template = {file: '', band: '', weff: 0.0, exptime: 0.0, factor: 0.0, $
      sblimit: 0.0, mode: 0.0, sigma: 0.0, skew: 0.0, nsky: 0L}

; wrap on each cluster    
    for ic = 14, 14 do begin
;   for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       splog, 'Working on cluster '+cluster

       factor = 1D12*10D^(-0.4*(allfiltinfo.zpt-kl*sample[ic].ebv))
       clusterpath = getenv('CLASH_ARCHIVE')+'/'+strtrim(sample[ic].dirname,2)+$
         '/HST/images/mosaicdrizzle_image_pipeline/'

       for ib = 0, nfilt-1 do begin
          drzfiles = file_search(clusterpath+['wfc3par?/'+cluster+'-wfc3par?_mosaic_065mas_wfc3ir_'+$
            strtrim(allfiltinfo[ib].short,2),'acspar?/'+cluster+'-acspar?_mosaic_065mas_acs_'+$
            strtrim(allfiltinfo[ib].short,2)]+'_drz.fits.gz',count=npar)
          whtfiles = repstr(drzfiles,'drz','wht')
          for ip = 0, npar-1 do begin
             skyinfo1 = skyinfo_template
             skyinfo1.file = file_basename(drzfiles[ip])
             skyinfo1.band = allfiltinfo[ib].short
             skyinfo1.weff = allfiltinfo[ib].weff
             skyinfo1.factor = factor[ib]

             image = mrdfits(drzfiles[ip],0,hdr,/silent)
             ivar = mrdfits(whtfiles[ip],0,ivarhdr,/silent)
             skyinfo1.exptime = sxpar(hdr,'EXPTIME')

             splog, 'Computing sky statistics for '+skyinfo1.band
             dobjects, image, plim=20.0, fobject=fobj, dpsf=dpsf
             good = where(ivar gt 0 and fobj eq -1,ngood)
             mmm, image[good], mode1, sigma1, skew1, nsky=nsky1, /debug;, /silent
             skyinfo1.mode = mode1*factor[ib]
             skyinfo1.sigma = sigma1*factor[ib]
             skyinfo1.skew = skew1*factor[ib]
             skyinfo1.nsky = nsky1

             if n_elements(skyinfo) eq 0 then skyinfo = skyinfo1 else skyinfo = [skyinfo,skyinfo1]
          endfor
       endfor 
       
       skyinfo.sblimit = -2.5*alog10(skyinfo.sigma)+5*alog10(pixscale)-2.5*alog10(1D-12)
       struct_print, skyinfo

       im_mwrfits, skyinfo, skyinfopath+'parallel-skyinfo-'+cluster+'.fits', /clobber
       delvarx, skyinfo
    endfor                       ; close cluster loop

return
end
    
