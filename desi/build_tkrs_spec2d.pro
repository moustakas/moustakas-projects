function build_spec2d, rotcurve, cat=cat, info=info, cwave=cwave, $
  cflux=cflux, header=hdr
; wavelength parameters    
    light = im_light(/kms)

; output wavelength vector    
    coeff0 = alog10(3600D)       ; starting rest wavelength [log-10 Angstrom]
    coeff1 = 1D-4                ; [dispersion in log-10 Angstroms]
    velsz = (10D^coeff1-1)*light ; [km/s]

    maxwave = 6700D
    nwavepix = long((alog10(maxwave)-coeff0)/coeff1+1)
    logwave = coeff0 + dindgen(nwavepix)*coeff1

; interpolate the continuum spectrum
    continuum = interpolate(cflux,findex(alog10(cwave),logwave)) ; [erg/s/cm^2/A]
    
; specify the lines we care about; for [OII] use a 1:0.71 flux ratio
; (f_3726 = 0.71 f_3728.8) following Weiner, and for the other lines
; adopt Mostek's values
    restwave = [3726.032D,3728.815D,4861.33D,4958.91D,5006.84D,6562.8D]
    lineratio = [0.415,0.585,0.4577,0.094,0.2821,2.002]
    linewave = alog10(restwave)
    nline = n_elements(linewave)

; specify the rest-frame emission-line EW (this will eventually come
; from a Monte Carlo simulation); then compute the total (rest-frame)
; [OII] flux
    ewoii = 20.0 ; [Angstrom]
    oiiwave = djs_mean([3726.032,3728.815])
    oiiflux = ewoii*interpolate(smooth(cflux,10),findex(cwave,oiiwave)) ; [erg/s/cm2] 

; build the spectrum    
    ny = 51                     ; [pixels]
    spec2d = fltarr(nwavepix,ny)

    ynorm = fltarr(ny) ; spatial normalization factor
    ynorm[rotcurve.row-1] = rotcurve.model_unblur_flux/max(rotcurve.model_unblur_flux)
    ynorm = ynorm/total(ynorm) ; normalize to unity
    
    for ii = 0, n_elements(rotcurve)-1 do begin
; add the continuum
       spec2d[*,rotcurve[ii].row-1] = ynorm[ii]*continuum
; add the emission lines (see Schlegel's linebackfit for the
; conversions to and from log-10 A)
       zshift = rotcurve[ii].model_unblur_vel/light           ; line-position due to the rotation curve
       sigma = rotcurve[ii].model_unblur_sigma/light/alog(10) ; line-width [log-10 Angstrom]
       oiiamplitude = oiiflux/alog(10)/oiiwave
       for jj = 0, nline-1 do begin
; emission-line amplitude          
          amplitude = ynorm[ii]*lineratio[jj]*oiiamplitude
;         splog, amplitude, lineratio[jj], ynorm[ii], 10^linewave[jj]
          linewave1 = linewave[jj]+alog10(1+zshift)
          lineflux = ynorm[ii]*0.5*lineratio[jj]*oiiflux*exp(-0.5*(logwave-linewave1)^2/sigma^2)
;         lineflux = amplitude*exp(-0.5*(logwave-linewave1)^2/sigma^2)/sqrt(2.0*!pi*sigma)*light
;         if jj eq 2 then begin
;            continuum = oiiflux*0.01
;            fluxerr = randomn(seed,nwavepix)*(continuum+lineflux)
;            rr = linebackfit(10^linewave[jj],logwave,continuum+lineflux,$
;              invvar=1/fluxerr^2,background=lineflux*0+oiiflux*0.01,yfit=yfit)
;            help, rr, /str
;            djs_plot, 10^logwave, continuum+lineflux, xr=[4840,4900], psym=10, ysty=3
;            djs_oplot, 10^logwave, yfit, color='orange', psym=10
;            stop
;         endif
          spec2d[*,rotcurve[ii].row-1] += lineflux
       endfor
;      print
    endfor

; finally redshift it all and make a header
    spec2d = spec2d/(1+cat.z)
    mkhdr, hdr, float(spec2d)
    sxdelpar, hdr, 'DATE'
    sxdelpar, hdr, 'COMMENT'
    sxaddpar, hdr, 'COEFF0', coeff0+alog10(1+cat.z)
    sxaddpar, hdr, 'COEFF1', coeff1
    
return, spec2d    
end

function read_rotcurve, rotfile
    readcol, rotfile, row, yoff, totflux, flux, vel, $
      errvel, sigma, errsigma, model_unblur_flux, $
      model_unblur_vel, model_unblur_sigma, $
      model_blur_flux, model_blur_vel, model_blur_sigma, $
      good, /silent, format='I,F,F,F,F,F,F,F,F,F,F,F,F,F,I'
    nrow = n_elements(row)
    curve = replicate({row: 0, $
      yoff: 0.0, $
      totflux: 0.0,$
      flux: 0.0, $
      vel: 0.0, $
      errvel: 0.0, $
      sigma: 0.0, $
      errsigma: 0.0, $
      model_unblur_flux: 0.0, $
      model_unblur_vel: 0.0, $
      model_unblur_sigma: 0.0, $
      model_blur_flux: 0.0, $
      model_blur_vel: 0.0, $
      model_blur_sigma: 0.0, $
      good: 0},nrow)

    curve.row = row
    curve.yoff = yoff
    curve.totflux = totflux
    curve.flux = flux
    curve.vel = vel
    curve.errvel = errvel
    curve.sigma = sigma
    curve.errsigma = errsigma
    curve.model_unblur_flux = model_unblur_flux
    curve.model_unblur_vel = model_unblur_vel
    curve.model_unblur_sigma = model_unblur_sigma
    curve.model_blur_flux = model_blur_flux
    curve.model_blur_vel = model_blur_vel
    curve.model_blur_sigma = model_blur_sigma
    curve.good = good
    
return, curve
end

pro build_tkrs_spec2d
; jm13jun15siena - build 2D spectra from Ben Weiner's rotation
; curve fits to the TKRS galaxies

    outpath = getenv('IM_PROJECTS_DIR')+'/desi/'
    datapath = outpath+'tkrs/'

; read the quality file
    qualfile = datapath+'list.all.id.rcqual.yrange'
    readcol, qualfile, id, quality, prob, minrow, maxrow, $
      sizepix, sizearcsec, format='L,I,I,X,F,F,F,F', /silent
    info = replicate({id: 0L, quality: 0, prob: 0, minrow: 0, $
      maxrow: 0, sizepix: 0, sizearcsec: 0.0},n_elements(id))
    info.id = id
    info.quality = quality
    info.prob = prob
    info.minrow = minrow
    info.maxrow = maxrow
    info.sizepix = sizepix
    info.sizearcsec = sizearcsec

; keep just Q>=3 objects    
    good = where(info.quality ge 3,ngood)
    
    cat = read_tkrs(/zcat)
    kcorr = read_tkrs(/kcorr)
    match, cat.id, info[good].id, m1, m2
    cat = cat[m1]
    kcorr = kcorr[m1]
    info = info[good[m2]]

; select galaxies at z~1    
    these = where(cat.z gt 0.95 and cat.z lt 1.05,ngal)
    cat = cat[these]
    kcorr = kcorr[these]
    info = info[these]
    im_mwrfits, cat, outpath+'sim_zcat.fits', /clobber
    
    prefix = string(cat.id,format='(I7.7)')
    rotfiles = file_search(datapath+prefix+'_*.rotcurv')
    fitfiles = file_search(datapath+prefix+'_*.fitdata')

; read all the fitdata files
    for ii = 0, ngal-1 do begin
       readcol, fitfiles[ii], id, flux, sigma, sigma_cor, $
         vsyst, radius, radius_err, vrot, vrot_err, $
         sigma2d, sigma2d_err, ngood, chi2, /silent
       if ii eq 0 then begin
          fitdata = replicate({$
            id: 0L,$
            flux: 0.0, $
            sigma: 0.0, $
            sigma_cor: 0.0, $
            vsyst: 0.0, $
            radius: 0.0, $
            radius_err: 0.0, $
            vrot: 0.0, $
            vrot_err: 0.0, $
            sigma2d: 0.0, $
            sigma2d_err: 0.0, $
            ngood: 0.0, $
            chi2: 0.0},ngal)
       endif
       fitdata[ii].id = id 
       fitdata[ii].flux = flux
       fitdata[ii].sigma = sigma
       fitdata[ii].sigma_cor = sigma_cor
       fitdata[ii].vsyst = vsyst
       fitdata[ii].radius = radius
       fitdata[ii].radius_err = radius_err
       fitdata[ii].vrot = vrot
       fitdata[ii].vrot_err = vrot_err
       fitdata[ii].sigma2d = sigma2d
       fitdata[ii].sigma2d_err = sigma2d_err
       fitdata[ii].ngood = ngood
       fitdata[ii].chi2 = chi2
    endfor

; build the 2D spectrum for each galaxy
    vname = 'default.nolines'
    k_load_vmatrix, vmatrix, lambda, vname=vname
    cwave = k_lambda_to_centers(lambda)
    for ii = 0, ngal-1 do begin
; rebuild the continuum model
;      cwave = wave*(1+cat[ii].z)
       cflux = vmatrix#kcorr[ii].coeffs;/(1+cat[ii].z)
       curve = read_rotcurve(rotfiles[ii])
       spec2d = build_spec2d(curve,cat=cat[ii],info=info[ii],$
         cwave=cwave,cflux=cflux,header=header)
       outfile = outpath+'sim_'+prefix[ii]+'.fits'
       im_mwrfits, spec2d, outfile, header, /clobber
stop
    endfor

    
    
stop    

return
end
