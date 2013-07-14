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
       oiiamplitude = oiiflux/alog(10)/oiiwave/(sqrt(2.0*!pi)*sigma)
;      oiiamplitude = oiiflux/sqrt(2.0*!pi)/sigma
       for jj = 0, nline-1 do begin
; emission-line amplitude          
          amplitude = ynorm[ii]*lineratio[jj]*oiiamplitude
;         splog, amplitude, lineratio[jj], ynorm[ii], 10^linewave[jj]
          linewave1 = linewave[jj]+alog10(1+zshift)
          lineflux = amplitude*exp(-0.5*(logwave-linewave1)^2/sigma^2)
          spec2d[*,rotcurve[ii].row-1] += lineflux
       endfor
;      print
    endfor
stop          

; finally redshift it all and make a header
    spec2d = spec2d/(1+cat.z)
    mkhdr, hdr, float(spec2d)
    sxdelpar, hdr, 'DATE'
    sxdelpar, hdr, 'COMMENT'
    sxaddpar, hdr, 'COEFF0', coeff0+alog10(1+cat.z)
    sxaddpar, hdr, 'COEFF1', coeff1
    
return, spec2d    
end

pro build_desi_tkrs_spec2d
; jm13jun15siena - build 2D spectra from Ben Weiner's rotation
; curve fits to the TKRS galaxies

    rootpath = getenv('IM_PROJECTS_DIR')+'/desi/'
    dataroot = 'tkrs'

; read the output of BUILD_DESI_TKRS_SAMPLE    
    info = mrdfits(rootpath+'tkrs_info.fits.gz',1)
    zcat = mrdfits(rootpath+'tkrs_zcat.fits.gz',1)
    ngal = n_elements(info)

; build the 2D spectrum for each galaxy
    vname = 'default.nolines'
    k_load_vmatrix, vmatrix, lambda, vname=vname
    cwave = k_lambda_to_centers(lambda)
    for ii = 0, ngal-1 do begin
; rebuild the continuum model
;      cwave = wave*(1+cat[ii].z)
       cflux = vmatrix#zcat[ii].coeffs;/(1+cat[ii].z)
       curve = read_tkrs_rotcurve(rootpath+info[ii].rotcurvefile)
;      plot, curve.row, curve.model_unblur_flux, psym=-8
       yfit = mpfitpeak(curve.row,curve.model_unblur_flux,aa)
       print, aa

;      im = dfakegal(sersicn=1.0,r50=10.0,ba=0.5,nx=100,ny=100,phi=90.0)

;      spec2d = build_spec2d(curve,cat=cat[ii],info=info[ii],$
;        cwave=cwave,cflux=cflux,header=header)
;      outfile = rootpath+'spec2d/model_'+prefix[ii]+'.fits'
;      im_mwrfits, spec2d, outfile, header, /clobber
    endfor

    
    
stop    

return
end
