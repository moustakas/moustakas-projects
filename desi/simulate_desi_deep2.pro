function get_emlines, ewoii, cflux=cflux, cwave=cwave, z=z, $
  newratio=newratio, linesigma=linesigma

    if n_elements(linesigma) eq 0 then linesigma = 100.0 ; [km/s]
    light = im_light(/kms)

    coeff0 = alog10(3600D)       ; starting rest wavelength [log-10 Angstrom]
    coeff1 = 1D-4                ; [dispersion in log-10 Angstroms]
    velsz = (10D^coeff1-1)*light ; [km/s]

    maxwave = 6700D
    nwavepix = long((alog10(maxwave)-coeff0)/coeff1+1)
    logwave = coeff0 + dindgen(nwavepix)*coeff1

; specify the lines we care about; for [OII] use a 1:0.71 flux ratio
; (f_3726 = 0.71 f_3728.8) following Weiner, and for the other lines
; adopt Mostek's values
    restwave = [3726.032D,3728.815D,4861.33D,4958.91D,5006.84D,6562.8D]
    lineratio = [0.415,0.585,0.4577,0.094,0.2821,2.002]
    if keyword_set(newratio) then lineratio = [0.415,0.585,0.9,0.04,0.12,1.75]
    linewave = alog10(restwave)
    nline = n_elements(linewave)

; specify the rest-frame emission-line EW (this will eventually come
; from a Monte Carlo simulation); then compute the total (rest-frame)
; [OII] flux
    oiiwave = djs_mean([3726.032,3728.815])
    oiiflux = ewoii*cflux ; [erg/s/cm2] 

; add the emission lines (see Schlegel's linebackfit for the
; conversions to and from log-10 A)
    zshift = 0.0
    sigma = linesigma/light/alog(10) ; line-width [log-10 Angstrom]
    oiiamplitude = oiiflux/alog(10)/oiiwave/(sqrt(2.0*!pi)*sigma)
    lineflux = logwave*0
    for jj = 0, nline-1 do begin
; emission-line amplitude          
       amplitude = lineratio[jj]*oiiamplitude
       linewave1 = linewave[jj]+alog10(1+zshift)
       lineflux += amplitude*exp(-0.5*(logwave-linewave1)^2/sigma^2)
    endfor

    lineflux = interpolate(lineflux/(1.0+z),findex((10^logwave)*(1+z),cwave),missing=0.0)
    
return, lineflux
end

pro simulate_desi_deep2
; jm13dec18siena - simulate DESI spectra using DEEP2

    prefix = 'desi_deep2'
    isedfit_dir = getenv('IM_PROJECTS_DIR')+'/desi/deep2/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

    ised = read_isedfit(isedfit_paramfile)
    ngal = n_elements(ised)

    ispec = mrdfits(isedfit_dir+'deep2_ispec.fits.gz',1)
    
    for ii = 0, ngal-1 do begin

; restore the best-fit iSEDfit model and resample to constant
; log-lambda 
       ised1 = read_isedfit(isedfit_paramfile,$
         /getmodels,index=ii,/flam)
       wave = ised1.wave
       cflux = ised1.flux

; get the continuum around [OII]; will eventually make this more
; sophisticated
       ewoii = ispec[ii].oii_3727_1_ew[0] + ispec[ii].oii_3727_2_ew[0]
       cflux_3727 = interpol(cflux,wave,3727.42)
       
       lineflux = get_emlines(ewoii,cflux=cflux_3727,$
         cwave=wave,z=ised1.z,linesigma=10.0)
       flux = cflux + lineflux

       djs_plot, wave/1D4, flux, xr=[0.7,0.75], xsty=3       
       djs_oplot, wave/1D4, cflux, color='orange'
       
       
stop       
    endfor
    
return
end
