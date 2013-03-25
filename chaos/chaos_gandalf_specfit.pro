pro chaos_rebin, spec1d, velscale=velscale, fluxscale=fluxscale, $
  rebin_wave=wave, rebin_flux=flux, rebin_ferr=ferr, snr=snr
; internal routine to logarithmically rebin the observed spectrum
    
; read the spectrum
    linearflux = spec1d.flux/fluxscale ; note!
    linearferr = spec1d.ferr/fluxscale
    linearwave = spec1d.wave

; rebin logarithmically in wavelength
    log_rebin, minmax(linearwave), linearflux, flux, $
      wave, velscale=velscale
    log_rebin, minmax(linearwave), linearferr^2, var, $
      velscale=velscale
    ferr = sqrt(abs(var))

    good = where(linearferr gt 0.0,ngood)
    snr = djs_median(linearflux[good]/linearferr[good])
    
return
end

pro chaos_zabs_vdisp, wave, flux, ferr, tempwave=tempwave, $
  tempflux=tempflux, velscale=velscale, zabs_guess=zabs_guess, $
  vdisp_guess=vdisp_guess, alllinefile=alllinefile, final_zabs=zabs, $
  final_vdisp=vdisp, err_final_zabs=zabs_err, err_final_vdisp=vdisp_err, $
  debug=debug
; internal support routine to fit the blue part of the spectrum to get
; ZABS and VDISP

    light = 2.99792458D5 ; speed of light [km/s]

; fitting parameters
    if (n_elements(zabs_guess) eq 0) then zabs_guess = 0.0
    if (n_elements(vdisp_guess) eq 0) then vdisp_guess = 150.0 ; [km/s]
    if (n_elements(vmaxshift) eq 0) then vmaxshift = 500D ; maximum velocity shift [km/s]
    if (n_elements(sigmamax) eq 0) then sigmamax = 500D   ; maximum velocity dispersion [km/s]
    npix = n_elements(wave)
    
    min_fitwave1 = 3700.0 ; rest-frame
    max_fitwave1 = 4850.0 ; 4200.0

    linearwave = exp(wave)
    restwave = wave - alog(zabs_guess+1.0D)
    restflux = flux*(zabs_guess+1.0D)
    restferr = ferr*(zabs_guess+1.0D)

; initial guesses: compute the *rest-frame* velocity offset between
; the templates and the data and restrict the wavelength range of the
; templates
    keep = where((exp(tempwave) gt exp(min(restwave))-10.0) and $
      (exp(tempwave) lt exp(max(restwave))+10.0),nkeep)
    if (nkeep eq 0) then message, 'Problem here!'
    fit_tempwave = tempwave[keep]
    fit_tempflux = tempflux[keep,*]

    vsys = alog(zabs_guess+1.0D)*light ; systemic velocity
    voffset = -(min(restwave)-min(fit_tempwave))*light
    start = [vsys*0.0+voffset,vdisp_guess]

; fit the continuum over a restricted wavelength range (around CaII)
; and mask emission lines
    linepars = read_gandalf_elinelist(alllinefile)
    goodpixels = gandalf_mask_emission_lines(npix,vsys,$
      linepars,velscale,wave[0],velscale/light,$
      l_rf_range=[min_fitwave1,max_fitwave1],sigma=300.0)

    ages_ppxf, fit_tempflux, restflux, restferr, velscale, $
      start, sol, goodpixels=goodpixels, plot=doplot, $
      moments=2, degree=degree, error=err, bestfit=bestfit, $
      /clean, /quiet, vmaxshift=vmaxshift, sigmamax=sigmamax
    chi2 = sol[6]
    err = err*sqrt(chi2)        ; scale by chi^2/dof

    zabs = exp((vsys+sol[0]-voffset)/light)-1.0D ; new absorption-line redshift
    zabs_err = zabs*(err[0]/light)
    vdisp = sol[1]
    vdisp_err = err[1]
    
    if keyword_set(debug) then begin
       djs_plot, exp(restwave), restflux, ps=10, xsty=3, ysty=3, xr=[3650,4400]
       djs_oplot, exp(restwave), bestfit, ps=10, color='yellow'
       splog, zabs_guess, zabs, zabs_err, vdisp, vdisp_err, chi2
       splog, 'Waiting for keystroke...' & cc = get_kbrd(1)
    endif
    
end
    
function chaos_fit_continuum, wave, flux, ferr, tempwave=tempwave, $
  tempflux=tempflux, velscale=velscale, zabs=zabs, vdisp=vdisp, $
  ebv_guess=ebv_guess, alllinefile=alllinefile, goodpixels=goodpixels, $
  final_ebv=ebv, chi2=chi2, continuum_coeff=continuum_coeff, debug=debug
; internal support routine to fit the continuum at a *fixed* ZABS and
; VDISP 

    light = 2.99792458D5 ; speed of light [km/s]
    if (n_elements(ebv_guess) eq 0) then ebv_guess = 0.1
    if (n_elements(vmaxshift) eq 0) then vmaxshift = 300D ; maximum velocity shift [km/s]
    if (n_elements(sigmamax) eq 0) then sigmamax = 500D   ; maximum velocity dispersion [km/s]
;   if (n_elements(vdisp) eq 0) then vdisp = 150.0        ; velocity dispersion guess [km/s]
    npix = n_elements(wave)
    
    restwave = wave - alog(zabs+1.0D)
    restflux = flux*(zabs+1.0D)
    restferr = ferr*(zabs+1.0D)
    
; initial guesses: compute the velocity offset between the templates
; and the data and restrict the wavelength range of the templates
    vsys = alog(zabs+1.0D)*light ; systemic velocity
    voffset = -(min(restwave)-min(tempwave))*light
    start = [vsys*0.0+voffset,vdisp]

    linepars = read_gandalf_elinelist(alllinefile)
    goodpixels = gandalf_mask_emission_lines(npix,vsys,$
      linepars,velscale,wave[0],velscale/light,$
      sigma=500.0)
    goodpixels = cmset_op(goodpixels,'and',$
      where(restferr lt 1E5))

    ebv = ebv_guess
    ages_ppxf, tempflux, restflux, restferr, velscale, start, $
      sol, goodpixels=goodpixels, plot=doplot, moments=0, degree=-1, $
      weights=continuum_coeff, bestfit=continuum, /clean, /quiet, $
      vmaxshift=vmaxshift, sigmamax=sigmamax, error=err, $
      lambda=exp(restwave), reddening=ebv
    chi2 = sol[6]

    if keyword_set(debug) then begin
       djs_plot, exp(restwave), restflux, ps=10, xsty=3, ysty=3, xr=[3650,4400]
       djs_oplot, exp(restwave), continuum, ps=10, color='blue'
 ;     djs_oplot, exp(restwave[goodpixels]), restflux[goodpixels], psym=4, color='red'
;      djs_oplot, exp(restwave[goodlinepixels]), restflux[goodlinepixels], ps=4, color='blue'
       splog, 'Waiting for keystroke...' & cc = get_kbrd(1)
    endif
    
return, continuum
end
    
function chaos_fit_lines, wave, flux, ferr, continuum, smooth_continuum, $
  zabs=zabs, vdisp=vdisp, velscale=velscale, linefile=linefile, $
  line_inst_vdisp=line_inst_vdisp, tempwave=tempwave, broad=broad, $
  linepars=linepars, esol=esol, echi2=echi2, etemplates=etemplates, $
  debug=debug
; construct the pure emission-line spectrum          

    light = 2.99792458D5 ; speed of light [km/s]

    if keyword_set(broad) then begin
       vmaxshift_iter1 = 500D
       vmaxshift_iter2 = 200D
       sigmamax_iter1 = 5000D
       sigmamax_iter2 = 5000D
       sigma_smooth = 1500D
    endif else begin
       vmaxshift_iter1 = 200D
       vmaxshift_iter2 = 100D
       sigmamax_iter1 = 300D
       sigmamax_iter2 = 300D
       sigma_smooth = 1000D
    endelse

    restwave = wave - alog(zabs+1.0D)
    restflux = flux*(zabs+1.0D)
    restferr = ferr*(zabs+1.0D)
    npix = n_elements(wave)
  
    restlineflux = restflux - continuum - smooth_continuum 

    vsys = alog(zabs+1.0D)*light ; systemic velocity
    voffset = -(min(restwave)-min(tempwave))*light
    kinematics  = [vsys*0.0+voffset,vdisp,0.0,0.0,0.0,0.0]    ; for AGES_GANDALF

; now read the emission-line parameter file of lines we *want* to fit 
    alllinepars = read_gandalf_elinelist(linefile,$
      fitlines=linepars,/actionfit)
    junk = gandalf_mask_emission_lines(npix,vsys,linepars,$
      velscale,wave[0],velscale/light,sigma=sigma_smooth)
    goodpixels = gandalf_mask_emission_lines(npix,vsys,alllinepars,$
      velscale,wave[0],velscale/light,sigma=sigma_smooth)
    
    goodlinepixels = lindgen(npix)
    remove, goodpixels, goodlinepixels
    goodlinepixels = cmset_op(goodlinepixels,'and',$
      where(restferr lt 1E5))
          
; now call GANDALF! for the continuum templates pass an array of
; zeros, and just fit for the emission lines
    sol = kinematics
    im_gandalf, tempwave*0.0, restlineflux, restferr, velscale, sol, $
      linepars, restwave[0], velscale/light, goodpixels=goodlinepixels, $
      int_disp=line_inst_vdisp, bestfit=bestlinefit, weights=lweights, $
      emission_templates=etemplates, error=esol, l0_templ=tempwave[0], $
      /for_errors, plot=plot, quiet=1, chi2=echi2, vmaxshift_iter1=vmaxshift_iter1, $
      vmaxshift_iter2=vmaxshift_iter2, sigmamax_iter1=sigmamax_iter1, $
      sigmamax_iter2=sigmamax_iter2

; parse the output
    if keyword_set(debug) then begin
       djs_plot, exp(restwave), restlineflux, ps=10, xr=[6400,6750] ; 
;      djs_plot, exp(restwave), restlineflux, ps=10;, xr=[3650,4100] ; xr=[6400,6750] ; 
       djs_oplot, exp(restwave), bestlinefit, ps=10, color='green'
;      djs_oplot, exp(restwave[goodpixels]), restlineflux[goodpixels], ps=4, color='blue'
;      djs_oplot, exp(restwave[goodlinepixels]), restlineflux[goodlinepixels], ps=4, color='red'
       splog, 'Waiting for keystroke...' & cc = get_kbrd(1)
    endif
       
;   bestfit = continuum + bestlinefit + smooth_continuum
;   mask = restferr lt 1E5
;   new_sol = gandalf_clean_emission(restwave,restflux,$
;     bestfit,mask,etemplates,linepars,sol,esol,$
;     snrcut=snrcut_line,velscale=velscale,debug=0,$
;     new_etemplates=new_etemplates,new_linepars=new_linepars)
;   elinefit = gandalf_parse_elinefit(new_sol,esol,$
;     new_linepars,vsys=vsys,fluxscale=1.0)

return, sol
end

function chaos_smooth_residuals, wave, flux, continuum, $
  zabs=zabs, velscale=velscale, linefile=linefile, $
  broad=broad, debug=debug
; internal support routine to smooth the residuals

    light = 2.99792458D5 ; speed of light [km/s]
    vsys = alog(zabs+1.0D)*light ; systemic velocity
    if keyword_set(broad) then $
      smooth_sigma = 500.0 else $
        smooth_sigma = 500.0
    
    restwave = wave - alog(zabs+1.0D)
    restflux = flux*(zabs+1.0D)
    npix = n_elements(wave)

    residuals = restflux - continuum
    smooth_continuum = residuals*0.0

    emask = iemission_mask(exp(wave),z=zabs,vdisp=250.0,$
      /sky,/nebular,/telluric,/qso)
    good = where((emask ne 0.0),ngood)
    if (ngood ne 0L) then begin
       smooth1 = medsmooth(residuals[good],31)
       linterp, wave[good], smooth1, wave, smooth1
       smooth_continuum = smooth(smooth1,51,/edge_truncate)
    endif

    if keyword_set(debug) then begin
       djs_plot, exp(restwave), residuals, xsty=3, ysty=3, ps=10, xr=[5500,6800] ; xr=[4700,5050]
       djs_oplot, exp(restwave[good]), residuals[good], ps=6, sym=0.4, color='orange'
       djs_oplot, exp(restwave), smooth1, ps=10, color='green'
       djs_oplot, exp(restwave), smooth_continuum, ps=10, color='red'
       splog, 'Waiting for keystroke...' & cc = get_kbrd(1)
    endif

return, smooth_continuum
end

pro chaos_gandalf_specfit, debug=debug, solar=solar
; jm10dec06ucsd - fit some Chaosspec spectra for Papovich's postdoc 
    
; path names and emission-line file names
    specfitpath = getenv('CHAOS_DATA')+'/'
    alllinefile = specfitpath+'gandalf_elinelist_all.dat'
    linefile = specfitpath+'gandalf_elinelist_v2.0.dat'

; read the templates (see BUILD_AGES_PPXF_TEMPLATES); resample and
; convolve to AGES pixel size and instrumental resolution
    velscale = ages_ppxf_velscale()
    inst_vdisp = ages_ppxf_instvdisp()
    line_inst_vdisp = 0.0 ; inst_vdisp/3.0 ; [km/s]
    if (n_elements(fluxscale) eq 0) then fluxscale = 1D-17

    tempflux = read_chaos_ppxf_templates(tempwave,ntemp=ntemp,$
      velscale=velscale,inst_vdisp=inst_vdisp,tempinfo=info,$
      solar=solar)
    metal = info.metallicity
    nmetal = n_elements(metal)

; read the data, redshifts, and specify the output file name 
    spec1d = read_chaos_sample()
    nobj = n_elements(spec1d)

    if keyword_set(solar) then suffix = '_solar' else suffix = '_subsolar'
    specdatafile = specfitpath+'chaos_specdata_raw'+suffix+'.fits'

; read the emission-line file here to get the total number of lines
; that we are going to be fitting
    junk = read_gandalf_elinelist(linefile,nlines=nlines,$
      fitlines=linepars,/actionfit)

; initialize the output data structure
    nfinalpix = 5000
    specdata_template = {$
      id:         0, $
      z:        0.0, $
      zabs:     0.0, $
      zabs_err:     0.0, $
      vdisp:        0.0, $
      vdisp_err:    0.0, $
      inst_vdisp:           inst_vdisp,$
      line_inst_vdisp: line_inst_vdisp,$
      continuum_snr:               0.0,$
      continuum_coeff:     fltarr(ntemp), $
      continuum_ebv:       0.0,$
      continuum_Z:         0.0,$
      continuum_chi2:      0.0,$ ; reduced chi^2
      emission_line_chi2:  0.0, $

      wave:          dblarr(nfinalpix), $ ; double!
      flux:          fltarr(nfinalpix), $
      ferr:          fltarr(nfinalpix), $
      continuum:     fltarr(nfinalpix), $
      smooth_continuum:        fltarr(nfinalpix),$
      etemplates:       fltarr(nfinalpix,nlines),$
      fitlines:                 intarr(nlines)-1,$
      sol:                      fltarr(4,nlines),$
      esol:                     fltarr(4,nlines)}

; fit each object using GANDALF/PPXF
    t0 = systime(1)
;   for iobj = 11, 11 do begin
    for iobj = 0, nobj-1 do begin
       print, format='("Fitting object ",I0,"/",I0,A10,$)', $
         iobj, nobj-1, string(13b)

; logarithmically rebin the spectrum
       chaos_rebin, spec1d[iobj], velscale=velscale, fluxscale=fluxscale, $
         rebin_wave=wave, rebin_flux=flux, rebin_ferr=ferr, snr=continuum_snr

; loop on each metallicity       
       for imetal = 0, nmetal-1 do begin
          tempflux1 = tempflux[*,*,imetal]
          
;; get the absorption-line redshift and velocity dispersion; if the
;; velocity dispersion is not well-measured then put a floor 
;          chaos_zabs_vdisp, wave, flux, ferr, tempwave=tempwave, $
;            tempflux=tempflux1, velscale=velscale, zabs_guess=spec1d[iobj].z,$
;            vdisp_guess=vdisp_guess, alllinefile=alllinefile, final_zabs=zabs1,$
;            final_vdisp=vdisp1, err_final_zabs=zabs_err1, err_final_vdisp=vdisp_err1, $
;            debug=debug
;          
;          vdisp_snr = vdisp1/(vdisp_err1+(vdisp_err1 eq 0))*(vdisp_err1 ne 0)
;          if (vdisp_snr lt 1.0) then begin
;             vdisp1 = 150.0
;             vdisp_err1 = -1.0
;          endif

;         splog, 'Using default velocity dispersion and redshift!'
          zabs1 = spec1d[iobj].z
          zabs_err1 = -1.0
          vdisp1 = 150.0
          vdisp_err1 = -1.0
          
; fit the full continuum once using the solar-metallicity templates to
; get a better estimate of the absorption-line redshift
          continuum1 = chaos_fit_continuum(wave,flux,ferr,tempwave=tempwave, $
            tempflux=tempflux1,velscale=velscale,zabs=zabs1,vdisp=vdisp1,$
            ebv_guess=ebv_guess,alllinefile=alllinefile,goodpixels=goodpixels,$
            final_ebv=ebv1,chi2=chi21,continuum_coeff=coeff1,debug=debug)

          if (imetal eq 0) then begin
             allzabs = zabs1
             allvdisp = vdisp1
             allzabs_err = zabs_err1
             allvdisp_err = vdisp_err1
             allebv = ebv1
             allchi2 = chi21
             allcoeff = coeff1
             continuum = continuum1
          endif else begin
             allzabs = [allzabs,zabs1]
             allvdisp = [allvdisp,vdisp1]
             allzabs_err = [allzabs,zabs_err1]
             allvdisp_err = [allvdisp_err,vdisp_err1]
             allebv = [allebv,ebv1]
             allchi2 = [allchi2,chi21]
             allcoeff = [[allcoeff],[coeff1]]
             continuum = [[continuum],[continuum1]]
          endelse
       endfor
; find the metallicity that has the minimum chi^2
       continuum_chi2 = min(allchi2,zbest)
       zabs = allzabs[zbest]
       vdisp = allvdisp[zbest]
       zabs_err = allzabs_err[zbest]
       vdisp_err = allvdisp_err[zbest]
       ebv = allebv[zbest]
       continuum = continuum[*,zbest]
       continuum_coeff = allcoeff[*,zbest]
       continuum_Z = metal[zbest]

; smooth the residuals before before fitting the lines
       smooth_continuum = chaos_smooth_residuals(wave,flux,$
         continuum,zabs=zabs,velscale=velscale,broad=broad,$
         linefile=linefile,debug=debug)

; finally fit the emission lines
       sol = chaos_fit_lines(wave,flux,ferr,continuum,smooth_continuum,$
         zabs=zabs,vdisp=vdisp,velscale=velscale,linefile=linefile,$
         line_inst_vdisp=line_inst_vdisp,tempwave=tempwave,$
         broad=broad,linepars=linepars,esol=esol,etemplates=etemplates,$
         echi2=echi2,debug=debug)

; pack everything into a structure          
       specdata1 = specdata_template
       restwave = wave - alog(zabs+1.0D)
       restflux = flux*(zabs+1.0D)
       restferr = ferr*(zabs+1.0D)
       npix = n_elements(wave)
       npix = n_elements(wave)

       specdata1.wave[0:npix-1] = double(restwave)
       specdata1.flux[0:npix-1] = float(restflux*fluxscale)
       specdata1.ferr[0:npix-1] = float(restferr*fluxscale)
       specdata1.continuum[0:npix-1] = float(continuum*fluxscale)
       specdata1.smooth_continuum[0:npix-1] = float(smooth_continuum*fluxscale)

       fitlines = where((linepars.action ne 'i') and $
        (linepars.kind eq 'l'),nfitlines)
       if (nfitlines ne 0) then begin
          specdata1.fitlines[0:nfitlines-1] = fitlines
          specdata1.sol[*,0:nfitlines-1] = sol
          specdata1.esol[*,0:nfitlines-1] = esol
          specdata1.etemplates[0:npix-1,0:nfitlines-1] = float(etemplates*fluxscale)
       endif
       
       specdata1.continuum_snr = continuum_snr
       specdata1.continuum_Z = continuum_Z
       specdata1.continuum_chi2 = continuum_chi2
       specdata1.emission_line_chi2 = echi2
       specdata1.continuum_coeff = continuum_coeff[0:ntemp-1]*fluxscale
       specdata1.continuum_ebv = ebv
       specdata1.id = spec1d[iobj].id
       specdata1.z = spec1d[iobj].z
       specdata1.zabs = zabs
       specdata1.zabs_err = zabs_err
       specdata1.vdisp = vdisp
       specdata1.vdisp_err = vdisp_err
       
       if (n_elements(specdata) eq 0) then specdata = specdata1 else $
         specdata = [temporary(specdata),specdata1]
       
    endfor
    im_mwrfits, specdata, specdatafile, /clobber
    splog, 'Total time = ', (systime(1)-t0)/60.0, ' minutes'
    
return
end






