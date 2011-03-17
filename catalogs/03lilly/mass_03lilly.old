pro mass_03lilly, phot, result, band_shift=band_shift, debug=debug, write=write
; jm06nov20nyu - based on MASS_05SAVAGLIO
    
; define the cosmology and some constants

    if (n_elements(band_shift) eq 0L) then band_shift = 0.0

    red, omega0=0.3, omegalambda=0.7, h100=0.70
    omega0 = redomega0() & omegal = redomegal() & h100 = redh100()

    imf_chabrier_to_salpeter = +0.25
    light = 2.99792458D18       ; speed of light [A/s]

    vname = 'default.nolines'   ; NOTE!
    
; read the photometry    
    
    path = getenv('CATALOGS_DIR')+'/03lilly/'

    if (n_elements(phot) eq 0L) then begin

       phot1 = rsex(path+'cfrs_phot_from_savaglio.dat')

       phot = {galaxy: '', z: 0.0, $
         phot_b: -999.0, phot_b_err: -999.0, $
         phot_v: -999.0, phot_v_err: -999.0, $
         phot_i: -999.0, phot_i_err: -999.0, $
         phot_k: -999.0, phot_k_err: -999.0}
       phot = replicate(phot,n_elements(phot1))

       phot.galaxy = phot1.galaxy
       phot.z      = phot1.z

; the cut on one magnitude of error excludes the K-band magnitude for
; CFRS_03.1375 (z=0.635)       
       
       good = where((phot1.i_iso lt 99.99) and (phot1.b3 lt 99.99) and (phot1.b3_err lt 1.00),ngood)
       if (ngood ne 0L) then begin
          phot[good].phot_b = phot1[good].i_iso + (phot1[good].b3-phot1[good].i3)
          phot[good].phot_b_err = phot1[good].b3_err
       endif
       
       good = where((phot1.i_iso lt 99.99) and (phot1.v3 lt 99.99) and (phot1.v3_err lt 1.00),ngood)
       if (ngood ne 0L) then begin
          phot[good].phot_v = phot1[good].i_iso + (phot1[good].v3-phot1[good].i3)
          phot[good].phot_v_err = phot1[good].v3_err
       endif
       
       good = where((phot1.i_iso lt 99.99) and (phot1.i3 lt 99.99) and (phot1.i3_err lt 1.00),ngood)
       if (ngood ne 0L) then begin
          phot[good].phot_i = phot1[good].i_iso + (phot1[good].i3-phot1[good].i3)
          phot[good].phot_i_err = phot1[good].i3_err
       endif
       
       good = where((phot1.i_iso lt 99.99) and (phot1.k3 lt 99.99) and (phot1.k3_err lt 1.00),ngood)
       if (ngood ne 0L) then begin
          phot[good].phot_k = phot1[good].i_iso + (phot1[good].k3-phot1[good].i3)
          phot[good].phot_k_err = phot1[good].k3_err
       endif
       
    endif
    ngalaxy = n_elements(phot)

; define the rest-frame filters of interest

    splog, 'Loading rest-frame filters:'
;   restfilters = ['bessell_'+['B','V','R'],'sdss_'+['u0','g0','r0'],'twomass_'+['J']]+'.par'
;   restbandpasses = ['B','V','R','sdss_u','sdss_g','sdss_r','J']
;   restvega = [1,1,1,0,0,0,1]
    restfilters = ['bessell_'+['B','V','R'],'twomass_'+['J']]+'.par'
    restbandpasses = ['B','V','R','J']
    restvega = [1,1,1,1]
    nrestfilter = n_elements(restfilters)

    restvega2ab = k_vega2ab(filterlist=restfilters,/kurucz)
    notvega = where((restvega eq 0L),nnotvega)
    if (nnotvega ne 0L) then restvega2ab[notvega] = 0.0
    
; define the filters (NDWFS/SDSS/FLAMEX) for which we have
; observed-frame photometry 

    splog, 'Defining filter functions with observed photometry:'
    obsfilters = ['bessell_'+['B','V','I'],'twomass_Ks']+'.par'
    obsbands = 'phot_'+['b','v','i','k']
    obsvega2ab = k_vega2ab(filterlist=obsfilters,/kurucz)
    nobsfilter = n_elements(obsfilters)

; convert the observed photometry to maggies and store

    splog, 'Converting observed photometry to maggies.'
    obsmaggies = convert_to_maggies(phot,obsbands,obsvega2ab,$
      maggies_invvar=obsmaggies_ivar,mags=mobs_ab,err_mags=mobs_ab_err)
    setzero = where(obsmaggies_ivar eq 1D-32,nsetzero)
    obsmaggies_err = 1.0/sqrt(obsmaggies_ivar)
    if (nsetzero ne 0L) then obsmaggies_ivar[setzero] = 0.0

    redshift = phot.z

; compute the k-corrections      

    splog, 'Fitting k-correct templates.'
    t0 = systime(1)
    kcorrect, obsmaggies, obsmaggies_ivar, redshift, first_kcorrect, $
      band_shift=band_shift, filterlist=obsfilters, filterpath=filters_path, $
      rmatrix=rmatrix, zvals=zvals, lambda=lambda, vmatrix=vmatrix, $
      coeffs=coeffs, rmaggies=obsmaggies_synth, chi2=chi2, maxiter=maxiter, omega0=omega0, $
      omegal0=omegal, mass=mass, intsfh=intsfh, mets=mets, b300=b300, $
      b1000=b1000, mtol=mtol, vname=vname, /silent
    splog, 'Total time = '+string((systime(1)-t0)/60.0,format='(G0.0)')+' minutes.'

; generate the maggies lookup table; only works to z=2 (ie, not for
; all quasars!) and then compute k-corrections
    
    splog, 'Computing k-corrections.'

    k_projection_table, rmatrix, vmatrix, lambda, zvals, restfilters, $ 
      zmin=zmin, zmax=zmax, nz=nz, filterpath=filters_path, /silent

    k_reconstruct_maggies, coeffs, replicate(band_shift,ngalaxy), $
      restmaggies_synth, rmatrix=rmatrix, zvals=zvals
    restmaggies_synth = restmaggies_synth/(1.0+band_shift)

; figure out the bandpass that minimizes the k-correction 
    
    obands = lindgen(nrestfilter) # replicate(1L,ngalaxy) ; index array

    lambda_in = k_lambda_eff(filterlist=obsfilters)
    lambda_out = k_lambda_eff(filterlist=restfilters,band_shift=band_shift)

    for i = 0L, ngalaxy-1L do begin
       for j = 0L, nrestfilter-1L do begin
          dmin = min(abs(lambda_in/(1.0+redshift[i])-lambda_out[j]),imin)
          obands[j,i]= imin
       endfor
    endfor

    kcorrect = fltarr(nrestfilter,ngalaxy)
    for i = 0L, ngalaxy-1L do for j = 0L, nrestfilter-1L do $
      kcorrect[j,i] = restmaggies_synth[j,i]/obsmaggies_synth[obands[j,i],i]
    kcorrect = 2.5*alog10(kcorrect)

; calculate rest-frame magnitudes using the observed photometry when
; available, otherwise the synthesized photometry (should be identical
; when the fit is good); alternatively, we could find the next nearest
; bandpass to compute the k-correction

    mrest_ab = -2.5*alog10(restmaggies_synth)
    mrest_ab_ivar = mrest_ab*0.0

;   absmag = mrest_ab*0.0
;   for i = 0L, nrestfilter-1L do absmag[i,*] = -2.5*alog10(restmaggies_synth[i,*]) - $ ; test
;     lf_distmod(redshift,omega0=omega0,omegal0=omegal0)-5.0*alog10(1.0/h100) 
    
    for j = 0L, ngalaxy-1L do begin

       igood = where((obsmaggies_ivar[obands[*,j],j] gt 0.0) and (obsmaggies[obands[*,j],j] gt 0.0),ngood)

       if (ngood gt 0L) then begin
          mrest_ab[igood,j] = -2.5*alog10(obsmaggies[obands[igood,j],j]) - kcorrect[igood,j]
;         absmag[igood,j] = -2.5*alog10(obsmaggies[obands[igood,j],j])-lf_distmod(redshift[j],$ ; test
;           omega0=omega0,omegal0=omegal0)-5.0*alog10(1.0/h100)-kcorrect[igood,j]
          mrest_ab_ivar[igood,j] = obsmaggies[obands[igood,j],j]^2.0 * $
            obsmaggies_ivar[obands[igood,j],j]*(0.4*alog(10.0))^2.0
       endif

    endfor

; initialize the output data structure

    result_template = {$
      galaxy:                                     '', $
      z:                                      -999.0, $
      kcorr_mass:                             -999.0, $
      kcorr_mets:                             -999.0, $
      kcorr_intsfh:                           -999.0, $
      kcorr_coeffs:                        fltarr(5), $  ; eigentemplate coefficients
      kcorr_chi2:                             -999.0, $ ; reduced chi2
      kcorr_zzmax:                            -999.0, $ ; compute z/zmax (see below)
      kcorr_kcorr:         fltarr(nrestfilter)-999.0, $ ; k-correction
      kcorr_mobs_ab:        fltarr(nobsfilter)-999.0, $ ; "actual observed photometry"
      kcorr_mobs_ab_err:    fltarr(nobsfilter)-999.0, $ ; error in MOBS_AB
      kcorr_mobs_ab_synth:        fltarr(nobsfilter), $ ; "synthesized observed photometry"
      kcorr_mrest:         fltarr(nrestfilter)-999.0, $ ; k-corrected, rest frame, using "actual observed photometry"
      kcorr_mrest_err:     fltarr(nrestfilter)-999.0, $ ; error in MREST
      kcorr_mrest_synth:         fltarr(nrestfilter)}   ; synthesized, rest frame

    for i = 0L, nrestfilter-1L do result_template = create_struct(result_template,'M_'+restbandpasses[i],-999.0)

    result = replicate(result_template,ngalaxy)
    
; store all the results; convert the stellar mass to the Salpeter IMF
; and to h=0.7

    result.galaxy            = phot.galaxy
    result.z                 = phot.z
    result.kcorr_mobs_ab     = mobs_ab
    result.kcorr_mobs_ab_err = mobs_ab_err

    good = where((mass gt 1E6) and finite(mass),ngood)
    if (ngood ne 0L) then begin

       result[good].kcorr_mass = alog10(mass[good]) - alog10(h100) + imf_chabrier_to_salpeter
       result[good].kcorr_mets = mets[good]
       result[good].kcorr_intsfh = alog10(intsfh[good])
       result[good].kcorr_coeffs = coeffs[*,good]
       result[good].kcorr_chi2 = chi2

       result[good].kcorr_kcorr = kcorrect[*,good]
       result[good].kcorr_mobs_ab_synth = -2.5*alog10(obsmaggies_synth[*,good])

       result[good].kcorr_mrest = mrest_ab[*,good] - rebin(restvega2ab,nrestfilter,ngood)
       result[good].kcorr_mrest_synth = -2.5*alog10(restmaggies_synth[*,good]) - rebin(restvega2ab,nrestfilter,ngood)

; compute absolute magnitudes       
       
       for iband = 0L, nrestfilter-1L do begin
          true = tag_exist(result,'M_'+restbandpasses[iband],index=bandindx)
          result[good].(bandindx) = result[good].kcorr_mrest[iband] - lf_distmod(redshift[good],$
            omega0=omega0,omegal0=omegal0) - 5.0*alog10(1.0/h100) 
       endfor
       
    endif 

; debugging plot    
    
    if keyword_set(debug) then begin

       plotsym, 0, 2.0, /fill, color=djs_icolor('red')
       im_window, 0, xratio=0.5, yratio=0.4

       obsfiltinfo = im_filterspecs(filterlist=obsfilters)
       k_load_vmatrix, vmatrix, restwave, vfile=vfile, vpath=vpath, lfile=lfile, vname=vname

       for i = 0L, ngalaxy-1L do begin

          good = where((result[i].kcorr_mobs_ab gt 0.0),ngood)

          if (ngood ne 0L) then begin
             
             wave = restwave*(1.0+phot[i].z)
             restflux = vmatrix#double(result[i].kcorr_coeffs) ; f_lambda
             flux = restflux/(1.0+phot[i].z)
             flux_fnu = flux*wave^2.0/light*10.0^(0.4*48.6) ; f_nu
             flux_ab = -2.5*alog10(flux_fnu>1D-50) ; AB
             
             xrange = [min(obsfiltinfo[good].weff-1.5*obsfiltinfo[good].fwhm),$
               max(obsfiltinfo[good].weff+1.5*obsfiltinfo[good].fwhm)] ;/(1.0+zobj)
             xrange[0] = xrange[0]>90.0
             get_element, wave, xrange, xx

             if keyword_set(fnu) then begin
                
                yrange = [min(obsmaggies[good,i]-obsmaggies_err[good,i]),$
                  max(obsmaggies[good,i]+obsmaggies_err[good,i])]
                yrange[0] = yrange[0]<min(flux_fnu[xx[0]:xx[1]])
                yrange[1] = yrange[1]>max(flux_fnu[xx[0]:xx[1]])
                
                plot, [0], [0], /nodata, xsty=3, ysty=3, charsize=1.5, charthick=2.0, xthick=2.0, $
                  ythick=2.0, xtitle='Observed Wavelength ['+angstrom()+']', ytitle=textoidl('f_{\nu}'), yrange=yrange, $
                  xrange=xrange, title=string(phot[i].galaxy,format='(A0)')+' z = '+$
                  strtrim(string(phot[i].z,format='(F12.3)'),2)
                oploterror, obsfiltinfo[good].weff, obsmaggies[good,i], obsfiltinfo[good].fwhm, $
                  obsmaggies_err[good,i], ps=8
                oplot, wave, flux_fnu, line=0.1

             endif else begin

                yrange = fltarr(2)
                yrange[0] = max(-2.5*alog10(obsmaggies[good,i]))>max(flux_ab[xx[0]:xx[1]]*1.00)
                yrange[1] = min(-2.5*alog10(obsmaggies[good,i]))<min(flux_ab[xx[0]:xx[1]]*0.98)

                plot, [0], [0], /nodata, xsty=3, ysty=3, charsize=1.5, charthick=2.0, xthick=2.0, $
                  ythick=2.0, xtitle='Observed Wavelength ['+angstrom()+']', ytitle=textoidl('m_{AB}'), $
                  yrange=yrange, xrange=xrange, title=string(phot[i].galaxy,format='(A0)')+' z = '+$
                  strtrim(string(phot[i].z,format='(F12.3)'),2)
                oploterror, obsfiltinfo[good].weff, -2.5*alog10(obsmaggies[good,i]), obsfiltinfo[good].fwhm, $
                  2.5*obsmaggies_err[good,i]/obsmaggies[good,i]/alog(10.0), ps=8
                oplot, wave, flux_ab, line=0.1

             endelse
             
             label = textoidl('\chi^{2}_{\nu} = '+strtrim(string(result[i].kcorr_chi2,format='(F12.2)'),2))
             legend, label, /right, /bottom, box=0, charsize=1.5, charthick=2.0

             cc = get_kbrd(1)

          endif

       endfor

    endif
    
; write out    
    
    if keyword_set(write) then begin

       splog, 'Writing '+path+'mass_03lilly_kcorr.fits.gz'
       mwrfits, result, path+'mass_03lilly_kcorr.fits', /create
       spawn, ['gzip -f '+path+'mass_03lilly_kcorr.fits'], /sh
       
    endif

return
end
