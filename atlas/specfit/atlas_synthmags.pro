;+
; NAME:
;       ATLAS_SYNTHMAGS
;
; PURPOSE:
;       Synthesize a variety of magnitudes from the integrated
;       spectral atlas.  Use the best-fitting Bruzual & Charlot models
;       to extend the spectra to bluer and redder wavelengths.
;
; CALLING SEQUENCE:
;       atlas_synthmags, outfile=, /doplot, /write
;
; INPUTS:
;       None.
;
; OPTIONAL INPUTS:
;       outfile - output file name
;
; KEYWORD PARAMETERS:
;       doplot - plot the individual SED's and filter curves
;       write  - write out the synthesized magnitudes structure 
;
; OUTPUTS:
;       See WRITE keyword.
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Feb 02, U of A - written
;       jm04mar10uofa - added 2MASS magnitudes and general cleanup;
;                       documented
;       jm04mar14uofa - added OUTFILE optional input
;       jm04jul28uofa - extensive updates, utilizing Blanton's
;                       k-correct filter routines; also use
;                       Monte-Carlo based error estimation
;-

pro atlas_synthmags, atlasinfo1, atlas1, outfile=outfile, doplot=doplot, write=write

    if n_elements(outfile) eq 0L then outfile = 'integrated_atlas_synthmags.fits'
    
    datapath = atlas_path(/atlas1d)
    outpath = atlas_path(/specfit)

; read the structures and sort by RA    
    
    if (n_elements(atlasinfo1) eq 0L) then atlasinfo1 = atlas_read_info()
    if (n_elements(atlas1) eq 0L) then atlas1 = read_integrated()

    atlasinfo = atlasinfo1[sort(atlasinfo1.ra)]
    atlas = atlas1[sort(atlas1.ra)]
    natlas = n_elements(atlas)
    natlasinfo = n_elements(atlasinfo)

; only synthesize magnitudes for objects that have been fitted

;   skipindx = cmset_op(atlasinfo.atlas_id,'AND',/not2,atlas.atlas_id,/index)
;   fitindx = cmset_op(atlasinfo.atlas_id,'AND',atlas.atlas_id,/index)
    
;   skipindx = cmset_op(atlasinfo.atlas_id,'AND',/not2,atlas.atlas_id,/index)
;   niceprint, atlasfit.galaxy, atlas[skipindx].galaxy
    
    eigendir = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='templates')

    nmonte = 100L
    
    filterlist = [$
      'bessell_U',$
      'bessell_B',$
      'bessell_V',$
      'bessell_R',$
      'bessell_I',$
      'twomass_J',$
      'twomass_H',$
      'twomass_Ks',$
      'sdss_u0',$
      'sdss_g0',$
      'sdss_r0',$
      'sdss_i0',$
      'sdss_z0' $
      ]+'.par'
    nfilters = n_elements(filterlist)

    filt = im_filterspecs(filterlist=filterlist,/verbose)
    vega2ab = filt.vega2ab
    
    sdss = where(filt.sdssflag,nsdss)
    if (nsdss ne 0L) then vega2ab[sdss] = 0.0 ; keep the SDSS magnitudes on the AB system
    
    Umsun = filt[0].solarmags
    Bmsun = filt[1].solarmags
    Vmsun = filt[2].solarmags
    Rmsun = filt[3].solarmags
    Imsun = filt[4].solarmags

    Jmsun = filt[5].solarmags
    Hmsun = filt[6].solarmags
    Ksmsun = filt[7].solarmags

    umsun_sdss = filt[8].solarmags
    gmsun_sdss = filt[9].solarmags
    rmsun_sdss = filt[10].solarmags
    imsun_sdss = filt[11].solarmags
    zmsun_sdss = filt[12].solarmags

    k_load_filters, filterlist, filter_nlambda, filter_lambda, filter_pass

    ufilt = {$
      filtw: filter_lambda[0L:filter_nlambda[0]-1L,0], $
      filtf: filter_pass[0L:filter_nlambda[0]-1L,0]}
    bfilt = {$
      filtw: filter_lambda[0L:filter_nlambda[1]-1L,1], $
      filtf: filter_pass[0L:filter_nlambda[1]-1L,1]}
    vfilt = {$
      filtw: filter_lambda[0L:filter_nlambda[2]-1L,2], $
      filtf: filter_pass[0L:filter_nlambda[2]-1L,2]}
    rfilt = {$
      filtw: filter_lambda[0L:filter_nlambda[3]-1L,3], $
      filtf: filter_pass[0L:filter_nlambda[3]-1L,3]}
    ifilt = {$
      filtw: filter_lambda[0L:filter_nlambda[4]-1L,4], $
      filtf: filter_pass[0L:filter_nlambda[4]-1L,4]}
    
    jfilt = {$
      filtw: filter_lambda[0L:filter_nlambda[5]-1L,5], $
      filtf: filter_pass[0L:filter_nlambda[5]-1L,5]}
    hfilt = {$
      filtw: filter_lambda[0L:filter_nlambda[6]-1L,6], $
      filtf: filter_pass[0L:filter_nlambda[6]-1L,6]}
    ksfilt = {$
      filtw: filter_lambda[0L:filter_nlambda[7]-1L,7], $
      filtf: filter_pass[0L:filter_nlambda[7]-1L,7]}
    
    sdss_u = {$
      filtw: filter_lambda[0L:filter_nlambda[8]-1L,8], $
      filtf: filter_pass[0L:filter_nlambda[8]-1L,8]}
    sdss_g = {$
      filtw: filter_lambda[0L:filter_nlambda[9]-1L,9], $
      filtf: filter_pass[0L:filter_nlambda[9]-1L,9]}
    sdss_r = {$
      filtw: filter_lambda[0L:filter_nlambda[10]-1L,10], $
      filtf: filter_pass[0L:filter_nlambda[10]-1L,10]}
    sdss_i = {$
      filtw: filter_lambda[0L:filter_nlambda[11]-1L,11], $
      filtf: filter_pass[0L:filter_nlambda[11]-1L,11]}
    sdss_z = {$
      filtw: filter_lambda[0L:filter_nlambda[12]-1L,12], $
      filtf: filter_pass[0L:filter_nlambda[12]-1L,12]}

    if keyword_set(doplot) then im_window, 0, xratio=0.9, yratio=0.7
    
; initialize the output data structure
    
    result = {$
      galaxy:             ' ', $
      distance:           0.0, $
      distance_err:       0.0, $

      synth_U:         -999.0, $ ; emission-line corrected magnitude
      synth_U_err:     -999.0, $
      synth_U_obs:     -999.0, $ ; observed magnitude (including emission lines) 
      synth_U_obs_err: -999.0, $
      synth_U_lum:     -999.0, $
      synth_U_lum_err: -999.0, $
      synth_M_U:       -999.0, $
      synth_M_U_err:   -999.0, $

      synth_B:         -999.0, $
      synth_B_err:     -999.0, $
      synth_B_obs:     -999.0, $
      synth_B_obs_err: -999.0, $
      synth_B_lum:     -999.0, $
      synth_B_lum_err: -999.0, $
      synth_M_B:       -999.0, $
      synth_M_B_err:   -999.0, $

      synth_V:         -999.0, $
      synth_V_err:     -999.0, $
      synth_V_obs:     -999.0, $
      synth_V_obs_err: -999.0, $
      synth_V_lum:     -999.0, $
      synth_V_lum_err: -999.0, $
      synth_M_V:       -999.0, $
      synth_M_V_err:   -999.0, $

      synth_R:         -999.0, $
      synth_R_err:     -999.0, $
      synth_R_obs:     -999.0, $
      synth_R_obs_err: -999.0, $
      synth_R_lum:     -999.0, $
      synth_R_lum_err: -999.0, $
      synth_M_r:       -999.0, $
      synth_M_r_err:   -999.0, $

      synth_I:         -999.0, $
      synth_I_err:     -999.0, $
      synth_I_obs:     -999.0, $
      synth_I_obs_err: -999.0, $
      synth_I_lum:     -999.0, $
      synth_I_lum_err: -999.0, $
      synth_M_I:       -999.0, $
      synth_M_I_err:   -999.0, $

      synth_J:         -999.0, $
      synth_J_err:     -999.0, $
      synth_J_obs:     -999.0, $
      synth_J_obs_err: -999.0, $
      synth_J_lum:     -999.0, $
      synth_J_lum_err: -999.0, $
      synth_M_J:       -999.0, $
      synth_M_J_err:   -999.0, $

      synth_H:         -999.0, $
      synth_H_err:     -999.0, $
      synth_H_obs:     -999.0, $
      synth_H_obs_err: -999.0, $
      synth_H_lum:     -999.0, $
      synth_H_lum_err: -999.0, $
      synth_M_H:       -999.0, $
      synth_M_H_err:   -999.0, $

      synth_Ks:         -999.0, $
      synth_Ks_err:     -999.0, $
      synth_Ks_obs:     -999.0, $
      synth_Ks_obs_err: -999.0, $
      synth_Ks_lum:     -999.0, $
      synth_Ks_lum_err: -999.0, $
      synth_M_Ks:       -999.0, $
      synth_M_Ks_err:   -999.0, $

      synth_sdss_u:            -999.0, $
      synth_sdss_u_err:        -999.0, $
      synth_sdss_u_obs:        -999.0, $
      synth_sdss_u_obs_err:    -999.0, $
      synth_sdss_u_lum:        -999.0, $
      synth_sdss_u_lum_err:    -999.0, $
      synth_sdss_M_u:          -999.0, $
      synth_sdss_M_u_err:      -999.0, $

      synth_sdss_g:            -999.0, $
      synth_sdss_g_err:        -999.0, $
      synth_sdss_g_obs:        -999.0, $
      synth_sdss_g_obs_err:    -999.0, $
      synth_sdss_g_lum:        -999.0, $
      synth_sdss_g_lum_err:    -999.0, $
      synth_sdss_M_g:          -999.0, $
      synth_sdss_M_g_err:      -999.0, $

      synth_sdss_r:            -999.0, $
      synth_sdss_r_err:        -999.0, $
      synth_sdss_r_obs:        -999.0, $
      synth_sdss_r_obs_err:    -999.0, $
      synth_sdss_r_lum:        -999.0, $
      synth_sdss_r_lum_err:    -999.0, $
      synth_sdss_M_r:          -999.0, $
      synth_sdss_M_r_err:      -999.0, $

      synth_sdss_i:            -999.0, $
      synth_sdss_i_err:        -999.0, $
      synth_sdss_i_obs:        -999.0, $
      synth_sdss_i_obs_err:    -999.0, $
      synth_sdss_i_lum:        -999.0, $
      synth_sdss_i_lum_err:    -999.0, $
      synth_sdss_M_i:          -999.0, $
      synth_sdss_M_i_err:      -999.0, $

      synth_sdss_z:            -999.0, $
      synth_sdss_z_err:        -999.0, $
      synth_sdss_z_obs:        -999.0, $
      synth_sdss_z_obs_err:    -999.0, $
      synth_sdss_z_lum:        -999.0, $
      synth_sdss_z_lum_err:    -999.0, $
      synth_sdss_M_z:          -999.0, $
      synth_sdss_M_z_err:      -999.0}

    bigresult = replicate(result,natlasinfo)
    result = replicate(result,natlas)

    result.galaxy       = atlas.galaxy
    result.distance     = atlas.distance
    result.distance_err = atlas.distance_err

    bigresult.galaxy       = atlasinfo.galaxy
    bigresult.distance     = atlasinfo.distance
    bigresult.distance_err = atlasinfo.distance_err

; read the Bruzual & Charlot models 

    eigenflux_Z004 = mrdfits(eigendir+'BC03_Z004_salpeter_templates.fits',2,eigenhead)
    eigenflux_Z02 = mrdfits(eigendir+'BC03_Z02_salpeter_templates.fits',2)
    eigenflux_Z05 = mrdfits(eigendir+'BC03_Z05_salpeter_templates.fits',2)

    eigenwave = make_wave(eigenhead)
    neigen = (size(eigenflux_Z004,/dimension))[1]

    bigwave = eigenwave
    bestspec = bigwave*0.0
    redcurve = k_lambda(eigenwave,/charlot)

    t0 = systime(1)
;   for i = 0L, 9L do begin
    for i = 0L, natlas-1L do begin

       print, format='("Synthesizing magnitudes for object ",I0,"/",I0,".",A1,$)', $
         i+1, natlas, string(13b)

       galaxy = strtrim(atlas[i].galaxy,2)
       nicegalaxy = atlas[i].nice_galaxy
       z = atlas[i].z_abs
       medsnr = atlas[i].continuum_snr ; median S/N in the continuum
       
; select the appropriate BC03 model       

       if strmatch(atlas[i].templates,'*Z004*',/fold) then eigenflux = eigenflux_Z004
       if strmatch(atlas[i].templates,'*Z02*',/fold) then eigenflux = eigenflux_Z02
       if strmatch(atlas[i].templates,'*Z05*',/fold) then eigenflux = eigenflux_Z05
       
; read the rest-frame spectrum and emission-line spectrum; we have to
; read the error spectrum from the original FITS file because we don't
; store it in the ISPECLINEFIT(); redshift the best-fitting spectrum
       
       specdata = read_atlas_specfit(galaxy,/silent) ; this will need to be modified!
       refwave = reform(specdata[*,0])*(1+z)
       flux = reform(specdata[*,1])/(1+z)
       testflux = mrdfits(datapath+strtrim(atlas[i].drift_file,2),0,/silent)
       ferr = mrdfits(datapath+strtrim(atlas[i].drift_file,2),1,/silent)
       eflux = reform(specdata[*,3])/(1+z)
       npix = n_elements(flux)
       
       flux_noelines = flux - eflux ; spectrum with no emission lines
       
       get_element, refwave, djs_mean(refwave), snrindx

;      spec = rd1dspec(atlas[i].driftfile,/silent,datapath=datapath)
;      flux = spec.spec
;      ferr = spec.sigspec
;      refwave = spec.wave

; compute the wavelength vector at the reference redshift; interpolate
; the spectrum and error spectrum onto BIGWAVE, replacing extrapolated
; values with 1.0
       
       linterp, refwave, flux, bigwave, bigflux, missing=1.0
       linterp, refwave, ferr, bigwave, bigfluxerr, missing=1.0
       linterp, refwave, flux_noelines, bigwave, bigflux_noelines, missing=1.0

       nbigpix = n_elements(bigflux)

; mask telluric bands

       tmask = telluric_mask(bigwave) & these = where(tmask eq 0,nthese)
       if (nthese ne 0L) then begin
          bigflux[these] = 1.0
          bigfluxerr[these] = 1.0
          bigflux_noelines[these] = 1.0
       endif
       
;      plot, bigwave, bigflux, yr=minmax(flux), xsty=3, xr=minmax(refwave)

; reconstruct the best-fitting continuum spectrum, including reddening
; and extrapolate the integrated spectrum to longer and shorter
; wavelengths using the BC03 models

       backcoeff = [atlas[i].continuum_coeff,atlas[i].continuum_ebv,$
         atlas[i].continuum_ebv*0.0+1.0]
       bestspec = backmodel(eigenflux,backcoeff,redcurve=redcurve,$
         nstar=atlas[i].ntemplate,nback=0,ndust=atlas[i].ntemplate)
       linterp, eigenwave, bestspec, bigwave, bestspec

; replace the observed galaxy spectrum where there are no data with
; the best-fitting BC03 model

       nodata = where(bigflux eq 1.0,nno)
       if (nno ne 0L) then bigflux[nodata] = bestspec[nodata]
       if (nno ne 0L) then bigflux_noelines[nodata] = bestspec[nodata]
       if (nno ne 0L) then bigfluxerr[nodata] = bestspec[nodata]/sqrt(bestspec[nodata]/bestspec[snrindx])/medsnr

;      ploterror, bigwave, bigflux, bigfluxerr, ps=10, xsty=3, ysty=3

; synthesize magnitudes! 

       bigwave_edges = k_lambda_to_edges(bigwave)

       bigflux_monte = (bigflux # replicate(1.0,nmonte)) + $
         randomn(seed,nbigpix,nmonte)*(bigfluxerr # replicate(1.0,nmonte))

       maggies = reform(k_project_filters(bigwave_edges,bigflux,filterlist=filterlist))
       maggies_noelines = reform(k_project_filters(bigwave_edges,bigflux_noelines,filterlist=filterlist))
       maggies_monte = reform(k_project_filters(bigwave_edges,bigflux_monte,filterlist=filterlist))

       vegamags = -2.5*alog10(maggies) - vega2ab
       vegamags_noelines = -2.5*alog10(maggies_noelines) - vega2ab

       vegamags_err = vegamags*0.0
       for ifilter = 0L, nfilters-1L do vegamags_err[ifilter] = stddev(-2.5*alog10(maggies_monte[*,ifilter]))

       result[i].synth_U          = vegamags_noelines[0]
       result[i].synth_U_err      = vegamags_err[0]
       result[i].synth_U_obs      = vegamags[0]
       result[i].synth_U_obs_err  = vegamags_err[0]

       result[i].synth_B          = vegamags_noelines[1]
       result[i].synth_B_err      = vegamags_err[1]
       result[i].synth_B_obs      = vegamags[1]
       result[i].synth_B_obs_err  = vegamags_err[1]

       result[i].synth_V          = vegamags_noelines[2]
       result[i].synth_V_err      = vegamags_err[2]
       result[i].synth_V_obs      = vegamags[2]
       result[i].synth_V_obs_err  = vegamags_err[2]

       result[i].synth_R          = vegamags_noelines[3]
       result[i].synth_R_err      = vegamags_err[3]
       result[i].synth_R_obs      = vegamags[3]
       result[i].synth_R_obs_err  = vegamags_err[3]

       result[i].synth_I          = vegamags_noelines[4]
       result[i].synth_I_err      = vegamags_err[4]
       result[i].synth_I_obs      = vegamags[4]
       result[i].synth_I_obs_err  = vegamags_err[4]

       result[i].synth_J          = vegamags_noelines[5]
       result[i].synth_J_err      = vegamags_err[5]
       result[i].synth_J_obs      = vegamags[5]
       result[i].synth_J_obs_err  = vegamags_err[5]

       result[i].synth_H          = vegamags_noelines[6]
       result[i].synth_H_err      = vegamags_err[6]
       result[i].synth_H_obs      = vegamags[6]
       result[i].synth_H_obs_err  = vegamags_err[6]

       result[i].synth_Ks         = vegamags_noelines[7]
       result[i].synth_Ks_err     = vegamags_err[7]
       result[i].synth_Ks_obs     = vegamags[7]
       result[i].synth_Ks_obs_err = vegamags_err[7]

       result[i].synth_sdss_u         = vegamags_noelines[8]
       result[i].synth_sdss_u_err     = vegamags_err[8]
       result[i].synth_sdss_u_obs     = vegamags[8]
       result[i].synth_sdss_u_obs_err = vegamags_err[8]

       result[i].synth_sdss_g         = vegamags_noelines[9]
       result[i].synth_sdss_g_err     = vegamags_err[9]
       result[i].synth_sdss_g_obs     = vegamags[9]
       result[i].synth_sdss_g_obs_err = vegamags_err[9]

       result[i].synth_sdss_r         = vegamags_noelines[10]
       result[i].synth_sdss_r_err     = vegamags_err[10]
       result[i].synth_sdss_r_obs     = vegamags[10]
       result[i].synth_sdss_r_obs_err = vegamags_err[10]

       result[i].synth_sdss_i         = vegamags_noelines[11]
       result[i].synth_sdss_i_err     = vegamags_err[11]
       result[i].synth_sdss_i_obs     = vegamags[11]
       result[i].synth_sdss_i_obs_err = vegamags_err[11]

       result[i].synth_sdss_z         = vegamags_noelines[12]
       result[i].synth_sdss_z_err     = vegamags_err[12]
       result[i].synth_sdss_z_obs     = vegamags[12]
       result[i].synth_sdss_z_obs_err = vegamags_err[12]

       if keyword_set(doplot) then begin

;         xrange = [6600.0,7200]
          xrange = [2500.0,11000.0]
          get_element, bigwave, xrange, ww

          xwave = bigwave[ww[0]:ww[1]]
          yflux = bigflux[ww[0]:ww[1]]
          yflux_noelines = bigflux_noelines[ww[0]:ww[1]]
          norm = max(yflux)
          yrange = [min(yflux)>0,max(yflux)]/norm
          
          djs_plot, xwave, yflux/norm, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
            xthick=2.0, ythick=2.0, charsize=2.0, charthick=2.0, xtitle='Wavelength', $
            ytitle='Normalized Flux', title='Johnson/Cousins Synthetic Magnitudes', $
            color='orange', thick=1.0, ps=10

          djs_oplot, xwave, yflux_noelines/norm, thick=1.0, ps=10
          djs_oplot, bigwave[nodata], bigflux[nodata]/norm, ps=3, color='cyan'
          legend, nicegalaxy, /right, /top, box=0, charsize=1.5, charthick=2.0

;         djs_oplot, refwave, testflux/norm, thick=3.0, color='yellow'

          djs_oplot, ufilt.filtw, ufilt.filtf/max(ufilt.filtf), line=0, color='purple'
          djs_oplot, bfilt.filtw, bfilt.filtf/max(bfilt.filtf), line=0, color='blue'
          djs_oplot, vfilt.filtw, vfilt.filtf/max(vfilt.filtf), line=0, color='yellow'
          djs_oplot, rfilt.filtw, rfilt.filtf/max(rfilt.filtf), line=0, color='green'
          djs_oplot, ifilt.filtw, ifilt.filtf/max(rfilt.filtf), line=0, color='red'

          cc = get_kbrd(1)
          
; ---------------------------------------------------------------------------
; SDSS AB magnitudes
; ---------------------------------------------------------------------------

          xrange = [2500.0,11000.0]
          get_element, bigwave, xrange, ww
          xwave = bigwave[ww[0]:ww[1]]
          yflux = bigflux[ww[0]:ww[1]]
          norm = max(yflux)
          yrange = [min(yflux)>0,max(yflux)]/norm
  
          djs_plot, xwave, yflux/norm, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
            xthick=2.0, ythick=2.0, charsize=2.0, charthick=2.0, xtitle='Wavelength', $
            ytitle='Normalized Flux', title='SDSS Synthetic Magnitudes', color='orange', thick=1.0
          djs_oplot, xwave, yflux_noelines/norm, thick=1.0
          legend, nicegalaxy, /right, /top, box=0, charsize=1.5, charthick=2.0
          djs_oplot, bigwave[nodata], bigflux[nodata]/norm, ps=3, color='cyan'

          djs_oplot, sdss_u.filtw, sdss_u.filtf/max(sdss_u.filtf), line=0, color='purple'
          djs_oplot, sdss_g.filtw, sdss_g.filtf/max(sdss_g.filtf), line=0, color='blue'
          djs_oplot, sdss_r.filtw, sdss_r.filtf/max(sdss_r.filtf), line=0, color='yellow'
          djs_oplot, sdss_i.filtw, sdss_i.filtf/max(sdss_i.filtf), line=0, color='green'
          djs_oplot, sdss_z.filtw, sdss_z.filtf/max(sdss_z.filtf), line=0, color='red'
  
          cc = get_kbrd(1)
          
; ---------------------------------------------------------------------------
; 2MASS Vega magnitudes
; ---------------------------------------------------------------------------

          xrange = [10000.0,24000.0]
          get_element, bigwave, xrange, ww

          xwave = bigwave[ww[0]:ww[1]]
          yflux = bigflux[ww[0]:ww[1]]
          yflux_noelines = bigflux_noelines[ww[0]:ww[1]]
          norm = max(yflux)
          yrange = [min(yflux)>0,max(yflux)]/norm
          
          djs_plot, xwave, yflux/norm, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
            xthick=2.0, ythick=2.0, charsize=2.0, charthick=2.0, xtitle='Wavelength', $
            ytitle='Normalized Flux', title='2MASS Synthetic Magnitudes', $
            color='orange', thick=1.0
          djs_oplot, xwave, yflux_noelines/norm, thick=1.0
          djs_oplot, bigwave[nodata], bigflux[nodata]/norm, ps=3, color='cyan'
          legend, nicegalaxy, /right, /top, box=0, charsize=1.5, charthick=2.0

          djs_oplot, jfilt.filtw, jfilt.filtf/max(jfilt.filtf), line=0, color='blue'
          djs_oplot, hfilt.filtw, hfilt.filtf/max(hfilt.filtf), line=0, color='yellow'
          djs_oplot, ksfilt.filtw, ksfilt.filtf/max(ksfilt.filtf), line=0, color='red'

          cc = get_kbrd(1)

       endif

       print, atlas[i].rc3_b, result[i].synth_b, atlas[i].rc3_v, result[i].synth_v
       
    endfor
    splog, 'Total time = '+string((systime(1)-t0)/60.0,format='(G0.0)')+' minutes.'

; compute luminosities and absolute magnitudes
    
    splog, 'Computing absolute magnitudes and luminosities.'

    band = ['U','B','V','R','I','twomass_J','twomass_H','twomass_Ks','sdss_u','sdss_g','sdss_r','sdss_i','sdss_z']

    tags = 'synth_'+['u','b','v','r','i','j','h','ks','sdss_u','sdss_g','sdss_r','sdss_i','sdss_z']
    tags_err = tags+'_err'

    abstags = 'synth_'+['m_u','m_b','m_v','m_r','m_i','m_j','m_h','m_ks','sdss_m_u','sdss_m_g','sdss_m_r','sdss_m_i','sdss_m_z']
    abstags_err = abstags+'_err'
    
    lumtags = 'synth_'+['u','b','v','r','i','j','h','ks','sdss_'+['u','g','r','i','z']]+'_lum'
;   lumtags = 'synth_'+[['u','b','v','r','i','j','h','ks']+'_lum','sdss_lum_'+['u','g','r','i','z']]
    lumtags_err = lumtags+'_err'

    for iband = 0L, n_elements(tags)-1L do begin

       true = tag_exist(result,tags[iband],index=tagsindx)
       true = tag_exist(result,tags_err[iband],index=tagsindx_err)
       true = tag_exist(result,abstags[iband],index=abstagsindx)
       true = tag_exist(result,abstags_err[iband],index=abstagsindx_err)
       true = tag_exist(result,lumtags[iband],index=lumtagsindx)
       true = tag_exist(result,lumtags_err[iband],index=lumtagsindx_err)

       good = where((result.distance gt -900.0) and (result.(tagsindx) gt -900.0),ngood)
       if (ngood ne 0L) then begin

          mags = im_absolute_magnitudes(band[iband],result[good].(tagsindx),$
            result[good].distance,mag_err=result[good].(tagsindx_err),$
            distance_err=result[good].distance_err)

          result[good].(abstagsindx)     = mags.absmag
          result[good].(abstagsindx_err) = mags.absmag_err
          
          result[good].(lumtagsindx)     = mags.lum
          result[good].(lumtagsindx_err) = mags.lum_err
          
       endif
       
    endfor

; copy "result" into "bigresult"

    doit = match_string(atlas.galaxy,atlasinfo.galaxy,/exact,findex=index)
    niceprint, atlas.galaxy, atlasinfo[index].galaxy

    bigresult[index] = result

; write out
    
    if keyword_set(write) then begin

       splog, 'Writing '+outpath+outfile+'.'
       mwrfits, bigresult, outpath+outfile, /create
       spawn, ['gzip -f '+outpath+outfile], /sh

    endif

stop

return
end
    
