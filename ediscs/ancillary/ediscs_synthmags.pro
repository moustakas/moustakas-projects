pro ediscs_synthmags, result, model=model, specfit=specfit, debug=debug, write=write
; jm07apr01nyu - synthesize magnitudes for the EDISCS spectra 

; read the data    
    datapath = ediscs_path(/specfit)
    spec1dpath = ediscs_path(/spec1d)
    eigendir = ediscs_path(/specfit)
    version = ediscs_version(/specfit)

    ediscs = read_ediscs(/specfit)
    ancillary = read_ediscs(/ancillary)
;   ww = where(strmatch(ediscs.galaxy,'*noname*',/fold))
;   ediscs = ediscs[ww]
;   ancillary = ancillary[ww]
    ngal = n_elements(ediscs)

; read the spectral fitting results
    if (n_elements(specfit) eq 0L) then begin
       splog, 'Reading '+datapath+'ediscs_specfit_'+version+'.fits.gz'
       specfit = mrdfits(datapath+'ediscs_specfit_'+version+'.fits.gz',0,/silent)
    endif
    
; reconstruct the best-fitting BC03 models
    if (n_elements(model) eq 0) then begin
       splog, 'Reconstructing the BC03 model fits'
       model = irestore_speclinefit(ediscs,templatedir=eigendir)
    endif

; load the observed and rest-frame filters 
    ubvri_band_shift = 0.0 ; z=0.0
    ubvri_filters = bessell_filterlist()

;   restfilters = ['bessell_'+['U','B','V','R'],'sdss_'+['u0','g0','r0']]+'.par'
;   restbandpasses = ['U','B','V','R','sdss_u','sdss_g','sdss_r']
;   restvega = [1,1,1,1,0,0,0]
;   nrestfilter = n_elements(restfilters)
;   restvega2ab = k_vega2ab(filterlist=restfilters,/kurucz,/silent)
;   notvega = where((restvega eq 0L),nnotvega)
;   if (nnotvega ne 0L) then restvega2ab[notvega] = 0.0

    obsfilters = ['FORS_B','FORS_V','FORS2_R','FORS_I']+'_ccd_atm.par'
    obsbandpasses = ['B','V','R','I']
    obsvega2ab = k_vega2ab(filterlist=obsfilters,/kurucz,/silent)
    nobsfilter = n_elements(obsfilters)

    k_load_filters, obsfilters, filt_nlam, filt_lam, filt_resp

; output structure
    result = {$
      galaxy:                                  '', $
      z:                                   -999.0, $
      synth_obsfilters:                obsfilters, $
      synth_obsmaggies:  fltarr(nobsfilter)-999.0, $
      synth_ubvri_absmag:         fltarr(5)-999.0, $
      synth_ubvri_kcorrect:       fltarr(5)-999.0, $
      synth_model_ubvri_absmag:   fltarr(5)-999.0, $
      synth_model_ubvri_kcorrect: fltarr(5)-999.0}
    result = replicate(result,ngal)

    result.galaxy = ediscs.galaxy
    result.z = ediscs.z_abs
;   result.z = ediscs.z_obj
    
; now do it    
    if keyword_set(write) then begin
       psfile = datapath+'ediscs_synthmags_'+version+'.ps'
       im_plotconfig, 8, pos, psfile=psfile
    endif else im_window, 0, xratio=0.7, yratio=0.6

    plotscale = 1D18

    t0 = systime(1)
;   for igal = 1654, ngal-1L do begin
    for igal = 0L, ngal-1L do begin

       print, format='("Synthesizing magnitudes for galaxy ",I0,"/",I0,".",A10,$)', $
         igal+1, ngal, string(13b)

       zobj = result[igal].z

; restore the original spectra and best-fitting models
       good = where((specfit[*,0,igal] ne 0.0),npix)
       gal_restwave = specfit[good,0,igal]
       gal_restflux = specfit[good,1,igal]
       gal_restcflux = specfit[good,2,igal] ; continuum
       gal_restlflux = specfit[good,3,igal] ; emission lines
       gal_restivar = specfit[good,4,igal]
       gal_wave = gal_restwave*(1.0+zobj)
       gal_flux = gal_restflux/(1.0+zobj)
       gal_ivar = gal_restivar*(1.0+zobj)^2.0

; interpolate over highly discrepant pixels
       djs_iterstat, gal_restflux-(gal_restcflux+gal_restlflux), sigrej=3.0, mask=rejmask
       gal_restflux = djs_maskinterp(gal_restflux,(rejmask eq 0),gal_restwave,/const)
       
; interpolate the data onto the wavelength grid of the BC03 model,
; extrapolating redward and blueward
       restwave = model[igal].restwave
       linterp, gal_restwave, gal_restflux, restwave, restflux, missing=-1.0

       nodata = where((restflux eq -1.0),nno)
       if (nno ne 0L) then restflux[nodata] = model[igal].restflux[nodata]

       wave = restwave*(1.0+zobj)
       flux = restflux/(1.0+zobj)
       
; synthesize observed magnitudes both with and without the data
       obsmaggies = k_project_filters(k_lambda_to_edges(wave),$
         flux,filterlist=obsfilters,/silent)
       obsmaggies_model = k_project_filters(k_lambda_to_edges(model[igal].wave),$
         model[igal].flux,filterlist=obsfilters,/silent)
       obsmaggies_ivar = 1.0/(obsmaggies/20.0)^2.0
       result[igal].synth_obsmaggies = obsmaggies

; if OBSMAGGIES is negative then it means the spectrum is of too low
; S/N to be reliable; replace those magnitudes with the model
; magnitudes
       neg = where(obsmaggies lt 0.0,nneg)
       if (nneg ne 0) then obsmaggies[neg] = obsmaggies_model[neg]

; compute K-corrections and rest-frame quantities; note that
; BESTMAGGIES is identical to OBSMAGGIES, by construction; there is
; one object with S/N=0 in the continuum for which these K-corrections
; are undefined
       if (ediscs[igal].continuum_snr gt 0.0) then begin
          kcorr = im_simple_kcorrect(zobj,obsmaggies,obsmaggies_ivar,$
            obsfilters,ubvri_filters,restwave,restflux,absmag=absmag)
          kcorr_model = im_simple_kcorrect(zobj,obsmaggies,obsmaggies_ivar,$
            obsfilters,ubvri_filters,model[igal].restwave,$
            model[igal].restflux,absmag=absmag_model)

          result[igal].synth_ubvri_absmag = absmag
          result[igal].synth_ubvri_kcorrect = kcorr
          result[igal].synth_model_ubvri_absmag = absmag_model
          result[igal].synth_model_ubvri_kcorrect = kcorr_model
       endif
          
; QAplot       
       if keyword_set(debug) or keyword_set(write) then begin
          stats = im_stats(gal_flux,sigrej=2.5)
          xrange = [3500,9200]
;         xrange = minmax(gal_wave)+1000.0*[-1,1]
;         yrange = [min(gal_flux),stats.mean_rej+5.0*stats.sigma_rej]*plotscale
;         yrange = [stats.mean_rej+5.0*stats.sigma_rej]*plotscale
          yrange = [-0.15,1.2]*stats.maxrej

          plot, [0], [0], /nodata, xrange=xrange, yrange=plotscale*yrange, $
            xsty=3, ysty=3, charsize=1.8, xtitle=textoidl('Observed Wavelength (\AA)'), $
            ytitle=textoidl('Flux (10^{-18} '+flam_units()+')'), $
            title=repstr(strtrim(ancillary[igal].specfile,2),'.fits','')+$
            ' (z = '+strtrim(string(zobj,format='(F12.4)'),2)+')'
          djs_oplot, gal_wave, plotscale*gal_flux, color='blue'
          djs_oplot, wave, plotscale*model[igal].flux ;, color='dark green'
          for ifilt = 0L, nobsfilter-1L do djs_oplot, filt_lam[0L:filt_nlam[ifilt]-1L,ifilt], $
            plotscale*0.85*yrange[1]*filt_resp[0L:filt_nlam[ifilt]-1L,ifilt]/max(filt_resp[0L:filt_nlam[ifilt]-1L,ifilt]), $
            line=1, color='red'
;         djs_oplot, !x.crange, [0,0], line=0
          if keyword_set(debug) then cc = get_kbrd(1)
       endif
    endfor
    splog, 'Total time = '+string((systime(1)-t0)/60.0,format='(G0.0)')+' minutes.'

; write out the results
    if keyword_set(write) then begin
       im_plotconfig, psfile=psfile, /psclose

       outfile = 'ediscs_synthmags_'+version+'.fits'
       im_mwrfits, result, outfile
    endif

return
end
    
