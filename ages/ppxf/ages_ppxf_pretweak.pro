pro ages_ppxf_pretweak, pass1, firstpass=firstpass, lastpass=lastpass, $
  test=test, doplot=doplot
; jm09nov13ucsd - see the README
    
    light = 2.99792458D5        ; speed of light [km/s]

; path names and emission-line file name, for masking
    version = ages_version(/ppxf_specfit)
    zabsvdisppath = ages_path(/spec1d)+'fluxed/zabs_vdisp/'+version+'/'
    specfitpath = ages_path(/spec1d)+'fluxed/pretweak/'+version+'/' 
    linefile = ages_path(/ppxf)+'gandalf_elinelist_all.dat'

; figure out which passes we are going to fit
    if (n_elements(pass1) eq 0) then $
      allpass = ages_allpasses(/fluxed) else $
      allpass = string(pass1,format='(I3.3)')
    npass = n_elements(allpass)

    if (n_elements(firstpass) eq 0) then firstpass = 0
    if (n_elements(lastpass) eq 0) then lastpass = npass-1

; read the templates (see BUILD_AGES_PPXF_TEMPLATES); resample and
; convolve to AGES pixel size and instrumental resolution
    velscale = ages_ppxf_velscale()
    inst_vdisp = ages_ppxf_instvdisp()

    tempflux = read_ages_ppxf_templates(tempwave,ntemp=ntemp,$
      velscale=velscale,inst_vdisp=inst_vdisp,/solar)

; initial fitting parameters    
    ebv_guess = 0.1
    degree = -1 ; do not use the additive polynomials
    
; initialize the output data structure
    nfinalpix = 5000
    result1 = {$
      ages_id:      0L, $
      ra:         0.0D, $
      dec:        0.0D, $
      pass:          0, $
      aper:          0, $
      z:           0.0, $
      zabs:        0.0, $
      zabs_raw:    0.0, $
      vdisp:       0.0, $
      vdisp_raw:   0.0, $
      chi2:        0.0, $ ; reduced chi^2
      ebv:         0.0, $
;     weights:     fltarr(ntemp), $
;     polyweights: fltarr(1+degree), $
      wave:        dblarr(nfinalpix), $ ; double!
      flux:        fltarr(nfinalpix), $
      ferr:        fltarr(nfinalpix), $
      redleak:     fltarr(nfinalpix), $
      bestfit:     fltarr(nfinalpix)}
    
; fit each plate separately
    t0 = systime(1)
    for ipass = firstpass, lastpass do begin

; read the 1D spectra and the output from AGES_GET_ZABS_VDISP 
       suffix = string(allpass[ipass],format='(I0)')       
       if keyword_set(test) then suffix = suffix+'_test'
       outfile = specfitpath+'ppxf_'+suffix+'.fits'

; we want to avoide spectrophotometric errors in the blue and the red
; leak in the, well, red, so mask appropriately       
       max_fitwave = 8400.0 ; observed frame
       case string(allpass[ipass],format='(I0)') of
          '104': min_fitwave = 5000.0
          '105': min_fitwave = 5000.0
          '113': min_fitwave = 5200.0
          '115': min_fitwave = 5000.0
          '201': min_fitwave = 5000.0
          '203': min_fitwave = 5500.0
          '204': min_fitwave = 4200.0
          '205': min_fitwave = 5200.0
          '206': min_fitwave = 4600.0
          '208': min_fitwave = 5700.0
          '212': min_fitwave = 4400.0
          '213': min_fitwave = 5100.0
          '214': min_fitwave = 4500.0
          '301': min_fitwave = 5400.0
          '302': min_fitwave = 4400.0
          '303': min_fitwave = 4700.0
          '304': min_fitwave = 5400.0
          '305': min_fitwave = 4300.0
          '306': min_fitwave = 4300.0
          '307': min_fitwave = 5400.0
          '308': min_fitwave = 5400.0
          '309': min_fitwave = 5300.0
          '314': min_fitwave = 4700.0
          '409': min_fitwave = 4300.0
          '410': min_fitwave = 4300.0
          '411': min_fitwave = 4500.0
          '422': min_fitwave = 4500.0
          else: min_fitwave = 4000.0
       endcase
       
; read the output from AGES_GET_ZABS_VDISP
       zabsvdisp = read_ages_zabs_vdisp(allpass[ipass],$
         zabsvdisppath=zabsvdisppath)
       
; read all the spectra on this plate and subscript to the subset of
; objects fitted in AGES_GET_ZABS_VDISP
       spec1d = read_ages_spec1d(allpass[ipass],/fix)
       index = where_array(zabsvdisp.ages_id,spec1d.ages_id)
       nobj = n_elements(index)
       if (nobj eq 0) then message, 'Problem here!'

; pack the structure       
       result = replicate(result1,nobj)
       result.ages_id = spec1d.ages_id[index]
       result.ra = spec1d.ra[index]
       result.dec = spec1d.dec[index]
       result.pass = spec1d.pass
       result.aper = spec1d.aper[index]
       result.z = spec1d.z[index]

; copy over the redshift and vdisp results
       result.zabs = zabsvdisp.zabs
       result.zabs_raw = zabsvdisp.zabs_raw
       result.vdisp = zabsvdisp.vdisp
       result.vdisp_raw = zabsvdisp.vdisp_raw
       
; fit each object using PPXF
       t0 = systime(1)
       if (keyword_set(doplot) eq 0) then $
         psfile = specfitpath+'qaplot_ppxf_'+suffix+'.ps'
       im_plotconfig, 8, pos, psfile=psfile, charsize=1.6, $
         thick=1, ymargin=[0.6,1.1]
;      for iobj = 50, 60 do begin
       for iobj = 0, nobj-1 do begin
          splog, allpass[ipass], iobj
;         print, format='("Fitting object ",I0,"/",I0,A10,$)', $
;           iobj, nobj-1, string(13b)

; see AGES_GET_ZABS_VDISP
          zabs = result[iobj].zabs
          vdisp = result[iobj].vdisp
          vsys = alog(zabs+1.0D)*light ; systemic velocity
          
; interpolate the galaxy spectrum logarithmically in wavelength
          good = where((spec1d.wave ge spec1d.minwave[index[iobj]]) and $
            (spec1d.wave le spec1d.maxwave[index[iobj]]),ngood)
          wave = spec1d.wave[good]
          flux = spec1d.flux[good,index[iobj]]
          ferr = spec1d.ferr[good,index[iobj]]

          log_rebin, minmax(wave), flux, lnflux, lnwave, velscale=velscale
          log_rebin, minmax(wave), ferr^2, lnvar, velscale=velscale
          lnferr = sqrt(abs(lnvar))

          lnrestwave = lnwave - alog(1.0+zabs)
          lnrestflux = lnflux*(1.0+zabs)
          lnrestferr = lnferr*(1.0+zabs)
          npix = n_elements(lnrestwave)

; compute the *rest-frame* velocity offset between the templates and
; the data and restrict the wavelength range of the templates
          keep = where((exp(tempwave) gt exp(min(lnrestwave))-10.0) and $
            (exp(tempwave) lt exp(max(lnrestwave))+10.0),nkeep)
          if (nkeep eq 0) then message, 'Problem here!'
          fit_tempwave = tempwave[keep]
          fit_tempflux = tempflux[keep,*]
          voffset = -(min(lnrestwave)-min(fit_tempwave))*light
          start = [vsys*0.0+voffset,vdisp] ; initial/final solution

; mask emission lines and other crummy regions
          linepars = read_gandalf_elinelist(linefile)
          goodpixels = gandalf_mask_emission_lines(npix,vsys,$
            linepars,velscale,alog(wave[0]),velscale/light,$
            l_rf_range=[min_fitwave,max_fitwave]/(zabs+1.0D),$
            sigma=500.0)
;         djs_plot, exp(restwave), restflux, xsty=3, ysty=3
;         djs_oplot, exp(restwave[goodpixels]), restflux[goodpixels], psym=4, color='red'

; call PPXF! fix the velocity offset and velocity dispersion by
; setting MOMENTS=0, and just solve for the continuum, with reddening,
; but *no* additive polynomial freedom
          ebv = ebv_guess
          ages_ppxf, fit_tempflux, lnrestflux, lnrestferr, velscale, start, $
            sol, goodpixels=goodpixels, plot=doplot, moments=0, $
            degree=degree, weights=weights, bestfit=bestfit, quiet=1, $
            clean=1, reddening=ebv, lambda=exp(lnrestwave)

          result[iobj].chi2 = sol[6]
          result[iobj].ebv = ebv
;         result[iobj].weights = weights
;         result[iobj].polyweights = polyweights
          
; smooth the residuals heavily so that we can correct for the red leak
          smooth1 = medsmooth(lnrestflux-bestfit,151)
          redleak = smooth(smooth1,51,/edge_truncate)

; store the final spectrum in the *observed* frame (ln-binned)
          result[iobj].wave[0:npix-1] = double(lnrestwave + alog(zabs+1.0D))
          result[iobj].flux[0:npix-1] = float(lnrestflux/(zabs+1.0D))
          result[iobj].ferr[0:npix-1] = float(lnrestferr/(zabs+1.0D))
          result[iobj].bestfit[0:npix-1] = float(bestfit/(zabs+1.0D))
          result[iobj].redleak[0:npix-1] = float(redleak/(zabs+1.0D))

; make a QAplot
          stats = im_stats(lnrestflux,sigrej=5.0)
          yrange = 1E17*[-0.15,1.2]*stats.maxrej
          xrange1 = minmax(exp(lnrestwave))
          djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=3, ysty=3, $
            xrange=xrange1, yrange=yrange, xtitle='Rest Wavelength (\AA)', $
            ytitle='Flux (10^{-17} '+flam_units()+')', title='Pass '+$
            repstr(suffix,'_',' ')
          djs_oplot, exp(lnrestwave), 1E17*(lnrestflux-bestfit), psym=10, color='light grey'
          djs_oplot, exp(lnrestwave), 1E17*redleak, psym=10, color='orange'
          djs_oplot, exp(lnrestwave), 1E17*lnrestflux, psym=10
          djs_oplot, exp(lnrestwave), 1E17*bestfit, psym=10, color='red'
          djs_oplot, exp(lnrestwave[goodpixels]), 1E17*bestfit[goodpixels], psym=6, $
            symsize=0.1, color='blue'
          djs_oplot, !x.crange, [0,0], line=0, thick=1

;         test = restore_ages_ppxf_bestfit(result[iobj].weights,$
;           ebv=result[iobj].ebv,bestwave=testwave,velscale=velscale)
;         djs_oplot, exp(testwave), 1E17*test, color='yellow'

          im_legend, [$
            'Fiber '+string(result[iobj].aper,format='(I3.3)'),$
            'z_{AGES}='+strtrim(string(result[iobj].z,format='(F12.5)'),2),$
            'z_{PPXF}='+strtrim(string(result[iobj].zabs,format='(F12.5)'),2),$
            '\sigma='+strtrim(string(result[iobj].vdisp,format='(F12.1)'),2)+' km s^{-1}', $
            'E(B-V)='+strtrim(string(result[iobj].ebv,format='(F12.3)'),2),$
            '\chi^{2}_{\nu}='+strtrim(string(result[iobj].chi2,format='(F12.3)'),2)], $
            /left, /top, box=0, charsize=1.3, margin=0
       endfor
       if (keyword_set(doplot) eq 0) then $
         im_plotconfig, psfile=psfile, /gzip, /psclose
       splog, 'Total time = ', (systime(1)-t0)/60.0, ' minutes'
       im_mwrfits, result, outfile, /clobber
       
    endfor 

return
end
    
    
