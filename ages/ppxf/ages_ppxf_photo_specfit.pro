pro ages_gandalf_specfit, pass1, firstpass=firstpass, lastpass=lastpass, $
  test=test, doplot=doplot, index=index1
; jm09nov13ucsd - in this first pass     

INCOMPLETE CODE - THIS CODE FITS PHOTOMETRY!
    
    common ages_catalogs, allcodes, allphot
    
    light1 = 2.99792458D5        ; speed of light1 [km/s]
    light2 = 2.99792458D18       ; speed of light1 [A/s]

; path names and emission-line file name, for masking    
    version = ages_version(/ppxf_specfit)
    base_specfitpath = ages_path(/ppxf)
    spec1dpath = base_specfitpath+'fluxed/pretweak/'+version+'/' 
    specfitpath = base_specfitpath+'fluxed/tweak/'+version+'/' 
    linefile = base_specfitpath+'gandalf_elinelist_'+version+'.dat'

; figure out which passes we are going to fit
    if (n_elements(pass1) eq 0) then allpass = ages_allpasses() else $
      allpass = string(pass1,format='(I3.3)')
    npass = n_elements(allpass)

    if (n_elements(firstpass) eq 0) then firstpass = 0
    if (n_elements(lastpass) eq 0) then lastpass = npass-1

; read the templates (see BUILD_AGES_PPXF_TEMPLATES); resample and
; convolve to AGES pixel size and instrumental resolution
    velscale = ages_ppxf_velscale()
    inst_vdisp = ages_ppxf_instvdisp()

    tempflux = read_ages_ppxf_templates(tempwave,ntemp=ntemp,$
      velscale=velscale,inst_vdisp=inst_vdisp,linear_tempwave=linear_tempwave,$
      linear_tempflux=linear_tempflux)

; initial fitting parameters    
    vdisp_guess = 100.0D
    ebv_guess = 0.1
    degree = 3 ; additive Legendre polynomials
    
; initialize the output data structure
    result1 = {$
      ages_id:   0L, $
      ra:      0.0D, $
      dec:     0.0D, $
      pass:       0, $
      aper:       0, $
      z:        0.0, $
      zabs:     0.0, $
      vdisp:    0.0, $
      chi2:     0.0, $ ; reduced chi^2
      ebv:      0.0, $
      weights:     fltarr(ntemp), $
      polyweights: fltarr(degree+1)}
    
; read the AGES codes and photometric catalogs
    if (n_elements(allcodes) eq 0) then allcodes = mrdfits($
      ages_path(/analysis)+'catalog.codes.fits.gz',1)
    if (n_elements(allphot) eq 0) then allphot = mrdfits($
      ages_path(/analysis)+'ages_photometry_'+$
      ages_version(/phot)+'.fits.gz',1)

; convert the photometry to maggies and compute modelmaggies 
    ages_to_maggies, allphot, allmaggies, allivarmaggies, $
      filterlist=filterlist
    allerrmaggies = allivarmaggies*0.0+1E16
    notzero = where(allivarmaggies gt 0.0)
    allerrmaggies[notzero] = 1.0/sqrt(allivarmaggies[notzero])
    weff = k_lambda_eff(filterlist=filterlist)
    refband = 3 ; reference bandpass (R)

    modelmaggies = fltarr(n_elements(weff),ntemp)
    for ii = 0, ntemp-1 do modelmaggies[*,ii] = k_project_filters($
      k_lambda_to_edges(linear_tempwave),linear_tempflux[*,ii],$
      filterlist=filterlist,/silent)

; fit each plate separately
    t0 = systime(1)
    for ipass = firstpass, lastpass do begin

; read the output from AGES_PPXF_PRETWEAK
       suffix = string(allpass[ipass],format='(I0)')       
       thisfile = spec1dpath+'ages_ppxf_'+suffix+'.fits.gz'
       splog, 'Reading '+thisfile
       spec1d = mrdfits(thisfile,1)
       if keyword_set(test) then suffix = suffix+'_test'
       
       codes = allcodes[spec1d.ages_id]
       if (n_elements(index1) ne 0) then index = index1 else begin
;         index = lindgen(n_elements(spec1d))
          index = where((codes.qshort eq 0) and $
            (spec1d.z gt 0.0) and (spec1d.z lt 1.0))
       endelse 
       nobj = n_elements(index)
       if (nobj eq 0) then message, 'Problem here!'

; photometry
       thesemaggies = allmaggies[*,spec1d.ages_id]
       theseerrmaggies = allerrmaggies[*,spec1d.ages_id]
       maggies = thesemaggies[*,index]
       errmaggies = theseerrmaggies[*,index]

; choose the spectrophotometric tweak curve
       tweak1 = ages_choose_specphot_tweak(allpass[ipass],$
         tweakwave=tweakwave1)
       
       outfile = specfitpath+'ages_gandalf_'+suffix+'.fits'
       result = replicate(result1,nobj)
       result.ages_id = spec1d[index].ages_id
       result.ra = spec1d[index].ra
       result.dec = spec1d[index].dec
       result.pass = spec1d[index].pass
       result.aper = spec1d[index].aper
       result.z = spec1d[index].z
       
; fit each object using GANDALF/PPXF
       t0 = systime(1)
       if (keyword_set(doplot) eq 0) then $
         psfile = specfitpath+'qaplot_ages_gandalf_'+suffix+'.ps'
       im_plotconfig, 6, pos, psfile=psfile, yspace=0.8, $
         charsize=1.6
       for iobj = 0, nobj-1 do begin
;      for iobj = 0, nobj-1 do begin
          splog, format='("Fitting object ",I0,"/",I0)', iobj+1, nobj

          zabs = spec1d[index[iobj]].zabs
          vdisp = spec1d[index[iobj]].vdisp
          vsys = alog(1+zabs)*light1 ; systemic velocity

; the spectra written out by AGES_PPXF_PRETWEAK are logarithmically
; binned in wavelength and in the *observed* frame
          notzero = where((spec1d[index[iobj]].wave gt 0.0),nnotzero)
          if (nnotzero eq 0) then message, 'Problem here!'
          wave = spec1d[index[iobj]].wave[notzero]
          redleak = spec1d[index[iobj]].redleak[notzero]
          flux1 = spec1d[index[iobj]].flux[notzero]
          ferr1 = spec1d[index[iobj]].ferr[notzero]
          npix = n_elements(wave)

; correct for the red leak and tweak the fluxing
;         djs_plot, exp(wave), flux1, xsty=3, ysty=3
          isred = where((exp(wave) ge 8500.0),nisred)
          if (nisred ne 0) then flux1[isred] = flux1[isred] - redleak[isred]
;         djs_oplot, exp(wave), flux, color='red'

          linterp, tweakwave1, tweak1, exp(wave), tweak, missing=1.0
          flux = flux1*tweak
          ferr = ferr1*tweak
;         djs_oplot, exp(wave), flux1, color='cyan'

; now convert the data and the photometry to the rest-frame
          restwave = wave - alog(1.0+zabs)
          restflux = flux*(1.0+zabs)
          restferr = ferr*(1.0+zabs)

          restweff = weff/(1.0+zabs)
          restmaggies = maggies[*,iobj]*(1.0+zabs)
          resterrmaggies = errmaggies[*,iobj]*(1.0+zabs)

          scale = (k_project_filters(exp(restwave),restflux,$
            filterlist=filterlist[refband]))[0]/restmaggies[refband]
          restmaggies = restmaggies*scale
          resterrmaggies = resterrmaggies*scale
          
          photo = {weff: restweff, maggies: restmaggies, $
            errmaggies: resterrmaggies, modelmaggies: modelmaggies}

; compute the *rest-frame* velocity offset between the templates and
; the data and restrict the wavelength range of the templates
;         keep = where((exp(tempwave) gt min(exp(wave)/(1+zabs))-10.0) and $
;           (exp(tempwave) lt max(exp(wave)/(1+zabs))+10.0),nkeep)
          keep = where((exp(tempwave) gt min(exp(restwave))-10.0) and $
            (exp(tempwave) lt max(exp(restwave))+10.0),nkeep)
          if (nkeep eq 0) then message, 'Problem here!'
          fit_tempwave = tempwave[keep]
          fit_tempflux = tempflux[keep,*]
;         voffset = -(min(wave)-min(fit_tempwave))*light1
          voffset = -(min(restwave)-min(fit_tempwave))*light1
;         kinematics = [voffset,vdisp,0.0,0.0,0.0,0.0]
          start  = [alog(1.0+zabs*0.0)*light1+voffset,vdisp]

; mask all the emission lines and then call PPFX
          linepars = read_gandalf_elinelist(linefile)
          goodpixels = gandalf_mask_emission_lines(npix,0.0,$
            linepars,velscale,restwave[0],velscale/light1,$
            sigma=300.0)

          ebv = ebv_guess*0.0 ; NOTE!!
          im_ppxf, fit_tempflux, restflux, restferr, velscale, start, sol,$
            goodpixels=goodpixels, plot=doplot, moments=2, $
            degree=-1, bestfit=bestfit, weights=weights, $
;           reddening=ebv, lambda=exp(restwave), $
            clean=0, quiet=1, polyweights=polyweights, $
            photo=photo, outphoto=outphoto
          zabs = exp((vsys+sol[0]-voffset)/light1)-1.0 ; new absorption-line redshift

          splog, outphoto.photo_chi2
          niceprint, restmaggies/outphoto.bestmaggies

;         dfpsclose
;         djs_plot, exp(restwave), bestfit
;         jj = restore_ages_ppxf_bestfit(weights,bestwave=ww,velscale=velscale)
;         djs_oplot, exp(ww), jj, color='red'
          
;; finally read the emission-line parameter file and use it to mask out
;; crummy regions
;          alllinepars = read_gandalf_elinelist(linefile,$
;            fitlines=linepars,/actionfit)
;          goodpixels = gandalf_mask_emission_lines(npix,vsys,$
;            alllinepars,velscale,wave[0],velscale/light1,$
;            sigma=300.0)
;          
;; we're finally ready to call GANDALF!
;          sol = kinematics
;          gandalf, fit_tempflux, flux, ferr, velscale, sol, linepars, $
;            wave[0], velscale/light1, goodpixels=goodpixels, $
;            int_disp=inst_vdisp/3.0, bestfit=bestfit, $
;            emission_templates=emission_templates, $
;            weights=weights, /plot, l0_templ=fit_tempwave[0], $
;            for_errors=for_errors, error=esol, reddening=0.0

          result[iobj].zabs = zabs
          result[iobj].vdisp = sol[1]
          result[iobj].chi2 = sol[6]
          result[iobj].weights = weights
;         result[iobj].polyweights = polyweights
          result[iobj].ebv = ebv

; shift the observed spectrum to the revised absorption-line redshift

;         restweff = weff/(1.0+zabs)
;         restmaggies = maggies*(1.0+zabs)*scale
;         resterrmaggies = errmaggies*(1.0+zabs)*scale
          
; make a QAplot
          stats = im_stats(restflux,sigrej=5.0)
          yrange = 1E17*[-0.15,1.2]*stats.maxrej
          xrange1 = minmax(exp(restwave))
          djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=3, ysty=3, $
            xrange=xrange1, yrange=yrange, xtitle='Rest Wavelength (\AA)', $
            ytitle='Flux (10^{-17} '+flam_units()+')', title='Pass '+$
            repstr(suffix,'_',' ')
;         djs_oplot, exp(restwave), 1E17*(restflux-bestfit), psym=10, color='light grey'
          djs_oplot, exp(restwave), 1E17*restflux, psym=10
          djs_oplot, exp(restwave), 1E17*bestfit, psym=10, color='red'
          djs_oplot, exp(restwave[goodpixels]), 1E17*bestfit[goodpixels], psym=6, $
            symsize=0.1, color='blue'
          djs_oplot, !x.crange, [0,0], line=0, thick=1

          im_legend, [$
            'Fiber '+string(result[iobj].aper,format='(I3.3)'),$
            'z_{AGES}='+strtrim(string(result[iobj].z,format='(F12.5)'),2),$
            'z_{PPXF}='+strtrim(string(result[iobj].zabs,format='(F12.5)'),2),$
            '\sigma='+strtrim(string(result[iobj].vdisp,format='(F12.1)'),2)+' km s^{-1}', $
            'E(B-V)='+strtrim(string(result[iobj].ebv,format='(F12.3)'),2),$
            '\chi^{2}_{\nu}='+strtrim(string(result[iobj].chi2,format='(F12.3)'),2)], $
            /left, /top, box=0, charsize=1.3, margin=0

; photo plot
          bestflux1 = restore_ages_ppxf_bestfit(weights,$
            ebv=ebv,bestwave=bestwave1)

;         plot, bestwave1, bestflux1, xrange=[1000,3E4], /xlog

          xrange1 = [min(restweff)-100.0,max(restweff)-100.0]
          inrange = where((bestwave1 gt xrange1[0]) and $
            (bestwave1 lt xrange1[1]))
          bestwave = bestwave1[inrange]
          bestflux = bestflux1[inrange]*bestwave1[inrange]^2.0/light2
          bestflux = -2.5*alog10(bestflux)-48.6

          bestmab = -2.5*alog10(outphoto.bestmaggies)
          good = where(restmaggies gt 0.0)
          mab = -2.5*alog10(restmaggies[good])

          yrange1 = [max(mab)>max(bestmab),min(mab)<min(bestmab)]
          yrange1 = yrange1+[+0.3,-0.3]
;         yrange1 = [max(bestflux),min(bestflux)]
          djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=3, ysty=3, $
            xrange=xrange1, yrange=yrange1, xtitle='Rest Wavelength (\AA)', $
            ytitle='m_{AB}', /xlog
          djs_oplot, bestwave, bestflux, psym=10, color='light grey'
          djs_oplot, restweff, bestmab, psym=symcat(6,thick=3), $
            color='blue', symsize=2.0
          djs_oplot, restweff[good], mab, psym=symcat(16), $
            color='red', symsize=2.0
          im_legend, '\chi^{2}_{\nu}='+strtrim(string(outphoto.photo_chi2,format='(F12.3)'),2), $
            /left, /top, box=0, charsize=1.3, margin=0
       endfor
       if (keyword_set(doplot) eq 0) then $
         im_plotconfig, psfile=psfile, /gzip, /psclose
       splog, 'Total time = ', (systime(1)-t0)/60.0, ' minutes'
       im_mwrfits, result, outfile, /clobber

stop       
       
    endfor 

return
end
    
    
