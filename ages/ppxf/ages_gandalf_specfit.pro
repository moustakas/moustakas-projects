;+
; NAME:
;   AGES_GANDALF_SPECFIT
;
; PURPOSE:
;   Fit the AGES spectra using PPXF/GANDALF.
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;   pass - reduce this pass number (default is to process all the
;     plates) 
;
; KEYWORD PARAMETERS: 
;   test - test output files
;   doplot - render the GANDALF QAplot
;
; OUTPUTS: 
;   See the code.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Nov 13, UCSD
;
; Copyright (C) 2009, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

function ages_fit_continuum, wave, flux, ferr, tempwave=tempwave, $
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

; compute the velocity offset between the templates and the data and
; restrict the wavelength range of the templates
    keep = where((exp(tempwave) gt exp(min(restwave))-10.0) and $
      (exp(tempwave) lt exp(max(restwave))+10.0),nkeep)
    if (nkeep eq 0) then message, 'Problem here!'
    fit_tempwave = tempwave[keep]
    fit_tempflux = tempflux[keep,*]

    vsys = alog(zabs+1.0D)*light ; systemic velocity
    voffset = -(min(restwave)-min(fit_tempwave))*light
    start = [vsys*0.0+voffset,vdisp]

; mask all emission lines, then fit for the continuum, the reddening,
; and the polynomial templates, fixing the velocity dispersion and
; absorption-line redshift (i.e., MOMENTS=0); also remove pixels with
; very large errors
    linepars = read_gandalf_elinelist(alllinefile)
    goodpixels = gandalf_mask_emission_lines(npix,vsys,$
      linepars,velscale,wave[0],velscale/light,$
      sigma=500.0)
    goodpixels = cmset_op(goodpixels,'and',where(restferr lt 1E5))

    ebv = ebv_guess
    ages_ppxf, fit_tempflux, restflux, restferr, velscale, start, $
      sol, goodpixels=goodpixels, plot=doplot, moments=0, degree=-1, $
      weights=continuum_coeff, bestfit=continuum, /clean, /quiet, $
      vmaxshift=vmaxshift, sigmamax=sigmamax, error=err, $
      lambda=exp(restwave), reddening=ebv
    chi2 = sol[6]

; special case of no continuum fitted    
    if (total(continuum_coeff,/double) eq 0.0) then ebv = 0.0 ; note!
    
    if keyword_set(debug) then begin
       djs_plot, exp(restwave), restflux, ps=10, xsty=3, ysty=3, xr=[3650,4400]
       djs_oplot, exp(restwave), continuum, ps=10, color='blue'
       djs_oplot, exp(restwave[goodpixels]), restflux[goodpixels], psym=6, color='red', symsize=0.2
;      djs_oplot, exp(restwave[goodlinepixels]), restflux[goodlinepixels], ps=6, color='blue', symsize=0.2
       splog, 'Waiting for keystroke...' & cc = get_kbrd(1)
    endif

return, continuum
end
    
function ages_smooth_residuals, wave, flux, continuum, zabs=zabs, $
  velscale=velscale, linefile=linefile, broad=broad, debug=debug
; internal support routine to smooth the residuals

    light = 2.99792458D5 ; speed of light [km/s]
    vsys = alog(zabs+1.0D)*light ; systemic velocity
    if keyword_set(broad) then begin
       smooth_sigma = 1500.0
       sigrej = 1.0
    endif else begin
       smooth_sigma = 500.0
       sigrej = 1.0
    endelse
    
    restwave = wave - alog(zabs+1.0D)
    restflux = flux*(zabs+1.0D)
    residuals = restflux - continuum
    npix = n_elements(wave)

; find all pixels that *could* be affected by emission lines    
    linepars = read_gandalf_elinelist(linefile)
    good1 = gandalf_mask_emission_lines(npix,vsys,linepars,$
      velscale,wave[0],velscale/light,sigma=smooth_sigma)
    bad1 = lindgen(npix)
    remove, good1, bad1

; within the affected regions, mask all deviant pixels
    msk = residuals*0.0+1.0 ; all are good
    msk[bad1] = 0.0
    msk[bad1] = abs(residuals[bad1]) lt sigrej*djsig(residuals[good1]) ; keep non-deviant pixels
    good = where(msk,ngood,comp=bad,ncomp=nbad)

; replace the line-affected pixels with noise to make the smoothing
; well-behaved, below
    smooth1 = residuals
    if (nbad ne 0) then smooth1[bad] = djs_median(residuals[good]) + $
      randomn(seed,nbad)*djsig(residuals[good])

    smooth2 = medsmooth(smooth1,151)
    smooth_continuum = smooth(smooth2,51,/edge)
    
    if keyword_set(debug) then begin
       djs_plot, exp(restwave), residuals, xsty=3, ysty=3, ps=10;, xr=[5500,6800] ; xr=[4700,5050]
       djs_oplot, exp(restwave[good]), residuals[good], ps=6, sym=0.2, color='orange'
       djs_oplot, exp(restwave[bad]), residuals[bad], ps=6, sym=0.2, color='blue'
;      djs_oplot, exp(restwave), smooth1, ps=10, color='green'
;      djs_oplot, exp(restwave), smooth2, ps=10, color='blue'
       djs_oplot, exp(restwave), smooth_continuum, ps=10, color='red'
       splog, 'Waiting for keystroke...' & cc = get_kbrd(1)
    endif

return, smooth_continuum
end

function ages_fit_lines, wave, flux, ferr, continuum, smooth_continuum, $
  zabs=zabs, vdisp=vdisp, velscale=velscale, linefile=linefile, $
  line_inst_vdisp=line_inst_vdisp, tempwave=tempwave, broad=broad, $
  linepars=linepars, esol=esol, echi2=echi2, etemplates=etemplates, $
  debug=debug, vmaxshift_iter1=vmaxshift_iter1, vmaxshift_iter2=vmaxshift_iter2, $
  sigmamax_iter1=sigmamax_iter1, sigmamax_iter2=sigmamax_iter2
; construct the pure emission-line spectrum          

    light = 2.99792458D5 ; speed of light [km/s]

    if keyword_set(broad) then begin
       if (n_elements(vmaxshift_iter1) eq 0) then vmaxshift_iter1 = 2000.0D
       if (n_elements(vmaxshift_iter2) eq 0) then vmaxshift_iter2 = 2000.0D
       if (n_elements(sigmamax_iter1) eq 0) then sigmamax_iter1 = 1D4 ; 2300D ; 1800D
       if (n_elements(sigmamax_iter2) eq 0) then sigmamax_iter2 = 1D4 ; 2300D ; 1800D
       sigma_smooth = 600.0
    endif else begin
       if (n_elements(vmaxshift_iter1) eq 0) then vmaxshift_iter1 = 500D
       if (n_elements(vmaxshift_iter2) eq 0) then vmaxshift_iter2 = 310.0D
       if (n_elements(sigmamax_iter1) eq 0) then sigmamax_iter1 = 520.0D
       if (n_elements(sigmamax_iter2) eq 0) then sigmamax_iter2 = 520.0D
       sigma_smooth = 400.0
    endelse

    restwave = wave - alog(zabs+1.0D)
    restflux = flux*(zabs+1.0D)
    restferr = ferr*(zabs+1.0D)
    npix = n_elements(wave)
  
    restlineflux = restflux - continuum - smooth_continuum 

    keep = where((exp(tempwave) gt exp(min(restwave))-10.0) and $
      (exp(tempwave) lt exp(max(restwave))+10.0),nkeep)
    if (nkeep eq 0) then message, 'Problem here!'
    fit_tempwave = tempwave[keep]
    
    vsys = alog(zabs+1.0D)*light ; systemic velocity
    voffset = -(min(restwave)-min(fit_tempwave))*light
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
    im_gandalf, fit_tempwave*0.0, restlineflux, restferr, velscale, sol, $
;   im_gandalf, reform(fit_tempwave,nkeep,1)*0.0, restlineflux, restferr, velscale, sol, $
      linepars, restwave[0], velscale/light, goodpixels=goodlinepixels, $
      int_disp=line_inst_vdisp, bestfit=bestlinefit, weights=lweights, $
      emission_templates=etemplates, error=esol, l0_templ=fit_tempwave[0], $
      /for_errors, plot=plot, quiet=1, chi2=echi2, vmaxshift_iter1=vmaxshift_iter1, $
      vmaxshift_iter2=vmaxshift_iter2, sigmamax_iter1=sigmamax_iter1, $
      sigmamax_iter2=sigmamax_iter2

; parse the output
    if keyword_set(debug) then begin
       djs_plot, exp(restwave), restlineflux, ps=10, xr=[3650,4100]
;      djs_plot, exp(restwave), restlineflux, ps=10, xr=[6400,6750];, xr=[3650,4100] ; xr=[6400,6750] ; 
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

pro ages_gandalf_specfit, pass1, debug=debug, broad=broad, solar=solar
    
    light = 2.99792458D5 ; speed of light [km/s]
    
; path names and emission-line file name
    version = ages_version(/ppxf_specfit)
    zabsvdisppath = ages_path(/spec1d)+'fluxed/zabs_vdisp/'+version+'/'
    spec1dpath = ages_path(/spec1d)+'fluxed/pretweak/'+version+'/' 
    specfitpath = ages_path(/spec1d)+'fluxed/tweak/'+version+'/' 
    alllinefile = ages_path(/ppxf)+'gandalf_elinelist_all.dat'
    linefile = ages_path(/ppxf)+'gandalf_elinelist_'+version+'.dat'
    if keyword_set(broad) then begin
       alllinefile = repstr(alllinefile,'elinelist_','elinelist_broad_')
       linefile = repstr(linefile,'elinelist_','elinelist_broad_')
    endif

; figure out which passes we are going to fit; exclude the unfluxed
; plates! 
    if (n_elements(pass1) eq 0) then begin
       if keyword_set(broad) then $
         allpass = ages_balmer_isbroad(/justpass) else $
           allpass = ages_allpasses(/fluxed) 
    endif else begin
       allpass = string(pass1,format='(I3.3)')
    endelse
    npass = n_elements(allpass)

; read the templates (see BUILD_AGES_PPXF_TEMPLATES); resample and
; convolve to AGES pixel size and instrumental resolution
    velscale = ages_ppxf_velscale()
    fluxscale = ages_ppxf_fluxscale()
    inst_vdisp = ages_ppxf_instvdisp()
    line_inst_vdisp = 5.0 ; inst_vdisp/3.0 ; [km/s]

    tempflux = read_ages_ppxf_templates(tempwave,ntemp=ntemp,$
      velscale=velscale,inst_vdisp=inst_vdisp,tempinfo=info,$
      solar=solar)
    metal = info.metallicity
    nmetal = n_elements(metal)

; read the emission-line file here to get the total number of lines
; that we are going to be fitting
    junk = read_gandalf_elinelist(linefile,nlines=nlines,$
      fitlines=linepars,/actionfit)
    
; initialize the output data structure
    nfinalpix = 5000
    specdata_template = {$
      ages_id:   0L, $
;     ra:      0.0D, $
;     dec:     0.0D, $
      pass:       0, $
      aper:       0, $
      z:        0.0, $
      zabs:     0.0, $
      zabs_raw:    0.0, $
      zabs_err:     0.0, $
      vdisp:        0.0, $
      vdisp_raw:    0.0, $
      vdisp_err:    0.0, $
      inst_vdisp: inst_vdisp, $
      line_inst_vdisp: line_inst_vdisp, $
      continuum_snr:   0.0, $
      continuum_coeff:     fltarr(ntemp), $
      continuum_ebv:       0.0, $
      continuum_Z:         0.0,$
      continuum_chi2:  0.0, $ ; reduced chi^2
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

; fit each plate separately
    t1 = systime(1)
    for ipass = 0, npass-1 do begin
       splog, 'Pass ', allpass[ipass]
       
; output file names
       suffix = string(allpass[ipass],format='(I0)')       
       if keyword_set(solar) then suffix = suffix+'_solar'
       specdatafile = specfitpath+'specdata_raw_'+suffix+'.fits'
       if keyword_set(broad) then specdatafile = repstr(specdatafile,'raw_','raw_broad_')

; read the outputs from AGES_PPXF_PRETWEAK and AGES_GET_ZABS_VDISP
       spec1dfile = spec1dpath+'ppxf_'+string(allpass[ipass],format='(I0)')+'.fits.gz'
       splog, 'Reading '+spec1dfile
       spec1d = mrdfits(spec1dfile,1)
       nobj = n_elements(spec1d)

       zabsvdisp = read_ages_zabs_vdisp(allpass[ipass],$
         zabsvdisppath=zabsvdisppath)
       if (nobj ne n_elements(zabsvdisp)) then message, 'Problem here'

; restrict to the set of broad-line AGN, if desired
       if keyword_set(broad) then begin
          isbroad = ages_balmer_isbroad(spec1d)
          if (isbroad[0] eq -1) then begin
             splog, 'No broad-line AGN in PASS '+$
               strtrim(allpass[ipass])
             continue
          endif else begin
             spec1d = spec1d[isbroad]
             zabsvdisp = zabsvdisp[isbroad]
             nobj = n_elements(spec1d)
          endelse
       endif
       
; choose the spectrophotometric tweak curve
       tweak1 = ages_choose_specphot_tweak(allpass[ipass],tweakwave=tweakwave1)

; fit each object using GANDALF/PPXF
       t0 = systime(1)
;      for iobj = 3, 3 do begin
;      for iobj = 164, 164 do begin
       for iobj = 0, nobj-1 do begin
          splog, strtrim(allpass[ipass],2)+': '+string(iobj+1,format='(I3.3)')+'/'+$
            string(nobj,format='(I3.3)')
;         print, format='("Fitting object ",I0,"/",I0,A10,$)', $
;           iobj, nobj-1, string(13b)

; the spectra written out by AGES_PPXF_PRETWEAK are logarithmically
; binned in wavelength and in the *observed* frame
          notzero = where((spec1d[iobj].wave gt 0.0),nnotzero)
          if (nnotzero eq 0) then message, 'Problem here!'
          wave = spec1d[iobj].wave[notzero]
          redleak = spec1d[iobj].redleak[notzero]
          flux1 = spec1d[iobj].flux[notzero]
          ferr1 = spec1d[iobj].ferr[notzero]
          npix = n_elements(wave)

; correct for the red leak and tweak the fluxing
          isred = where((exp(wave) ge 8500.0),nisred)
          if (nisred ne 0) then flux1[isred] = flux1[isred] - redleak[isred]

          linterp, tweakwave1, tweak1, exp(wave), tweak, missing=1.0
          flux = flux1*tweak/fluxscale
          ferr = ferr1*tweak/fluxscale

; see AGES_GET_ZABS_VDISP
          zabs1 = zabsvdisp[iobj].zabs ; zabsvdisp[iobj].z
          vdisp1 = zabsvdisp[iobj].vdisp

; loop on each metallicity       
          for imetal = 0, nmetal-1 do begin
             tempflux1 = tempflux[*,*,imetal]
          
; fit the full continuum once using the solar-metallicity templates to
; get a better estimate of the absorption-line redshift
             continuum1 = ages_fit_continuum(wave,flux,ferr,tempwave=tempwave, $
               tempflux=tempflux1,velscale=velscale,zabs=zabs1,vdisp=vdisp1,$
               ebv_guess=ebv_guess,alllinefile=alllinefile,goodpixels=goodpixels,$
               final_ebv=ebv1,chi2=chi21,continuum_coeff=coeff1,$
               debug=debug)
             
             if (imetal eq 0) then begin
                allzabs = zabs1
                allvdisp = vdisp1
                allebv = ebv1
                allchi2 = chi21
                allcoeff = coeff1
                continuum = continuum1
             endif else begin
                allzabs = [allzabs,zabs1]
                allvdisp = [allvdisp,vdisp1]
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
          ebv = allebv[zbest]
          continuum = continuum[*,zbest]
          continuum_coeff = allcoeff[*,zbest]
          continuum_Z = metal[zbest]
          
; smooth the residuals before before fitting the lines
          smooth_continuum = ages_smooth_residuals(wave,flux,$
            continuum,zabs=zabs,velscale=velscale,broad=broad,$
            linefile=linefile,debug=debug)

; finally fit the emission lines, including special cases
          delvarx, vmaxshift_iter1
          if (allpass[ipass] eq 506) and (zabsvdisp[iobj].aper eq 57) then vmaxshift_iter1 = 1000D
;         if (allpass[ipass] eq 304) and (zabsvdisp[iobj].aper eq 290) then begin
;            sigmamax_iter1 = 1000D
;         endif
          sol = ages_fit_lines(wave,flux,ferr,continuum,smooth_continuum,$
            zabs=zabs,vdisp=vdisp,velscale=velscale,linefile=linefile,$
            line_inst_vdisp=line_inst_vdisp,tempwave=tempwave,$
            broad=broad,linepars=linepars,esol=esol,etemplates=etemplates,$
            echi2=echi2,debug=debug,vmaxshift_iter1=vmaxshift_iter1)

; pack everything into a structure          
          specdata1 = specdata_template
          struct_assign, spec1d[iobj], specdata1, /nozero
          struct_assign, zabsvdisp[iobj], specdata1, /nozero

          restwave = wave - alog(zabs+1.0D)
          restflux = flux*(zabs+1.0D)
          restferr = ferr*(zabs+1.0D)
          specdata1.wave[0:npix-1] = double(restwave); + alog(1.0+zabs)
          specdata1.flux[0:npix-1] = float(restflux*fluxscale);/(1.0+zabs)
          specdata1.ferr[0:npix-1] = float(restferr*fluxscale);/(1.0+zabs)
          specdata1.continuum[0:npix-1] = float(continuum*fluxscale);/(1.0+zabs)
          specdata1.smooth_continuum[0:npix-1] = float(smooth_continuum*fluxscale);/(1.0+zabs)

; flag the indices of the lines fitted
          fitlines = where((linepars.action ne 'i') and $
            (linepars.kind eq 'l'),nfitlines)
          if (nfitlines ne 0) then begin
             specdata1.fitlines[0:nfitlines-1] = fitlines
             specdata1.sol[*,0:nfitlines-1] = sol
             specdata1.esol[*,0:nfitlines-1] = esol
             specdata1.etemplates[0:npix-1,0:nfitlines-1] = float(etemplates*fluxscale)
          endif

          specdata1.continuum_snr = zabsvdisp[iobj].snr
          specdata1.continuum_Z = continuum_Z
          specdata1.continuum_chi2 = continuum_chi2
          specdata1.continuum_coeff = continuum_coeff[0:ntemp-1]*fluxscale
          specdata1.continuum_ebv = ebv
          specdata1.emission_line_chi2 = echi2

          if (n_elements(specdata) eq 0) then specdata = specdata1 else $
            specdata = [temporary(specdata),specdata1]
       endfor 
       im_mwrfits, specdata, specdatafile, /clobber
       splog, 'Total time = ', (systime(1)-t0)/60.0, ' minutes'
       delvarx, specdata ; delete memory!!
    endfor 
    splog, 'Total time for all plates = ', (systime(1)-t1)/60.0, ' minutes'

return
end
