;+
; NAME:
;   DEEP2_GANDALF_SPECFIT_DR4
;
; PURPOSE:
;   Fit the DEEP2 spectra using PPXF/GANDALF.
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
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
;   J. Moustakas, 2013 Dec 20 - based on AGES_GANDALF_SPECFIT_DR4
;
; Copyright (C) 2013, John Moustakas
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

function deep2_fit_unfluxed_continuum, wave, flux, ferr, zabs=zabs, $
  velscale=velscale, linefile=linefile, broad=broad, debug=debug

    light = 2.99792458D5 ; speed of light [km/s]
    vsys = alog(zabs+1.0D)*light ; systemic velocity
    smooth_sigma = 500.0
    sigrej = 1.0
    
    restwave = wave - alog(zabs+1.0D)
    restflux = flux*(zabs+1.0D)
    restferr = ferr*(zabs+1.0D)
    npix = n_elements(wave)

; find all pixels that *could* be affected by emission lines; within
; the affected regions, mask all deviant pixels 
    msk = restflux*0.0+1.0 ; all are good

    linepars = read_gandalf_elinelist(linefile)
    good1 = gandalf_mask_emission_lines(npix,vsys,linepars,$
      velscale,wave[0],velscale/light,sigma=smooth_sigma)
    if n_elements(good1) lt npix then begin
       bad1 = lindgen(npix)
       remove, good1, bad1

       msk[bad1] = 0.0
       msk[bad1] = abs(restflux[bad1]) lt sigrej*djsig(restflux[good1]) ; keep non-deviant pixels
       good = where(msk,ngood,comp=bad,ncomp=nbad)
    endif else begin
       ngood = npix
       good = lindgen(ngood)
       bad = -1
       nbad = 0
    endelse

; iterate twice, on the second iteration replacing the line-affected
; pixels with noise to make the fitting better-behaved
    tempflux = restflux
    invvar = msk/restferr^2.0
    for ii = 0, 1 do begin
       medwidth = 0.05
       sset = bspline_iterfit(restwave,tempflux,invvar=invvar,$
         bkspace=medwidth,nord=4.0,lower=1.5,upper=1.5,/silent,$
         outmask=outmsk)
       continuum = bspline_valu(restwave,sset)

       if (ii eq 0) then if (nbad ne 0) then tempflux[bad] = $
         interpol(continuum,restwave,restwave[bad]) + $
         randomn(seed,nbad)*djsig(restflux[good]-continuum[good])
       invvar = 1.0/restferr^2.0
       
;      djs_plot, restwave, restflux
;      djs_oplot, restwave, continuum, color='blue'
;      djs_oplot, restwave, tempflux, color='red'
    endfor

; simple smoothing    
;   smooth2 = medsmooth(smooth1,151)
;   continuum = smooth(smooth2,51,/edge)
    
    if keyword_set(debug) then begin
       rej = where(outmsk eq 0)
       djs_plot, exp(restwave), restflux, xsty=3, ysty=3, ps=10
;      djs_plot, exp(restwave), restflux, xsty=3, ysty=3, ps=10, xr=[4800,5050]
;      djs_plot, exp(restwave), restflux, xsty=3, ysty=3, ps=10;, xr=[5500,6800] ; xr=[4700,5050]
       djs_oplot, exp(restwave[rej]), restflux[rej], ps=6, sym=0.2, color='blue'
       djs_oplot, exp(restwave), continuum, ps=10, color='red'
       splog, 'Waiting for keystroke...' & cc = get_kbrd(1)
    endif

return, continuum
end

function deep2_fit_unfluxed_lines, wave, flux, ferr, continuum, $
  zabs=zabs, velscale=velscale, linefile=linefile, broad=broad, $
  line_inst_vdisp=line_inst_vdisp, linepars=linepars, esol=esol, $
  echi2=echi2, etemplates=etemplates, debug=debug, $
  vmaxshift_iter1=vmaxshift_iter1
; construct the pure emission-line spectrum          

    light = 2.99792458D5 ; speed of light [km/s]

    if keyword_set(broad) then begin
       vmaxshift_iter1 = 2000.0D
       vmaxshift_iter2 = 2000.0D
       sigmamax_iter1 = 1D4 ; 2300D ; 1800D
       sigmamax_iter2 = 1D4 ; 2300D ; 1800D
       sigma_smooth = 600.0
    endif else begin
       if (n_elements(vmaxshift_iter1) eq 0) then vmaxshift_iter1 = 500D
       vmaxshift_iter2 = 310.0D
       sigmamax_iter1 = 520.0D
       sigmamax_iter2 = 520.0D
       sigma_smooth = 300.0
    endelse

    restwave = wave - alog(zabs+1.0D)
    restflux = flux*(zabs+1.0D)
    restferr = ferr*(zabs+1.0D)
    npix = n_elements(wave)
  
    restlineflux = restflux - continuum

    vsys = alog(zabs+1.0D)*light ; systemic velocity
    kinematics  = [0.0,150.0,0.0,0.0,0.0,0.0]    ; for IM_GANDALF
    
; now read the emission-line parameter file of lines we *want* to fit 
    alllinepars = read_gandalf_elinelist(linefile,$
      fitlines=linepars,/actionfit)
    junk = gandalf_mask_emission_lines(npix,vsys,linepars,$
      velscale,wave[0],velscale/light,sigma=sigma_smooth)
    goodpixels = gandalf_mask_emission_lines(npix,vsys,alllinepars,$
      velscale,wave[0],velscale/light,sigma=sigma_smooth)
    
    if n_elements(goodpixels) lt npix then begin
       goodlinepixels = lindgen(npix)
       remove, goodpixels, goodlinepixels
       goodlinepixels = cmset_op(goodlinepixels,'and',$
         where(restferr lt 1E5))
       if goodlinepixels[0] eq -1 then goodlinepixels = where(restferr lt 1E5)
    endif else goodlinepixels = where(restferr lt 1E5)

; now call GANDALF! for the continuum templates pass an array of
; zeros, and just fit for the emission lines
    sol = kinematics
    im_gandalf, restlineflux*0.0, restlineflux, restferr, velscale, sol, $
;   im_gandalf, reform(fit_tempwave,nkeep,1)*0.0, restlineflux, restferr, velscale, sol, $
      linepars, restwave[0], velscale/light, goodpixels=goodlinepixels, $
      int_disp=line_inst_vdisp, bestfit=bestlinefit, weights=lweights, $
      emission_templates=etemplates, error=esol, $
      /for_errors, plot=plot, quiet=1, chi2=echi2, vmaxshift_iter1=vmaxshift_iter1, $
      vmaxshift_iter2=vmaxshift_iter2, sigmamax_iter1=sigmamax_iter1, $
      sigmamax_iter2=sigmamax_iter2

; parse the output
    if keyword_set(debug) then begin
;      djs_plot, exp(restwave), restlineflux, ps=10, xr=[4800,5100]
;      djs_plot, exp(restwave), restlineflux, ps=10;, xr=[3650,4100]
       djs_plot, exp(restwave), restlineflux, ps=10, xr=[6000,6950]
       djs_oplot, exp(restwave), bestlinefit, ps=10, color='green'
;      djs_oplot, exp(restwave[goodpixels]), restlineflux[goodpixels], ps=4, color='blue'
;      djs_oplot, exp(restwave[goodlinepixels]), restlineflux[goodlinepixels], ps=4, color='red'
       splog, 'Waiting for keystroke...' & cc = get_kbrd(1)
    endif
    
return, sol
end

pro deep2_gandalf_specfit_dr4, zcat, debug=debug, broad=broad, $
  thismask=thismask, firstmask=firstmask, lastmask=lastmask
    
    light = 2.99792458D5 ; speed of light [km/s]
    
; path names and emission-line file name
    version = deep2_version(/ppxf)
    specfitpath = deep2_path(/ppxf,/dr4)
    spec1dpath = deep2_path(/dr4)
    linefile = specfitpath+'gandalf_elinelist_'+version+'.dat'
    if keyword_set(broad) then $
      linefile = repstr(linefile,'elinelist_','elinelist_broad_')

    if n_elements(zcat) eq 0L then zcat = read_deep2_zcat() ; Q34 sample

; figure out which masks we are going to fit
    if n_elements(thismask) eq 0 then $
      thismask = fix(zcat[uniq(zcat.mask,$
      sort(zcat.mask))].mask)
    nmask = n_elements(thismask)

    if (n_elements(firstmask) eq 0L) then firstmask = 0L
    if (n_elements(lastmask) eq 0L) then lastmask = nmask-1L

; fit each mask separately; spectral parameters are from J. Newman+13,
; Table 2
    velscale = deep2_ppxf_velscale()
    line_inst_vdisp = 20.0 ; Gaussian-sigma resolution at z=1 [km/s]
;   line_inst_vdisp = 24.0 ; Gaussian-sigma resolution at z=1 [km/s]
    
; read the emission-line file here to get the total number of lines
; that we are going to be fitting
    junk = read_gandalf_elinelist(linefile,nlines=nlines,$
      fitlines=linepars,/actionfit)
    
; initialize the output data structure
    nfinalpix = 9000L
;   nfinalpix = 2*4096L+500L
    specdata_template = {$
      z:        0.0, $
      line_inst_vdisp: line_inst_vdisp, $
      continuum_snr:   0.0, $
      emission_line_chi2:  0.0, $

      wave:          dblarr(nfinalpix), $ ; double!
      flux:          fltarr(nfinalpix), $
      ferr:          fltarr(nfinalpix), $
      continuum:     fltarr(nfinalpix), $
      etemplates:       fltarr(nfinalpix,nlines),$
      fitlines:                 intarr(nlines)-1,$
      sol:                      fltarr(4,nlines),$
      esol:                     fltarr(4,nlines)}

; fit each mask separately
    t1 = systime(1)
    for imask = firstmask, lastmask do begin
       splog, 'Mask '+strtrim(thismask[imask],2)

; output file names
       suffix = string(thismask[imask],format='(I0)')       
       specdatafile = specfitpath+'specdata_raw_'+suffix+'.fits'
       if keyword_set(broad) then specdatafile = repstr(specdatafile,'raw_','raw_broad_')

; read all the spectra off this mask so we don't have to access
       these = where(thismask[imask] eq fix(zcat.mask),nobj)
       if (nobj eq 0L) then begin
          splog, 'No spectra for mask '+strtrim(thismask[imask],2)
          continue
       endif
       zcat_mask = zcat[these]

; fit each object using GANDALF/PPXF
       t0 = systime(1)
       for iobj = 4, 4 do begin
;      for iobj = 0, nobj-1 do begin
;         print, iobj
          print, format='("Object ",I0,"/",I0,A10,$)', iobj, nobj, string(13b)

          zabs = zcat_mask[iobj].zbest
          vsys = alog(zabs+1.0D)*light ; systemic velocity

; possibly interpolate where the spectra do not overlap and set the
; inverse variance there equal to zero; also possibly normalize the
; spectra to be the same where they (nearly) overlap 
          spec1 = mrdfits(strtrim(spec1dpath+zcat_mask[iobj].file,2)+'.gz',1,/silent)
          spec2 = mrdfits(strtrim(spec1dpath+zcat_mask[iobj].file,2)+'.gz',2,/silent)

          obswave = [spec1.lambda,spec2.lambda]
          obsflux = [spec1.spec,spec2.spec]
          obsinvvar = [spec1.ivar,spec2.ivar]

          srt = sort(obswave)
          obswave = obswave[srt]
          obsflux = obsflux[srt]
          obsinvvar = obsinvvar[srt]          
          obsferr = obsinvvar*0.0+1E16
          good = where(obsinvvar gt 0.0,ngood)
          if ngood eq 0L then message, 'Bad bad!'
          obsferr[good] = 1.0/sqrt(obsinvvar[good])

;         djs_plot, [0], [0], xr=minmax(obswave), yr=minmax(obsflux)
;         djs_oplot, spec2.lambda, spec2.spec, color='red', ps=10
;         djs_oplot, spec1.lambda, spec1.spec, color='blue', ps=10
;         cc = get_kbrd(1)
          
; interpolate the galaxy spectrum logarithmically in wavelength;
; special case for mask 3252 which extends much further into the blue
; than any of the other 400 masks; the redshifts are typically z~0.7,
; so without too much loss of information set the minimum wavelength
; to 6350.0
          if thismask[imask] eq 3252 then begin
             good = where((obswave ge 6500.0) and $
               (obswave le 9200.0),ngood)
          endif else begin
             good = where((obswave ge zcat_mask[iobj].minwave) and $
               (obswave le zcat_mask[iobj].maxwave),ngood)
          endelse
          obswave = obswave[good]
          obsflux = obsflux[good]
          obsferr = obsferr[good]
          obsinvvar = obsinvvar[good]

          log_rebin, minmax(obswave), obsflux, flux, wave, velscale=velscale
          log_rebin, minmax(obswave), obsferr^2, var, velscale=velscale
          ferr = sqrt(abs(var))
          npix = n_elements(wave)

; get the mean S/N
          inrange = where((obsflux gt 0.0) and (obsferr gt 0.0))
;         inrange = where((obsflux gt 0.0) and (obswave gt 7900.0) and $
;           (obswave lt 7700.0) and (obsferr gt 0.0))
          snr = djs_median(obsflux[inrange]/obsferr[inrange])
          
; fit the continuum using a very similiar algorithm to
; AGES_SMOOTH_RESIDUALS (see ages_gandalf_specfit)
          continuum = deep2_fit_unfluxed_continuum(wave,flux,ferr,zabs=zabs,$
            velscale=velscale,broad=broad,linefile=linefile,debug=debug)

; now fit the emission lines, including special cases
          delvarx, vmaxshift_iter1
;         vmaxshift_iter1 = 1000D
          sol = deep2_fit_unfluxed_lines(wave,flux,ferr,continuum,zabs=zabs,$
            velscale=velscale,linefile=linefile,line_inst_vdisp=line_inst_vdisp,$
            broad=broad,linepars=linepars,esol=esol,etemplates=etemplates,$
            echi2=echi2,debug=debug,vmaxshift_iter1=vmaxshift_iter1)
stop
; pack everything into a structure          
          specdata1 = struct_addtags(struct_trimtags(zcat_mask[iobj],$
            select=['galaxy','file','zcatindx','minwave','maxwave','objno',$
            'ra','dec','mask','slit']),specdata_template)
          specdata1.z = zabs

          restwave = wave - alog(zabs+1.0D)
          restflux = flux*(zabs+1.0D)
          restferr = ferr*(zabs+1.0D)
          specdata1.wave[0:npix-1] = double(restwave); + alog(1.0+zabs)
          specdata1.flux[0:npix-1] = float(restflux);/(1.0+zabs)
          specdata1.ferr[0:npix-1] = float(restferr);/(1.0+zabs)
          specdata1.continuum[0:npix-1] = float(continuum);/(1.0+zabs)

; flag the indices of the lines fitted
          fitlines = where((linepars.action ne 'i') and $
            (linepars.kind eq 'l'),nfitlines)
          if (nfitlines ne 0) then begin
             specdata1.fitlines[0:nfitlines-1] = fitlines
             specdata1.sol[*,0:nfitlines-1] = sol
             specdata1.esol[*,0:nfitlines-1] = esol
             specdata1.etemplates[0:npix-1,0:nfitlines-1] = float(etemplates)
             specdata1.emission_line_chi2 = echi2
          endif

          specdata1.continuum_snr = snr

          if (n_elements(specdata) eq 0) then specdata = specdata1 else $
            specdata = [temporary(specdata),specdata1]
       endfor 
       im_mwrfits, specdata, specdatafile, /clobber
       splog, 'Total time = ', (systime(1)-t0)/60.0, ' minutes'
       delvarx, specdata ; delete memory!!
    endfor 
    splog, 'Total time for all masks = ', (systime(1)-t1)/60.0, ' minutes'

return
end
