;+
; NAME:
;   AGES_GANDALF_SPECFIT_UNFLUXED
;
; PURPOSE:
;   Fit the *unfluxed* AGES spectra using PPXF/GANDALF.
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
;   J. Moustakas, 2010 Dec 08, UCSD
;
; Copyright (C) 2010, John Moustakas
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

function select_spec1d, old, index=index
; select certain unfluxed spec1d spectra
    new = {$
      ages_id: old.ages_id[index],$
      galaxy:  old.galaxy[index],$
      ra:      old.ra[index],$
      dec:     old.dec[index],$
      pass:    old.pass,$
      aper:    old.aper[index],$
      class:   old.class[index],$
      z:       old.z[index],$
      spzbest: old.spzbest[index],$
      zmerge_multiple: old.zmerge_multiple[index],$
      minwave: old.minwave[index],$
      maxwave: old.maxwave[index],$
      wave:    old.wave,$
      flux:    old.flux[*,index],$
      ferr:    old.ferr[*,index],$
      skyflux: old.skyflux,$
      skyferr: old.skyferr}
    
return, new
end
    
function ages_fit_unfluxed_continuum, wave, flux, ferr, zabs=zabs, $
  velscale=velscale, linefile=linefile, broad=broad, debug=debug

    light = 2.99792458D5 ; speed of light [km/s]
    vsys = alog(zabs+1.0D)*light ; systemic velocity
    smooth_sigma = 500.0
    sigrej = 1.0
    
    restwave = wave - alog(zabs+1.0D)
    restflux = flux*(zabs+1.0D)
    restferr = ferr*(zabs+1.0D)
    npix = n_elements(wave)

; find all pixels that *could* be affected by emission lines    
    linepars = read_gandalf_elinelist(linefile)
    good1 = gandalf_mask_emission_lines(npix,vsys,linepars,$
      velscale,wave[0],velscale/light,sigma=smooth_sigma)
    bad1 = lindgen(npix)
    remove, good1, bad1

; within the affected regions, mask all deviant pixels
    msk = restflux*0.0+1.0 ; all are good
    msk[bad1] = 0.0
    msk[bad1] = abs(restflux[bad1]) lt sigrej*djsig(restflux[good1]) ; keep non-deviant pixels
    good = where(msk,ngood,comp=bad,ncomp=nbad)

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

function ages_fit_unfluxed_lines, wave, flux, ferr, continuum, $
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
       sigma_smooth = 400.0
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

    goodlinepixels = lindgen(npix)
    remove, goodpixels, goodlinepixels
    goodlinepixels = cmset_op(goodlinepixels,'and',$
      where(restferr lt 1E5))

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

pro ages_gandalf_specfit_unfluxed, pass1, debug=debug, broad=broad
    
    light = 2.99792458D5 ; speed of light [km/s]
    
; path names and emission-line file name
    version = ages_version(/ppxf_specfit)

    specfitpath = ages_path(/spec1d)+'unfluxed/ppxf/'+version+'/' 
    linefile = ages_path(/ppxf)+'gandalf_elinelist_'+version+'.dat'
    if keyword_set(broad) then $
      linefile = repstr(linefile,'elinelist_','elinelist_broad_')

; AGES codes    
    allcodes = mrdfits(ages_path(/cat)+'catalog.codes.fits.gz',1)    
    
; figure out which passes we are going to fit
    if (n_elements(pass1) eq 0) then begin
;      allpass = ages_allpasses()
       if keyword_set(broad) then $
         allpass = ages_balmer_isbroad(/justpass) else $
           allpass = ages_allpasses()
    endif else begin
       allpass = string(pass1,format='(I3.3)')
    endelse
    npass = n_elements(allpass)

    velscale = ages_ppxf_velscale()
    fluxscale = ages_ppxf_fluxscale()
    line_inst_vdisp = 5.0 ; inst_vdisp/3.0 ; [km/s]

; read the emission-line file here to get the total number of lines
; that we are going to be fitting
    junk = read_gandalf_elinelist(linefile,nlines=nlines,$
      fitlines=linepars,/actionfit)
    
; initialize the output data structure
    nfinalpix = 5000
    specdata_template = {$
      ages_id:   0L, $
      pass:       0, $
      aper:       0, $
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

; fit each plate separately
    t1 = systime(1)
    for ipass = 0, npass-1 do begin
       splog, 'Pass ', allpass[ipass]
       
; output file names
       suffix = string(allpass[ipass],format='(I0)')       
       specdatafile = specfitpath+'specdata_raw_'+suffix+'.fits'
       if keyword_set(broad) then specdatafile = repstr(specdatafile,'raw_','raw_broad_')

; read the data
       spec1d = read_ages_spec1d(allpass[ipass],/fix,/unfluxed)

; choose the objects we're going to fit; for simplicity just fit the
; subset of galaxies fitted by AGES_GET_ZABS_VDISP
       codes = allcodes[spec1d.ages_id]
       index = where((codes.gshort gt 0) and $
         (spec1d.z gt 0.001) and (spec1d.z lt 1.0))
       spec1d = select_spec1d(spec1d,index=index)

; restrict to the broad-line AGN, if desired       
       if keyword_set(broad) then begin
          isbroad = ages_balmer_isbroad(spec1d)
          if (isbroad[0] eq -1) then begin
             splog, 'No broad-line AGN in PASS '+$
               strtrim(allpass[ipass])
             continue
          endif else begin
             spec1d = select_spec1d(spec1d,index=isbroad)
          endelse
       endif
       nobj = n_elements(spec1d.ages_id)

; fit each object using GANDALF/PPXF
       t0 = systime(1)
;      for iobj = 3, 3 do begin
;      for iobj = 41, 41 do begin
       for iobj = 0, nobj-1 do begin
          splog, strtrim(allpass[ipass],2)+': '+string(iobj+1,format='(I3.3)')+'/'+$
            string(nobj,format='(I3.3)')

          zabs = spec1d.z[iobj]
          vsys = alog(zabs+1.0D)*light ; systemic velocity

; interpolate the galaxy spectrum logarithmically in wavelength
          good = where((spec1d.wave ge spec1d.minwave[iobj]) and $
            (spec1d.wave le spec1d.maxwave[iobj]),ngood)
          obswave = spec1d.wave[good]
          obsflux = spec1d.flux[good,iobj]
          obsferr = spec1d.ferr[good,iobj]

          log_rebin, minmax(obswave), obsflux, flux, wave, velscale=velscale
          log_rebin, minmax(obswave), obsferr^2, var, velscale=velscale
          ferr = sqrt(abs(var))

          restwave = wave - alog(zabs+1.0D)
          restflux = flux*(zabs+1.0D)
          restferr = ferr*(zabs+1.0D)
          npix = n_elements(restwave)
          
; get the mean S/N
          inrange = where((obsflux gt 0.0) and (obswave gt 6400.0) and $
            (obswave lt 6600.0) and (obsferr gt 0.0))
          snr = djs_median(obsflux[inrange]/obsferr[inrange])
          
; fit the continuum using a very similiar algorithm to
; AGES_SMOOTH_RESIDUALS (see ages_gandalf_specfit)
          continuum = ages_fit_unfluxed_continuum(wave,flux,ferr,zabs=zabs,$
            velscale=velscale,broad=broad,linefile=linefile,debug=debug)

; now fit the emission lines, including special cases
          delvarx, vmaxshift_iter1
          if (allpass[ipass] eq 506) and (spec1d.aper[iobj] eq 57) then vmaxshift_iter1 = 1000D
          sol = ages_fit_unfluxed_lines(wave,flux,ferr,continuum,zabs=zabs,$
            velscale=velscale,linefile=linefile,line_inst_vdisp=line_inst_vdisp,$
            broad=broad,linepars=linepars,esol=esol,etemplates=etemplates,$
            echi2=echi2,debug=debug,vmaxshift_iter1=vmaxshift_iter1)

; pack everything into a structure          
          specdata1 = specdata_template
          specdata1.pass = spec1d.pass
          specdata1.ages_id = spec1d.ages_id[iobj]
          specdata1.aper = spec1d.aper[iobj]
          specdata1.z = spec1d.z[iobj]

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
             specdata1.etemplates[0:npix-1,0:nfitlines-1] = float(etemplates*fluxscale)
          endif

          specdata1.continuum_snr = snr
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
