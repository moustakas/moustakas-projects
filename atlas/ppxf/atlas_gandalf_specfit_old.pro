;+
; NAME:
;   ATLAS_GANDALF_SPECFIT
;
; PURPOSE:
;   Fit the ATLAS and SINGS spectra using PPXF/GANDALF.
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
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Dec 17, UCSD
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

pro atlas_gandalf_specfit, thisfile, test=test, $
  doplot=doplot, broad=broad
    
    light = 2.99792458D5 ; speed of light [km/s]
    fluxscale = 1D-17
    
; path names and emission-line file name
    version = atlas_version(/ppxf_specfit)
    spec1dpath = getenv('IM_DATA_DIR')+'/atlas/spec1d/'
    specfitpath = atlas_path(/ppxf)
    alllinefile = specfitpath+'gandalf_elinelist_all.dat'
    linefile = specfitpath+'gandalf_elinelist_'+version+'.dat'
    if keyword_set(broad) then begin
       alllinefile = repstr(alllinefile,'elinelist_','elinelist_broad_')
       linefile = repstr(linefile,'elinelist_','elinelist_broad_')
    endif

; read the emission-line file here to get the total number of lines
; that we are going to be fitting
    junk = read_gandalf_elinelist(linefile,nlines=nlines,$
      fitlines=linepars,/actionfit)
    
; read the templates (see BUILD_AGES_PPXF_TEMPLATES); resample and
; convolve to AGES pixel size and instrumental resolution
    velscale = atlas_ppxf_velscale()
    inst_vdisp = atlas_ppxf_instvdisp()
    line_inst_vdisp = 5.0 ; inst_vdisp/3.0 ; [km/s]

    tempflux = read_atlas_ppxf_templates(tempwave,ntemp=ntemp,$
      velscale=velscale,inst_vdisp=inst_vdisp)

; initial fitting parameters    
    vdisp_guess = 150.0
    ebv_guess = 0.1
    degree = -1 ; additive Legendre polynomials
    vmaxshift = 2000D   ; maximum velocity shift [km/s]
    sigmamax = 500D ; maximum velocity dispersion [km/s]

    if keyword_set(broad) then begin
       vmaxshift_iter1 = 2000.0D
       vmaxshift_iter2 = 2000.0D
       sigmamax_iter1 = 1D4 ; 2300D ; 1800D
       sigmamax_iter2 = 1D4 ; 2300D ; 1800D
       smooth_sigma = 600.0
    endif else begin
       vmaxshift_iter1 = 500D
       vmaxshift_iter2 = 300D
       sigmamax_iter1 = 500D
       sigmamax_iter2 = 500D
       smooth_sigma = 400.0
    endelse

; generate the file list to be fitted
    objfiles = file_search(spec1dpath+'*.fits.gz',count=nobj)
    if (n_elements(thisfile) eq 1) then begin
       keep = where(strmatch(objfiles,'*'+thisfile*'*'),nkeep)
       if (nkeep eq 0) then begin
          splog, 'File '+thisfile+' not found!'
          return
       endif
       objfiles = objfiles[keep]
       nobj = nkeep
    endif

;; initialize the output data structure
;    nfinalpix = 5000
;    specdata_template = {$
;      ages_id:   0L, $
;;     ra:      0.0D, $
;;     dec:     0.0D, $
;      pass:       0, $
;      aper:       0, $
;      z:        0.0, $
;      zabs:     0.0, $
;      zabs_raw:    0.0, $
;      zabs_err:     0.0, $
;      vdisp:        0.0, $
;      vdisp_raw:    0.0, $
;      vdisp_err:    0.0, $
;      inst_vdisp: inst_vdisp, $
;      line_inst_vdisp: line_inst_vdisp, $
;      continuum_snr:   0.0, $
;      continuum_coeff:     fltarr(ntemp), $
;;     continuum_polycoeff: fltarr(degree+1), $
;      continuum_ebv:       0.0, $
;      continuum_chi2:  0.0, $ ; reduced chi^2
;      emission_line_chi2:  0.0, $
;
;      wave:          dblarr(nfinalpix), $ ; double!
;      flux:          fltarr(nfinalpix), $
;      ferr:          fltarr(nfinalpix), $
;      continuum:     fltarr(nfinalpix), $
;      smooth_continuum:        fltarr(nfinalpix),$
;      etemplates:       fltarr(nfinalpix,nlines),$
;      sol:                      fltarr(4,nlines),$
;      esol:                     fltarr(4,nlines)}

; fit each object using GANDALF/PPXF
    t0 = systime(1)
;   for iobj = 5, 5 do begin
    for iobj = 0, nobj-1 do begin
;      print, format='("Fitting object ",I0,"/",I0,A10,$)', $
;        iobj, nobj-1, string(13b)
       splog, file_basename(objfiles[iobj])+': '+$
         string(iobj+1,format='(I3.3)')+'/'+$
         string(nobj,format='(I3.3)')

; read the spectrum
       linearflux = mrdfits(objfiles[iobj],0,hdr,/silent)/fluxscale ; note!
       linearferr = mrdfits(objfiles[iobj],1,/silent)/fluxscale
       linearwave = make_wave(hdr)
       zabs_guess = sxpar(hdr,'ZNED') ; from NED
;      zabs_guess = 0.006

; rebin logarithmically in wavelength
       log_rebin, minmax(linearwave), linearflux, $
         flux, wave, velscale=velscale
       log_rebin, minmax(linearwave), linearferr^2, $
         var, velscale=velscale
       ferr = sqrt(abs(var))
       npix = n_elements(wave)

; fit the continuum, iterating twice
       zabs_guess1 = zabs_guess
       vdisp_guess1 = vdisp_guess
       ebv_guess1 = ebv_guess
       for iter = 0, 0 do begin
;      for iter = 0, 1 do begin
; shift to the rest frame
          restwave = wave - alog(zabs_guess1+1.0D)
          restflux = flux*(zabs_guess1+1.0D)
          restferr = ferr*(zabs_guess1+1.0D)
       
; initial guesses: compute the velocity offset between the templates
; and the data and restrict the wavelength range of the templates
          vsys = alog(zabs_guess1+1.0D)*light ; systemic velocity
          voffset = -(min(restwave)-min(tempwave))*light
          start = [vsys*0.0+voffset,vdisp_guess1]

          linepars = read_gandalf_elinelist(alllinefile)
          goodpixels = gandalf_mask_emission_lines(npix,vsys*0.0,$
            linepars,velscale,wave[0],velscale/light,$
            sigma=500.0)
          goodpixels = cmset_op(goodpixels,'and',$
            where(restferr lt 1E5))

          ages_ppxf, tempflux, restflux, restferr, velscale, start, $
            sol, goodpixels=goodpixels, plot=doplot, moments=2, degree=degree, $
            weights=weights, polyweights=polyweights, bestfit=bestcontinuum, /clean, $
            /quiet, vmaxshift=vmaxshift, sigmamax=sigmamax, error=err, $
            lambda=exp(restwave), reddening=ebv_guess1
          err = err*sqrt(sol[6]) ; scale by chi^2/dof

          vdisp = sol[1]
          ebv = ebv_guess1
          zabs = exp((vsys+sol[0]-voffset)/light)-1.0D ; new absorption-line redshift
          zabs_err = zabs*(err[0]/light)
          splog, zabs, zabs_err, vdisp, ebv

          zabs_guess1 = zabs
          vdisp_guess1 = vdisp
          
          djs_plot, exp(restwave), restflux, ps=10, xr=[3700,4200]
          djs_oplot, exp(restwave), bestcontinuum, ps=10, color='blue'
          djs_oplot, exp(restwave[goodpixels]), restflux[goodpixels], psym=4, color='red'
;         djs_oplot, exp(restwave[goodlinepixels]), restflux[goodlinepixels], ps=4, color='blue'
          cc = get_kbrd(1)
       endfor

; shift to the rest frame one last time
       restwave = wave - alog(zabs_guess1+1.0D)
       restflux = flux*(zabs_guess1+1.0D)
       restferr = ferr*(zabs_guess1+1.0D)

       vsys = alog(zabs+1.0D)*light ; systemic velocity
       voffset = -(min(restwave)-min(tempwave))*light
       kinematics  = [vsys*0.0+voffset,vdisp,0.0,0.0,0.0,0.0]    ; for AGES_GANDALF

; smooth the residuals (semi-intelligently) before fitting the lines 
       goodpixels1 = gandalf_mask_emission_lines(npix,vsys,$
         linepars,velscale,wave[0],velscale/light,$
         sigma=smooth_sigma)
       badpixels = lindgen(npix)
       remove, goodpixels1, badpixels

       residuals = restflux - bestcontinuum
       msk = restflux*0.0+1.0   ; all are good
       if keyword_set(broad) then sigrej = 1.0 else sigrej = 1.5
       djs_iterstat, residuals, sigrej=sigrej, mask=msk1
       msk1[goodpixels1] = 1    ; keep regions not affected by lines
       msk = msk and msk1

       smooth1 = medsmooth(residuals,151)
       smooth2 = smooth(smooth1,51,/edge)
       linterp, restwave[where(msk)], smooth2[where(msk)], restwave, smooth3
       smooth_continuum = smooth(smooth3,51,/edge)

       djs_plot, exp(restwave), residuals, ysty=3 ;, xr=[2500,3000] ; 4800,5050]
       djs_oplot, exp(restwave[where(msk eq 1)]), residuals[where(msk eq 1)], ps=6, sym=0.7, color='blue'
       djs_oplot, exp(restwave[badpixels]), residuals[badpixels], ps=6, sym=0.2, color='green'
;      djs_oplot, exp(restwave[goodpixels1]), residuals[goodpixels1], ps=6, sym=0.2, color='orange'
       djs_oplot, exp(restwave), smooth_continuum, color='red'

; construct the pure emission-line spectrum          
       restlineflux = restflux - bestcontinuum - smooth_continuum 

; now read the emission-line parameter file of lines we *want* to fit 
       junk = read_gandalf_elinelist(linefile,$
         fitlines=linepars,/actionfit)
       junk = gandalf_mask_emission_lines(npix,vsys,$
         linepars,velscale,wave[0],velscale/light,$
         sigma=500.0)

       goodlinepixels = lindgen(npix)
       remove, goodpixels, goodlinepixels
       goodlinepixels = cmset_op(goodlinepixels,'and',$
         where(restferr lt 1E5))
          
       djs_plot, exp(restwave), restflux, xsty=3, ysty=3, xr=[3700,4200] ; 6500,6700]
       djs_oplot, exp(restwave[goodpixels]), restflux[goodpixels], psym=4, color='red'
       djs_oplot, exp(restwave[goodlinepixels]), restflux[goodlinepixels], ps=4, color='blue'

; now call GANDALF! for the continuum templates pass an array of
; zeros, and just fit for the emission lines
       sol = kinematics
       ages_gandalf, tempwave*0.0, restlineflux, restferr, velscale, sol, $
         linepars, restwave[0], velscale/light, goodpixels=goodlinepixels, $
         int_disp=line_inst_vdisp, bestfit=bestlinefit, weights=lweights, $
         emission_templates=etemplates, error=esol, l0_templ=tempwave[0], $
         /for_errors, plot=doplot, quiet=1, chi2=echi2, vmaxshift_iter1=vmaxshift_iter1, $
         vmaxshift_iter2=vmaxshift_iter2, sigmamax_iter1=sigmamax_iter1, $
         sigmamax_iter2=sigmamax_iter2

; parse the output       
       djs_plot, exp(restwave), restlineflux, ps=10, xr=[6400,6750]
       djs_oplot, exp(restwave), bestlinefit, ps=10, color='green'
       djs_oplot, exp(restwave[goodpixels]), restlineflux[goodpixels], ps=4, color='blue'
       djs_oplot, exp(restwave[goodlinepixels]), restlineflux[goodlinepixels], ps=4, color='red'
       
       bestfit = bestcontinuum + bestlinefit + smooth_continuum
       mask = restferr lt 1E5
       new_sol = gandalf_clean_emission(restwave,restflux,$
         bestfit,mask,etemplates,linepars,sol,esol,$
         snrcut=snrcut_line,velscale=velscale,debug=0,$
         new_etemplates=new_etemplates,new_linepars=new_linepars)
       elinefit = gandalf_parse_elinefit(new_sol,esol,$
         new_linepars,vsys=vsys,fluxscale=1.0)

;;; pack everything into a structure          
;;       specdata1 = specdata_template
;;       struct_assign, spec1d[iobj], specdata1, /nozero
;;       struct_assign, zabsvdisp[iobj], specdata1, /nozero
;;       
;;       fitlines = where((linepars.action ne 'i') and $
;;         (linepars.kind eq 'l'),nfitlines)
;;       if (nfitlines ne 0) then begin
;;          specdata1.sol[*,0:nfitlines-1] = sol
;;          specdata1.esol[*,0:nfitlines-1] = esol
;;       endif
;;       
;;       specdata1.continuum_snr = zabsvdisp[iobj].snr
;;       specdata1.continuum_chi2 = continuum_chi2
;;       specdata1.emission_line_chi2 = echi2
;;       specdata1.continuum_ebv = ebv
;;       specdata1.continuum_coeff = cweights[0:ntemp-1]*fluxscale
;;       
;;       specdata1.wave[0:npix-1] = double(restwave)                                 ; + alog(1.0+zabs)
;;       specdata1.flux[0:npix-1] = float(restflux*fluxscale)                        ;/(1.0+zabs)
;;       specdata1.ferr[0:npix-1] = float(restferr*fluxscale)                        ;/(1.0+zabs)
;;       specdata1.continuum[0:npix-1] = float(bestcontinuum*fluxscale)              ;/(1.0+zabs)
;;       specdata1.smooth_continuum[0:npix-1] = float(smooth_continuum*fluxscale)    ;/(1.0+zabs)
;;       specdata1.etemplates[0:npix-1,0:nfitlines-1] = float(etemplates*fluxscale)
;;       
;;       if (n_elements(specdata) eq 0) then specdata = specdata1 else $
;;         specdata = [temporary(specdata),specdata1]

stop       
       
    endfor
    im_mwrfits, specdata, specdatafile, /clobber
    splog, 'Total time = ', (systime(1)-t0)/60.0, ' minutes'
    
stop    

return
end
