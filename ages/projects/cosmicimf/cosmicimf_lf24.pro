;+
; NAME:
;   COSMICIMF_LF24
;
; PURPOSE:
;   Build the L(24) luminosity functions in AGES in several redshift
;   bins from z=0.05-0.75. 
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Mar 25, UCSD
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

function get_rhosfr, lf24, c_24_lo, c_24_hi, $
  c_24_lo_err, c_24_hi_err, rhosfr_err=rhosfr_err

    dl24 = 0.02
    l24axis = im_array(6.0,15.0,dl24)
    phi = lf_double_powerlaw(l24axis,lf24)

    ff = l24axis*0.0
    lo = where((10^l24axis lt 1.3D10),comp=hi)
    ff[lo] = c_24_lo*10^l24axis[lo]
    ff[hi] = c_24_lo*10^l24axis[hi]*(c_24_hi*10^l24axis[hi])^0.048

    rhosfr = total(alog(10)*dl24*ff*phi) ; [M_sun/Mpc^3]

; get the uncertainty by Monte Carlo; include all sources of error;
; THE ERRORS ARE WAY TOO BIG!
    if arg_present(rhosfr_err) then begin
       nmonte = 500
       delvarx, seed1, seed2, seed3
       alpha_monte = lf24.alpha + randomn(seed1,nmonte)*lf24.alpha_err
       beta_monte = lf24.beta + randomn(seed1,nmonte)*lf24.beta_err
       lstar_monte = lf24.lstar + randomn(seed1,nmonte)*lf24.lstar_err
       phistar_monte = lf24.phistar + randomn(seed1,nmonte)*$
         sqrt(lf24.phistar_err^2+lf24.phistar_lss_err^2)
       c_24_lo_monte = c_24_lo + randomn(seed1,nmonte)*c_24_lo_err
       c_24_hi_monte = c_24_hi + randomn(seed1,nmonte)*c_24_hi_err

       rhosfr_monte = fltarr(nmonte)
       lf24_monte = lf24
       for ii = 0, nmonte-1 do begin
          lf24_monte.phistar = phistar_monte[ii]
          lf24_monte.alpha = alpha_monte[ii]
          lf24_monte.beta = beta_monte[ii]
          lf24_monte.lstar = lstar_monte[ii]
          rhosfr_monte[ii] = get_rhosfr(lf24_monte,c_24_lo_monte[ii],$
            c_24_hi_monte[ii])
       endfor
       rhosfr_err = djsig(rhosfr_monte)
    endif
    
return, rhosfr
end

pro cosmicimf_lf24, debug=debug

    cosmicimfpath = ages_path(/projects)+'cosmicimf/'
    pegpath = cosmicimfpath+'pegase/'

    area = ages_survey_area()*!dtor^2
    h100 = 0.7
    dl24 = 0.02
    l24axis = im_array(6.0,15.0,dl24)

; read the sample, the SFRs, and the Pegase model parameters 
    allsample = read_cosmicimf_sample()
    allsfrs = read_cosmicimf_sample(/sfrs)
    peg = read_cosmicimf_sample(/pegase)
    nimf = n_elements(peg)
    
; compute the SFR density for each IMF using Wiphu's results at
; 0<z<1.2 and Shupe+98's local luminosity function and then
; write out; the errors derived using Monte Carlo are way too big;
; instead, use Wiphu's errors for all IMFs
    wiphu = lf24_wiphu(nimf=nimf)
    wiphu_rhosfr_err = alog(10)*[0.08,0.07,0.08,0.08,0.07,0.11,0.14]*$
      10^[-1.83,-1.60,-1.43,-1.26,-1.07,-0.96,-0.63]
    for ii = 0, nimf-1 do begin
       for jj = 0, n_elements(wiphu)-1 do begin
          wiphu[jj].rhosfr[ii] = get_rhosfr(wiphu[jj],peg[ii].c_24_lo,$
            peg[ii].c_24_hi,peg[ii].c_24_lo_err,peg[ii].c_24_hi_err);,$
;           rhosfr_err=rhosfr_err)
;         wiphu[jj].rhosfr_err[ii] = rhosfr_err
          wiphu[jj].rhosfr_err[ii] = wiphu_rhosfr_err[jj]
       endfor
    endfor

; use Wiphu's quoted error
    shupe = lf24_shupe(nimf=nimf)
    for ii = 0, nimf-1 do begin
       shupe.rhosfr[ii] = get_rhosfr(shupe,peg[ii].c_24_lo,$
         peg[ii].c_24_hi,peg[ii].c_24_lo_err,peg[ii].c_24_hi_err)
       shupe.rhosfr_err[ii] = alog(10)*0.04*10^(-1.71)
    endfor

    wiphufile = cosmicimfpath+'lf24_wiphu.fits'
    shupefile = cosmicimfpath+'lf24_shupe.fits'
    im_mwrfits, wiphu, wiphufile, /clobber
    im_mwrfits, shupe, shupefile, /clobber
    
; now move onto my results: read the redshift bins and the binsize 
    zbins = cosmicimf_zbins(nzbins,/lf24)
    binsize = cosmicimf_binsize(histmin=histmin,histmax=histmax,/lf24)
    limits = mrdfits(cosmicimfpath+'limits_lf24.fits.gz',1,/silent)

; initialize the output data structure and the output filenames; fix
; alpha, beta, and phistar at the local values from Wiphu's paper
    parinfo = init_lf24_parinfo(alpha=wiphu[0].alpha,beta=wiphu[0].beta,$
      phistar=wiphu[0].phistar)
    init_lf24_results, nzbins, nimf, lf24_fit=lf24_fit, lf24_data=lf24_data

; store the local SFRD (for each IMF) as the weighted average of
; Shupe+98 and Wiphu's point at z=0
    for ii = 0, nimf-1 do begin
       lf24_fit.rhosfr_local[ii] = im_weighted_mean([wiphu[0].rhosfr[ii],$
         shupe.rhosfr[ii]],wsigma=rhosfr_local_err)
       lf24_fit.rhosfr_local_err[ii] = rhosfr_local_err
    endfor
    
    fitfile = cosmicimfpath+'lf24_fit.fits'
    datafile = cosmicimfpath+'lf24_data.fits'
    
    for jj = 0, nzbins-1 do begin
       these = where((allsample.z gt zbins[jj].zlo) and $
         (allsample.z lt zbins[jj].zup) and (allsfrs.mips eq 1) and $
         (allsfrs.agn eq 0) and (allsfrs.zmax_24 gt allsample.z),ngal) ; wacky!
       sample = allsample[these]
       sfrs = allsfrs[these]

; compute Vmax       
       zmin = sample.zmin_noevol>zbins[jj].zlo
       zmax = (sample.zmax_noevol<sfrs.zmax_24)<zbins[jj].zup
       vmax = (area/3.0)*(lf_comvol(zmax)-lf_comvol(zmin))/h100^3.0 ; h=0.7
       vol = (area/3.0)*((lf_comvol(zbins[jj].zup)-$
         lf_comvol(zbins[jj].zlo)))[0]/h100^3.0 ; h=0.7

; build the Vmax-weighted 24-micron LF and then fit for L*
       im_lf_vmax, sfrs.l24, sample.final_weight/vmax/binsize, $
         binsize=binsize, histmin=histmin, histmax=histmax, $
         binlum=binlum, phi=phi, errphi=phierr, debug=0, $
         minlum=limits[jj].l24_lim, fullbin=fullbin, number=number

       full = where(fullbin)
       lf_fit_double_powerlaw, binlum[full], phi[full], $
         phierr[full], fit, parinfo=parinfo, quiet=1

; store the results       
       fill_lf24_results, jj, ngal, lf24_fit, lf24_data, fit, $
         binlum, phi, phierr, fullbin, number

; QAplot
       if keyword_set(debug) then begin
          ploterror, binlum, phi, phierr, ps=10, xsty=3, ysty=1, $
            xrange=[8.8,12], yrange=[1E-6,0.1], /ylog
          djs_oplot, l24axis, lf_double_powerlaw(l24axis,fit), color='blue'
          help, fit, /str
          cc = get_kbrd(1)
       endif

; jackknife to get the cosmic variance errors on the parameters; for
; the jackknife error equation see
; http://www.physics.utah.edu/~detar/phycs6730/handouts/jackknife/jackknife
       allfield = sample.field
       field = allfield[uniq(allfield,sort(allfield))]
       nfield = n_elements(field)
       init_lf24_results, nfield, nimf, lf24_fit=lf24_jfit, lf24_data=lf24_jdata
       for ff = 0, nfield-1 do begin
          keep = where(allfield ne field[ff],ngal)
          im_lf_vmax, sfrs[keep].l24, sample[keep].final_weight/vmax[keep]/binsize, $
            binsize=binsize, histmin=histmin, histmax=histmax, $
            binlum=binlum, phi=phi, errphi=phierr, debug=0, $
            minlum=limits[jj].l24_lim, fullbin=fullbin, number=number
          full = where(fullbin)
          lf_fit_double_powerlaw, binlum[full], phi[full], $
            phierr[full], jfit1, parinfo=parinfo, quiet=1
          fill_lf24_results, ff, ngal, lf24_jfit, lf24_jdata, $
            jfit1, binlum, phi, phierr, fullbin, number
       endfor
; jacknife errors on the *parameters*
       jfactor = sqrt((nfield-1)/float(nfield))
       lf24_fit[jj].lstar_cv_err = jfactor*sqrt(total((lf24_jfit.lstar-fit.lstar)^2))
    endfor ; close ZBINS

; finally loop on each redshift bin and IMF and compute the integrated
; SFR density
    for ii = 0, nimf-1 do begin
       for jj = 0, nzbins-1 do begin
          lf24_fit[jj].rhosfr[ii] = get_rhosfr(lf24_fit[jj],peg[ii].c_24_lo,$
            peg[ii].c_24_hi,peg[ii].c_24_lo_err,peg[ii].c_24_hi_err)
       endfor ; close ZBINS
    endfor ; close IMF

; write out
    im_mwrfits, lf24_fit, fitfile, /clobber
    im_mwrfits, lf24_data, datafile, /clobber

return
end

