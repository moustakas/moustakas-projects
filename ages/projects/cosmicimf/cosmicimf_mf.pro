;+
; NAME:
;   COSMICIMF_MF
;
; PURPOSE:
;   Build the stellar mass functions (MFs) in AGES in several redshift
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
;   J. Moustakas, 2010 Mar 23, UCSD
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

pro cosmicimf_mf, debug=debug, fitslope=fitslope

    cosmicimfpath = ages_path(/projects)+'cosmicimf/'

    area = ages_survey_area()*!dtor^2
    massaxis = im_array(7.0,12.0,0.02)
    h100 = 0.7
    nmonte = 300

    fitfile = cosmicimfpath+'mf_fit.fits'
    datafile = cosmicimfpath+'mf_data.fits'
    
; read the sample, the SFRs, and the Pegase model parameters 
    allsample = read_cosmicimf_sample()
    allsfrs = read_cosmicimf_sample(/sfrs)
    allmass = read_cosmicimf_sample(/mass)
    peg = read_cosmicimf_sample(/pegase)
    imf = strtrim(peg.imf,2)
    nimf = n_elements(imf)
    
; read the redshift bins and the binsize
    zbins = cosmicimf_zbins(nzbins)
    binsize = cosmicimf_binsize(histmin=histmin,histmax=histmax)

; initialize the output data structure; fix the faint-end slope and
; uncertainty unless /FITSLOPE
    alpha = -1.15D
    alpha_err = 0.1D
    parinfo = cosmicimf_init_mf_parinfo(alpha=alpha,fitslope=fitslope)
    init_cosmicimf_mf_results, nzbins, mf_fit=mf_fit, mf_data=mf_data

; loop on each IMF+Salpeter and measure the MF in each redshift bin
    for ii = 0, nimf-1 do begin
       splog, 'IMF '+imf[ii]
       limits = mrdfits(cosmicimfpath+'limits_mf.fits.gz',ii+1)
;      for jj = 5, nzbins-1 do begin
       for jj = 0, nzbins-1 do begin
          these = where((allsample.z gt zbins[jj].zlo) and $
            (allsample.z lt zbins[jj].zup) and $
            (allsample.zmax_evol gt allsample.z) and $
            (allsfrs.agn eq 0),ngal) ; wacky!
          sample = allsample[these]
;         mass = allmass[these].sfhgrid02[ii]
          mass = allmass[these].sfhgrid03[ii]

; compute Vmax       
;         zmin = sample.zmin_noevol>zbins[jj].zlo
;         zmax = sample.zmax_noevol<zbins[jj].zup
          zmin = sample.zmin_evol>zbins[jj].zlo
          zmax = sample.zmax_evol<zbins[jj].zup
          vmax = (area/3.0)*(lf_comvol(zmax)-lf_comvol(zmin))/h100^3.0 ; h=0.7
          vol = (area/3.0)*((lf_comvol(zbins[jj].zup)-$
            lf_comvol(zbins[jj].zlo)))[0]/h100^3.0 ; h=0.7

; build the Vmax-weighted MF and then fit it
          mf_vmax, mass, sample.final_weight/vmax/binsize, $
            binsize=binsize, histmin=histmin, histmax=histmax, $
            binmass=binmass, phi=phi, errphi=phierr, $
            minmass=limits[jj].mass_lim, fullbin=fullbin, $
            number=number, debug=0

          full = where(fullbin)
          mf_fit_schechter, binmass[full], phi[full], $
            phierr[full], fit, parinfo=parinfo, /quiet;, /norhoerr
          
; get the parameter errors using Monte Carlo; we need to do this here
; (rather than using the formal MPFIT errors) because we fix the
; faint-end slope *alpha* and there is considerable covariance between
; M* and alpha; as a result, the Monte Carlo uncertainties are factors
; of several times higher here than the formal errors
          parinfo_monte = parinfo
          delvarx, seed1, seed2, seed3
          phistar_monte = fit.phistar + randomn(seed1,nmonte)*fit.phistar_err
          mstar_monte = fit.mstar + randomn(seed2,nmonte)*fit.mstar_err
          alpha_monte = alpha + randomn(seed3,nmonte)*alpha_err
          for imonte = 0, nmonte-1 do begin
             parinfo_monte[0].value = phistar_monte[imonte]
             parinfo_monte[1].value = mstar_monte[imonte]
             parinfo_monte[2].value = alpha_monte[imonte]
             mf_fit_schechter, binmass[full], phi[full], $
               phierr[full], mfit1, parinfo=parinfo_monte, $
               /norhoerr, /quiet
             if (imonte eq 0) then mfit = mfit1 else mfit = [mfit,mfit1]
          endfor
          fit.phistar_err = djsig(mfit.phistar)
          fit.mstar_err = djsig(mfit.mstar)
          fit.alpha_err = alpha_err
          fit.rho_err = djsig(mfit.rho)

; store the results       
          fill_cosmicimf_mf_results, jj, ngal, mf_fit, mf_data, fit, $
            binmass, phi, phierr, fullbin, number
;         help, mf_fit[jj], /st

; QAplot
          if keyword_set(debug) then begin
             ploterror, binmass, phi, phierr, ps=10, xsty=3, ysty=1, $
               xrange=minmax(massaxis), yrange=[1E-6,0.1], /ylog
             djs_oplot, massaxis, mf_schechter(massaxis,fit), color='blue'
             oplot_bgd08_mf, massaxis, params=params, log=0, color='red', $
               line=5, /salpeter
             help, fit, /str
             cc = get_kbrd(1)
          endif

; jackknife the MF to get the cosmic variance errors on the
; parameters; for the jackknife error equation see
; http://www.physics.utah.edu/~detar/phycs6730/handouts/jackknife/jackknife
          allfield = sample.field
          field = allfield[uniq(allfield,sort(allfield))]
          nfield = n_elements(field)
          init_cosmicimf_mf_results, nfield, mf_fit=mf_jfit, mf_data=mf_jdata
          for ff = 0, nfield-1 do begin
             keep = where(allfield ne field[ff],ngal)
;            mf_vmax, allmass[these[keep]].sfhgrid02[ii], $
             mf_vmax, allmass[these[keep]].sfhgrid03[ii], $
               sample[keep].final_weight/vmax[keep]/binsize, $
               binsize=binsize, histmin=histmin, histmax=histmax, $
               binmass=binmass, phi=phi, errphi=phierr, $
               minmass=limits[jj].mass_lim, fullbin=fullbin, $
               number=number
             full = where(fullbin)
             mf_fit_schechter, binmass[full], phi[full], $
               phierr[full], jfit1, parinfo=parinfo, quiet=1, $
               alpha_err=alpha_err, /norhoerr
             fill_cosmicimf_mf_results, ff, ngal, mf_jfit, mf_jdata, $
               jfit1, binmass, phi, phierr, fullbin, number, /nolss
          endfor
; jacknife errors on the *parameters*
          jfactor = sqrt((nfield-1)/float(nfield))
          mf_fit[jj].phistar_cv_err = jfactor*sqrt(total((mf_jfit.phistar-fit.phistar)^2))
          mf_fit[jj].mstar_cv_err = jfactor*sqrt(total((mf_jfit.mstar-fit.mstar)^2))
          mf_fit[jj].alpha_cv_err = jfactor*sqrt(total((mf_jfit.alpha-fit.alpha)^2))
          mf_fit[jj].rho_cv_err = jfactor*sqrt(total((mf_jfit.rho-fit.rho)^2))
       endfor ; close ZBINS
; write out       
       if (ii eq nimf-1) then begin
          splog, 'Writing '+fitfile
          splog, 'Writing '+datafile
       endif
       mwrfits, mf_fit, fitfile, create=(ii eq 0)
       mwrfits, mf_data, datafile, create=(ii eq 0)
    endfor ; close IMF
    spawn, 'gzip -f '+fitfile, /sh
    spawn, 'gzip -f '+datafile, /sh
    
return
end

