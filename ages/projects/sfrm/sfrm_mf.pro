;+
; NAME:
;   SFRM_MF
;
; PURPOSE:
;   Build the stellar mass functions (MFs) in AGES for various 
;   subsamples in six redshift bins from z=0.05-0.75.
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;   quiescent - just the quiescent galaxies (see SELECT_QUIESCENT)
;   active - just the actively star-forming galaxies (complement of
;     QUIESCENT) 
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Feb 05, UCSD - 
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

pro sfrm_mf, quiescent=quiescent, active=active, $
  double_schechter=double_schechter, debug=debug

    sfrmpath = ages_path(/projects)+'sfrm/'

    area = ages_survey_area()*!dtor^2
    maxis1 = im_array(7.5,12.5,0.01)
    h100 = 0.7

    suffix = 'all'
    if keyword_set(quiescent) then suffix = 'quiescent'
    if keyword_set(active) then suffix = 'active'
    if keyword_set(double_schechter) then $
      dsuffix = '_double' else dsuffix = ''

    fitfile = sfrmpath+'mf_fit_'+suffix+dsuffix+'.fits'
    datafile = sfrmpath+'mf_data_'+suffix+'.fits'
    
; read the sample, the predefined redshift bins, and the stellar mass
; limits, and initialize the output data structure
    parent = read_sfrm_sample()
    binsize = sfrm_binsize(histmin=histmin,histmax=histmax)

    zbins = sfrm_zbins(nzbins)
    limits = mrdfits(sfrmpath+'sfrm_limits_'+suffix+'.fits.gz',1,/silent)

    if keyword_set(double_schechter) then fixslope = 0 else fixslope = 1
;   fixslope = 0
    parinfo = init_sfrm_mf_parinfo(fixslope=fixslope,quiescent=quiescent,$
      active=active,double_schechter=double_schechter)
    init_sfrm_mf_results, nzbins, mf_fit=mf_fit, mf_data=mf_data, $
      double_schechter=double_schechter

; now measure the mass function in each redshift bin
    for ii = 0, nzbins-1 do begin
       these = where((parent.z gt zbins[ii].zlo) and $
         (parent.z lt zbins[ii].zup))
       sample = parent[these]

       nuvmr = sample.k_galex_absmag_01[1]-sample.k_ugriz_absmag_01[2]
;      rmj = sample.ugriz_absmag[2]-(sample.ubvrijhk_absmag[5]+jv2ab)
       qq = select_quiescent(nuvmr,active=aa)
;      qq = select_quiescent(nuvmr,rmj,active=aa)
       if keyword_set(quiescent) then sample = sample[qq]
       if keyword_set(active) then sample = sample[aa]
       ngal = n_elements(sample)

; compute Vmax       
       zmin = sample.zmin_noevol>zbins[ii].zlo
       zmax = sample.zmax_noevol<zbins[ii].zup
       vmax = (area/3.0)*(lf_comvol(zmax)-lf_comvol(zmin))/h100^3.0 ; h=0.7
       vol = (area/3.0)*((lf_comvol(zbins[ii].zup)-$
         lf_comvol(zbins[ii].zlo)))[0]/h100^3.0 ; h=0.7

; build the Vmax-weighted mass function and then fit it
       splog, 'Using K-correct mass!!!!!'
       mf_vmax, sample.k_mass, sample.final_weight/vmax/binsize, $
         binsize=binsize, histmin=histmin, histmax=histmax, $
         binmass=binmass, phi=phi, errphi=phierr, $
         minmass=limits.minmass[ii], fullbin=fullbin, number=number
       
       full = where(fullbin)
       if keyword_set(double_schechter) then begin
          mf_fit_schechter_plus, binmass[full], phi[full], $
            phierr[full], fit, parinfo=parinfo
       endif else begin
          mf_fit_schechter, binmass[full], phi[full], $
            phierr[full], fit, parinfo=parinfo
       endelse

; store the results       
       fill_sfrm_results, ii, ngal, mf_fit, mf_data, fit, $
         binmass, phi, phierr, fullbin, number
       
; QAplot
       if keyword_set(debug) then begin
          ploterror, binmass, phi, phierr, ps=10, xsty=3, ysty=1, $
            xrange=minmax(maxis1), yrange=[1E-6,0.1], /ylog
          if keyword_set(double_schechter) then $
            djs_oplot, maxis1, mf_schechter_plus(maxis1,fit), color='blue' else $
              djs_oplot, maxis1, mf_schechter(maxis1,fit), color='blue'
          oplot_bgd08_mf, maxis1, params=params, log=0, color='red', line=5
          cc = get_kbrd(1)
       endif
       
; jackknife the MF to get the cosmic variance errors on the
; parameters; for the jackknife error equation see
; http://www.physics.utah.edu/~detar/phycs6730/handouts/jackknife/jackknife
       allfield = sample.field
       field = allfield[uniq(allfield,sort(allfield))]
       nfield = n_elements(field)
       init_mf_results, nfield, mf_fit=mf_jfit, mf_data=mf_jdata
       for ff = 0, nfield-1 do begin
          keep = where(allfield ne field[ff],ngal)
          mf_vmax, sample[keep].k_mass, sample[keep].final_weight/vmax[keep]/binsize, $
            binsize=binsize, histmin=histmin, histmax=histmax, $
            binmass=binmass, phi=phi, errphi=phierr, $
            minmass=limits.minmass[ii], fullbin=fullbin, number=number
          full = where(fullbin)
          if keyword_set(double_schechter) then begin
             mf_fit_schechter_plus, binmass[full], phi[full], $
               phierr[full], jfit1, parinfo=parinfo
          endif else begin
             mf_fit_schechter, binmass[full], phi[full], $
               phierr[full], jfit1, parinfo=parinfo
          endelse
          fill_sfrm_results, ff, ngal, mf_jfit, mf_jdata, $
            jfit1, binmass, phi, phierr, fullbin, number
       endfor
; jacknife errors on the *parameters*
       jfactor = sqrt((nfield-1)/float(nfield))
       mf_fit[ii].phistar_cv_err = jfactor*sqrt(total((mf_jfit.phistar-fit.phistar)^2))
       mf_fit[ii].mstar_cv_err = jfactor*sqrt(total((mf_jfit.mstar-fit.mstar)^2))
       mf_fit[ii].alpha_cv_err = jfactor*sqrt(total((mf_jfit.alpha-fit.alpha)^2))
       if keyword_set(double_schechter) then begin
          mf_fit[ii].phiplus_cv_err = jfactor*sqrt(total((mf_jfit.phiplus-fit.phiplus)^2))
          mf_fit[ii].alphaplus_cv_err = jfactor*sqrt(total((mf_jfit.alphaplus-fit.alphaplus)^2))
       endif
; jacknife errors on the *MFs* themselves
       splog, 'Fix the errors when there are fewer than 3 points!!'
       for jj = 0, mf_jdata[0].nbins-1 do begin
          good = where(mf_jdata.fullbin[jj] eq 1,ngood)
;         if (ngood gt 0) and (ngood le 3) then stop ; splog, 'Need statistics!'
          if (ngood gt 0) then mf_data[ii].phierr_cv[jj] = djsig(mf_jdata.phi[jj])/sqrt(ngood)
       endfor

; vary the way we derive stellar masses to get the errors on mass
; fitting 
       nmodel = n_elements(sample[0].model_mass)
       init_mf_results, nmodel, mf_fit=mf_mfit, mf_data=mf_mdata
       for mm = 0, nmodel-1 do begin
          mf_vmax, sample.model_mass[mm], sample.final_weight/vmax/binsize, $
            binsize=binsize, histmin=histmin, histmax=histmax, $
            binmass=binmass, phi=phi, errphi=phierr, $
            minmass=limits.minmass[ii], fullbin=fullbin, number=number
          full = where(fullbin)
          if keyword_set(double_schechter) then begin
             mf_fit_schechter_plus, binmass[full], phi[full], $
               phierr[full], mfit1, parinfo=parinfo
          endif else begin
             mf_fit_schechter, binmass[full], phi[full], $
               phierr[full], mfit1, parinfo=parinfo
          endelse
          fill_sfrm_results, mm, ngal, mf_mfit, mf_mdata, $
            mfit1, binmass, phi, phierr, fullbin, number
       endfor
; formal model errors on the parameters
       mf_fit[ii].phistar_model_err = djsig(mf_mfit.phistar)
       mf_fit[ii].mstar_model_err = djsig(mf_mfit.mstar)
       mf_fit[ii].alpha_model_err = djsig(mf_mfit.alpha)
       if keyword_set(double_schechter) then begin
          mf_fit[ii].phiplus_model_err = djsig(mf_mfit.phiplus)
          mf_fit[ii].alphaplus_model_err = djsig(mf_mfit.alphaplus)
       endif
; model errors on the *MFs* themselves
       for jj = 0, mf_mdata[0].nbins-1 do begin
          good = where(mf_mdata.fullbin[jj] eq 1,ngood)
          if (ngood gt 0) and (ngood le 3) then splog, 'Need statistics!'
          if (ngood gt 0) then mf_data[ii].phierr_model[jj] = djsig(mf_mdata.phi[jj])/sqrt(ngood)
       endfor
    endfor ; close ZBINS

; write out
    im_mwrfits, mf_fit, fitfile, /clobber
    im_mwrfits, mf_data, datafile, /clobber

return
end

