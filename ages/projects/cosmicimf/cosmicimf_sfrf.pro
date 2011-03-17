;+
; NAME:
;   COSMICIMF_MF
;
; PURPOSE:
;   Build the SFR functions (SFRFs) in AGES in several redshift bins
;   from z=0.05-0.75. 
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
;   J. Moustakas, 2010 Mar 18, UCSD - 
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

pro cosmicimf_sfrf, debug=debug

    cosmicimfpath = ages_path(/projects)+'cosmicimf/'

    area = ages_survey_area()*!dtor^2
    sfraxis1 = im_array(-1.0,2.0,0.02)
    h100 = 0.7

    fitfile = cosmicimfpath+'sfrf_fit.fits'
    datafile = cosmicimfpath+'sfrf_data.fits'
    
; read the sample, the SFRs, and the Pegase model parameters 
    allsample = read_cosmicimf_sample()
    allsfrs = read_cosmicimf_sample(/sfrs)
    peg = read_cosmicimf_sample(/pegase)
    nimf = n_elements(peg)
    
; read the redshift bins and the binsize
    zbins = cosmicimf_zbins(nzbins,/sfrf)
    binsize = cosmicimf_binsize(histmin=histmin,histmax=histmax,/sfrf)
    limits = mrdfits(cosmicimfpath+'cosmicimf_limits_sfrf.fits.gz',1,/silent)

; initialize the output data structure    
    parinfo = init_sfrf_parinfo(fixslope=fixslope,schechter=schechter)
;    init_sfrf_results, nzbins, sfrf_fit=sfrf_fit, sfrf_data=sfrf_data, $
;      double_schechter=double_schechter

; loop on each IMF+Salpeter and measure the SFRF in each redshift bin
    for ii = -1, -1 do begin
;   for ii = -1, nimf-1 do begin
       for jj = 0, nzbins-1 do begin
          these = where((allsample.z gt zbins[jj].zlo) and $
            (allsample.z lt zbins[jj].zup) and (allsfrs.agn eq 0) and $
            (allsfrs.zmax_24 gt allsample.z),ngal) ; wacky!
          sample = allsample[these]
          sfrs = allsfrs[these]

; compute Vmax       
          zmin = sample.zmin_noevol>zbins[jj].zlo
          zmax = (sample.zmax_noevol<sfrs.zmax_24)<zbins[jj].zup
          vmax = (area/3.0)*(lf_comvol(zmax)-lf_comvol(zmin))/h100^3.0 ; h=0.7
          vol = (area/3.0)*((lf_comvol(zbins[jj].zup)-$
            lf_comvol(zbins[jj].zlo)))[0]/h100^3.0 ; h=0.7

; compute the SFR using this IMF and Kennicutt+98
          c_1500_k98 = -alog10(1.4D-28) - alog10(1500.0/im_light(/ang))
          if (ii eq -1) then c_1500 = peg[0].c_1500_salp else c_1500 = peg[ii].c_1500
          sfr = sfrs.l1500_cor + alog10(3.826D33) - c_1500 ; [M_sun/yr]
          sfr_uncor = sfrs.l1500 + alog10(3.826D33) - c_1500 ; [M_sun/yr]
          sfr_k98 =  sfrs.l1500_cor + alog10(3.826D33) - c_1500_k98

; build the Vmax-weighted SFRF and then fit it
          abovelimit = where(sfrs.l1500 gt limits[jj].l1500_lim)
          minsfr = min(sfr[abovelimit])
          
          sfrf_vmax, sfr, sample.final_weight/vmax/binsize, $
            binsize=binsize, histmin=histmin, histmax=histmax, $
            binsfr=binsfr, phi=phi, errphi=phierr, $
            minsfr=minsfr, fullbin=fullbin, number=number, debug=0
;         print, alog10(im_integral(10^binsfr,phi)) ; [M_sun/yr/Mpc^3]

          full = where(fullbin)
          sfrf_fit_double_powerlaw, 10^binsfr[full], phi[full], $
            phierr[full], fit, parinfo=parinfo
;         sfrf_fit_schechter, 10^binsfr[full], phi[full], $
;           phierr[full], fit, parinfo=parinfo

; store the results       
;         fill_cosmicimf_results, jj, ngal, sfrf_fit, sfrf_data, fit, $
;           binsfr, phi, phierr, fullbin, number
          
; QAplot
          if keyword_set(debug) then begin
             ploterror, binsfr, phi, phierr, ps=10, xsty=3, ysty=1, $
               xrange=minmax(sfraxis1), yrange=[1E-6,0.1], /ylog
             djs_oplot, sfraxis1, sfrf_double_powerlaw(10^sfraxis1,fit), color='blue'
;            djs_oplot, sfraxis1, sfrf_schechter(10^sfraxis1,fit), color='blue'
             oplot_bgd08_mf, sfraxis1, params=params, log=0, color='red', line=5
             cc = get_kbrd(1)
          endif

stop          
          
; jackknife the MF to get the cosmic variance errors on the
; parameters; for the jackknife error equation see
; http://www.physics.utah.edu/~detar/phycs6730/handouts/jackknife/jackknife
          allfield = sample.field
          field = allfield[uniq(allfield,sort(allfield))]
          nfield = n_elements(field)
          init_sfrf_results, nfield, sfrf_fit=sfrf_jfit, sfrf_data=sfrf_jdata
          for ff = 0, nfield-1 do begin
             keep = where(allfield ne field[ff],ngal)
             sfrf_vmax, sfr[keep], sample[keep].final_weight/vmax[keep]/binsize, $
               binsize=binsize, histmin=histmin, histmax=histmax, $
               binsfr=binsfr, phi=phi, errphi=phierr, $
               minsfr=minsfr, fullbin=fullbin, number=number
             full = where(fullbin)
             if keyword_set(double_schechter) then begin
                sfrf_fit_schechter_plus, 10^binsfr[full], phi[full], $
                  phierr[full], jfit1, parinfo=parinfo
             endif else begin
                sfrf_fit_schechter, 10^binsfr[full], phi[full], $
                  phierr[full], jfit1, parinfo=parinfo
             endelse
             fill_cosmicimf_results, ff, ngal, sfrf_jfit, sfrf_jdata, $
               jfit1, binsfr, phi, phierr, fullbin, number
          endfor
; jacknife errors on the *parameters*
          jfactor = sqrt((nfield-1)/float(nfield))
          sfrf_fit[jj].phistar_cv_err = jfactor*sqrt(total((sfrf_jfit.phistar-fit.phistar)^2))
          sfrf_fit[jj].mstar_cv_err = jfactor*sqrt(total((sfrf_jfit.mstar-fit.mstar)^2))
          sfrf_fit[jj].alpha_cv_err = jfactor*sqrt(total((sfrf_jfit.alpha-fit.alpha)^2))
          if keyword_set(double_schechter) then begin
             sfrf_fit[jj].phiplus_cv_err = jfactor*sqrt(total((sfrf_jfit.phiplus-fit.phiplus)^2))
             sfrf_fit[jj].alphaplus_cv_err = jfactor*sqrt(total((sfrf_jfit.alphaplus-fit.alphaplus)^2))
          endif
; jacknife errors on the *MFs* themselves
          for kk = 0, sfrf_jdata[0].nbins-1 do begin
             good = where(sfrf_jdata.fullbin[kk] eq 1,ngood)
             if (ngood gt 0) and (ngood le 3) then stop ; splog, 'Need statistics!'
             if (ngood gt 0) then sfrf_data[jj].phierr_cv[kk] = djsig(sfrf_jdata.phi[kk])/sqrt(ngood)
          endfor
       endfor ; close ZBINS

; write out
       im_mwrfits, sfrf_fit, fitfile, /clobber
       im_mwrfits, sfrf_data, datafile, /clobber
    endfor ; close IMF
       
return
end

