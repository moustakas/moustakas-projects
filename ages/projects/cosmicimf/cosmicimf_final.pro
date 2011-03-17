;+
; NAME:
;   COSMICIMF_FINAL
;
; PURPOSE:
;   Build the final plots for the COSMICIMF project.
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

function qpint1d_func, x, time=time, sfrd=sfrd
return, interpol(sfrd,time,x)
end

pro cosmicimf_final, out

    cosmicimfpath = ages_path(/projects)+'cosmicimf/'
    pegpath = cosmicimfpath+'pegase/'
    finalfile = cosmicimfpath+'final.fits'
    
; read the results from COSMICIMF_MF and COSMICIMF_LF24 
    zbins = cosmicimf_zbins(nzbins)
    lf24 = mrdfits(cosmicimfpath+'lf24_fit.fits.gz',1)
    dmsalp = mrdfits(cosmicimfpath+'dmsalp.fits.gz',1)
    peg = read_cosmicimf_sample(/pegase)
    nimf = n_elements(peg)

    for ii = 0, nimf-1 do begin
       mf1 = mrdfits(cosmicimfpath+'mf_fit.fits.gz',ii+1,/silent)
       if (ii eq 0) then mf = mf1 else mf = [[mf],[mf1]]
    endfor
    
; initialize the output data structures; for each MF store the
; *measured* stellar mass density (for each IMF) and the *predicted*
; mass density based on the measured SFH (averaged over a handful of
; tau models - see below)
    final = {$
      rhostar:               fltarr(nimf),$ ; stellar mass density for each IMF [M_sun/Mpc^3]
      rhostar_err:           fltarr(nimf),$
      rhostar_local:         fltarr(nimf),$
      rhostar_local_err:     fltarr(nimf),$

      sfh_rhostar1_min:      fltarr(nimf),$ ; zform=1.2
      sfh_rhostar1_max:      fltarr(nimf),$
      sfh_rhostar1_best:     fltarr(nimf),$
      sfh_rhostar5_min:      fltarr(nimf),$ ; zform=5.0
      sfh_rhostar5_max:      fltarr(nimf),$
      sfh_rhostar5_best:     fltarr(nimf)}

;     sfh_rhostar1:           fltarr(nimf),$ ; zform=1.2
;     sfh_rhostar1_err:       fltarr(nimf),$
;     sfh_rhostar5:           fltarr(nimf),$ ; zform=5.0
;     sfh_rhostar5_err:       fltarr(nimf)}
    final = replicate(final,nzbins)

; --------------------------------------------------
; first, copy over the measured stellar mass density to the output
; data structure 
    final.rhostar = transpose(mf.rho)
    final.rhostar_err = transpose(sqrt(mf.rho_err^2+mf.rho_cv_err^2+mf.rho_lss_err^2))

; --------------------------------------------------
; loop through each IMF and compute the local stellar mass density for
; each IMF; below, we compute the *predicted* mass density
    cole = mf_cole()
    bell = mf_bell()
    panter = mf_panter()
    eke = mf_eke()
    gallazzi = mf_gallazzi()
    for ii = 0, nimf-1 do begin
; for the local stellar mass density use the weighted average of
; Bell+03, Cole+01, and AGES at z~0.1 (note: the value is the same for
; all redshift bins)
       conv = 10^(-dmsalp[ii].dmsalp_sfhgrid03) ; salpeter --> IMF[ii]
       conv_err = alog(10)*dmsalp[ii].dmsalp_sfhgrid03_err*conv

       rcole = conv*cole.rho
       rbell = conv*bell.rho
       rpanter = conv*panter.rho
       reke = conv*eke.rho
       rgallazzi = conv*gallazzi.rho
       rages = mf[0,ii].rho
;      allrho = [rcole,rbell,rpanter,reke]
       allrho = [rcole,rbell,rpanter,reke,rgallazzi]
;      allrho = [rcole,rbell]
;      allrho = [rcole,rbell,rages]
       
       rcole_err = im_compute_error(conv,conv_err,cole.rho,cole.rho_err)
       rbell_err = im_compute_error(conv,conv_err,bell.rho,bell.rho_err)
       rpanter_err = im_compute_error(conv,conv_err,panter.rho,panter.rho_err)
       reke_err = im_compute_error(conv,conv_err,eke.rho,eke.rho_err)
       rgallazzi_err = im_compute_error(conv,conv_err,gallazzi.rho,gallazzi.rho_err)
       rages_err = mf[0,ii].rho_err
       allrho_err = [rcole_err,rbell_err,rpanter_err,reke_err,rgallazzi_err]
;      allrho_err = [rcole_err,rbell_err]
;      allrho_err = [rcole_err,rbell_err,rages_err]

       final.rhostar_local[ii] = djs_mean(allrho)
       final.rhostar_local_err[ii] = djsig(allrho)
;      final.rhostar_local[ii] = im_weighted_mean(allrho,$
;        allrho_err,wsigma=rhostar_local_err)
;      final.rhostar_local_err[ii] = rhostar_local_err
    endfor       
    
; --------------------------------------------------
; now deal with the models; read the grid of models that fit the
; observed SFRD evolution for two different formation redshifts (see
; build_cosmicimf_pegase, /sfrd); normalize each model to the "min",
; "max", and "best" SFRD fit outputted by get_sfrd_evolution()
    zform = [1.2,5.0]
    zcut = [1.2,1.0]
    finalfile_models = cosmicimfpath+'final_models_zform'+$
      string(zform,format='(F3.1)')+'.fits'
    for iz = 0, 1 do begin
       sfhfile = 'sfrd_zform'+string(zform[iz],format='(F3.1)')+$
         '_'+['best','min','max']+'.fits'
       nsfh = n_elements(sfhfile)

       sfrd = 10^get_sfrd_evolution(zcut=zcut[iz],$
         zform=zform[iz],dz=0.02,zaxis=zaxis)
       sfrdmin = 10^get_sfrd_evolution(zcut=zcut[iz],$
         zform=zform[iz],dz=0.02,/sfrdmin)
       sfrdmax = 10^get_sfrd_evolution(zcut=zcut[iz],$
         zform=zform[iz],dz=0.02,/sfrdmax)
    
       allsfrd = [[sfrd],[sfrdmin],[sfrdmax]]
       taxis = getage(zaxis)-getage(zform[iz])
       taxis_z0 = getage(0.0)-getage(zaxis) ; lookback time
       taxis = getage(zaxis)-getage(zform[iz])
       nzaxis = n_elements(zaxis)

       zformmodels = {$
         sfhfile:  '',$
         zform:   zform[iz],$
         taxis: taxis,$
         zaxis: zaxis,$
         rhosfr:         fltarr(nzaxis,nimf),$
         sfh_rhostar:    fltarr(nzaxis,nimf),$ 
         sfh_rhostar_z0: fltarr(nzaxis,nimf)} ; integrate from z=0
       zformmodels = replicate(zformmodels,nsfh)
       zformmodels.sfhfile = sfhfile

       for isfh = 0, nsfh-1 do begin
          for ii = 0, nimf-1 do begin
             imfpath = pegpath+strtrim(peg[ii].imf,2)+'/'
             thismodel = cosmicimf_read_peg(imfpath+sfhfile[isfh])
; scale the model output such that the SFR at t=0 matches the
; *measured* SFR density at z=0, for each IMF
             rhosfr_local = allsfrd[0,isfh]*(lf24[0].rhosfr_local[ii]/allsfrd[0,0]) ; z=0
             scale = rhosfr_local/interpol(thismodel.sfr,thismodel.age/1E3,max(taxis))
             zformmodels[isfh].rhosfr[*,ii] = scale*interpol(thismodel.sfr,thismodel.age/1E3,taxis)
             zformmodels[isfh].sfh_rhostar[*,ii] = scale*interpol(thismodel.mstar,thismodel.age/1E3,taxis)
; integrate the SFH from z=0 to z=zform and convert to a stellar mass
; density so that we're not susceptible to the uncertainty in
; the cosmic SFH at z>1; normalize to the *observed* stellar-mass
; density (for each IMF) at z=0
             intsfh_norm = qpint1d('qpint1d_func',min(taxis_z0),max(taxis_z0),$
               functargs={time:taxis_z0,sfrd:zformmodels[isfh].rhosfr[*,ii]})
             for jj = 0, nzaxis-1 do begin
                intsfh = qpint1d('qpint1d_func',min(taxis_z0),taxis_z0[jj],$
                  functargs={time:taxis_z0,sfrd:zformmodels[isfh].rhosfr[*,ii]})
                zformmodels[isfh].sfh_rhostar_z0[jj,ii] = final[0].rhostar_local[ii]*$
                  ((1-intsfh/intsfh_norm)>0.01) ; the 0.01 is to keep the last redshift bin from being zero
             endfor
;            splog, alog10(rhosfr_local), scale
;            plot, zaxis, alog10(zformmodels[isfh].rhosfr[*,ii]), psym=-6, xsty=3, ysty=3
;            plot, zaxis, alog10(zformmodels[isfh].rhosfr[*,0]), ps=6
          endfor 
       endfor

; write out
       im_mwrfits, zformmodels, finalfile_models[iz], /clobber
       if (iz eq 0) then models1 = zformmodels else models5 = zformmodels
    endfor ; close ZFORM

; --------------------------------------------------
; finally compute the mean *predicted* stellar mass density at the
; center of each redshift bin
    for ii = 0, nimf-1 do begin
       if (ii eq 0) then imfpath = pegpath+'Salpeter/' else $
         imfpath = pegpath+strtrim(peg[ii].imf,2)+'/'
       for jj = 0, nzbins-1 do begin
; zform=1.2
          inrange = where((models1[0].zaxis gt zbins[jj].zlo) and $
            (models1[0].zaxis lt zbins[jj].zup))
          final[jj].sfh_rhostar1_best[ii] = djs_mean(models1[0].sfh_rhostar[inrange,ii]) ; best
          final[jj].sfh_rhostar1_min[ii] = djs_mean(models1[1].sfh_rhostar[inrange,ii]) ; min
          final[jj].sfh_rhostar1_max[ii] = djs_mean(models1[2].sfh_rhostar[inrange,ii]) ; max
; zform=5.0
          inrange = where((models5[0].zaxis gt zbins[jj].zlo) and $
            (models5[0].zaxis lt zbins[jj].zup))
          final[jj].sfh_rhostar5_best[ii] = djs_mean(models5[0].sfh_rhostar[inrange,ii]) ; best
          final[jj].sfh_rhostar5_min[ii] = djs_mean(models5[1].sfh_rhostar[inrange,ii]) ; min
          final[jj].sfh_rhostar5_max[ii] = djs_mean(models5[2].sfh_rhostar[inrange,ii]) ; max
       endfor
       splog, peg[ii].imf, alog10(final[0].sfh_rhostar1_max[ii]/final[0].sfh_rhostar1_min[ii]), $
         alog10(final[0].sfh_rhostar5_max[ii]/final[0].sfh_rhostar5_min[ii])
    endfor

; write out    
    im_mwrfits, final, finalfile, /clobber
    
return
end
