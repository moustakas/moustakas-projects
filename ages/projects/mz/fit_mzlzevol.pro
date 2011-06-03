function mlfit_lzevol, absmag, oh, oh_err, weight, z, $
  q0=q0, qz0=qz0, lzlocal=lzlocal, quiet=quiet
; fit the LZ relation with redshift

    if (n_elements(q0) eq 0) then q0 = 1.5 ; luminosity evolution [mag/z]
    
    nparams = 4
    parinfo = {value: 1D, limited: [0,0], limits: [0D,0D], fixed: 0}
    parinfo = replicate(parinfo,nparams)

; fix the first four parameters at their local values
    parinfo[0:1].value = lzlocal.coeff
;   parinfo[0:1].fixed = 1
    parinfo[1].fixed = 1

; P[2]: (O/H)*(z) = (O/H)*(z=0.1) + S*(z-qz0)
    parinfo[2].value = -0.1 ; [dex/z]
;   parinfo[2].limited = 1
;   parinfo[2].limits = [-1.0,0.0]
    
; P[3]: MB*(z) = MB*(z=0.1) + Q*(z-qz0)
    parinfo[3].value = q0 ; [mag/z]
    parinfo[3].fixed = 1  

; do the fit, pack it in, and return
    functargs = {z: z, qz0: qz0, pivotmag: lz_pivotmag()}
    ohweight = weight/oh_err^2

    params = mpfitfun('lzevol_func',absmag,oh,weight=ohweight,$
      parinfo=parinfo,functargs=functargs,perror=perror,dof=dof,$
      covar=covar,status=mpstatus,quiet=quiet,bestnorm=chi2)
    fit = {coeff: params, coeff_err: perror, chi2: chi2, dof: dof}
    
return, fit
end

function get_lzevol_mock, lzlocal, qz0=qz0
; estimate the "false" rate of metallicity evolution due to the loss
; of faint galaxies at higher redshift

; doesn't seem to matter!!    

    zbins = mz_zbins(nz,zmin=zmin,zmax=zmax)
    qb = 1.5
    
    nmock = 5000
    zz = randomu(seed,nmock)*(zmax-zmin)+zmin
    mb = randomu(seed,nmock)*(lzlocal.faintmag-lzlocal.brightmag)+lzlocal.brightmag
    oh = poly(mb-lz_pivotmag(),lzlocal.coeff)+randomn(seed,nmock)*lzlocal.scatter
    mb_evol = mb - qb*(zz-qz0)

    dmod = dmodulus(zz)
    keep = where((mb_evol lt 20.4-dmod^1.01) and (mb_evol gt 15.0-dmod))
    zz = zz[keep]
    mb = mb[keep]
    oh = oh[keep]
    mb_evol = mb_evol[keep]

    plot, zz, mb_evol, psym=6, yrange=[-16,-25]

    lzevol = mlfit_lzevol(mb_evol,oh,oh*0+0.1,oh*0.0+1,$
      zz,q0=0.0,qz0=qz0,lzlocal=lzlocal,quiet=0)
    
    
stop
    
    
;   mbaxis = range(lzlocal.faintmag,lzlocal.brightmag,100) ; for the plot
;   mbcut = [-17.0,-18.0,-19.0,-19.5,-20.5,-21]

    qb = 1.5
    for iz = 0, nz-1 do begin
       keep = where(mb lt mbcut[iz],nkeep)
       zfit1 = randomu(seed,nkeep)*(zbins[iz].zup-zbins[iz].zlo)+zbins[iz].zlo
       mb1 = mb[keep] - qb*(zfit1-qz0)
       oh1 = oh[keep]
       if (iz eq 0) then begin
          mbfit = mb1
          ohfit = oh1
;         zfit = keep*0+zbins[iz].zbin
          zfit = zfit1 
       endif else begin
          mbfit = [mbfit,mb1]
          ohfit = [ohfit,oh1]
;         zfit = [zfit,keep*0+zbins[iz].zbin]
          zfit = [zfit,zfit1]
       endelse
;      djs_plot, mb, oh, psym=6, xsty=3, ysty=3
;      djs_oplot, mb1, oh1, psym=6, color='orange'
;      djs_oplot, mbaxis, poly(mbaxis-lz_pivotmag(),lzlocal.coeff), line=0
;      djs_oplot, mbaxis, lzevol_func(mbaxis,lzevol.params,z=zbins[iz].zbin,$
;        qz0=qz0,pivotmag=lz_pivotmag())
    endfor

; no apparent evolution!!    
    lzevol = mlfit_lzevol(mbfit,ohfit,ohfit*0+0.1,$
      ohfit*0.0+1,zfit,q0=0.0,qz0=qz0,lzlocal=lzlocal,quiet=0)
stop

return, mock
end

function mlfit_mzevol, mass, oh, oh_err, weight, z, $
  p0=p0, r0=r0, qz0=qz0, mzlocal=mzlocal, quiet=quiet
; fit the MZ relation with redshift

    nparams = 5 ; see MLFIT_MZEVOL
    functargs = {z: z, qz0: qz0}

; closed-box model with linear evolution in log (O/H)* and log M* (so
; five total parameters)
    parinfo = {value: 1D, limited: [0,0], limits: [0D,0D], fixed: 0}
    parinfo = replicate(parinfo,nparams)

; fix the first three parameters at their local values
    parinfo[0:2].value = mzlocal.coeff
    parinfo[0:2].fixed = 1

; P[3]: log M*(z) = log M*(z=0.1) + R*(z-qz0)
    if (n_elements(r0) eq 0) then begin ; M* evolution [dex/z]
       parinfo[3].value = +1.0          ; [dex/z]
       parinfo[3].fixed = 0
;      parinfo[3].limited[0] = 1
    endif else begin
       parinfo[3].value = r0            ; [dex/z]
       parinfo[3].fixed = 1
    endelse

; P[4]: (O/H)*(z) = (O/H)*(z=0.1) + P*(z-qz0)
    if (n_elements(p0) eq 0) then begin
       parinfo[4].value = -0.1  ; [dex/z]
       parinfo[4].fixed = 0
;      parinfo[4].limited[1] = 1
    endif else begin
       parinfo[4].value = p0
       parinfo[4].fixed = 1
    endelse
    struct_print, parinfo
       
; do the fit, pack it in, and return
    ohweight = weight/oh_err^2
    params = mpfitfun('mzevol_func',mass,oh,weight=ohweight,$
      parinfo=parinfo,functargs=functargs,perror=perror,dof=dof,$
      covar=covar,status=mpstatus,quiet=quiet,bestnorm=chi2,yfit=yfit)
    fit = {params: params, perror: perror, $
      chi2: chi2, dof: dof, covar: covar}

;   djs_plot, mass, oh, psym=6, ysty=3
;   djs_oplot, mass, yfit, psym=6, color='red'
    
;; Monte Carlo to get the errors on the parameters
;    nmonte = 30
;    params_monte = fltarr(nparams,nmonte)
;    params_err = fltarr(nparams)
;    parinfo_monte = parinfo
;    for mm = 0, nmonte-1 do begin
;       parinfo_monte[0:2].value = mzlocal.coeff+randomn(seed,n_elements(mzlocal.coeff))*mzlocal.coeff_err
;       newoh = oh + randomn(seed,n_elements(oh))*oh_err
;       params_monte[*,mm] = mpfitfun('mzevol_func',mass,newoh,weight=ohweight,$
;         parinfo=parinfo_monte,functargs=functargs,quiet=1)
;    endfor
;    for mm = 0, nparams-1 do params_err[mm] = djsig(params_monte[mm,*])
;    
;    fit = {params: params, perror: params_err, $ ; perror, $
;      chi2: chi2, dof: dof, covar: covar}

return, fit
end

function get_mzevol, ohdust, ancillary, mass, calib=calib, $
  t04=t04, m91=m91, kk04=kk04, qz0=qz0, mzlocal=mzlocal, $
  sdssohdust=sdssohdust, sdssancillary=sdssancillary, $
  sdssmass=sdssmass
; internal support routine to measure the MZ evolution

    zbins = mz_zbins(nz)
    sdss_zbins = mz_zbins(/sdss)
    massbins = mz_massbins(nmassbins)
    
    mzevol = {calib: calib, ngal: n_elements(ohdust), qz0: qz0, $
; AGES #########################
; ngal, median mass, and median redshift in each mass & redshift bin
      ngal_bymass: lonarr(nmassbins,nz), medmass_bymass: fltarr(nmassbins,nz)-999, $
      medz_bymass: fltarr(nmassbins,nz)-999, $
; average metallicity in bins of mass and redshift
      ohmean_bymass: fltarr(nmassbins,nz)-999.0, ohmean_bymass_err: fltarr(nmassbins,nz)-999.0, $
      ohmean_bymass_sigma: fltarr(nmassbins,nz)-999.0, $
; change in metallicity with respect to the metallicity at z=0.1
      dlogoh_bymass: fltarr(nmassbins,nz)-999.0, dlogoh_bymass_err: fltarr(nmassbins,nz)-999.0, $
; SDSS #########################
      sdss_ngal_bymass: lonarr(nmassbins), sdss_medmass_bymass: fltarr(nmassbins)-999, $
      sdss_medz_bymass: fltarr(nmassbins)-999, $
      sdss_ohmean_bymass: fltarr(nmassbins)-999.0, sdss_ohmean_bymass_err: fltarr(nmassbins)-999.0, $
      sdss_ohmean_bymass_sigma: fltarr(nmassbins)-999.0, $
      sdss_dlogoh_bymass: fltarr(nmassbins)-999.0, sdss_dlogoh_bymass_err: fltarr(nmassbins)-999.0, $
; fit the dlogoh evolution      
      coeffs_bymass: fltarr(2,nmassbins)-999.0, coeffs_bymass_err: fltarr(2,nmassbins)-999.0}

; compute the mean metallicity offset from the local MZ relation in
; bins of redshift and stellar mass
    for iz = 0, nz-1 do begin
       zinfo = mzlz_grab_info(ohdust,ancillary,mass,t04=t04,m91=m91,$
         kk04=kk04,zmin=zbins[iz].zlo,zmax=zbins[iz].zup,/nolimit)

; --------------------------------------------------
; the following code block computes dlogoh wrt the local MZ relation;
; it works fine, but we adopt a different procedure, below     

; delta-(O/H) for galaxies with M>10^9; use Monte Carlo to account for
; the uncertainty in the local MZ relation
;          these = where(zinfo.mass gt 9.5,nobj)
;          ohlocal = mz_closedbox(zinfo.mass[these],mzlocal.coeff)
;          mzevol.dlogoh[iz] = im_weighted_mean(zinfo.oh[these]-ohlocal,$
;            weight=zinfo.weight[these]/zinfo.oh_err[these]^2,$
;            wmean_err=wmean_err,wsigma=wsigma)
;          mzevol.dlogoh_sigma[iz] = wsigma
;;         mzevol.dlogoh_err[iz] = wmean_err
;
;          nmonte = 30
;          dlogoh_monte = fltarr(nmonte)
;          for mm = 0, nmonte-1 do begin
;             ohlocal_monte = mz_closedbox(zinfo.mass,mzlocal.coeff+$
;               randomn(seed,n_elements(mzlocal.coeff))*mzlocal.coeff_err)
;             ohmonte = zinfo.oh + randomn(seed,n_elements(zinfo.oh))*zinfo.oh_err
;             dlogoh_monte[mm] = im_weighted_mean(ohmonte-ohlocal_monte,$
;               weight=zinfo.weight/zinfo.oh_err^2)
;          endfor
;          mzevol.dlogoh_err[iz] = sqrt(wmean_err^2+djsig(dlogoh_monte)^2)
; --------------------------------------------------

; mean metallicity in bins of mass
       for jj = 0, nmassbins-1 do begin
          these = where((zinfo.mass ge massbins[jj].lomass) and $
            (zinfo.mass lt massbins[jj].himass),nobj)
          mzevol.ngal_bymass[jj,iz] = nobj
          if (massbins[jj].lomass ge 11.0) and (zbins[iz].zbin gt 0.6) then mingal = 11 else mingal = 9 ; special case 
          if (nobj ge mingal) then begin
             mzevol.medmass_bymass[jj,iz] = djs_median(zinfo.mass[these])
             mzevol.medz_bymass[jj,iz] = djs_median(zinfo.z[these])

             mzevol.ohmean_bymass[jj,iz] = im_weighted_mean(zinfo.oh[these],$
               weight=zinfo.weight[these]/zinfo.oh_err[these]^2,$
               wmean_err=wmean_err,wsigma=wsigma)
             mzevol.ohmean_bymass_err[jj,iz] = wmean_err
             mzevol.ohmean_bymass_sigma[jj,iz] = wsigma
;            print, mzevol.ohmean_ngal[jj,iz], mzevol.ohmean[jj,iz], $
;              mzevol.ohmean_err[jj,iz], mzevol.ohmean_sigma[jj,iz]

; --------------------------------------------------
; the following code block computes dlogoh in bins of mass wrt the
; local MZ relation; it works fine, but we adopt a different
; procedure, below
             
;; delta-(O/H) in bins of mass
;            ohlocal = mz_closedbox(zinfo.mass[these],mzlocal.coeff)
;            mzevol.dlogoh_bymass[jj,iz] = im_weighted_mean(zinfo.oh[these]-ohlocal,$
;              weight=zinfo.weight[these]/zinfo.oh_err[these]^2,$
;              wmean_err=wmean_err,wsigma=wsigma)
;            mzevol.dlogoh_bymass_err[jj,iz] = wmean_err
;            mzevol.dlogoh_bymass_sigma[jj,iz] = wsigma
;;           ohevol = mzevol_func(zinfo.mass[these],mz.params,z=zbins[iz].zbin,qz0=qz0)
;;           ploterror, zinfo.mass[these], zinfo.oh[these], zinfo.oh_err[these], $
;;             psym=6, xr=[9,12], yr=[8,9.5], xsty=3, ysty=3
;;           djs_oplot, zinfo.mass[these], ohevol, line=0, color='blue'
;;           djs_oplot, zinfo.mass[these], ohlocal, line=0, color='red'
; --------------------------------------------------
          endif
       endfor
    endfor

; do the same for the SDSS sample    
    zinfo = mzlz_grab_info(sdssohdust,sdssancillary,sdssmass,$
      t04=t04,m91=m91,kk04=kk04,/nolimit)
    for jj = 0, nmassbins-1 do begin
       these = where((zinfo.mass ge massbins[jj].lomass) and $
         (zinfo.mass lt massbins[jj].himass),nobj)
       mzevol.sdss_ngal_bymass[jj] = nobj
       mzevol.sdss_medmass_bymass[jj] = djs_median(zinfo.mass[these])
       mzevol.sdss_medz_bymass[jj] = djs_median(zinfo.z[these])
       
       mzevol.sdss_ohmean_bymass[jj] = im_weighted_mean(zinfo.oh[these],$
         weight=zinfo.weight[these]/zinfo.oh_err[these]^2,$
         wmean_err=wmean_err,wsigma=wsigma)
       mzevol.sdss_ohmean_bymass_err[jj] = wmean_err
       mzevol.sdss_ohmean_bymass_sigma[jj] = wsigma
    endfor

; now loop back through and fit a linear model to the mean metallicity
; in each mass bin; compute dlogoh as the difference between the
; metallicities at each redshift and the metallicity of the *model*
; fit at z=0.1
    zbin = [sdss_zbins.zbin,zbins.zbin]
    for mm = 0, nmassbins-1 do begin
       ohmean = [mzevol.sdss_ohmean_bymass[mm],reform(mzevol.ohmean_bymass[mm,*])]
       ohmean_err = [mzevol.sdss_ohmean_bymass_err[mm],reform(mzevol.ohmean_bymass_err[mm,*])]
       medz = [mzevol.sdss_medz_bymass[mm],reform(mzevol.medz_bymass[mm,*])]
; ignore the highest redshift point in the log(M)>11 subsample
       gd = where(ohmean gt -900.0,ngd)
       if (massbins[mm].massbin gt 11.0) then $
         fitgd = where((ohmean gt -900.0) and (zbin lt 0.6),nfitgd) else fitgd = gd

       if (ngd ne 0) then begin
          mzevol.coeffs_bymass[*,mm] = linfit(medz[fitgd]-qz0,$
            ohmean[fitgd],measure_err=ohmean_err[fitgd],sigma=sig)
          mzevol.coeffs_bymass_err[*,mm] = sig

; AGES          
          agd = where(mzevol.ohmean_bymass[mm,*] gt -900.0,nagd)
          if (nagd ne 0) then mzevol.dlogoh_bymass[mm,agd] = mzevol.ohmean_bymass[mm,agd]-mzevol.coeffs_bymass[0,mm]
; SDSS
          if (mzevol.sdss_ohmean_bymass[mm] gt -900.0) then mzevol.sdss_dlogoh_bymass[mm] = $
            mzevol.sdss_ohmean_bymass[mm]-mzevol.coeffs_bymass[0,mm]
          
;         zaxis = range(0.0,0.8,50)
;         ploterror, mzevol.medz_bymass[mm,gd], mzevol.ohmean_bymass[mm,gd], $
;           mzevol.ohmean_bymass_err[mm,gd], psym=-6, sym=3, xr=[0,0.7], yr=[8.6,9.1]
;         djs_oplot, zaxis, poly(zaxis-qz0,mzevol.coeffs_bymass[*,mm]), line=0
;         kbrd = get_kbrd(1)
       endif
    endfor
    
return, mzevol
end

pro fit_mzlzevol, mzavg, clobber=clobber
; jm10may14ucsd - maximum likelihood fit of the metallicity evolution

    mzpath = mz_path()

; read the data    
    agesancillary = read_mz_sample(/mzhii_ancillary)
    agesmass = read_mz_sample(/mzhii_mass)
    agesohdust = read_mz_sample(/mzhii_log12oh)

    sdssancillary = read_mz_sample(/mzhii_ancillary,/sdss)
    sdssmass = read_mz_sample(/mzhii_mass,/sdss)
    sdssohdust = read_mz_sample(/mzhii_log12oh,/sdss)

; grid of evolutionary parameters
    qz0 = 0.1 ; reference redshift
    ncalib = 3
    zbins = mz_zbins(nz)
    massbins = mz_massbins(nmassbins)
    
; ##################################################
; fit each calibration separately 
    for ii = 0, ncalib-1 do begin
       case ii of
          0: begin
             t04 = 1 & m91 = 0 & kk04 = 0
             calib = 't04'
          end
          1: begin
             t04 = 0 & m91 = 1 & kk04 = 0
             calib = 'm91'
          end
          2: begin
             t04 = 0 & m91 = 0 & kk04 = 1
             calib = 'kk04'
          end
       endcase

; LZ relation evolution; note that the mock code is ok, but the
; results are not needed 
       lzlocal = mrdfits(mzpath+'lzlocal_sdss_ews_'+calib+'.fits.gz',1)
;      lzevol_mock = get_lzevol_mock(lzlocal,qz0=qz0) 
       
       info = mzlz_grab_info(agesohdust,agesancillary,agesmass,$
         t04=t04,m91=m91,kk04=kk04,/nolimit)
       lzevol = mlfit_lzevol(info.mb_ab,info.oh,info.oh_err,$
         info.weight,info.z,q0=0.0,qz0=qz0,lzlocal=lzlocal,$
         quiet=quiet)
       lzevol = struct_addtags({calib: calib, ngal: n_elements(info.oh), qz0: qz0},lzevol)
       if (ii eq 0) then alllzevol = lzevol else alllzevol = [alllzevol,lzevol]

       lzfile = mzpath+'lzevol_'+calib+'.fits'
       im_mwrfits, lzevol, lzfile, clobber=clobber

; MZ relation evolution
       mzlocal = mrdfits(mzpath+'mzlocal_sdss_ews_'+calib+'.fits.gz',1)
       mzevol = get_mzevol(agesohdust,agesancillary,agesmass,calib=calib,$
         t04=t04,m91=m91,kk04=kk04,qz0=qz0,mzlocal=mzlocal,$
         sdssohdust=sdssohdust,sdssancillary=sdssancillary,sdssmass=sdssmass)
       if (ii eq 0) then allmzevol = mzevol else allmzevol = [allmzevol,mzevol]

       mzfile = mzpath+'mzevol_'+calib+'.fits'
       im_mwrfits, mzevol, mzfile, clobber=clobber
    endfor           

; ##################################################
; average over the three calibrations
    mzavg = {qz0: qz0, coeffs: fltarr(2), coeffs_err: fltarr(2), $
      coeffs_bymass: fltarr(2,nmassbins), coeffs_bymass_err: fltarr(2,nmassbins), $
; AGES
      medz: fltarr(nz), dlogoh: fltarr(nz), dlogoh_err: fltarr(nz), $
      medz_bymass: fltarr(nmassbins,nz)-999, medmass_bymass: fltarr(nmassbins,nz)-999, $
      dlogoh_bymass: fltarr(nmassbins,nz)-999, dlogoh_bymass_err: fltarr(nmassbins,nz)-999, $
; SDSS
      sdss_medz: 0.0, sdss_dlogoh: 0.0, sdss_dlogoh_err: 0.0, $
      sdss_medz_bymass: fltarr(nmassbins)-999, sdss_medmass_bymass: fltarr(nmassbins)-999, $
      sdss_dlogoh_bymass: fltarr(nmassbins)-999, sdss_dlogoh_bymass_err: fltarr(nmassbins)-999, $
; metallicity evolution rate vs mass coefficients
      dlogohdz_medmass: fltarr(nmassbins)-999, dlogohdz_normmass: 0.0, $
      dlogohdz_coeff: fltarr(2), dlogohdz_coeff_err: fltarr(2), $
; evolutionary coefficients for both P and R
      p0: fltarr(ncalib), r0: fltarr(ncalib), p0r0_chi2: fltarr(ncalib), $
;     p0_err: fltarr(ncalib), r0_err: fltarr(ncalib), 
      p0avg: 0.0, p0avg_err: 0.0, r0avg: 0.0, r0avg_err: 0.0, p0r0avg_chi2: 0.0, $
      r0zero_p0avg: 0.0, r0zero_p0avg_err: 0.0, r0zero_chi2: 0.0, $
      p0zero_r0avg: 0.0, p0zero_r0avg_err: 0.0, p0zero_chi2: 0.0, $
; B-band metallicity and luminosity evolution
      lz_slope: 0.0, lz_slope_err: 0.0, lz_dlogohdz: 0.0, lz_dlogohdz_err: 0.0, $
      dmb_bymass: fltarr(nmassbins,nz)-999.0, dmb_bymass_err: fltarr(nmassbins,nz)-999.0}

; average the coefficients and dlog over the three calibrations at
; fixed stellar mass
    for mm = 0, nmassbins-1 do begin
       mzavg.coeffs_bymass[0,mm] = djs_mean(allmzevol.coeffs_bymass[0,mm])
       mzavg.coeffs_bymass[1,mm] = djs_mean(allmzevol.coeffs_bymass[1,mm])
       mzavg.coeffs_bymass_err[0,mm] = djsig(allmzevol.coeffs_bymass[0,mm])
       mzavg.coeffs_bymass_err[1,mm] = djsig(allmzevol.coeffs_bymass[1,mm])

; SDSS: average over all the calibrations at fixed stellar mass
       mzavg.sdss_medmass_bymass[mm] = djs_mean(allmzevol.sdss_medmass_bymass[mm])
       mzavg.sdss_medz_bymass[mm] = djs_mean(allmzevol.sdss_medz_bymass[mm])
       gd = where(allmzevol.sdss_dlogoh_bymass[mm] gt -900.0,ngd)
       if (ngd ne 0) then begin
          mzavg.sdss_dlogoh_bymass[mm] = djs_mean(allmzevol[gd].sdss_dlogoh_bymass[mm])
          mzavg.sdss_dlogoh_bymass_err[mm] = djsig(allmzevol[gd].sdss_dlogoh_bymass[mm])
       endif
       
; AGES: average over all the calibrations at fixed stellar mass and redshift
       for iz = 0, nz-1 do begin
          mzavg.medmass_bymass[mm,iz] = djs_mean(allmzevol.medmass_bymass[mm,iz])
          mzavg.medz_bymass[mm,iz] = djs_mean(allmzevol.medz_bymass[mm,iz])

          gd = where(allmzevol.dlogoh_bymass[mm,iz] gt -900.0,ngd)
          if (ngd ne 0) then begin
             mzavg.dlogoh_bymass[mm,iz] = djs_mean(allmzevol[gd].dlogoh_bymass[mm,iz])
             mzavg.dlogoh_bymass_err[mm,iz] = djsig(allmzevol[gd].dlogoh_bymass[mm,iz])
          endif
       endfor
    endfor

; now average over all the stellar mass bins
    mzavg.sdss_medz = djs_mean(mzavg.sdss_medz_bymass)
    gd = where(mzavg.sdss_dlogoh_bymass[*] gt -900.0,ngd)
    if (ngd ne 0) then begin
       mzavg.sdss_dlogoh = djs_mean(mzavg.sdss_dlogoh_bymass[gd])
       mzavg.sdss_dlogoh_err = djsig(mzavg.sdss_dlogoh_bymass[gd])
    endif

    for iz = 0, nz-1 do begin
       mzavg.medz[iz] = djs_mean(mzavg.medz_bymass[*,iz])
       gd = where(mzavg.dlogoh_bymass[*,iz] gt -900.0,ngd)
       if (ngd ne 0) then begin
          mzavg.dlogoh[iz] = djs_mean(mzavg.dlogoh_bymass[gd,iz])
          mzavg.dlogoh_err[iz] = djsig(mzavg.dlogoh_bymass[gd,iz])
       endif
    endfor
    
; average over all masses and calibrations; ignore the lowest mass
; bin; note that the intercept is defined to be zero at QZ0
    mzavg.coeffs[1] = djs_mean(mzavg.coeffs_bymass[1,0:nmassbins-2])
    mzavg.coeffs_err[1] = djsig(mzavg.coeffs_bymass[1,0:nmassbins-2])

; there is a strong correlation between the rate of metallicity
; evolution and stellar mass; fit it here, ignoring the lowest mass
; bin 
    mzavg.dlogohdz_normmass = 10.5
    for mm = 0, nmassbins-1 do begin       
       gd = where(mzavg.medmass_bymass[mm,*] gt -900.0,ngd)
       if (ngd ne 0) then mzavg.dlogohdz_medmass[mm] = $
         djs_median(mzavg.medmass_bymass[mm,gd])
    endfor
;   slopemass = massbins[0:nmassbins-2].massbin
    slopemass = mzavg.dlogohdz_medmass[0:nmassbins-2]
    slope = mzavg.coeffs_bymass[1,0:nmassbins-2]
    slopeerr = mzavg.coeffs_bymass_err[1,0:nmassbins-2]

    mzavg.dlogohdz_coeff = linfit(slopemass-mzavg.dlogohdz_normmass,slope,$
      measure_err=slopererr,sigma=coeff_err,covar=covar,chisq=chi2)
    mzavg.dlogohdz_coeff_err = coeff_err

; the coefficients above give us the shape of the MZ relation as a
; function of redshift; fit that evolution here with a simple model
; that allows for evolution in M* and (O/H)*; basically, we want to
; solve for P and R, averaged over all three calibrations
    nfake = 2000
    zmin = 0.1 & zmax = 0.9
    minmass = 9.3 & maxmass = 11.5
    zval = randomu(seed,nfake)*(zmax-zmin)+zmin
    mass = randomu(seed,nfake)*(maxmass-minmass)+minmass

    calib = ['kk04','t04','m91']
    for ii = 0, ncalib-1 do begin
       mzlocal = mrdfits(mzpath+'mzlocal_sdss_ews_'+calib[ii]+'.fits.gz',1)

       ohmodel = mz_closedbox(mass,mzlocal.coeff) + (zval-qz0)*$
         poly(mass-mzavg.dlogohdz_normmass,mzavg.dlogohdz_coeff)
       ohmodel_err = ohmodel*0.0+0.01

       mlfit1 = mlfit_mzevol(mass,ohmodel,ohmodel_err,$
         ohmodel*0+1,zval,qz0=qz0,mzlocal=mzlocal,quiet=0)
       mlfit_p0zero1 = mlfit_mzevol(mass,ohmodel,ohmodel_err,p0=0.0,$
         ohmodel*0+1,zval,qz0=qz0,mzlocal=mzlocal,quiet=0)
       mlfit_r0zero1 = mlfit_mzevol(mass,ohmodel,ohmodel_err,r0=0.0,$
         ohmodel*0+1,zval,qz0=qz0,mzlocal=mzlocal,quiet=0)

       if (ii eq 0) then begin
          mlfit = mlfit1
          mlfit_p0zero = mlfit_p0zero1
          mlfit_r0zero = mlfit_r0zero1
       endif else begin
          mlfit = [mlfit,mlfit1]
          mlfit_p0zero = [mlfit_p0zero,mlfit_p0zero1]
          mlfit_r0zero = [mlfit_r0zero,mlfit_r0zero1]
       endelse
    endfor

; each individual calibration; the errors are meaningless 
    mzavg.p0r0_chi2 = mlfit.chi2/float(mlfit.dof)
    mzavg.p0 = mlfit.params[4]
    mzavg.r0 = mlfit.params[3]
;   mzavg.p0_err = mlfit.perror[4]*sqrt(mzavg.p0r0_chi2)
;   mzavg.r0_err = mlfit.perror[3]*sqrt(mzavg.p0r0_chi2)

; average over all calibrations    
    mzavg.p0avg = djs_mean(mlfit.params[4])
    mzavg.r0avg = djs_mean(mlfit.params[3])
    mzavg.p0avg_err = djsig(mlfit.params[4])
    mzavg.r0avg_err = djsig(mlfit.params[3])
    mzavg.p0r0avg_chi2 = djs_mean(mlfit.chi2/float(mlfit.dof))

    mzavg.r0zero_p0avg = djs_mean(mlfit_r0zero.params[4])
    mzavg.r0zero_p0avg_err = djsig(mlfit_r0zero.params[4])
    mzavg.r0zero_chi2 = djs_mean(mlfit_r0zero.chi2/float(mlfit_r0zero.dof))

    mzavg.p0zero_r0avg = djs_mean(mlfit_p0zero.params[3])
    mzavg.p0zero_r0avg_err = djsig(mlfit_p0zero.params[3])
    mzavg.p0zero_chi2 = djs_mean(mlfit_p0zero.chi2/float(mlfit_p0zero.dof))

; average metallicity evolution rate via the LZ relation (assuming no
; luminosity evolution); also compute the average LZ slope
    mzavg.lz_slope = djs_mean(alllzevol.coeff[1])
    mzavg.lz_slope_err = djsig(alllzevol.coeff[1])
    mzavg.lz_dlogohdz = djs_mean(alllzevol.coeff[2]) ; dex per redshift
    mzavg.lz_dlogohdz_err = djsig(alllzevol.coeff[2])

; compare the amount of metallicity evolution from the MZ and LZ
; relations and attribute the differences to luminosity evolution
    for iz = 0, nz-1 do begin
       dlogoh_lum = mzavg.lz_dlogohdz*(zbins[iz].zbin-qz0)
       dlogoh_lum_err = mzavg.lz_dlogohdz_err*(zbins[iz].zbin-qz0)

       dlogoh_mass = mzavg.dlogoh_bymass[*,iz]
       dlogoh_mass_err = mzavg.dlogoh_bymass_err[*,iz]
       gd = where(dlogoh_mass gt -900.0,ngd)

       numer = dlogoh_mass[gd]-dlogoh_lum
       numer_err = dlogoh_mass_err[gd]
;      numer_err = sqrt(dlogoh_mass_err[gd]^2+dlogoh_lum_err^2)

       denom = abs(mzavg.lz_slope)
       denom_err = mzavg.lz_slope_err
       
       mzavg.dmb_bymass[gd,iz] = numer/denom ; [mag/z]
       mzavg.dmb_bymass_err[gd,iz] = im_compute_error(numer,numer_err,$
         denom,denom_err,/quotient)
    endfor

; now go write the paper!    
    mzavgfile = mzpath+'mzevol_avg.fits'
    im_mwrfits, mzavg, mzavgfile, clobber=clobber

return
end
