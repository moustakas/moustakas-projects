function mlfit_lzevol, absmag, oh, oh_err, weight, z, $
  q0=q0, qz0=qz0, lzlocal=lzlocal, quiet=quiet
; fit the LZ relation with redshift

    if (n_elements(q0) eq 0) then q0 = -1.5 ; luminosity evolution [mag/z]
    
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
    fit = {coeff: params, coeff_err: perror, chi2: chi2/dof}
    
return, fit
end

function get_lzevol_mock, lzlocal, qz0=qz0
; estimate the "false" rate of metallicity evolution due to the loss
; of faint galaxies at higher redshift

; doesn't seem to matter!!    

    mzpath = mz_path()
    limits = mrdfits(mzpath+'mz_limits.fits.gz',1)
    zbins = mz_zbins(nz,zmin=zmin,zmax=zmax)

    qb = 0.0 ; 1.5
    qb = 3.5
    
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
    
    mbaxis = range(lzlocal.faintmag,lzlocal.brightmag,100) ; for the plot
    mbcut = limits.mblimit_50

    qb = 3.5
;   qb = 0.0
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
       djs_plot, mb, oh, psym=6, xsty=3, ysty=3
       djs_oplot, mb1, oh1, psym=6, color='orange'
       djs_oplot, mbaxis, poly(mbaxis-lz_pivotmag(),lzlocal.coeff), line=0
       djs_oplot, mbaxis, lzevol_func(mbaxis,lzevol.coeff,z=zbins[iz].zbin,$
         qz0=qz0,pivotmag=lz_pivotmag())
       cc = get_kbrd(1)
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
       parinfo[3].value = +1.0         ; [dex/z]
       parinfo[3].fixed = 0
       parinfo[3].limited = 1
       parinfo[3].limits = [-3.0,+3.0]
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
      chi2: chi2/dof, covar: covar}
;   djs_plot, mass, oh, psym=6, ysty=3
;   djs_oplot, mass, yfit, psym=7, color='red'

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

    mzpath = mz_path()
    limits = mrdfits(mzpath+'mz_limits.fits.gz',1)

    zbins = mz_zbins(nz)
    sdss_zbins = mz_zbins(/sdss)
    allmassbins = mz_massbins(nallmassbins)
    
    mzevol = {calib: calib, ngal: n_elements(ohdust), qz0: qz0, $
;     ngal: lonarr(nz), medz: fltarr(nz)-999, $
; AGES #########################
; ngal, median mass, and median redshift in each mass & redshift bin
      ngal_bymass: lonarr(nallmassbins,nz), medmass_bymass: fltarr(nallmassbins,nz)-999, $
      medz_bymass: fltarr(nallmassbins,nz)-999, complete_bymass: intarr(nallmassbins,nz), $
; average metallicity in bins of mass and redshift
      ohmean_bymass: fltarr(nallmassbins,nz)-999.0, ohmean_bymass_err: fltarr(nallmassbins,nz)-999.0, $
      ohmean_bymass_sigma: fltarr(nallmassbins,nz)-999.0, $
; change in metallicity with respect to the metallicity at z=0.1
      dlogoh_bymass: fltarr(nallmassbins,nz)-999.0, dlogoh_bymass_err: fltarr(nallmassbins,nz)-999.0, $
; SDSS #########################
      sdss_ngal_bymass: lonarr(nallmassbins), sdss_medmass_bymass: fltarr(nallmassbins)-999, $
      sdss_medz_bymass: fltarr(nallmassbins)-999, $
      sdss_ohmean_bymass: fltarr(nallmassbins)-999.0, sdss_ohmean_bymass_err: fltarr(nallmassbins)-999.0, $
      sdss_ohmean_bymass_sigma: fltarr(nallmassbins)-999.0, $
      sdss_dlogoh_bymass: fltarr(nallmassbins)-999.0, sdss_dlogoh_bymass_err: fltarr(nallmassbins)-999.0, $
; fit the dlogoh evolution      
      coeffs_bymass: fltarr(2,nallmassbins)-999.0, coeffs_bymass_err: fltarr(2,nallmassbins)-999.0}

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
       ohminerr = 0.0
       for mm = 0, nallmassbins-1 do begin
          these = where((zinfo.mass ge allmassbins[mm].lomass) and $
            (zinfo.mass lt allmassbins[mm].himass),nobj)
          mzevol.ngal_bymass[mm,iz] = nobj
;         if (allmassbins[mm].lomass ge 11.0) and (zbins[iz].zbin gt 0.6) then mingal = 11 else mingal = 9 ; special case 
          mingal = 10
          if (nobj ge mingal) then begin
             mzevol.complete_bymass[mm,iz] = allmassbins[mm].lomass gt limits.masslimit_50[iz]

             mzevol.medmass_bymass[mm,iz] = im_weighted_mean(zinfo.mass[these],weight=zinfo.weight[these])
             mzevol.medz_bymass[mm,iz] = im_weighted_mean(zinfo.z[these],weight=zinfo.weight[these])

             oherr = sqrt(ohminerr^2+zinfo.oh_err[these]^2)
             mzevol.ohmean_bymass[mm,iz] = im_weighted_mean(zinfo.oh[these],$
               weight=zinfo.weight[these]/oherr^2,$
               wmean_err=wmean_err,wsigma=wsigma)
;            plot, zinfo.oh[these], zinfo.oh_err[these], psym=6, xsty=3, ysty=3, yr=[0,0.3], xr=[8.5,9.2]
;            cc = get_kbrd(1)
             
             mzevol.ohmean_bymass_err[mm,iz] = wmean_err
             mzevol.ohmean_bymass_sigma[mm,iz] = wsigma

; median 
;            quant = [1.0-gauss_pdf(2.0),0.5,gauss_pdf(2.0)]
;            qquant = weighted_quantile(zinfo.oh[these],$
;              zinfo.weight[these]/zinfo.oh_err[these]^2,quant=quant)
;            mzevol.ohmean_bymass[mm,iz] = qquant[1]
;            mzevol.ohmean_bymass_sigma[mm,iz] = (qquant[2]-qquant[0])/4.0
;            mzevol.ohmean_bymass_err[mm,iz] = mzevol.ohmean_bymass_sigma[mm,iz]/sqrt(nobj)
;            mzevol.medmass_bymass[mm,iz] = djs_median(zinfo.mass[these])
;            mzevol.medz_bymass[mm,iz] = djs_median(zinfo.z[these])

;            mzevol.ohmean_bymass_err[mm,iz] = wmean_err
;            mzevol.ohmean_bymass_sigma[mm,iz] = wsigma
;            print, mzevol.ohmean_ngal[mm,iz], mzevol.ohmean[mm,iz], $
;              mzevol.ohmean_err[mm,iz], mzevol.ohmean_sigma[mm,iz]

; --------------------------------------------------
; the following code block computes dlogoh in bins of mass wrt the
; local MZ relation; it works fine, but we adopt a different
; procedure, below
             
;; delta-(O/H) in bins of mass
;            ohlocal = mz_closedbox(zinfo.mass[these],mzlocal.coeff)
;            mzevol.dlogoh_bymass[mm,iz] = im_weighted_mean(zinfo.oh[these]-ohlocal,$
;              weight=zinfo.weight[these]/zinfo.oh_err[these]^2,$
;              wmean_err=wmean_err,wsigma=wsigma)
;            mzevol.dlogoh_bymass_err[mm,iz] = wmean_err
;            mzevol.dlogoh_bymass_sigma[mm,iz] = wsigma
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
    for mm = 0, nallmassbins-1 do begin
       these = where((zinfo.mass ge allmassbins[mm].lomass) and $
         (zinfo.mass lt allmassbins[mm].himass),nobj)
       mzevol.sdss_ngal_bymass[mm] = nobj

       mzevol.sdss_medmass_bymass[mm] = im_weighted_mean(zinfo.mass[these],weight=zinfo.weight[these])
       mzevol.sdss_medz_bymass[mm] = im_weighted_mean(zinfo.z[these],weight=zinfo.weight[these])
       
       mzevol.sdss_ohmean_bymass[mm] = im_weighted_mean(zinfo.oh[these],$
         weight=zinfo.weight[these]/zinfo.oh_err[these]^2,$
         wmean_err=wmean_err,wsigma=wsigma)
       mzevol.sdss_ohmean_bymass_err[mm] = wmean_err
       mzevol.sdss_ohmean_bymass_sigma[mm] = wsigma

;      mzevol.sdss_medmass_bymass[mm] = djs_median(zinfo.mass[these])
;      mzevol.sdss_medz_bymass[mm] = djs_median(zinfo.z[these])
    endfor

; now loop back through and fit a linear model to the mean metallicity
; in each mass bin; compute dlogoh as the difference between the
; metallicities at each redshift and the metallicity of the *model*
; fit at z=0.1
    zbin = [sdss_zbins.zbin,zbins.zbin]
    for mm = 0, nallmassbins-1 do begin
       complete = [1,reform(mzevol.complete_bymass[mm,*])]
       ohmean = [mzevol.sdss_ohmean_bymass[mm],reform(mzevol.ohmean_bymass[mm,*])]
       ohmean_err = [mzevol.sdss_ohmean_bymass_err[mm],reform(mzevol.ohmean_bymass_err[mm,*])]
       medz = [mzevol.sdss_medz_bymass[mm],reform(mzevol.medz_bymass[mm,*])]

       fitgd = where(complete and ohmean gt -900.0,ngd)
;      if (allmassbins[mm].massbin gt 11.0) then $
;        fitgd = where((ohmean gt -900.0) and (zbin lt 0.6),nfitgd) else fitgd = gd

       if (ngd gt 1) then begin
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
    
; ##################################################
; fit each calibration separately 
    allcalib = ['kk04','t04','m91']
    ncalib = n_elements(allcalib)
    for ii = 0, ncalib-1 do begin
       case ii of
          0: begin
             t04 = 0 & m91 = 0 & kk04 = 1
          end
          1: begin
             t04 = 1 & m91 = 0 & kk04 = 0
          end
          2: begin
             t04 = 0 & m91 = 1 & kk04 = 0
          end
       endcase

; LZ relation evolution; note that the mock code is ok, but the
; results are not needed 
       lzlocal = mrdfits(mzpath+'lzlocal_sdss_ews_'+allcalib[ii]+'.fits.gz',1)
;      lzevol_mock = get_lzevol_mock(lzlocal,qz0=qz0) 
       
       info = mzlz_grab_info(agesohdust,agesancillary,agesmass,$
         t04=t04,m91=m91,kk04=kk04,/nolimit)
       lzevol = mlfit_lzevol(info.mb_ab,info.oh,info.oh_err,$
         info.weight,info.z,q0=0.0,qz0=qz0,lzlocal=lzlocal,$
         quiet=quiet)
       lzevol = struct_addtags({calib: allcalib[ii], ngal: n_elements(info.oh), qz0: qz0},lzevol)
       if (ii eq 0) then alllzevol = lzevol else alllzevol = [alllzevol,lzevol]
       
       lzfile = mzpath+'lzevol_'+allcalib[ii]+'.fits'
       im_mwrfits, lzevol, lzfile, clobber=clobber

;      jj = mrdfits('matchedages_'+calib+'_q2.50-qz0.1.fits.gz',1)
;      jj = mrdfits('matchedages_'+calib+'_q1.50-qz0.1.fits.gz',1)
;      jj = mrdfits('matchedages_'+calib+'_q0.00-qz0.1.fits.gz',1)
;      bb = mlfit_lzevol(jj.mb,jj.oh,jj.oh_err,jj.weight,$
;        jj.z,q0=jj[0].q0,qz0=qz0,lzlocal=lzlocal,$
;        quiet=quiet)
       
; MZ relation evolution
       mzlocal = mrdfits(mzpath+'mzlocal_sdss_ews_'+allcalib[ii]+'.fits.gz',1)
       mzevol = get_mzevol(agesohdust,agesancillary,agesmass,calib=allcalib[ii],$
         t04=t04,m91=m91,kk04=kk04,qz0=qz0,mzlocal=mzlocal,$
         sdssohdust=sdssohdust,sdssancillary=sdssancillary,sdssmass=sdssmass)
       if (ii eq 0) then allmzevol = mzevol else allmzevol = [allmzevol,mzevol]

       mzfile = mzpath+'mzevol_'+allcalib[ii]+'.fits'
       im_mwrfits, mzevol, mzfile, clobber=clobber
    endfor           

; ##################################################
; average M91 and T04, but keep KK04 separate
    zbins = mz_zbins(nz)
    allmassbins = mz_massbins(nallmassbins)

    masscut = 9.6D ; 9.5D
    massbins = mz_massbins(/rev,masscut=masscut,masskeep=masskeep)
    nmassbins = n_elements(massbins)
    struct_print, massbins

; give the average metallicity offset wrt KK04, for the paper
    ohmean = fltarr(ncalib)
    for ii = 0, ncalib-1 do begin
       ww = where(allmzevol[ii].ohmean_bymass gt -900)
       ohmean[ii] = djs_mean((allmzevol[ii].ohmean_bymass)[ww])
    endfor
    niceprint, ohmean, ohmean[0]-ohmean, 10^(ohmean[0]-ohmean)-1
    
    mzavg = {qz0: qz0, dlogohdz_normmass: 10.5, calib: allcalib, $
; linear metallicity evolution
      coeffs: fltarr(2), coeffs_err: fltarr(2), $
      coeffs_all: fltarr(2,ncalib), coeffs_err_all: fltarr(2,ncalib), $

      coeffs_bymass: fltarr(2,nmassbins), coeffs_bymass_err: fltarr(2,nmassbins), $
      coeffs_bymass_all: fltarr(2,nmassbins,ncalib), coeffs_bymass_err_all: fltarr(2,nmassbins,ncalib), $
; change in metallicity wrt z=0.1; combine SDSS+AGES
      medz: fltarr(nz+1), dlogoh: fltarr(nz+1), dlogoh_err: fltarr(nz+1), $
      medz_all: fltarr(nz+1,ncalib), dlogoh_all: fltarr(nz+1,ncalib), dlogoh_err_all: fltarr(nz+1,ncalib), $

      medz_bymass: fltarr(nmassbins,nz+1)-999, medmass_bymass: fltarr(nmassbins,nz+1)-999, $
      medz_bymass_all: fltarr(nmassbins,nz+1,ncalib)-999, medmass_bymass_all: fltarr(nmassbins,nz+1,ncalib)-999, $

      dlogoh_bymass: fltarr(nmassbins,nz+1)-999, dlogoh_bymass_err: fltarr(nmassbins,nz+1)-999, $
      dlogoh_bymass_all: fltarr(nmassbins,nz+1,ncalib)-999, dlogoh_bymass_err_all: fltarr(nmassbins,nz+1,ncalib)-999, $

; metallicity evolution rate vs mass coefficients
      dlogohdz_medmass: fltarr(nmassbins)-999, $
      dlogohdz_medmass_all: fltarr(nmassbins,ncalib)-999, $

      dlogohdz_coeff: fltarr(2), dlogohdz_coeff_err: fltarr(2), dlogohdz_covar: fltarr(2,2), $
      dlogohdz_coeff_all: fltarr(2,ncalib), dlogohdz_coeff_err_all: fltarr(2,ncalib), $

; evolutionary coefficients for both P and R
      mzevol_coeffs: fltarr(5,3), mzevol_coeffs_err: fltarr(5,3), mzevol_chi2: fltarr(3), $
      mzevol_coeffs_r0zero: fltarr(5,3), mzevol_coeffs_err_r0zero: fltarr(5,3), mzevol_chi2_r0zero: fltarr(3), $
      mzevol_coeffs_avg: fltarr(5), mzevol_coeffs_avg_err: fltarr(5), $
      mzevol_coeffs_avg_r0zero: fltarr(5), mzevol_coeffs_avg_err_r0zero: fltarr(5), $
      p0avg: 0.0, p0avg_err: 0.0, r0avg: 0.0, r0avg_err: 0.0, p0r0_chi2: 0.0, $
      p0r0zero_avg: 0.0, p0r0zero_avg_err: 0.0, p0r0zero_chi2: 0.0, $
      
; B-band metallicity and luminosity evolution
      lzevol_coeffs_p0zero: fltarr(4,3), lzevol_coeffs_err_p0zero: fltarr(4,3), $
      lzevol_coeffs_p0fixed: fltarr(4,3), lzevol_coeffs_err_p0fixed: fltarr(4,3)}

;      lz_slope: fltarr(ncalib), lz_slope_err: fltarr(ncalib), $
;;     lz_slope_avg: 0.0, lz_slope_avg_err: 0.0, $
;      lz_s0: fltarr(ncalib), lz_s0_err: fltarr(ncalib), $
;      lz_s0_avg: 0.0, lz_s0_avg_err: 0.0, $
;      dmb_bymass: fltarr(nmassbins,nz)-999.0, dmb_bymass_err: fltarr(nmassbins,nz)-999.0}

; compute the *weighted* average of the coefficients and dlogoh over
; the three calibrations at fixed stellar mass
    for mm = 0, nmassbins-1 do begin
;      mzavg.coeffs_bymass[0,mm] = djs_mean(allmzevol.coeffs_bymass[0,masskeep[mm]])
;      mzavg.coeffs_bymass[1,mm] = djs_mean(allmzevol.coeffs_bymass[1,masskeep[mm]])
;      mzavg.coeffs_bymass_err[0,mm] = djsig(allmzevol.coeffs_bymass[0,masskeep[mm]])/sqrt(3)
;      mzavg.coeffs_bymass_err[1,mm] = djsig(allmzevol.coeffs_bymass[1,masskeep[mm]])/sqrt(3)

       mzavg.coeffs_bymass[0,mm] = im_weighted_mean(allmzevol.coeffs_bymass[0,masskeep[mm]],$
         errors=allmzevol.coeffs_bymass_err[0,masskeep[mm]],wmean_err=wmean0_err,wsigma=wsigma0)
       mzavg.coeffs_bymass[1,mm] = im_weighted_mean(allmzevol.coeffs_bymass[1,masskeep[mm]],$
         errors=allmzevol.coeffs_bymass_err[1,masskeep[mm]],wmean_err=wmean1_err,wsigma=wsigma1)
       mzavg.coeffs_bymass_err[0,mm] = wmean0_err ; wsigma0
       mzavg.coeffs_bymass_err[1,mm] = wmean1_err ; wsigma1
       
       mzavg.coeffs_bymass_all[0,mm,*] = allmzevol.coeffs_bymass[0,masskeep[mm]]
       mzavg.coeffs_bymass_err_all[0,mm,*] = allmzevol.coeffs_bymass_err[0,masskeep[mm]]
       mzavg.coeffs_bymass_all[1,mm,*] = allmzevol.coeffs_bymass[1,masskeep[mm]]
       mzavg.coeffs_bymass_err_all[1,mm,*] = allmzevol.coeffs_bymass_err[1,masskeep[mm]]

; SDSS+AGES: average over the calibrations at fixed stellar mass
       mzavg.medmass_bymass[mm,0] = djs_mean(allmzevol.sdss_medmass_bymass[masskeep[mm]])
       mzavg.medz_bymass[mm,0] = djs_mean(allmzevol.sdss_medz_bymass[masskeep[mm]])

       gd = where(allmzevol.sdss_dlogoh_bymass[masskeep[mm]] gt -900.0,ngd)
       if (ngd ne 0) then begin
          mzavg.dlogoh_bymass[mm,0] = djs_mean(allmzevol[gd].sdss_dlogoh_bymass[masskeep[mm]])
          mzavg.dlogoh_bymass_err[mm,0] = djsig(allmzevol[gd].sdss_dlogoh_bymass[masskeep[mm]])/sqrt(ngd)
       endif

       for iz = 1, nz do begin
          mzavg.medmass_bymass[mm,iz] = djs_mean(allmzevol.medmass_bymass[masskeep[mm],iz-1])
          mzavg.medz_bymass[mm,iz] = djs_mean(allmzevol.medz_bymass[masskeep[mm],iz-1])

          gd = where(allmzevol.dlogoh_bymass[masskeep[mm],iz-1] gt -900.0,ngd)
          if (ngd ne 0) then begin
             mzavg.dlogoh_bymass[mm,iz] = djs_mean(allmzevol[gd].dlogoh_bymass[masskeep[mm],iz-1])
             mzavg.dlogoh_bymass_err[mm,iz] = djsig(allmzevol[gd].dlogoh_bymass[masskeep[mm],iz-1])/sqrt(ngd)
          endif
       endfor 

; combine SDSS+AGES but keep the calibrations separate
       mzavg.medmass_bymass_all[mm,0,*] = allmzevol.sdss_medmass_bymass[masskeep[mm]]
       mzavg.medz_bymass_all[mm,0,*] = allmzevol.sdss_medz_bymass[masskeep[mm]]

       mzavg.dlogoh_bymass_all[mm,0,*] = allmzevol.sdss_dlogoh_bymass[masskeep[mm]]
       mzavg.dlogoh_bymass_err_all[mm,0,*] = allmzevol.sdss_dlogoh_bymass[masskeep[mm]]

       for iz = 1, nz do begin
          mzavg.medmass_bymass_all[mm,iz,*] = allmzevol.medmass_bymass[masskeep[mm],iz-1]
          mzavg.medz_bymass_all[mm,iz,*] = allmzevol.medz_bymass[masskeep[mm],iz-1]

          mzavg.dlogoh_bymass_all[mm,iz,*] = allmzevol.dlogoh_bymass[masskeep[mm],iz-1]
          mzavg.dlogoh_bymass_err_all[mm,iz,*] = allmzevol.dlogoh_bymass[masskeep[mm],iz-1]
       endfor 
    endfor

; now average over all the stellar mass bins at fixed redshift 
    for iz = 0, nz do begin
       gd = where(mzavg.medz_bymass[*,iz] gt -900.0,ngd)
       if (ngd ne 0) then begin
          mzavg.medz[iz] = djs_mean(mzavg.medz_bymass[gd,iz])
          mzavg.dlogoh[iz] = djs_mean(mzavg.dlogoh_bymass[gd,iz])
          mzavg.dlogoh_err[iz] = djsig(mzavg.dlogoh_bymass[gd,iz])
       endif

       for ii = 0, ncalib-1 do begin
          gd = where(mzavg.medz_bymass_all[*,iz,ii] gt -900.0,ngd)
          if (ngd ne 0) then begin
             mzavg.medz_all[iz,ii] = djs_mean(mzavg.medz_bymass_all[gd,iz,ii])
             mzavg.dlogoh_all[iz,ii] = djs_mean(mzavg.dlogoh_bymass_all[gd,iz,ii])
             mzavg.dlogoh_err_all[iz,ii] = djsig(mzavg.dlogoh_bymass_all[gd,iz,ii])
          endif
       endfor
    endfor 

; average over all masses and calibrations; note that the intercept is
; defined to be zero at QZ0
    mzavg.coeffs[1] = djs_mean(mzavg.coeffs_bymass[1,*])
    mzavg.coeffs_err[1] = djsig(mzavg.coeffs_bymass[1,*])

    for ii = 0, ncalib-1 do begin
       mzavg.coeffs_all[1,ii] = djs_mean(mzavg.coeffs_bymass_all[1,*,ii])
       mzavg.coeffs_err_all[1,ii] = djsig(mzavg.coeffs_bymass_all[1,*,ii])
    endfor

; fit the correlation between the rate of metallicity evolution and
; stellar mass
    for mm = 0, nmassbins-1 do begin       
       gd = where(mzavg.medmass_bymass[mm,*] gt -900.0,ngd)
       if (ngd ne 0) then mzavg.dlogohdz_medmass[mm] = $
         djs_mean(mzavg.medmass_bymass[mm,gd])

       for ii = 0, ncalib-1 do begin
          gd = where(mzavg.medmass_bymass_all[mm,*,ii] gt -900.0,ngd)
          if (ngd ne 0) then mzavg.dlogohdz_medmass_all[mm,ii] = $
            djs_mean(mzavg.medmass_bymass_all[mm,gd,ii])
       endfor
    endfor

    mzavg.dlogohdz_coeff = linfit(mzavg.dlogohdz_medmass-mzavg.dlogohdz_normmass,$
      mzavg.coeffs_bymass[1,*],measure_err=mzavg.coeffs_bymass_err[1,*],$
      sigma=coeff_err,covar=covar,chisq=chi2)
    mzavg.dlogohdz_coeff_err = coeff_err
    mzavg.dlogohdz_covar = covar

    for ii = 0, ncalib-1 do begin
       mzavg.dlogohdz_coeff_all[*,ii] = linfit(mzavg.dlogohdz_medmass_all[*,ii]-mzavg.dlogohdz_normmass,$
         mzavg.coeffs_bymass_all[1,*,ii],measure_err=mzavg.coeffs_bymass_err_all[1,*,ii],$
         sigma=coeff_err_all,covar=covar,chisq=chi2)
       mzavg.dlogohdz_coeff_err_all[*,ii] = coeff_err_all
    endfor
    
;    maxis = range(8,12,50)
;    ploterror, mzavg.dlogohdz_medmass, mzavg.coeffs_bymass[1,*], $
;      mzavg.coeffs_bymass_err[1,*], psym=6, yr=[-0.4,-0.05], xsty=3, $
;      ysty=3, symsize=3, /trad, xr=[9.8,11.1]
;    djs_oplot, maxis, poly(maxis-mzavg.dlogohdz_normmass,mzavg.dlogohdz_coeff)
;
;    for ii = 0, ncalib-1 do begin
;       oploterror, mzavg.dlogohdz_medmass_all[*,ii], mzavg.coeffs_bymass_all[1,*,ii], $
;         mzavg.coeffs_bymass_err_all[1,*,ii], psym=7, symsize=3, /trad
;       djs_oplot, maxis, poly(maxis-mzavg.dlogohdz_normmass,mzavg.dlogohdz_coeff_all[*,ii]), line=5
;    endfor
;   niceprint, mzavg.dlogohdz_coeff, mzavg.dlogohdz_coeff_err & print
;   print, mzavg.dlogohdz_coeff_all
    
; the coefficients above give us the shape of the MZ relation as a
; function of redshift; fit that evolution here with a simple model
; that allows for evolution in M* and (O/H)*; basically, we want to
; solve for P and R, both individually and averaged over all three
; calibrations 
    nfake = 2000
    zmin = 0.1 & zmax = 0.9
    minmass = 9.3 & maxmass = 11.5
    zval = randomu(seed,nfake)*(zmax-zmin)+zmin
    mass = randomu(seed,nfake)*(maxmass-minmass)+minmass

    for ii = 0, ncalib-1 do begin
       mzlocal = mrdfits(mzpath+'mzlocal_sdss_ews_'+allcalib[ii]+'.fits.gz',1)

       mass = [[reform(allmzevol[ii].sdss_medmass_bymass,nallmassbins,1)],[allmzevol[ii].medmass_bymass]]
       redshift = [[reform(allmzevol[ii].sdss_medz_bymass,nallmassbins,1)],[allmzevol[ii].medz_bymass]]
       complete = [[fltarr(nallmassbins,1)+1],[allmzevol[ii].complete_bymass]]
       ohmean = [[allmzevol[ii].sdss_ohmean_bymass],[allmzevol[ii].ohmean_bymass]]
       ohmeanerr = [[allmzevol[ii].sdss_ohmean_bymass_err],[allmzevol[ii].ohmean_bymass_err]]

; scale by the reduced chi^2       
       fitgd = where(complete and ohmean gt -900.0,ngd)
       mlfit1 = mlfit_mzevol(mass[fitgd],ohmean[fitgd],ohmeanerr[fitgd],$
         redshift[fitgd]*0+1,redshift[fitgd],qz0=qz0,mzlocal=mzlocal,quiet=0,$
         r0=0.0)
       mzavg.mzevol_coeffs_r0zero[*,ii] = mlfit1.params
       mzavg.mzevol_coeffs_err_r0zero[*,ii] = mlfit1.perror*sqrt(mlfit1.chi2)
       mzavg.mzevol_chi2_r0zero[ii] = mlfit1.chi2

       mlfit1 = mlfit_mzevol(mass[fitgd],ohmean[fitgd],ohmeanerr[fitgd],$
         redshift[fitgd]*0+1,redshift[fitgd],qz0=qz0,mzlocal=mzlocal,quiet=0)
       mzavg.mzevol_coeffs[*,ii] = mlfit1.params
       mzavg.mzevol_coeffs_err[*,ii] = mlfit1.perror*sqrt(mlfit1.chi2)
       mzavg.mzevol_chi2[ii] = mlfit1.chi2
       
;      mlfit1 = mlfit_mzevol(mass[fitgd],ohmean[fitgd],ohmeanerr[fitgd],$
;        redshift[fitgd]*0+1,redshift[fitgd],qz0=qz0,mzlocal=mzlocal,quiet=0)

;     maxis = range(9.5,11.3,50)
;     ploterror, mass[fitgd], ohmean[fitgd], ohmeanerr[fitgd], psym=6, ysty=3, /trad, $
;       xrange=[9.5,11.3], yrange=[8.8,9.2], xsty=3
;     for iz = 0, nz do djs_oplot, maxis, mzevol_func(maxis,mlfit1.params,z=mzavg.medz[iz],qz0=mzavg.qz0)
       
;       mlfit_p0zero1 = mlfit_mzevol(mass,ohmodel,ohmodel_err,p0=0.0,$
;         ohmodel*0+1,zval,qz0=qz0,mzlocal=mzlocal,quiet=0)
;       mlfit_r0zero1 = mlfit_mzevol(mass,ohmodel,ohmodel_err,r0=0.0,$
;         ohmodel*0+1,zval,qz0=qz0,mzlocal=mzlocal,quiet=0)

       
; now get QB from the evolution of the LZ relation       
       lzevol = mrdfits(mzpath+'lzevol_'+allcalib[ii]+'.fits.gz',1)
       mzavg.lzevol_coeffs_p0zero[*,ii] = lzevol.coeff
       mzavg.lzevol_coeffs_err_p0zero[*,ii] = lzevol.coeff_err

; Q = (P-S)/c1       
       mzavg.lzevol_coeffs_p0fixed[0:1,ii] = lzevol.coeff[0:1] ; local LZ
       mzavg.lzevol_coeffs_p0fixed[2,ii] = mzavg.mzevol_coeffs_r0zero[4,ii] ; from MZ
       mzavg.lzevol_coeffs_p0fixed[3,ii] = (mzavg.mzevol_coeffs_r0zero[4,ii]-lzevol.coeff[2])/lzevol.coeff[1] ; =QB
    endfor
    
; average over all calibrations    
    for pp = 0, 4 do begin
       mzavg.mzevol_coeffs_avg[pp] = djs_mean(mzavg.mzevol_coeffs[pp,*])
       mzavg.mzevol_coeffs_avg_err[pp] = djsig(mzavg.mzevol_coeffs[pp,*])
       mzavg.mzevol_coeffs_avg_r0zero[pp] = djs_mean(mzavg.mzevol_coeffs_r0zero[pp,*])
       mzavg.mzevol_coeffs_avg_err_r0zero[pp] = djsig(mzavg.mzevol_coeffs_r0zero[pp,*])
    endfor
    
    mzavg.r0avg = djs_mean(mzavg.mzevol_coeffs[3,*])
    mzavg.r0avg_err = djsig(mzavg.mzevol_coeffs[3,*])
    mzavg.p0avg = djs_mean(mzavg.mzevol_coeffs[4,*])
    mzavg.p0avg_err = djsig(mzavg.mzevol_coeffs[4,*])
    mzavg.p0r0_chi2 = djs_mean(mzavg.mzevol_chi2)

    mzavg.p0r0zero_avg = djs_mean(mzavg.mzevol_coeffs_r0zero[4,*])
    mzavg.p0r0zero_avg_err = djsig(mzavg.mzevol_coeffs_r0zero[4,*])
    mzavg.p0r0zero_chi2 = djs_mean(mzavg.mzevol_chi2_r0zero)

; now deal with the LZ relation
    
;; average metallicity evolution rate via the LZ relation (assuming no
;; luminosity evolution); also compute the average LZ slope
;    mzavg.lz_slope = alllzevol.coeff[1]
;    mzavg.lz_slope_err = alllzevol.coeff_err[1]
;    mzavg.lz_s0 = alllzevol.coeff[2] ; dex per redshift
;    mzavg.lz_s0_err = alllzevol.coeff_err[2]
;
;;   mzavg.lz_slope_avg = djs_mean(mzavg.lz_slope)
;;   mzavg.lz_slope_avg_err = djsig(mzavg.lz_slope)
;    mzavg.lz_s0_avg = djs_mean(mzavg.lz_s0)
;    mzavg.lz_s0_avg_err = djsig(mzavg.lz_s0)

;; compare the amount of metallicity evolution from the MZ and LZ
;; relations and attribute the differences to luminosity evolution
;    for iz = 0, nz-1 do begin
;       dlogoh_lum = mzavg.lz_dlogohdz*(zbins[iz].zbin-qz0)
;       dlogoh_lum_err = mzavg.lz_dlogohdz_err*(zbins[iz].zbin-qz0)
;
;       dlogoh_mass = mzavg.dlogoh_bymass[*,iz]
;       dlogoh_mass_err = mzavg.dlogoh_bymass_err[*,iz]
;       gd = where(dlogoh_mass gt -900.0,ngd)
;
;       numer = dlogoh_mass[gd]-dlogoh_lum
;       numer_err = dlogoh_mass_err[gd]
;;      numer_err = sqrt(dlogoh_mass_err[gd]^2+dlogoh_lum_err^2)
;
;       denom = abs(mzavg.lz_slope)
;       denom_err = mzavg.lz_slope_err
;       
;       mzavg.dmb_bymass[gd,iz] = numer/denom ; [mag/z]
;       mzavg.dmb_bymass_err[gd,iz] = im_compute_error(numer,numer_err,$
;         denom,denom_err,/quotient)
;    endfor

; now go write the paper!    
    mzavgfile = mzpath+'mzevol_avg.fits'
    im_mwrfits, mzavg, mzavgfile, clobber=clobber

return
end
