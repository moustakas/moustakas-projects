function mlfit_lzevol, absmag, oh, oh_err, weight, z, $
  q0=q0, qz0=qz0, lzlocal=lzlocal
; fit the LZ relation with redshift

    if (n_elements(q0) eq 0) then q0 = 1.5 ; luminosity evolution [mag/z]
    
    nparams = 4
    parinfo = {value: 1D, limited: [0,0], limits: [0D,0D], fixed: 0}
    parinfo = replicate(parinfo,nparams)

; fix the first four parameters at their local values
    parinfo[0:1].value = lzlocal.coeff
    parinfo[1].fixed = 1
;   parinfo[0:1].fixed = 1

; P[2]: (O/H)*(z) = (O/H)*(z=0.1) + S*(z-qz0)
    parinfo[2].value = -0.1 ; [dex/z]
;   parinfo[2].limited = 1
;   parinfo[2].limits = [-1.0,0.0]
    
; P[3]: MB*(z) = MB*(z=0.1) + Q*(z-qz0)
    parinfo[3].value = q0 ; [mag/z]
    parinfo[3].fixed = 1  

; do the fit, pack it in, and return
    functargs = {z: z, qz0: qz0, pivotmag: lz_pivotmag()}
    params = mpfitfun('lzevol_func',absmag,oh,oh_err/sqrt(weight),$
      parinfo=parinfo,functargs=functargs,perror=perror,dof=dof,$
      covar=covar,status=mpstatus,quiet=1,bestnorm=chi2)
    fit = {params: params, perror: perror, chi2: chi2, $
      dof: dof};, covar: covar}
    
return, fit
end

function mlfit_mzevol, mass, oh, oh_err, weight, z, $
  r0=r0, qz0=qz0, mzlocal=mzlocal
; fit the MZ relation with redshift

    if (n_elements(r0) eq 0) then r0 = 0.0 ; mass evolution [dex/z]

; closed-box model with linear evolution in (O/H)* and M* (so five
; total parameters)
    nparams = 5
    parinfo = {value: 1D, limited: [0,0], limits: [0D,0D], fixed: 0}
    parinfo = replicate(parinfo,nparams)

; fix the first three parameters at their local values
    parinfo[0:2].value = mzlocal.coeff
    parinfo[0:2].fixed = 1

; P[4]: (O/H)*(z) = (O/H)*(z=0.1) + P*(z-qz0)
;   parinfo[3].value = -0.0 ; [dex/z]
;   parinfo[3].fixed = 1
    parinfo[3].value = -0.1 ; [dex/z]
    parinfo[3].fixed = 0

;   parinfo[3].value = -0.1     ; [dex/z]
;   parinfo[4].limited = 1
;   parinfo[4].limits = [-1.0,0.0]
    
; P[5]: M*(z) = M*(z=0.1) + R*(z-qz0)
    parinfo[4].value = r0 ; [dex/z]
    parinfo[4].fixed = 1    
;   parinfo[4].fixed = 0
;   parinfo[4].limited[1] = 1
;   parinfo[4].limits[1] = 0.0

; do the fit, pack it in, and return
    functargs = {z: z, qz0: qz0}
;   params = mpfitfun('mzevol_func',mass,oh,weight=weight*0.0+1.0,$
;     parinfo=parinfo,functargs=functargs,perror=perror,dof=dof,$
;     covar=covar,status=mpstatus,quiet=0,bestnorm=chi2,yfit=yfit)

;   oh_err1 = sqrt((oh_err/sqrt(weight))^2+0.01^2)
    oh_err1 = oh_err/sqrt(weight)
    params = mpfitfun('mzevol_func',mass,oh,oh_err1,$
      parinfo=parinfo,functargs=functargs,perror=perror,dof=dof,$
      covar=covar,status=mpstatus,quiet=0,bestnorm=chi2,yfit=yfit)
    fit = {params: params, perror: perror, chi2: chi2, $
      dof: dof};, covar: covar}
;   print
    
return, fit
end

pro fit_mzlzevol, clobber=clobber
; jm10may14ucsd - maximum likelihood fit of the metallicity evolution

    qz0 = 0.1 ; reference redshift

; read the local MZ/LZ relations (see FIT_MZLZLOCAL)
    mzpath = ages_path(/projects)+'mz/'
    mzlocal = mrdfits(mzpath+'mzlocal_sdss_closedbox.fits.gz',1)
;   mzlocal = mrdfits(mzpath+'mzlocal_sdss_brokenpl.fits.gz',1)
    lzlocal = mrdfits(mzpath+'lzlocal_sdss.fits.gz',1)
    ncalib = n_elements(lzlocal)

    mzfile = mzpath+'mzevol.fits'
    lzfile = mzpath+'lzevol_B.fits'

    zbins = mz_zbins(nz)
    
; read the data    
    agesancillary = read_mz_sample(/mzhii_ancillary)
    agesmass = read_mz_sample(/mzhii_mass)
    agesohdust = read_mz_sample(/mzhii_log12oh)

; grid of evolutionary parameters
    nmzparams = 5 ; see MLFIT_MZEVOL
    nlzparams = 4 ; see MLFIT_LZEVOL

    nr0grid = 3 ; 9
    nq0grid = 4 ; 6
    r0grid = range(0.0,0.5,nr0grid)
    q0grid = [0.0,range(1.2,2.0,nq0grid-1)]

    mzevol = replicate({calib: '', ngal: 0, qz0: qz0, $
      ohstar: fltarr(nz,nr0grid), ohstar_err: fltarr(nz,nr0grid), $
      dlogoh: fltarr(nz,nr0grid), $ ; dlogoh_err: fltarr(nz,nr0grid), $
      dlogoh_avg: fltarr(nz,nr0grid), dlogoh_avg_err: fltarr(nz,nr0grid), $
      coeffs: fltarr(nmzparams,nr0grid), $
      pavg: fltarr(nr0grid), pavg_err: fltarr(nr0grid)},ncalib)
    lzevol = replicate({calib: '', ngal: 0, qz0: qz0, $
      ohstar: fltarr(nz,nq0grid), ohstar_err: fltarr(nz,nq0grid), $
      dlogoh: fltarr(nz,nq0grid), $ ; dlogoh_err: fltarr(nz,nq0grid), $
      dlogoh_avg: fltarr(nz,nq0grid), dlogoh_avg_err: fltarr(nz,nq0grid), $
      coeffs: fltarr(nlzparams,nq0grid), $
      savg: fltarr(nq0grid), savg_err: fltarr(nq0grid)},ncalib)

    mzevol.calib = mzlocal.calib
    lzevol.calib = lzlocal.calib
    
; fit each calibration separately
;splog, 'Testing!!'
;   for ii = 2, ncalib-1 do begin
    for ii = 0, ncalib-1 do begin
       t04 = 0 & m91 = 0 & kk04 = 0
       case ii of
          0: t04 = 1
          1: m91 = 1
          2: kk04 = 1
       endcase
       if keyword_set(t04) then calib = 't04'
       if keyword_set(m91) then calib = 'm91'
       if keyword_set(kk04) then calib = 'kk04'

       mzlocal1 = mzlocal[where(strtrim(mzlocal.calib,2) eq calib)]
       lzlocal1 = lzlocal[where(strtrim(lzlocal.calib,2) eq calib)]
       info = mzlz_grab_info(agesohdust,agesancillary,agesmass,$
         flux=flux,t04=t04,m91=m91,kk04=kk04);,zmin=0.15)
       mzevol[ii].ngal = info.ngal
       lzevol[ii].ngal = info.ngal
       
; MZ relation
       for jj = 0, nr0grid-1 do begin
          mz = mlfit_mzevol(info.mass,info.oh,info.oh_err,$
            info.weight,info.z,r0=r0grid[jj],qz0=qz0,mzlocal=mzlocal1)
          mzevol[ii].coeffs[*,jj] = mz.params
          noevol_params = mz.params*[1,1,1,1,0.0,1]
; compute the mean offset of the data from the evolving MZ model
          for iz = 0, nz-1 do begin
             zinfo = mzlz_grab_info(agesohdust,agesancillary,agesmass,$
               flux=flux,t04=t04,m91=m91,kk04=kk04,zmin=zbins[iz].zlo,$
               zmax=zbins[iz].zup)
; (O/H) at log(M/Msun)=10.5
             ohmodel_noevol = mzevol_func(zinfo.mass,noevol_params,z=zinfo.z,qz0=qz0)
             mzevol[ii].dlogoh[iz,jj] = im_weighted_mean(zinfo.oh-ohmodel_noevol,$
               errors=zinfo.oh_err/sqrt(zinfo.weight),wmean_err=wmean_err)
             mzevol[ii].ohstar[iz,jj] = mzevol[ii].dlogoh[iz,jj] + $
               mzevol_func(mz_pivotmass(),noevol_params,z=zbins[iz].zbin,qz0=qz0) ; O/H @10.5
             mzevol[ii].ohstar_err[iz,jj] = wmean_err
          endfor
; adjust dlogoh to be relative to the first redshift bin
          mzevol[ii].dlogoh[*,jj] = mzevol[ii].dlogoh[*,jj]-mzevol[ii].dlogoh[0,jj]
       endfor

;      mzevol[ii].ohstar[iz,jj] = weighted_quantile(zinfo.oh-ohmodel_noevol,zinfo.weight)+ $
;        mzevol_func(mz_pivotmass(),noevol_params,z=zbins[iz].zbin,qz0=qz0) ; O/H @10.5

; LZ relation       
       for jj = 0, nq0grid-1 do begin
          lz = mlfit_lzevol(info.mb_ab,info.oh,info.oh_err,$
            info.weight,info.z,q0=q0grid[jj],qz0=qz0,lzlocal=lzlocal1)
          lzevol[ii].coeffs[*,jj] = lz.params
          noevol_params = lz.params*[1,1,0,1]
;         splog, lz.params
; compute the mean offset of the data from the evolving LZ model
          for iz = 0, nz-1 do begin
             zinfo = mzlz_grab_info(agesohdust,agesancillary,agesmass,$
               flux=flux,t04=t04,m91=m91,kk04=kk04,zmin=zbins[iz].zlo,$
               zmax=zbins[iz].zup)
             ohmodel_noevol = lzevol_func(zinfo.mb_ab,noevol_params,z=zinfo.z,$
               qz0=qz0,pivotmag=lz_pivotmag())
             lzevol[ii].dlogoh[iz,jj] = im_weighted_mean(zinfo.oh-ohmodel_noevol,$
               errors=zinfo.oh_err/sqrt(zinfo.weight),wmean_err=wmean_err)
             lzevol[ii].ohstar[iz,jj] = lzevol[ii].dlogoh[iz,jj] + $
               lzevol_func(lz_pivotmag(),noevol_params,z=zbins[iz].zbin,$
               qz0=qz0,pivotmag=lz_pivotmag()) ; O/H @-20.5
             lzevol[ii].ohstar_err[iz,jj] = wmean_err
; if (jj eq 2) and (iz eq 4) then stop
;            print, q0grid[jj], mean(zinfo.oh-ohmodel_noevol), $
;              lzevol_func(lz_pivotmag(),noevol_params,z=zbins[iz].zbin,$
;              qz0=qz0,pivotmag=lz_pivotmag()), mean(zinfo.oh-ohmodel_noevol)+$
;              lzevol_func(lz_pivotmag(),noevol_params,z=zbins[iz].zbin,$
;              qz0=qz0,pivotmag=lz_pivotmag())
          endfor
; adjust dlogoh to be relative to the first redshift bin
          lzevol[ii].dlogoh[*,jj] = lzevol[ii].dlogoh[*,jj]-lzevol[ii].dlogoh[0,jj]
       endfor 
    endfor          

; compute the average dlogoh across all three calibrations;
; unfortunately we'll have three copies of the same numbers, but it's
; easier than making a whole new structure
    for ii = 0, ncalib-1 do for jj = 0, nr0grid-1 do for iz = 0, nz-1 do begin
       mzevol[ii].dlogoh_avg[iz,jj] = djs_mean(mzevol.dlogoh[iz,jj])
       mzevol[ii].dlogoh_avg_err[iz,jj] = djsig(mzevol.dlogoh[iz,jj])
    endfor
    for ii = 0, ncalib-1 do for jj = 0, nq0grid-1 do for iz = 0, nz-1 do begin
       lzevol[ii].dlogoh_avg[iz,jj] = djs_mean(lzevol.dlogoh[iz,jj])
       lzevol[ii].dlogoh_avg_err[iz,jj] = djsig(lzevol.dlogoh[iz,jj])
    endfor

; finally also compute the mean evolution *rate* based on the MZ and
; LZ relations, averaging across the three calibrations
    for ii = 0, ncalib-1 do for jj = 0, nr0grid-1 do begin
       mzevol[ii].pavg[jj] = djs_mean(mzevol.coeffs[4,jj]) ; mean metallicity evolution
       mzevol[ii].pavg_err[jj] = djsig(mzevol.coeffs[4,jj])
    endfor
    for ii = 0, ncalib-1 do for jj = 0, nq0grid-1 do begin
       lzevol[ii].savg[jj] = djs_mean(lzevol.coeffs[2,jj]) ; mean metallicity evolution
       lzevol[ii].savg_err[jj] = djsig(lzevol.coeffs[2,jj])
    endfor    
    
; write out    
    im_mwrfits, mzevol, mzfile, clobber=clobber
    im_mwrfits, lzevol, lzfile, clobber=clobber

return
end
