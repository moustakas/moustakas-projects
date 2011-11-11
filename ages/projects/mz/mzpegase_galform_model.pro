function mzpegase_galform_model, zobs, alpha=alpha, beta=beta, mass0=mass0, $
  zform0=zform0, tau0=tau0
; jm10nov05ucsd - generate a very simple galaxy formation model,
; whereby tau is a simple power-law function of the total baryonic
; mass and the formation redshift is a free parameter

    if (n_elements(zobs) eq 0) then zobs = 0.1
    
; galaxy formation model:
;   tau = tau0*(M/M0)^alpha
;   (1+zf) = (1+zf,0)(M/M0)^beta

; fiducial normalization values    
    if (n_elements(zform0) eq 0) then zform0 = 2.0
    if (n_elements(mass0) eq 0) then mass0 = 10.0^11 ; [Msun]
    if (n_elements(tau0) eq 0) then tau0 = 3.0 ; [Gyr]

; fiducial parameter values    
    if (n_elements(alpha) eq 0) then alpha = -0.5
    if (n_elements(beta) eq 0) then beta = +0.2
;   if (beta lt 0.2) then splog, 'Beware of beta<0.2!'

; read the output from mzpegase_build_models()
    pegmodels = mzpegase_read_models()
    nage = n_elements(pegmodels[0].age)

; generate a uniform grid in tau, but do not extrapolate the
; precomputed model grid
    bigmbaryon = 10D^range(9.5,11.5,50)
    bigtau = tau0*(bigmbaryon/mass0)^alpha

    ntau = 50 ; number of output models
    mintau = min(bigtau>min(pegmodels.tau))
    maxtau = max(bigtau<max(pegmodels.tau))
    tau = range(mintau,maxtau,ntau,/log)
    mbaryon = mass0*(tau/tau0)^(1.0/alpha)
    
    zform = zform0*(mbaryon/mass0)^beta
    flag = where(zform le zobs)
    if flag[0] ne -1 then message, 'Too young - zform<zobs!'

; now pack the output structure    
    model = replicate({alpha: alpha, beta: beta, mass0: mass0, zform0: zform0, $
      tau0: tau0, tau: 0.0, zform: 0.0, zobs: zobs, age: 0.0, $
      mbaryon: 0.0, mstar: 0.0, mgas: 0.0, sfr: 0.0, sfrm: 0.0, $
      zgas: 0.0, log12oh: 0.0},ntau)
    model.tau = tau
    model.zform = zform
    model.mbaryon = alog10(mbaryon)
    model.age = getage(zobs)-getage(model.zform) ; age at z=zobs

    tauindx = findex(pegmodels.tau,model.tau)
    ageindx = findex(pegmodels[0].age,model.age)

; scale by the total baryon mass
    for ii = 0, ntau-1 do begin
       model[ii].zgas = interpolate(interpolate(pegmodels.zgas,tauindx[ii]),ageindx[ii])
       model[ii].log12oh = interpolate(interpolate(pegmodels.log12oh,tauindx[ii]),ageindx[ii])
       model[ii].mstar = interpolate(interpolate(pegmodels.mstar,tauindx[ii]),ageindx[ii])
       model[ii].mgas = interpolate(interpolate(pegmodels.mgas,tauindx[ii]),ageindx[ii])
       model[ii].sfr = interpolate(interpolate(pegmodels.sfr,tauindx[ii]),ageindx[ii])
    endfor

    mscale = 10.0^model.mbaryon
    model.mstar = alog10(model.mstar*mscale)
    model.mgas = alog10(model.mgas*mscale)
    model.sfr = alog10(model.sfr*mscale)
    
return, model
end

