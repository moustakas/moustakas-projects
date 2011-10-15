function mzpegase_galform_model, zobs, alpha=alpha, beta=beta, mass0=mass0, zform0=zform0, tau0=tau0
; jm10nov05ucsd - generate a very simple galaxy formation model,
; whereby tau is a simple power-law function of the total baryonic
; mass and the formation redshift is a free parameter

    if (n_elements(zobs) eq 0) then zobs = 0.1
    
; galaxy formation model:
;   tau = tau0*(M/M0)^alpha
;   (1+zf) = (1+zf,0)(M/M0)^beta

; fiducial normalization values    
    if (n_elements(zform0) eq 0) then zform0 = 2.0
    if (n_elements(mass0) eq 0) then mass0 = 10.0^10.5 ; [Msun]
    if (n_elements(tau0) eq 0) then tau0 = 10.0 ; [Gyr]

; fiducial parameter values    
    if (n_elements(alpha) eq 0) then alpha = -1.0
    if (n_elements(beta) eq 0) then beta = +0.3
    
; read the output from mzpegase_build_models()
    pegmodels = mzpegase_read_models()

; generate a uniform grid in mbaryon, but do not extrapolate the
; precomputed model grid
    mbaryon = 10D^range(8.0,13.0,25)
    tau = tau0*(mbaryon/mass0)^alpha
    zform = (1+zform0)*(mbaryon/mass0)^beta-1
    
    keep = where((tau ge min(pegmodels.tau)) and $
      (tau le max(pegmodels.tau)),ntau)
    model = replicate({zobs: zobs, alpha: alpha, beta: beta, $
      mass0: mass0, zform0: zform0, tau0: tau0, tau: 0.0, zform: 0.0, $
      mbaryon: 0.0, mstar: 0.0, mgas: 0.0, sfr: 0.0, sfrm: 0.0, $
      zgas: 0.0, log12oh: 0.0},ntau)
    model.tau = tau[keep]
    model.zform = zform[keep]
    model.mbaryon = alog10(mbaryon[keep])

; get the physical quantites of interest at ZOBS by applying the
; mass-dependent formation redshift
stop
    for ii = 0, ntau-1 do begin
       aindx = findex(model[0].age,getage(zaxis)-getage(zform))

       
       model[ii].age = interpolate(model[ii].age,aindx)
       model[ii].zgas = interpolate(model[ii].zgas,aindx)
       model[ii].log12oh = interpolate(model[ii].log12oh,aindx)
       model[ii].mstar = interpolate(model[ii].mstar,aindx)
       model[ii].mgas = interpolate(model[ii].mgas,aindx)
       model[ii].sfr = interpolate(model[ii].sfr,aindx)
    endfor
    

    
stop    
    
    zaxis = range(0.0,zform,n_elements(pegmodels[0].age))
    
    nage = n_elements(model[0].age)
    massscale = 10D^rebin(reform(model.mbaryon,1,ntau),nage,ntau)

    tindx = findex(pegmodels.tau,model.tau)
    model.age     = interpolate(pegmodels.age,tindx)
    model.zgas    = interpolate(pegmodels.zgas,tindx)
    model.log12oh = interpolate(pegmodels.log12oh,tindx)
    model.sfr     = interpolate(pegmodels.sfr,tindx)*massscale
    model.mstar   = alog10(interpolate(pegmodels.mstar,tindx)*massscale)
    model.mgas    = alog10(interpolate(pegmodels.mgas,tindx)*massscale)

;   plot, alog10(mbaryon), tau0*(mbaryon/mass0)^alpha, /ylog

; remap the ages onto redshift, given zform

    
return, model
end

