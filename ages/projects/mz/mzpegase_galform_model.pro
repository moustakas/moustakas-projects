function mzpegase_galform_model, zform, alpha=alpha, mass0=mass0, tau0=tau0
; jm10nov05ucsd - generate a very simple galaxy formation model,
; whereby tau is a simple power-law function of the total baryonic
; mass and the formation redshift is a free parameter

; galaxy formation model:    
;   tau = tau0*(M/M0)^alpha

    if (n_elements(zform) eq 0) then zform = 2.0
    if (n_elements(alpha) eq 0) then alpha = -0.8
    if (n_elements(mass0) eq 0) then mass0 = 10D^10.5
    if (n_elements(tau0) eq 0) then tau0 = 10.0

; read the output from mzpegase_build_models()
    pegmodels = mzpegase_read_models()
    zaxis = range(0.0,zform,n_elements(pegmodels[0].age))
    
; generate a uniform grid in mbaryon, but do not extrapolate the
; precomputed model grid
    mbaryon = 10D^range(8.0,13.0,25)
    tau = tau0*(mbaryon/mass0)^alpha
    
    keep = where((tau ge min(pegmodels.tau)) and $
      (tau le max(pegmodels.tau)),ntau)
    model = im_empty_structure(pegmodels,ncopies=ntau)
    model = struct_addtags(temporary(model),replicate($
      {zaxis: zaxis, zform: zform, alpha: alpha, mass0: mass0, $
      tau0: tau0},ntau))
    model.tau = tau[keep]
    model.mbaryon = alog10(mbaryon[keep])
    model.maxzform = pegmodels[0].maxzform

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
    aindx = findex(model[0].age,getage(zaxis)-getage(zform))
    for ii = 0, ntau-1 do begin
       model[ii].age = interpolate(model[ii].age,aindx)
       model[ii].zgas = interpolate(model[ii].zgas,aindx)
       model[ii].log12oh = interpolate(model[ii].log12oh,aindx)
       model[ii].mstar = interpolate(model[ii].mstar,aindx)
       model[ii].mgas = interpolate(model[ii].mgas,aindx)
       model[ii].sfr = interpolate(model[ii].sfr,aindx)
    endfor
    
return, model
end

