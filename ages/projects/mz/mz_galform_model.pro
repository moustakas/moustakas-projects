function mz_galform_model, in_info, alpha=alpha, mass0=mass0, tau0=tau0
; jm10nov05ucsd - generate a very simple galaxy formation model,
; whereby tau is a simple power-law function of the total baryonic
; mass; IN_INFO is the output from get_pegase_info()

; galaxy formation model:    
;   tau = tau0*(M/M0)^alpha

    if (n_elements(alpha) eq 0) then alpha = -0.8
    if (n_elements(mass0) eq 0) then mass0 = 10D^10.5
    if (n_elements(tau0) eq 0) then tau0 = 10.0

; generate a uniform grid in mbaryon, but do not extrapolate the
; precomputed model grid
    mbaryon = 10D^range(8.0,13.0,25)
    tau = tau0*(mbaryon/mass0)^alpha

    keep = where((tau ge min(in_info.tau)) and $
      (tau le max(in_info.tau)),ntau)
    info = im_empty_structure(in_info,ncopies=ntau)
    info = struct_addtags(temporary(info),replicate($
      {alpha: alpha, mass0: mass0, tau0: tau0},ntau))
    info.tau = tau[keep]
    info.mbaryon = alog10(mbaryon[keep])
    info.maxzform = in_info[0].maxzform

    nage = n_elements(info[0].age)
    massscale = 10D^rebin(reform(info.mbaryon,1,ntau),nage,ntau)

    tindx = findex(in_info.tau,info.tau)
    info.age     = interpolate(in_info.age,tindx)
    info.zgas    = interpolate(in_info.zgas,tindx)
    info.log12oh = interpolate(in_info.log12oh,tindx)
    info.sfr     = interpolate(in_info.sfr,tindx)*massscale
    info.mstar   = alog10(interpolate(in_info.mstar,tindx)*massscale)
    info.mgas    = alog10(interpolate(in_info.mgas,tindx)*massscale)

;   plot, alog10(mbaryon), tau0*(mbaryon/mass0)^alpha, /ylog
    
return, info
end

