pro priors

    common sirtf_simulations

    alpha = 2.0
    beta = 2.0

    lf = sirtf.lf
    nlf = lf.nlf
    zarray = *sirtf.redshift.zarray
    dz = sirtf.redshift.dz
    nz = n_elements(zarray)

    dvdz = dvcomoving(zarray)/1D18           ; [Mpc^3]
    dvdz[0] = 1.0

    temp = fltarr(nlf,nz)
    lumz = (transpose((1+zarray)^alpha # ((dblarr(nlf)+1)))) * (lf.lum # ((dblarr(nz)+1))) ; [nlf,nz]
    phiz = (transpose((1+zarray)^beta # ((dblarr(nlf)+1)))) * (lf.phi # ((dblarr(nz)+1)))  ; [nlf,nz]

    dlum = alog(10D)*lumz*lf.logbinsz/2.5D ; [nlf,nz]

    for i = 0L, nz-1L do temp[*,i] = phiz[*,i] * dlum[*,i] * dvdz[i] * dz

    plot, zarray, temp[0,*]/total(temp[0,*])
    oplot, zarray, temp[10,*]/total(temp[10,*]), line=2

    plot, lumz[*,0], temp[*,0], /ylog, /xlog, yr=[1E-10,1E6] ; LLF
    oplot, lumz[*,33], temp[*,33]             ; LF z=1

return
end
