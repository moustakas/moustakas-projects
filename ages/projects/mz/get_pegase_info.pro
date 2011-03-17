function get_pegase_info
; jm10nov05ucsd - gather the info we need from the precomputed Pegase
; grid (see build_mz_pegase_models)

    zsun = 0.0134 ; Asplund+09
    log12ohsun = 8.69
    maxzform = 4.0 ; maximum zf allowed
    
    pegpath = ages_path(/projects)+'mz/pegase/'
    taufile = file_search(pegpath+'kroupa_tau_*.fits.gz',count=ntau)

    nage = 100
    info = {$
      tau:              0.0,$
      maxzform:    maxzform,$
      mbaryon:          0.0,$ ; see qaplot_mz_pegase_models
      age:     fltarr(nage),$
      zgas:    fltarr(nage),$
      log12oh: fltarr(nage),$
      mstar:   fltarr(nage),$
      mgas:    fltarr(nage),$
      sfr:     fltarr(nage)}
    info = replicate(info,ntau)

    age = range(0.1,getage(0.0)-getage(maxzform),nage) ; [Gyr]
    for ii = 0, ntau-1 do begin
       peg = mrdfits(taufile[ii],1,/silent)
       findx = findex(peg.age/1D9,age)
       info[ii].tau = peg.tau
       info[ii].age = age
       info[ii].zgas = interpolate(peg.zgas,findx)
       info[ii].log12oh = alog10(info[ii].zgas/zsun)+log12ohsun
       info[ii].mstar = interpolate(peg.mstar,findx)
       info[ii].mgas = interpolate(peg.mgas,findx)
       info[ii].sfr = interpolate(peg.sfr,findx)
    endfor
    info = info[sort(info.tau)]
    info = info[where(info.tau gt 0)]

return, info
end
