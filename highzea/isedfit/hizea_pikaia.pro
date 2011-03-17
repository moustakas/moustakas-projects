function hizea_likelihood, n, p, _extra=extra

stop    
    
    mass = p[0]
    tauindx = p[1]
    ageindx = p[2]
;   tauindx = round(p[1])
;   ageindx = round(p[2])

    imodelmaggies = mass*interpolate(extra.modelmaggies,ageindx,tauindx)
;   imodelmaggies = mass*extra.modelmaggies[*,ageindx,tauindx]
    chi2 = total(extra.ivarmaggies*(extra.maggies-imodelmaggies)^2.0)
;   print, chi2, mass*1D10, extra.tau[tauindx], extra.age[ageindx]
    
return, exp(-0.5*chi2)
end

pro hizea_pikaia

    light = 2.99792458D18 ; speed of light [A/s]

    filterlist = ['sdss_'+['u0','g0','r0','i0','z0'],$
      'spitzer_irac_'+['ch1','ch2','ch3','ch4']]+'.par'
    filtinfo = im_filterspecs(filterlist=filterlist,/verbose)
    nfilt = n_elements(filterlist)

    path = '~/home/research/projects/hizea/'
    readcol, path+'SDSSJ231122.58-083923.7.flux', weff, flux, $
      ferr, filter, format='F,F,F,A', /silent, comment='#'
    zobj = 0.73

    flam2maggies = filtinfo.weff^2.0/light*10.0^(0.4*48.6)
    maggies = flux[2:10]*flam2maggies
    errmaggies = ferr[2:10]*flam2maggies
    ivarmaggies = 1.0/errmaggies^2.0

    out = mrdfits(path+'hizea_chi2grid.fits',1)    

; do the traditional minimization

    tau = out.tau
    age = out[0].age
    ntau = n_elements(tau)
    nage = n_elements(age)
    
    finalchi2 = min(out.chi2,chi2indx)
    indx = array_indices([nage,ntau],chi2indx,/dim)
    ageindx = indx[0] & tauindx = indx[1]
    mstar = out[tauindx].mass[ageindx]
    finalage = age[ageindx]
    finaltau = tau[tauindx]

;   niceprint, mstar*out[tauindx].modelmaggies[*,ageindx], maggies

; now minimize using the genetic algorithm

    maggiesnorm = 1D-9
    massnorm = 1D10
    nmaggies = maggies/maggiesnorm
    nivarmaggies = ivarmaggies*maggiesnorm^2.0
    nmodelmaggies = massnorm*out.modelmaggies/maggiesnorm

    nparam = 3
    prange = dblarr(2,nparam)
;   prange[*,0] = [1D5,1D13]/massnorm ; mass
;   prange[*,1] = [0.0,1.0] ; tau index number
;   prange[*,2] = [120,130] ; age index number
    prange[*,0] = [1D11,1D12]/massnorm ; mass
    prange[*,1] = [0.0,ntau-1] ; tau index number
    prange[*,2] = [0.0,nage-1] ; age index number

    functargs = {maggies: nmaggies, ivarmaggies: nivarmaggies, $
      modelmaggies: nmodelmaggies, tau: tau, age: age}

    pikaia, n, ctrl, x, f, status, fname='hizea_likelihood', $
      x=.4,y=.6

    print, x
    
stop    

return
end
