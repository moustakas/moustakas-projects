pro hizea_iterproc, func, p, iter, interrupt, functargs=functargs, $
  oref=oref, _extra=extra

    compile_opt idl2,hidden

    oref->get_property, best_parms=best
    best[0] = best[0]*1D10

    print, 'Iteration '+strtrim(iter,2), strjoin(best,' '), $
      strtrim(call_function(func,best,_extra=functargs),2)

return
end
    
function hizea_chi2, p, _extra=extra

    mass = p[0]
    tauindx = p[1]
    ageindx = p[2]
;   tauindx = round(p[1])
;   ageindx = round(p[2])

    imodelmaggies = mass*interpolate(extra.modelmaggies,ageindx,tauindx)
;   imodelmaggies = mass*extra.modelmaggies[*,ageindx,tauindx]
    chi2 = total(extra.ivarmaggies*(extra.maggies-imodelmaggies)^2.0)
;   print, chi2, mass*1D10, extra.tau[tauindx], extra.age[ageindx]
    
return, chi2
end

pro hizea_mass, dogrid=dogrid
; jm09feb02nyu 

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

    niceprint, mstar*out[tauindx].modelmaggies[*,ageindx], maggies

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

    func = 'hizea_chi2'
    ftol = 0.001
    npop = 500
    pmutate = 0.5
    pcross = 0.9
;   stretch_factor = 10.0
    gene_length = 20.0
    itmax = 20
;   iterproc = 'hizea_iterproc'

    p = rmd_ga(ftol,function_value=function_value,function_name=func,$
      prange=prange,ncalls=ncalls,quiet=0,pcross=pcross,gene_length=gene_length,$
      pmutate=pmutate,stretch_factor=stretch_factor,itmax=itmax,iterproc=iterproc,$
      iterargs=iterargs,npop=npop,functargs=functargs,/boltzman)

    print, 'mass', mstar, p[0]*massnorm
    print, 'tau', finaltau, tau[p[1]]
    print, 'age', finalage, age[p[2]]
    print, 'chi2', finalchi2, function_value

    print

;   print, p[0]*massnorm*interpolate(out.modelmaggies,p[2],p[1])
;   print, total(ivarmaggies*(maggies-p[0]*massnorm*interpolate(out.modelmaggies,p[2],p[1]))^2.0)
    
stop    
    
; minimize and make a plot

    chi2min = min(out.chi2,chi2indx)
    indx = array_indices([nage,ntau],chi2indx,/dim)
    ageindx = indx[0]
    tauindx = indx[1]

    mstar = out[tauindx].mass[ageindx]
    niceprint, mstar*out[tauindx].modelmaggies[*,ageindx], maggies

    allsed = mrdfits(taufiles[tauindx],1,/silent)
    sedwave = allsed.wave*(1+zobj)
    get_element, allsed.age, out[0].age, these

;   zfactor = 1.0
    zfactor = 1.0/((1+zobj)*(4.0*!dpi*dluminosity(zobj,/cm)^2))
    sedflam = mstar*lsun*allsed.flux[*,these[ageindx]]*zfactor ; [erg/s/cm2/A]
    
    sedfnu = sedflam*sedwave^2.0/light
    sedmab = -2.5*alog10(sedfnu)-48.6

    dfpsplot, path+'SDSSJ2311_mass.ps', /color, /landscape
    im_plotfaves, /post
    djs_plot, [0], [0], /nodata, xrange=[0.1,20.0], yrange=[21.5,18], $
      /xsty, /ysty, /xlog, ytitle='m_{AB}', xtitle='Observed Wavelength (\mu'+'m)', $
      charsize=1.8
    label = ['\chi^{2} = '+strtrim(string(chi2min,format='(F12.1)'),2),$
      'log (M_{*}/M_{\odot}) = '+string(alog10(mstar),format='(F5.2)'),$
      '\tau = '+string(out[tauindx].tau,format='(F4.1)')+' Gyr',$
      't = '+string(out[tauindx].age[ageindx]/1E6,format='(F5.1)')+' Myr',$
      'Z/Z_{\odot} = 1']
    im_legend, label, /left, /top, box=0, charsize=1.8
    djs_oplot, sedwave/1E4, sedmab, line=0, color='green'
    plotsym, 8, 2.0, fill=1
    oploterror, filtinfo.weff/1E4, mab, filtinfo.fwhm/1E4, errmab, $
      psym=8, errthick=4.0
    im_plotfaves
    dfpsclose

    
    
stop    

return
end
    
