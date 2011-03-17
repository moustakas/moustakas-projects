pro munoz_ubvri

    ss = rsex('table6.sex')
    munoz_to_maggies, ss, maggies, ivarmaggies, filterlist=filterlist
    nfilt = n_elements(filterlist)
    nobj = n_elements(ss)
    
    redshift = 70.0*replicate(15.0,nobj)/2.99D5

    kcorr = im_kcorrect(redshift,maggies,ivarmaggies,$
      filterlist,bessell_filterlist(),chi2=chi2,$
      coeffs=coeffs,mass=mass,bestmaggies=bestmaggies)

    res = replicate({z: 0.0, chi2: 0.0, mass: 0.0, $
      maggies: fltarr(nfilt), ivarmaggies: fltarr(nfilt), $
      bestmaggies: fltarr(nfilt), coeffs: fltarr(5)},nobj)
    res.maggies = maggies
    res.ivarmaggies = ivarmaggies
    res.z = redshift
    res.chi2 = chi2
    res.mass = alog10(mass)
    res.maggies = maggies
    res.ivarmaggies = ivarmaggies
    res.bestmaggies = bestmaggies
    res.coeffs = coeffs

    psfile = 'sings.ps'
    kcorrect_qaplot, res, filterlist, psfile=psfile, /clobber

stop    

    
return
end
