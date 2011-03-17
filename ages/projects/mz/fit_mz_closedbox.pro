function fit_mz_closedbox, mass, oh, weight, oh_err=oh_err, $
  binsize=binsize, minmass=minmass, maxmass=maxmass, $
  mingal=mingal, fit_minmass=fit_minmass, verbose=verbose
; jm10nov09ucsd - fit the MZ relation using a model motivated by
; closed-box chemical evolution
    
    ngal = n_elements(mass)
    if (n_elements(binsize) eq 0) then binsize = 0.1
    if (n_elements(minmass) eq 0) then minmass = min(mass)
    if (n_elements(maxmass) eq 0) then maxmass = max(mass)
    if (n_elements(mingal) eq 0) then mingal = 100
    if (n_elements(fit_minmass) eq 0) then fit_minmass = 9.0

; see MZ_CLOSEDBOX() for the model
    nparams = 3 ; 4
    parinfo = {value: 1D, limited: [0,0], limits: [0D,0D], fixed: 0}
    parinfo = replicate(parinfo,nparams)
    parinfo.value = [9.0,10.5,1.0] ; initial guess
;   parinfo.value = [9.0,10.5,1.0,1.0] ; initial guess

;   parinfo[1].limited[0] = 1 ; all positive
;   parinfo[1].limits[0] = 0.0
    
; bin the MZ relation; let IM_MEDXBIN() deal with WEIGHT undefined 
    mzbin = im_medxbin(mass,oh,binsize,weight=weight,$
      minx=minmass,maxx=maxmass,minpts=mingal,verbose=verbose)

; fit to the full dataset *and* to the binned points
;   splog, 'Fitting binned points'
    coeff = mpfitfun('mz_closedbox',mzbin.xbin,mzbin.medy,$
      weight=(mzbin.xbin gt fit_minmass)/mzbin.sigymean^2,$
      parinfo=parinfo,functargs=functargs,perror=coeff_err,$
      dof=dof,covar=covar,status=mpstatus,quiet=1,bestnorm=chi2,$
      yfit=yfit)
    chi2 = chi2/float(dof)
    scatter = djsig(oh-mz_closedbox(mass,coeff))

;   djs_plot, mzbin.xbin, mzbin.medy, psym=-6, ysty=3, yr=[8.4,9.5]
;   djs_oplot, mzbin.xbin, yfit, psym=6, color='red'
;   cc = get_kbrd(1)
    
;   splog, 'Fitting all points'    
    if (n_elements(oh_err) ne ngal) then message, 'Need uncertainties on O/H'
    if (n_elements(weight) eq 0L) then weight = mass*0.0+1.0
    weighted_oh_err = oh_err/sqrt(weight)
    coeff_all = mpfitfun('mz_closedbox',mass,oh,$
      weight=(mass gt fit_minmass)/weighted_oh_err^2,$
      parinfo=parinfo,functargs=functargs,perror=coeff_err_all,$
      dof=dof_all,covar=covar_all,status=mpstatus_all,quiet=1,$
      bestnorm=chi2_all,yfit=yfit_all)
    chi2_all = chi2_all/float(dof_all)
    scatter_all = djsig(oh-mz_closedbox(mass,coeff_all))

;    djs_plot, mass, oh, psym=3, ysty=3, yr=[8.4,9.3]
;    djs_oplot, mzbin.xbin, mzbin.medy, psym=6, color='red', sym=2
;    djs_oplot, mzbin.xbin, mzbin.meany, psym=6, color='dark green', sym=2
;    mm = range(8.0,11.2,30)
;    djs_oplot, mm, mz_closedbox(mm,coeff), color='orange', thick=4
;    djs_oplot, mm, mz_closedbox(mm,coeff_all), color='blue', thick=4
;    splog, coeff
;    splog, coeff_all
;;   cc = get_kbrd(1)

    fit = {$
      ngal:          ngal, $
      coeff:         coeff, $
      coeff_err:     coeff_err, $
      scatter:       scatter, $
      chi2:          chi2, $
      coeff_all:     coeff_all, $
      coeff_all_err: coeff_err_all, $
      scatter_all:   scatter_all, $
      chi2_all:      chi2_all, $
      bin_mass:      mzbin.xbin, $
      bin_oh:        mzbin.medy,$
      bin_oh_err:    mzbin.sigymean}
    
return, fit
end

