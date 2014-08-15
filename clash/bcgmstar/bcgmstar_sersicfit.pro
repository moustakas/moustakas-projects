pro bcgmstar_sersic_allbands, rr, sb, wave, scoeff, sb_ivar=sb_ivar, $
  init_params=init_params, fixed=fixed, sersicfit=sersicfit, $
  fixdevac=fixdevac, verbose=verbose, results=results, lambda_ref=lambda_ref, $
  alphabetazero=alphabetazero
; fit a single-Sersic function to all the bands simultaneously by
; allowing the half-light radius and Sersic n parameter to vary as a
; power-law function of wavelength, while allowing the surface
; brightness at re in each band to be free

; parse the data
    xx = rr
    xxsb = -2.5*alog10(sb)
    xxsb_ivar = sb_ivar*(alog(10)*sb/2.5)^2.0
;   xxsb = sb
;   xxsb_ivar = sb_ivar

; figure out how many bands we need to fit and then set up the
; parameters 
    uwave = wave[uniq(wave,reverse(sort(wave)))] ; you must reverse!
;   uwave = wave[uniq(wave,sort(wave))]
    nband = n_elements(uwave)

    get_element, uwave, lambda_ref, this
    wave_ref = uwave[this] ; =F160W
    nparam = nband+4

    parinfo = replicate({parname: '', value: 0D, fixed: 0, $
      limited: [0,0], limits: [0D,0D]},nparam)

; 1<n_ref<10
    parinfo[0].parname = 'n_ref'
    parinfo[0].limited = [1,1]
    parinfo[0].limits = [0.1D,10D]
    parinfo[0].value = 4D
;   parinfo[0].limits = alog10([0.1D,10D])
;   parinfo[0].value = alog10(4D)
; 0.01<re_ref<500
    parinfo[1].parname = 're_ref'
    parinfo[1].limited = [1,1]
    parinfo[1].limits = [0.01D,500D]
    parinfo[1].value = 20D ; kpc
;   parinfo[1].limits = alog10([0.01D,500D])
;   parinfo[1].value = alog10(20D) ; kpc
; alpha1
    parinfo[2].parname = 'alpha1'
    parinfo[2].limited = [1,1]
    parinfo[2].limits = [-2,1]
    parinfo[2].value = 0.1D
    if keyword_set(alphabetazero) then begin
       parinfo[2].value = 0D
       parinfo[2].fixed = 1
    endif
; beta1
    parinfo[3].parname = 'beta1'
    parinfo[3].limited = [1,1]
    parinfo[3].limits = [-2,1.0]
    parinfo[3].value = 0.1D
    if keyword_set(alphabetazero) then begin
       parinfo[3].value = 0D
       parinfo[3].fixed = 1
    endif
    
; sbe - surface brightness at re_ref
    for ib = 0, nband-1 do begin
       parinfo[4+ib].parname = 'sbe_'+strtrim(uwave[ib],2)
       parinfo[4+ib].limited = [1,1]
       parinfo[4+ib].limits = [10D,35D]
;      parinfo[4+ib].limited = [1,0]
;      parinfo[4+ib].limits = [1D-14,0D]
;      parinfo[4+ib].limits = [0D,0D]
       these = where(uwave[ib] eq wave)
       parinfo[4+ib].value = median(xxsb[these])
    endfor

    if n_elements(fixed) eq nparam then parinfo.fixed = fixed
    if n_elements(init_params) eq nparam then parinfo.value = init_params

; fix the Sersic model to be a de Vaucouleurs profile at all
; wavelengths
    if keyword_set(fixdevac) then begin
       parinfo[0].value = 4D ; n_ref=4 is fixed
;      parinfo[0].value = alog10(4D) ; n_ref=4 is fixed
       parinfo[2].value = 0D ; alpha1 is fixed (n=constant)
       parinfo[0].fixed = 1
       parinfo[2].fixed = 1
    endif
;   struct_print, parinfo

    params = mpfitfun('bcgmstar_sersic_allbands_func',xx,xxsb,parinfo=parinfo,$
      yfit=sersicfit,perror=perror,covar=covar,weights=xxsb_ivar,dof=dof,$
      bestnorm=chi2,status=status,niter=niter,quiet=keyword_set(verbose) eq 0,$
;     xtol=1D-12,$
      functargs={parinfo: parinfo, wave: wave, lambda_ref: lambda_ref})
    factor = sqrt(chi2/dof)

;   djs_plot, xx, xxsb, psym=8, /xlog, xsty=3, ysty=3
;   djs_oplot, xx, sersicfit, psym=6, color='orange'

; build two output structures: one that has the formal fitting results
; (alpha, beta, etc.)...
    results = {$
      sersic_niter:  niter,$
      sersic_status: status,$
      sersic_chi2:     chi2,$
      sersic_dof:       dof,$
      sersic_covar:   covar,$
      sersic_params:  params,$
      sersic_perror:  perror,$

      sersic_wave:     uwave,$
      sersic_wave_ref: wave_ref,$
      
      sersic_n_ref:  params[0],$
      sersic_re_ref: params[1],$
;     sersic_n_ref:  10D^params[0],$
;     sersic_re_ref: 10D^params[1],$
      sersic_alpha1:  params[2],$
      sersic_beta1:   params[3],$

      sersic_n_ref_err:  perror[0],$
      sersic_re_ref_err: perror[1],$
;     sersic_n_ref_err:  perror[0]*10D^params[0]*alog(10),$
;     sersic_re_ref_err: perror[1]*10D^params[1]*alog(10),$
      sersic_alpha1_err:  perror[2],$
      sersic_beta1_err:   perror[3],$

      sersic_sbe:       params[4:4+nband-1],$
      sersic_sbe_err:   perror[4:4+nband-1]}
    
; ...and a second one that has the model evaluated for each filter and
; which takes into account all the covariances in the parameters
    for ib = 0, nband-1 do begin
       n = results.sersic_n_ref*(uwave[ib]/wave_ref)^results.sersic_alpha1
       re = results.sersic_re_ref*(uwave[ib]/wave_ref)^results.sersic_beta1

; these uncertainties are only approximate because they ignore the
; covariance between the parameters (note that I can't use
; MPRANDOMN because the covariance matrix is not positive-definite)
       n_err = sqrt((results.sersic_n_ref_err/results.sersic_n_ref/alog(10))^2 + $
         (results.sersic_alpha1_err*(uwave[ib]/wave_ref))^2)*alog(10)*n
       re_err = sqrt((results.sersic_re_ref_err/results.sersic_re_ref/alog(10))^2 + $
         (results.sersic_beta1_err*(uwave[ib]/wave_ref))^2)*alog(10)*re
       
       sbe = params[4+ib]
       sbe_err = perror[4+ib]

       scoeff1 = {$
         wave: uwave[ib], $
         sersic_all_sbe:     sbe,$
         sersic_all_re:      re,$
         sersic_all_n:       n,$
         sersic_all_sbe_err: sbe_err,$
         sersic_all_re_err:  re_err,$
         sersic_all_n_err:   n_err}
       if ib eq 0 then scoeff = scoeff1 else scoeff = [scoeff,scoeff1]
    endfor

    struct_print, scoeff
;   help, results, /str

return
end

pro bcgmstar_sersic2_allbands, rr, sb, wave, scoeff, sb_ivar=sb_ivar, $
  init_params=init_params, fixed=fixed, sersicfit=sersicfit, $
  fixdevac=fixdevac, verbose=verbose, results=results, lambda_ref=lambda_ref, $
  alphabetazero=alphabetazero
; fit a double-Sersic function to all the bands simultaneously by
; allowing the half-light radius and Sersic n parameter to vary as a
; power-law function of wavelength, while allowing the surface
; brightness at re in each band to be free 

; parse the data
    xx = rr
    xxsb = -2.5*alog10(sb)
    xxsb_ivar = sb_ivar*(alog(10)*sb/2.5)^2.0
;   xx = rr
;   xxsb = sb
;   xxsb_ivar = sb_ivar

; figure out how many bands we need to fit and then set up the
; parameters 
    uwave = wave[uniq(wave,reverse(sort(wave)))] ; you must reverse!
;   uwave = wave[uniq(wave,sort(wave))]
    nband = n_elements(uwave)

    get_element, uwave, lambda_ref, this
    wave_ref = uwave[this] ; =F160W
    nparam = 2*nband+8

    parinfo = replicate({parname: '', value: 0D, fixed: 0, $
      limited: [0,0], limits: [0D,0D]},nparam)

; 1<n1_ref<10
    parinfo[0].parname = 'n1_ref'
    parinfo[0].limited = [1,1]
    parinfo[0].limits = [0.1D,10D]
    parinfo[0].value = 1D
; 1<n2_ref<10
    parinfo[1].parname = 'n2_ref'
    parinfo[1].limited = [1,1]
    parinfo[1].limits = [0.1D,10D]
    parinfo[1].value = 4D
; 0.01<re1_ref<500
    parinfo[2].parname = 're1_ref'
    parinfo[2].limited = [1,1]
;   parinfo[2].limits = [0.1D,20D]
    parinfo[2].limits = [0.1D,500D]
    parinfo[2].value = 10D ; kpc
; 0.01<re2_ref<500
    parinfo[3].parname = 're2_ref'
    parinfo[3].limited = [1,1]
    parinfo[3].limits = [0.1D,500D]
    parinfo[3].value = 100D ; kpc

; alpha1
    parinfo[4].parname = 'alpha1'
    parinfo[4].limited = [1,1]
    parinfo[4].limits = [-1,1]
    parinfo[4].value = 0.1D
    if keyword_set(alphabetazero) then begin
       parinfo[4].value = 0D
       parinfo[4].fixed = 1
    endif
; alpha2
    parinfo[5].parname = 'alpha2'
    parinfo[5].limited = [1,1]
    parinfo[5].limits = [-1,1]
    parinfo[5].value = 0.1D
    if keyword_set(alphabetazero) then begin
       parinfo[5].value = 0D
       parinfo[5].fixed = 1
    endif
; beta1
    parinfo[6].parname = 'beta1'
    parinfo[6].limited = [1,1]
    parinfo[6].limits = [-1,1]
    parinfo[6].value = 0.1D
    if keyword_set(alphabetazero) then begin
       parinfo[6].value = 0D
       parinfo[6].fixed = 1
    endif
; beta2
    parinfo[7].parname = 'beta2'
    parinfo[7].limited = [1,1]
    parinfo[7].limits = [-1,1]
    parinfo[7].value = 0.1D
    if keyword_set(alphabetazero) then begin
       parinfo[7].value = 0D
       parinfo[7].fixed = 1
    endif
    
; sbe1 - surface brightness at re1
    for ib = 0, nband-1 do begin
       parinfo[8+ib].parname = 'sbe1_'+strtrim(uwave[ib],2)
       parinfo[8+ib].limited = [1,1]
       parinfo[8+ib].limits = [10D,35D]
;      parinfo[8+ib].limited = [1,0]
;      parinfo[8+ib].limits = [1D-14,0D]
;      parinfo[8+ib].limits = [0D,0D]
       these = where(uwave[ib] eq wave)
       parinfo[8+ib].value = median(xxsb[these])
    endfor

; sbe2 - surface brightness at re2
    for ib = 0, nband-1 do begin
       parinfo[8+nband+ib].parname = 'sbe2_'+strtrim(uwave[ib],2)
       parinfo[8+nband+ib].limited = [1,1]
       parinfo[8+nband+ib].limits = [10D,35D]
;      parinfo[8+nband+ib].limited = [1,0]
;      parinfo[8+nband+ib].limits = [1D-15,0D]
;      parinfo[8+nband+ib].limits = [0D,0D]
       these = where(uwave[ib] eq wave)
       parinfo[8+nband+ib].value = median(xxsb[these])
    endfor

    if n_elements(fixed) eq nparam then parinfo.fixed = fixed
    if n_elements(init_params) eq nparam then parinfo.value = init_params

;; fix the two components to be deVac with no wavelength dependence 
;    parinfo[0].value = 1D
;    parinfo[1].value = 4D
;    parinfo[4].value = 0D
;    parinfo[5].value = 0D
;    parinfo[6].value = 0D
;    parinfo[7].value = 0D
;
;    parinfo[0].fixed = 1
;    parinfo[1].fixed = 1
;    parinfo[4].fixed = 1
;    parinfo[5].fixed = 1
;    parinfo[6].fixed = 1
;    parinfo[7].fixed = 1
    
    params = mpfitfun('bcgmstar_sersic2_allbands_func',xx,xxsb,parinfo=parinfo,$
      yfit=sersicfit,perror=perror,covar=covar,weights=xxsb_ivar,dof=dof,$
      bestnorm=chi2,status=status,niter=niter,quiet=keyword_set(verbose) eq 0,$
      xtol=1D-12,$
      functargs={parinfo: parinfo, wave: wave, lambda_ref: lambda_ref})
    factor = sqrt(chi2/dof)

;  djs_plot, xx, xxsb, psym=8, /xlog
;  djs_oplot, xx, sersicfit, psym=6, color='orange'
  
; build two output structures: one that has the formal fitting results
; (alpha, beta, etc.)...
    results = {$
      sersic2_niter:  niter,$
      sersic2_status: status,$
      sersic2_chi2:     chi2,$
      sersic2_dof:       dof,$
      sersic2_covar:   covar,$
      sersic2_params:  params,$
      sersic2_perror:  perror,$

      sersic2_wave:     uwave,$
      sersic2_wave_ref: wave_ref,$

      sersic2_n1_ref:  params[0],$
      sersic2_n2_ref:  params[1],$
      sersic2_re1_ref: params[2],$
      sersic2_re2_ref: params[3],$
      sersic2_alpha1:  params[4],$
      sersic2_alpha2:  params[5],$
      sersic2_beta1:   params[6],$
      sersic2_beta2:   params[7],$

      sersic2_n1_ref_err:  perror[0],$
      sersic2_n2_ref_err:  perror[1],$
      sersic2_re1_ref_err: perror[2],$
      sersic2_re2_ref_err: perror[3],$
;     sersic2_n1_ref_err:  perror[0]*10D^params[0]*alog(10),$
;     sersic2_n2_ref_err:  perror[1]*10D^params[1]*alog(10),$
;     sersic2_re1_ref_err: perror[2]*10D^params[2]*alog(10),$
;     sersic2_re2_ref_err: perror[3]*10D^params[3]*alog(10),$
      sersic2_alpha1_err:  perror[4],$
      sersic2_alpha2_err:  perror[5],$
      sersic2_beta1_err:   perror[6],$
      sersic2_beta2_err:   perror[7],$

      sersic2_sbe1:       params[8:8+nband-1],$
      sersic2_sbe2:       params[8+nband:8+2*nband-1],$
      sersic2_sbe1_err:   perror[8:8+nband-1],$
      sersic2_sbe2_err:   perror[8+nband:8+2*nband-1]}
    
; ...and a second one that has the model evaluated for each filter and
; which takes into account all the covariances in the parameters
    for ib = 0, nband-1 do begin
       n1 = results.sersic2_n1_ref*(uwave[ib]/wave_ref)^results.sersic2_alpha1
       n2 = results.sersic2_n2_ref*(uwave[ib]/wave_ref)^results.sersic2_alpha2
       re1 = results.sersic2_re1_ref*(uwave[ib]/wave_ref)^results.sersic2_beta1
       re2 = results.sersic2_re2_ref*(uwave[ib]/wave_ref)^results.sersic2_beta2

; these uncertainties are only approximate because they ignore the
; covariance between the parameters (note that I can't use
; MPRANDOMN because the covariance matrix is not positive-definite)
       n1_err = sqrt((results.sersic2_n1_ref_err/results.sersic2_n1_ref/alog(10))^2 + $
         (results.sersic2_alpha1_err*(uwave[ib]/wave_ref))^2)*alog(10)*n1
       n2_err = sqrt((results.sersic2_n2_ref_err/results.sersic2_n2_ref/alog(10))^2 + $
         (results.sersic2_alpha2_err*(uwave[ib]/wave_ref))^2)*alog(10)*n2

       re1_err = sqrt((results.sersic2_re1_ref_err/results.sersic2_re1_ref/alog(10))^2 + $
         (results.sersic2_beta1_err*(uwave[ib]/wave_ref))^2)*alog(10)*re1
       re2_err = sqrt((results.sersic2_re2_ref_err/results.sersic2_re2_ref/alog(10))^2 + $
         (results.sersic2_beta2_err*(uwave[ib]/wave_ref))^2)*alog(10)*re2
       
       sbe1 = params[8+ib]
       sbe2 = params[8+nband+ib]
       sbe1_err = perror[8+ib]
       sbe2_err = perror[8+nband+ib]

       scoeff1 = {$
;        wave: uwave[ib], $
         sersic2_all_sbe1:     sbe1,$
         sersic2_all_re1:      re1,$
         sersic2_all_n1:       n1,$
         sersic2_all_sbe1_err: sbe1_err,$
         sersic2_all_re1_err:  re1_err,$
         sersic2_all_n1_err:   n1_err,$

         sersic2_all_sbe2:     sbe2,$
         sersic2_all_re2:      re2,$
         sersic2_all_n2:       n2,$
         sersic2_all_sbe2_err: sbe2_err,$
         sersic2_all_re2_err:  re2_err,$
         sersic2_all_n2_err:   n2_err}
       if ib eq 0 then scoeff = scoeff1 else scoeff = [scoeff,scoeff1]
    endfor
;   struct_print, scoeff
;   help, results, scoeff, /str

return
end

pro bcgmstar_sersic2, rr, sb, scoeff, sb_ivar=sb_ivar, $
  init_params=init_params, fixed=fixed, sersicfit=sersicfit, $
  fixdevac=fixdevac, verbose=verbose
; fit a double Sersic function

; parse the data; take the log
    xx = rr
    xxsb = -2.5*alog10(sb)
    xxsb_ivar = sb_ivar*(alog(10)*sb/2.5)^2.0
;   good = where(sb_ivar gt 0.0 and finite(sb),nn)
;   xx = rr[good]
;   xxsb = sb[good]
;   xxsb_ivar = sb_ivar[good] 

    nparam = 6
    parinfo = replicate({value: 0D, fixed: 0, limited: [0,0], $
      limits: [0D,0D]},nparam)

; sbe1 > 0 
    parinfo[0].limited = [1,1]
    parinfo[0].limits = [10D,35D]
;   parinfo[0].limits = [1D-12,0D]
    parinfo[0].value = median(xxsb)
; re1>0
    parinfo[1].limited = [1,1]
    parinfo[1].limits = [0.5D,500D]
    parinfo[1].value = 1D
; 1<n1<10
    parinfo[2].limited = [1,1]
    parinfo[2].limits = [0.1D,10D]
    parinfo[2].value = 4D

; sbe2 > 0
    parinfo[3].limited = [1,1]
    parinfo[3].limits = [10D,35D]
    parinfo[3].value = median(xxsb)
; re2>0
    parinfo[4].limited = [1,1]
    parinfo[4].limits = [0.01D,500D]
    parinfo[4].value = 10D
; 1<n2<10
    parinfo[5].limited = [1,1]
    parinfo[5].limits = [0.1D,10D]        
    parinfo[5].value = 4D
    
    if n_elements(init_params) eq nparam then parinfo.value = init_params
    if n_elements(fixed) eq nparam then parinfo.fixed = fixed

    params = mpfitfun('bcgmstar_sersic2_func',xx,xxsb,parinfo=parinfo,$
      yfit=sersicfit,perror=perror,covar=covar,weights=xxsb_ivar,$
      dof=dof,bestnorm=chi2,status=status,quiet=keyword_set(verbose) eq 0,$
      functargs={parinfo: parinfo})

;   djs_plot, rr, -2.5*alog10(xxsb), /xlog, psym=8, xsty=3, ysty=3, $
;     yr=-2.5*alog10(minmax(xxsb))
;   djs_oplot, rr, -2.5*alog10(sersicfit), color='red'
;   djs_oplot, rr, bcgmstar_sersic_func(rr,[-2.5*alog10(params[0]),params[1],params[2]]), $
;     color='cyan', line=5
;   djs_oplot, rr, bcgmstar_sersic_func(rr,[-2.5*alog10(params[3]),params[4],params[5]]), $
;     color='orange', line=5
;   cc = get_kbrd(1)
    
    scoeff = {$
      sersic2_status: status,$
      sersic2_chi2:     chi2,$
      sersic2_dof:       dof,$
      sersic2_covar:   covar,$

      sersic2_sbe1: params[0], $
      sersic2_re1:  params[1], $
      sersic2_n1:   params[2], $
      sersic2_sbe1_err: perror[0], $
      sersic2_re1_err:  perror[1], $
      sersic2_n1_err:   perror[2], $

      sersic2_sbe2: params[3], $
      sersic2_re2:  params[4], $
      sersic2_n2:   params[5], $
      sersic2_sbe2_err: perror[3], $
      sersic2_re2_err:  perror[4], $
      sersic2_n2_err:   perror[5]}

;     sersic2_total1:         0.0, $
;     sersic2_total2: 0.0}
;   scoeff.sersic2_total1 = cumsersic_total(params[0:2])
;   scoeff.sersic2_total2 = cumsersic_total(params[3:5])
    
return
end

pro bcgmstar_sersic, rr, sb, scoeff, sb_ivar=sb_ivar, $
  init_params=init_params, fixed=fixed, sersicfit=sersicfit, $
  verbose=verbose
; fit a single Sersic function; RR should be in kpc; the model
; returned is in magnitudes 

; parse the data; take the log
    xx = rr
    xxsb = -2.5*alog10(sb)
    xxsb_ivar = sb_ivar*(alog(10)*sb/2.5)^2.0
;   xxsb_ivar = sb_ivar[good]*sb[good]^2.0
;   xxsb = sb[good]
;   xxsb_ivar = sb_ivar[good] 

    nparam = 3
    parinfo = replicate({value: 0D, fixed: 0, limited: [0,0], $
      limits: [0D,0D]},nparam)

; sbe - surface brightness at re
    parinfo[0].limited = [1,1]
    parinfo[0].limits = [10D,35D]
    parinfo[0].value = median(xxsb)
; re>0
    parinfo[1].limited = [1,1]
    parinfo[1].limits = [0.01D,500D]
    parinfo[1].value = 10D
; 0<n<10
    parinfo[2].limited = [1,1]
    parinfo[2].limits = [0.1D,10D]
    parinfo[2].value = 4D

    if n_elements(init_params) eq nparam then parinfo.value = init_params
    if n_elements(fixed) eq nparam then parinfo.fixed = fixed
    
    params = mpfitfun('bcgmstar_sersic_func',xx,xxsb,parinfo=parinfo,$
      perror=perror,covar=covar,weights=xxsb_ivar,bestnorm=chi2,dof=dof,$
      status=status,quiet=keyword_set(verbose) eq 0,yfit=sersicfit,$
      functargs={parinfo: parinfo})

    factor = sqrt(chi2/dof)
    scoeff = {$
      sersic_status: status,$
      sersic_chi2:     chi2,$
      sersic_dof:       dof,$
      sersic_covar:   covar,$

      sersic_sbe: params[0],$
      sersic_re:  params[1],$
      sersic_n:   params[2],$
      sersic_sbe_err: perror[0],$ ; *factor,$
      sersic_re_err:  perror[1],$ ; *factor,$
      sersic_n_err:   perror[2]}  ; *factor,$

;     sersic_total: cumsersic_total(params)}

return
end

pro bcgmstar_sersicfit, dofit_allbands=dofit_allbands, dofit_oneband=dofit_oneband, $
  qaplot_allbands=qaplot_allbands, qaplot_oneband=qaplot_oneband, $
  alphabetazero=alphabetazero, verbose=verbose, clobber=clobber
; jm13oct22siena - fit various Sersic models to the output of
; BCGMSTAR_ELLIPSE 

; dofit_allbands - fit a wavelength-dependent single and double Sersic
;   model to all bands simultaneously
; dofit_oneband - fit single and double Sersic models to every band
;   independently
; alphabetazero - to be used in tandem with /dofit_allbands; fixes
;   alpha and beta at zero (i.e., no wavelength-dependent fitting)

    if keyword_set(qaplot_allbands) and keyword_set(qaplot_oneband) then begin
       splog, 'qaplot_allbands and qaplot_oneband cannot be set simultaneously'
       return
    endif
    
    ellpath = bcgmstar_path(/ellipse)
    sersicpath = bcgmstar_path(/sersic)
    qapath = bcgmstar_path()+'qaplots-sersic/'

; read the sample
    sample = read_bcgmstar_sample()
    sample = sample[0]
    ncl = n_elements(sample)

    pixscale = 0.065D           ; [arcsec/pixel]
    errfloor = 0.0D ; 0.02      ; magnitude error floor on my SB measurements 
    lambda_ref = 8002.83D
    
; ##################################################
; fit a wavelength-dependent single and double Sersic model to all
; bands simultaneously  
    if keyword_set(dofit_allbands) then begin
; wrap on each cluster    
;      for ic = 0, 0 do begin
       for ic = 0, ncl-1 do begin
          arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
          
          cluster = strtrim(sample[ic].shortname,2)
          splog, 'Simultaneously Sersic fitting cluster '+cluster
          
; read the data and Marc's SB profiles to determine the last radius at
; which the models are reliable
          modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
          allband = strtrim(modphot.band,2)

          pp = read_bcg_profiles(cluster,these_filters=allband)
          nfilt = n_elements(modphot)

; output structure
          out = struct_addtags(struct_trimtags(modphot,$
            select=['file','band','weff','sblimit','ra','dec','mge_*']),$
            replicate({amin_pixels: 0.0, amax_pixels: 0.0, amin_kpc: 0.0, amax_kpc: 0.0, $
            rmin_kpc: 0.0, rmax_kpc: 0.0},nfilt))

          init_sb = fltarr(nfilt)
          for ib = 0, nfilt-1 do begin
             ppgood = where(pp[ib].sma gt -90)
             amax_fit = max(pp[ib].sma[ppgood],mxindx)
             amax_sblimit = interpol(pp[ib].sma[ppgood],pp[ib].mu[ppgood],$
               modphot[ib].sblimit-2.5*alog10(2)) ; 3-sigma
             amax_kpc = amax_fit < amax_sblimit

             modgood = where(modphot[ib].majora*pixscale*arcsec2kpc le amax_kpc and $
               modphot[ib].sb0fit gt 0 and modphot[ib].sb0fit_ivar gt 0,nmodgood)
             amin_kpc = min(modphot[ib].majora[modgood])*pixscale*arcsec2kpc
             
             sb = modphot[ib].sb0fit[modgood]*1D
             sb_var_floor = (sb*errfloor)^2.0
             sb_ivar = 1D/(1D/modphot[ib].sb0fit_ivar[modgood]+sb_var_floor)
             radius_kpc = modphot[ib].radius_kpc[modgood] ; [kpc]

             rmin_kpc = min(radius_kpc)
             rmax_kpc = max(radius_kpc)
                
             out[ib].amin_pixels = amin_kpc/pixscale/arcsec2kpc
             out[ib].amax_pixels = amax_kpc/pixscale/arcsec2kpc
             out[ib].amin_kpc = amin_kpc
             out[ib].amax_kpc = amax_kpc
             out[ib].rmin_kpc = rmin_kpc
             out[ib].rmax_kpc = rmax_kpc
             
             srt = sort(radius_kpc)
             radius_kpc = radius_kpc[srt]
             sb = sb[srt]
             sb_ivar = sb_ivar[srt]
             init_sb[ib] = djs_mean(sb)

             if ib eq 0 then begin
                fit_sb = sb
                fit_sb_ivar = sb_ivar
                fit_radius_kpc = radius_kpc
                fit_wave = replicate(modphot[ib].weff,nmodgood)
             endif else begin
                fit_sb = [fit_sb,sb]
                fit_sb_ivar = [fit_sb_ivar,sb_ivar]
                fit_radius_kpc = [fit_radius_kpc,radius_kpc]
                fit_wave = [fit_wave,replicate(modphot[ib].weff,nmodgood)]
             endelse
          endfor 

; the optimal model for A2261 and A209 is two Sersic functions with no
; wavelength dependence; set the initial conditions here based on
; looking at the qa_sersic_a209.pdf and qa_sersic_a2261.pdf files
          case cluster of
             'a2261': begin
                init_sersic2 = [0.6D,1.3D,6D,30D,$
                  0D,0D,0D,0D,-2.5*alog10(init_sb)-1.0,-2.5*alog10(init_sb)]
             end
             'a209': begin
                init_sersic2 = [0.5D,3D,0.6D,30D,$
                  0D,0D,0D,0D,-2.5*alog10(init_sb)-1.0,-2.5*alog10(init_sb)]
             end
             else: delvarx, init_sersic2
          endcase
          
          bcgmstar_sersic_allbands, fit_radius_kpc, fit_sb, fit_wave, $
            allsersic, sb_ivar=fit_sb_ivar, sersicfit=sersicfit, $
            fixdevac=fixdevac, verbose=verbose, results=results, $
            init_params=init_sersic, fixed=fixed, lambda_ref=lambda_ref, $
            alphabetazero=alphabetazero
          out = struct_addtags(out,allsersic)

          bcgmstar_sersic2_allbands, fit_radius_kpc, fit_sb, fit_wave, $
            allsersic2, sb_ivar=fit_sb_ivar, sersicfit=sersic2fit, $
            fixdevac=fixdevac, verbose=verbose, results=results2, $
            init_params=init_sersic2, fixed=fixed, lambda_ref=lambda_ref, $
            alphabetazero=alphabetazero
          out = struct_addtags(out,allsersic2)

          ploterror, fit_radius_kpc, -2.5*alog10(fit_sb), $
            1.0/sqrt(fit_sb_ivar*(alog(10)*fit_sb/2.5)^2.0), $
            /xlog, psym=8, /trad, yrange=[25,17]
          djs_oplot, fit_radius_kpc, sersicfit, color='green', psym=6
          djs_oplot, fit_radius_kpc, sersic2fit, color='orange', psym=7

; write out
          if keyword_set(alphabetazero) then suffix = '-alphabetazero' else suffix = ''
          im_mwrfits, out, sersicpath+cluster+'-allsersic'+suffix+'.fits', clobber=clobber
          im_mwrfits, struct_addtags(results,results2), sersicpath+cluster+$
            '-allsersic'+suffix+'-results.fits', clobber=clobber
       endfor
    endif 
    
; ##################################################
; fit single and double Sersic models to every band independently 
    if keyword_set(dofit_oneband) then begin
; wrap on each cluster    
;      for ic = 0, 0 do begin
       for ic = 0, ncl-1 do begin
          arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
          
          cluster = strtrim(sample[ic].shortname,2)
          splog, 'Sersic fitting cluster '+cluster
          
; read the data and Marc's SB profiles to determine the last radius at
; which the models are reliable
          modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
          allband = strtrim(modphot.band,2)

          pp = read_bcg_profiles(cluster,these_filters=allband)
          nfilt = n_elements(modphot)

; output structure
          out = struct_addtags(struct_trimtags(modphot,$
            select=['file','band','weff','sblimit','ra','dec','mge_*']),$
            replicate({amin_pixels: 0.0, amax_pixels: 0.0, amin_kpc: 0.0, amax_kpc: 0.0, $
            rmin_kpc: 0.0, rmax_kpc: 0.0},nfilt))

          for ib = 0, nfilt-1 do begin
; if AMAX_KPC occurs when the SB profile is below the formal 1-sigma
; surface brightness limit of the data, then cut it off
             ppgood = where(pp[ib].sma gt -90)
             amax_fit = max(pp[ib].sma[ppgood],mxindx)
             amax_sblimit = interpol(pp[ib].sma[ppgood],pp[ib].mu[ppgood],$
               modphot[ib].sblimit-2.5*alog10(2)) ; 3-sigma
             amax_kpc = amax_fit < amax_sblimit

;            amax_kpc = max(pp[ib].sma,mxindx) ; [kpc]
;            if pp[ib].mu[mxindx] gt modphot[ib].sblimit then $
;              amax_kpc = pp[ib].sma[mxindx-1]

             modgood = where(modphot[ib].majora*pixscale*arcsec2kpc le amax_kpc and $
               modphot[ib].sb0fit gt 0 and modphot[ib].sb0fit_ivar gt 0,nmodgood)
             amin_kpc = min(modphot[ib].majora[modgood])*pixscale*arcsec2kpc

; the equivalent radius needs to be sorted (note: the ellipse
; parameters in BCGMSTAR_ELLIPSE are monotonic in *semi-major axis*,
; not equivalent radius!)  also add a minimum error floor to the
; surface brightnesses
             sb = modphot[ib].sb0fit[modgood]*1D
             sb_var_floor = (sb*errfloor)^2.0
             sb_ivar = 1D/(1D/modphot[ib].sb0fit_ivar[modgood]+sb_var_floor)
             radius_kpc = modphot[ib].radius_kpc[modgood] ; [kpc]
             
             rmin_kpc = min(radius_kpc)
             rmax_kpc = max(radius_kpc)
             
             out[ib].amin_pixels = amin_kpc/pixscale/arcsec2kpc
             out[ib].amax_pixels = amax_kpc/pixscale/arcsec2kpc
             out[ib].amin_kpc = amin_kpc
             out[ib].amax_kpc = amax_kpc
             out[ib].rmin_kpc = rmin_kpc
             out[ib].rmax_kpc = rmax_kpc
             
             srt = sort(radius_kpc)
             radius_kpc = radius_kpc[srt]
             sb = sb[srt]
             sb_ivar = sb_ivar[srt]
             
; fit with a single-Sersic and then a double-Sersic
             bcgmstar_sersic, radius_kpc, sb, sersic, sb_ivar=sb_ivar, verbose=verbose
             bcgmstar_sersic2, radius_kpc, sb, sersic2, sb_ivar=sb_ivar, $
               fixdevac=fixdevac, verbose=verbose

; expand the output structure
             if ib eq 0 then begin
                out = struct_addtags(out,im_empty_structure(sersic[0],ncopies=nfilt))
                out = struct_addtags(out,im_empty_structure(sersic2[0],ncopies=nfilt))
             endif 
             
             out[ib] = im_struct_assign(sersic,out[ib],/nozero)
             out[ib] = im_struct_assign(sersic2,out[ib],/nozero)
          endfor 

; write out
          im_mwrfits, out, sersicpath+cluster+'-sersic.fits', clobber=clobber
       endfor
    endif
    
; ##################################################
; build QAplots
    ncol = 3 ; number of columns
    rr = [0,range(0.01,200,500,/log)] ; equivalent radius [kpc]

; -------------------------    
; QAplot: SB profiles and the Sersic fits
    if keyword_set(qaplot_allbands) or keyword_set(qaplot_oneband) then begin
       if keyword_set(alphabetazero) then begin
          suffix1 = 'alphabetazero_'
          suffix2 = '-alphabetazero'
       endif else begin
          suffix1 = ''
          suffix2 = ''
       endelse

       for ic = 0, ncl-1 do begin
          cluster = strtrim(sample[ic].shortname,2)
          splog, cluster
          if keyword_set(qaplot_allbands) then psfile = qapath+'qa_allsersic_'+suffix1+cluster+'.ps' else $
            psfile = qapath+'qa_sersic_'+cluster+'.ps'
          im_plotconfig, 0, pos, psfile=psfile, charsize=1.3

          arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]

          if keyword_set(qaplot_allbands) then begin
             sersic = mrdfits(sersicpath+cluster+'-allsersic'+suffix2+'.fits.gz',1,/silent)
             sersic_results = mrdfits(sersicpath+cluster+'-allsersic'+suffix2+'-results.fits.gz',1,/silent)
          endif else begin
             sersic = mrdfits(sersicpath+cluster+'-sersic.fits.gz',1,/silent)
          endelse
          modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
          nfilt = n_elements(modphot)

          allband = strtrim(modphot.band,2)
          pp = read_bcg_profiles(cluster,these_filters=allband)
          
          if keyword_set(qaplot_allbands) then nrow = ceil((nfilt+2)/float(ncol)) else $
            nrow = ceil(nfilt/float(ncol))
          pos = im_getposition(nx=ncol,ny=nrow,yspace=0.0,xspace=0.0,$
            xmargin=[0.9,0.4],width=2.4)
          count = 0
;         for ib = nfilt-1, 0, -1 do begin ; reverse order
          for ib = 0, nfilt-1 do begin
             band = strtrim(strupcase(modphot[ib].band),2)

             modgood = where(modphot[ib].majora*pixscale*arcsec2kpc le sersic[ib].amax_kpc and $
               modphot[ib].sb0fit gt 0 and modphot[ib].sb0fit_ivar gt 0)
             modbad = where(modphot[ib].majora*pixscale*arcsec2kpc gt sersic[ib].amax_kpc and $
               modphot[ib].sb0fit gt 0 and modphot[ib].sb0fit_ivar gt 0)
             
             sb = modphot[ib].sb0fit[modgood]*1D
             sb_var_floor = (sb*errfloor)^2.0
             sb_ivar = 1D/(1D/modphot[ib].sb0fit_ivar[modgood]+sb_var_floor)

             radius_kpc = modphot[ib].radius_kpc[modgood]  ; [kpc]

             if count eq 1 then title = strupcase(cluster) else delvarx, title
             if count ge nfilt-3 then begin
                delvarx, xtickname
             endif else begin
                xtickname = replicate(' ',10)
             endelse
             if (count mod 3) eq 0 then begin
                delvarx, ytickname
             endif else begin
                ytickname = replicate(' ',10)
             endelse
             
             djs_plot, radius_kpc, -2.5*alog10(sb), psym=symcat(16), /xlog, noerase=count gt 0, $
               xrange=[0.3,200], xsty=1, yrange=[30,16], position=pos[*,count], $
               xtickname=xtickname, ytickname=ytickname, title=title, $
               symsize=0.5, ytickinterval=4, ysty=1

             if keyword_set(qaplot_allbands) then begin
                label = [$
                  '\mu_{e}='+strtrim(string(sersic[ib].sersic_all_sbe,format='(F12.1)'),2)+','+$
                  'n='+strtrim(string(sersic[ib].sersic_all_n,format='(F12.2)'),2)+','+$
                  'r_{e}='+strtrim(string(sersic[ib].sersic_all_re,format='(F12.2)'),2)+' kpc',$
                  '\mu_{e1}='+strtrim(string(sersic[ib].sersic2_all_sbe1,format='(F12.1)'),2)+','+$
                  'n_{1}='+strtrim(string(sersic[ib].sersic2_all_n1,format='(F12.2)'),2)+','+$
                  'r_{e1}='+strtrim(string(sersic[ib].sersic2_all_re1,format='(F12.2)'),2)+' kpc',$
                  '\mu_{e2}='+strtrim(string(sersic[ib].sersic2_all_sbe2,format='(F12.1)'),2)+','+$
                  'n_{2}='+strtrim(string(sersic[ib].sersic2_all_n2,format='(F12.2)'),2)+','+$
                  'r_{e2}='+strtrim(string(sersic[ib].sersic2_all_re2,format='(F12.2)'),2)+' kpc']
             endif else begin
                label = [$
                  '\chi^{2}_{\nu, single, double}='+$
                  strtrim(string(sersic[ib].sersic_chi2/sersic[ib].sersic_dof,format='(F12.2)'),2)+', '+$
                  strtrim(string(sersic[ib].sersic2_chi2/sersic[ib].sersic2_dof,format='(F12.2)'),2),$
                  '\mu_{e}='+strtrim(string(sersic[ib].sersic_sbe,format='(F12.1)'),2)+','+$
                  'n='+strtrim(string(sersic[ib].sersic_n,format='(F12.2)'),2)+','+$
                  'r_{e}='+strtrim(string(sersic[ib].sersic_re,format='(F12.1)'),2)+' kpc',$
                  '\mu_{e1}='+strtrim(string(sersic[ib].sersic2_sbe1,format='(F12.1)'),2)+','+$
                  'n_{1}='+strtrim(string(sersic[ib].sersic2_n1,format='(F12.2)'),2)+','+$
                  'r_{e1}='+strtrim(string(sersic[ib].sersic2_re1,format='(F12.2)'),2)+' kpc',$
                  '\mu_{e2}='+strtrim(string(sersic[ib].sersic2_sbe2,format='(F12.1)'),2)+','+$
                  'n_{2}='+strtrim(string(sersic[ib].sersic2_n2,format='(F12.2)'),2)+','+$
                  'r_{e2}='+strtrim(string(sersic[ib].sersic2_re2,format='(F12.2)'),2)+' kpc']
             endelse 
             im_legend, label, /left, /bottom, box=0, margin=0, charsize=0.7, charthick=1.8

             if modbad[0] ne -1 then begin
                djs_oplot, modphot[ib].radius_kpc[modbad], -2.5*alog10(modphot[ib].sb0fit[modbad]), $
                  psym=symcat(9), color=cgcolor('medium grey'), symsize=0.5
             endif
             
             im_legend, band, /right, /top, box=0, margin=0, charsize=1.0

             if keyword_set(qaplot_allbands) then begin
;               djs_oplot, rr, -2.5*alog10(bcgmstar_sersic2_allbands_func(rr,$
;                 sersic_results.params,wave=rr*0+sersic[ib].wave)), color=cgcolor('forest green')
                djs_oplot, rr, bcgmstar_sersic2_func(rr,params=sersic[0],/allbands), $
                  color=cgcolor('forest green')
                djs_oplot, rr, bcgmstar_sersic2_func(rr,params=sersic[ib],/allbands), $
                  color=cgcolor('firebrick')
                djs_oplot, rr, bcgmstar_sersic_func(rr,[sersic[ib].sersic2_all_sbe1,$
                  sersic[ib].sersic2_all_re1,sersic[ib].sersic2_all_n1]), color=cgcolor('orange'), line=2
                djs_oplot, rr, bcgmstar_sersic_func(rr,[sersic[ib].sersic2_all_sbe2,$
                  sersic[ib].sersic2_all_re2,sersic[ib].sersic2_all_n2]), color=cgcolor('orange'), line=2

                djs_oplot, rr, bcgmstar_sersic_func(rr,params=sersic[ib],/allbands), $
                  color=cgcolor('purple'), thick=4
             endif else begin
                djs_oplot, rr, bcgmstar_sersic_func(rr,params=sersic[ib]), $
                  color=cgcolor('firebrick')
                if ib gt 0 then djs_oplot, rr, bcgmstar_sersic_func(rr,$
                  params=sersic[0]), color=cgcolor('forest green')
             
                djs_oplot, rr, bcgmstar_sersic2_func(rr,params=sersic[ib]), $
                  color=cgcolor('dodger blue')
                if sersic[ib].sersic2_sbe1 eq 0.0 or sersic[ib].sersic2_sbe2 eq 0.0 then begin
                   splog, '  '+band+': second Sersic dropped!'
                endif else begin
                   djs_oplot, rr, bcgmstar_sersic_func(rr,[sersic[ib].sersic2_sbe1,$
                     sersic[ib].sersic2_re1,sersic[ib].sersic2_n1]), color=cgcolor('orange'), line=2
                   djs_oplot, rr, bcgmstar_sersic_func(rr,[sersic[ib].sersic2_sbe2,$
                     sersic[ib].sersic2_re2,sersic[ib].sersic2_n2]), color=cgcolor('orange'), line=2
                endelse
             endelse

             rmax_fit = interpol(modphot[ib].radius_kpc[modgood],$
               modphot[ib].majora[modgood]*pixscale*arcsec2kpc,max(pp[ib].sma))
             
             djs_oplot, rmax_fit*[1,1], !y.crange, line=1, thick=1
             djs_oplot, 10^!x.crange, modphot[ib].sblimit*[1,1], line=1, thick=1
             djs_oplot, 10^!x.crange, (modphot[ib].sblimit-2.5*alog10(2))*[1,1], line=1, thick=1
             djs_oplot, 10^!x.crange, (modphot[ib].sblimit-2.5*alog10(3))*[1,1], line=1, thick=1
;            djs_oplot, [70.0,10^!x.crange[1]], modphot[ib].sblimit*[1,1], line=0
;            djs_oplot, [6.0,10^!x.crange[1]], modphot[ib].sblimit*[1,1], line=0
;            djs_oplot, sersic[ib].rmax_kpc*[1,1], [!y.crange[0]-0.2,!y.crange[0]-6], line=0
             count++
          endfor 

          if keyword_set(qaplot_allbands) then begin
             label = ['\chi^{2}_{\nu}='+strtrim(string(sersic_results.sersic_chi2/$
               sersic_results.sersic_dof,format='(F12.2)'),2),$
               '\nu='+strtrim(string(sersic_results.sersic_dof,format='(I0)'),2),$
               '\alpha_{1}='+strtrim(string(sersic_results.sersic_alpha1,format='(F12.3)'),2)+'\pm'+$
               strtrim(string(sersic_results.sersic_alpha1_err,format='(F12.3)'),2),$
               '\beta_{1}='+strtrim(string(sersic_results.sersic_beta1,format='(F12.3)'),2)+'\pm'+$
               strtrim(string(sersic_results.sersic_beta1_err,format='(F12.3)'),2)]

             label2 = ['\chi^{2}_{\nu}='+strtrim(string(sersic_results.sersic2_chi2/$
               sersic_results.sersic2_dof,format='(F12.2)'),2),$
               '\nu='+strtrim(string(sersic_results.sersic2_dof,format='(I0)'),2),$
               '\alpha_{1}='+strtrim(string(sersic_results.sersic2_alpha1,format='(F12.3)'),2)+'\pm'+$
               strtrim(string(sersic_results.sersic2_alpha1_err,format='(F12.3)'),2),$
               '\alpha_{2}='+strtrim(string(sersic_results.sersic2_alpha2,format='(F12.3)'),2)+'\pm'+$
               strtrim(string(sersic_results.sersic2_alpha2_err,format='(F12.3)'),2),$
               '\beta_{1}='+strtrim(string(sersic_results.sersic2_beta1,format='(F12.3)'),2)+'\pm'+$
               strtrim(string(sersic_results.sersic2_beta1_err,format='(F12.3)'),2),$
               '\beta_{2}='+strtrim(string(sersic_results.sersic2_beta2,format='(F12.3)'),2)+'\pm'+$
               strtrim(string(sersic_results.sersic2_beta2_err,format='(F12.3)'),2)]

             djs_plot, [0], [0], /nodata, /noerase, position=pos[*,count], $
               xsty=5, ysty=5
             im_legend, label, /left, /bottom, box=0, margin=0, charsize=1, charthick=1.8

             djs_plot, [0], [0], /nodata, /noerase, position=pos[*,count+1], $
               xsty=5, ysty=5
             im_legend, label2, /left, /bottom, box=0, margin=0, charsize=1, charthick=1.8
          endif
          
          xyouts, min(pos[0,*])-0.06, (max(pos[3,*])-min(pos[1,*]))/2.0+min(pos[1,*]), $
            textoidl('\mu (mag arcsec^{-2})'), orientation=90, align=0.5, charsize=1.4, /norm
          xyouts, (max(pos[2,*])-min(pos[0,*]))/2.0+min(pos[0,*]), min(pos[1,*])-0.06, $
            textoidl('Equivalent Radius (kpc)'), align=0.5, charsize=1.4, /norm
          im_plotconfig, psfile=psfile, /psclose, /pdf
       endfor 

    endif
       
return
end



;; this (working) code is to use the F160W to constrain the free
;; parameters of the Sersic model
;          if ib gt 0 then begin
;             bcgmstar_sersic, radius_kpc, sb, sersic, sb_ivar=sb_ivar, $
;               init_params=[out[0].sersic_sb0,out[0].sersic_k,$
;               out[0].sersic_n], fixed=[0,1,1]
;          endif else begin
;             bcgmstar_sersic, radius_kpc, sb, sersic, sb_ivar=sb_ivar
;          endelse




;; -------------------------    
;; QAplot: SEDs
;    if keyword_set(qaplot_seds) then begin
;       splog, 'This is a dumb QAplot'
;       psfile = qapath+'qa_seds.ps'
;       im_plotconfig, 0, pos, psfile=psfile, charsize=1.8, $
;         height=5.0
;    
;       xrange = [0.3,2.0]
;       yrange = [26,10]
;    
;       for ic = 0, ncl-1 do begin
;          cluster = strtrim(sample[ic].shortname,2)
;          phot = mrdfits(sersicpath+cluster+'-phot.fits.gz',1,/silent)
;          nrad = n_elements(phot[0].photradius_kpc)
;          nfilt = n_elements(phot)
;
;          djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;            xrange=xrange, yrange=yrange, xtitle='Wavelength \lambda (\mu'+'m)', $
;            ytitle='Magnitude (AB)', /xlog, title=strupcase(cluster)
;
;; integrated light       
;          mab = maggies2mag(phot.maggies_int,ivarmaggies=phot.ivarmaggies_int,$
;            magerr=maberr,lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper)
;          used = where(mab gt -90.0,nused)
;          upper = where(mab lt -90.0 and mabupper gt -90,nupper)
;          
;          oploterror, phot[used].weff/1D4, mab[used], mabhierr[used], $
;            psym=-symcat(16), symsize=1.5, color=cgcolor('firebrick'), $
;            /hibar, errcolor=cgcolor('firebrick')
;          oploterror, phot[used].weff/1D4, mab[used], mabloerr[used], psym=3, $
;            color=cgcolor('firebrick'), /lobar, errcolor=cgcolor('firebrick')
;
;; radial bins       
;          for ir = 0, nrad-1 do begin
;             mab = maggies2mag(phot.maggies[ir],ivarmaggies=phot.ivarmaggies[ir],$
;               magerr=maberr,lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper)
;             used = where(mab gt -90.0,nused)
;             upper = where(mab lt -90.0 and mabupper gt -90,nupper)
;             if (nused ne 0) then begin
;                oploterror, phot[used].weff/1D4, mab[used], mabhierr[used], $
;                  psym=-symcat(16), symsize=1.5, color=cgcolor('dodger blue'), $
;                  /hibar, errcolor=cgcolor('dodger blue')
;                oploterror, phot[used].weff/1D4, mab[used], mabloerr[used], psym=3, $
;                  color=cgcolor('dodger blue'), /lobar, errcolor=cgcolor('dodger blue')
;             endif
;             if (nupper ne 0) then begin
;                djs_oplot, [phot[upper].weff/1D4], [mabupper[upper]], $
;                  psym=symcat(11,thick=6), symsize=2.0, color=cgcolor('forest green')
;             endif
;          endfor 
;       endfor 
;       im_plotconfig, psfile=psfile, /psclose, /pdf
;    endif


;; ##################################################
;; fit single and double Sersic models to every band independently 
;    if keyword_set(dofit) then begin
;; wrap on each cluster    
;;      for ic = 0, 0 do begin
;       for ic = 0, ncl-1 do begin
;          arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
;          
;          cluster = strtrim(sample[ic].shortname,2)
;          splog, 'Sersic fitting cluster '+cluster
;          
;; read the data and Marc's SB profiles to determine the last radius at
;; which the models are reliable
;          modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
;          allband = strtrim(modphot.band,2)
;
;          pp = read_bcg_profiles(cluster,these_filters=allband)
;          nfilt = n_elements(modphot)
;
;; first fit the (rest-frame) near-IR bands simultaneously, then the
;; optical bands together, and finally the bluest bands; we do this
;; because to help constrain the low S/N bands
;          delvarx, nir, opt, blue, nir_multiband
;          nir_multiband = 0 ; by default fit the near-IR bands individually
;          case cluster of
;             'a209': begin
;                opt = -1
;                blue = where(allband eq 'f475w' or allband eq 'f435w' or allband eq 'f390w')
;                nir = lindgen(nfilt) & remove, [blue], nir ; everything else
;             end
;             'a383': begin
;                opt = where(allband eq 'f775w' or allband eq 'f814w')
;                blue = where(allband eq 'f475w' or allband eq 'f435w')
;                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
;             end
;             'macs0744': begin
;                opt = -1
;                blue = where(allband eq 'f475w' or allband eq 'f555w')
;                nir = lindgen(nfilt) & remove, [blue], nir ; everything else
;             end
;             'a611': begin
;                opt = where(allband eq 'f850lp' or allband eq 'f606w' or $
;                  allband eq 'f814w' or allband eq 'f775w')
;                blue = where(allband eq 'f475w' or allband eq 'f435w')
;                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
;                nir_multiband = 1 ; constrain!
;             end
;             'macs1149': begin
;                opt = where(allband eq 'f850lp' or allband eq 'f606w' or $
;                  allband eq 'f814w' or allband eq 'f775w')
;                blue = where(allband eq 'f475w' or allband eq 'f435w')
;                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
;                nir_multiband = 1 ; constrain!
;             end
;             'a1423': begin
;                opt = where(allband eq 'f850lp' or allband eq 'f814w' or $
;                  allband eq 'f775w' or allband eq 'f625w')
;                blue = where(allband eq 'f606w' or allband eq 'f475w' or $
;                  allband eq 'f435w' or allband eq 'f390w')
;                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
;             end
;             'macs1206': begin
;                opt = where(allband eq 'f850lp' or allband eq 'f814w' or $
;                  allband eq 'f775w' or allband eq 'f625w')
;                blue = where(allband eq 'f606w' or allband eq 'f475w' or $
;                  allband eq 'f435w' or allband eq 'f390w')
;                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
;                nir_multiband = 1 ; constrain!
;             end
;             'clj1226': begin
;                blue = -1
;                opt = where(allband eq 'f850lp' or allband eq 'f814w' or $
;                  allband eq 'f775w' or allband eq 'f625w' or allband eq 'f606w')
;                nir = lindgen(nfilt) & remove, [opt], nir ; everything else
;                nir_multiband = 1 ; constrain!
;             end
;             'macs1311': begin
;                opt = where(allband eq 'f775w' or allband eq 'f625w' or allband eq 'f606w')
;                blue = where(allband eq 'f475w' or allband eq 'f435w')
;                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
;             end
;             'macs1720': begin
;                opt = where(allband eq 'f850lp' or allband eq 'f814w' or $
;                  allband eq 'f775w' or allband eq 'f625w')
;                blue = where(allband eq 'f606w' or allband eq 'f475w')
;                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
;             end
;             'macs2129': begin
;                opt = where(allband eq 'f850lp' or allband eq 'f814w')
;                blue = where(allband eq 'f475w' or allband eq 'f435w')
;                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
;             end
;             'rxj2129': begin
;                opt = where(allband eq 'f850lp' or allband eq 'f814w' or $
;                  allband eq 'f775w' or allband eq 'f625w' or allband eq 'f606w')
;                blue = where(allband eq 'f475w' or allband eq 'f435w')
;                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
;                nir_multiband = 1 ; constrain!
;             end
;; the blue bands for ms2137 don't fit well!
;             'ms2137': begin
;                opt = where(allband eq 'f850lp' or allband eq 'f814w' or $
;                  allband eq 'f775w' or allband eq 'f625w')
;                blue = where(allband eq 'f475w' or allband eq 'f435w' or allband eq 'f390w')
;                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
;                nir_multiband = 1 ; constrain!
;             end
;             'rxj2248': begin
;                opt = where(allband eq 'f850lp' or allband eq 'f814w' or $
;                  allband eq 'f775w' or allband eq 'f625w' or allband eq 'f606w')
;                blue = where(allband eq 'f475w' or allband eq 'f435w')
;                nir = lindgen(nfilt) & remove, [opt,blue], nir ; everything else
;                nir_multiband = 1 ; constrain!
;             end
;             else: begin
;                opt = -1
;                blue = where(modphot.weff/1D4 lt 0.5,comp=nir)
;             end
;          endcase
;          nnir = n_elements(nir)
;          if opt[0] eq -1 then nopt = 0 else nopt = n_elements(opt)
;          if blue[0] eq -1 then nblue = 0 else nblue = n_elements(blue)
;
;;         nopt = 0
;;         nir = where(modphot.weff/(1+sample[ic].z)/1D4 ge 0.4,nnir,comp=blue,ncomp=nblue)
;;         splog, 'Tying blue bands: '+strjoin(modphot[blue].band,', ')
;          
;;         nir = where(modphot.weff/(1+sample[ic].z)/1D4 ge 0.6,nnir)
;;         opt = where(modphot.weff/(1+sample[ic].z)/1D4 lt 0.6 and $
;;           modphot.weff/(1+sample[ic].z)/1D4 gt 0.4,nopt)
;;         blue = where(modphot.weff/(1+sample[ic].z)/1D4 le 0.4,nblue)
;
;; output structure
;          out = struct_addtags(struct_trimtags(modphot,$
;            select=['file','band','weff','sblimit','ra','dec','mge_*']),$
;            replicate({amin_pixels: 0.0, amax_pixels: 0.0, amin_kpc: 0.0, amax_kpc: 0.0, $
;            rmin_kpc: 0.0, rmax_kpc: 0.0, $
;            sersic_covar_nir:  fltarr(2+nnir,2+nnir),$
;            sersic_covar_opt:  fltarr(2+nopt,2+nopt),$
;            sersic_covar_blue: fltarr(2+nblue,2+nblue)},nfilt))
;          if dosersic2 then begin
;             out = struct_addtags(out,replicate({$
;               sersic2_covar_nir:  fltarr(2*nnir+4,2*nnir+4),$
;               sersic2_covar_opt:  fltarr(2*nopt+4,2*nopt+4),$
;               sersic2_covar_blue: fltarr(2*nblue+4,2*nblue+4)},nfilt))
;          endif
;
;; near-IR, fitted individually
;          if nir_multiband eq 0 then begin
;             for ib = 0, nnir-1 do begin
;
;; if AMAX_KPC occurs when the SB profile is below the formal 1-sigma
;; surface brightness limit of the data, then cut it off
;                amax_kpc = max(pp[nir[ib]].sma,mxindx) ; [kpc]
;                if pp[nir[ib]].mu[mxindx] gt modphot[nir[ib]].sblimit then $
;                  amax_kpc = pp[nir[ib]].sma[mxindx-1]
;
;                modgood = where(modphot[nir[ib]].majora*pixscale*arcsec2kpc le amax_kpc and $
;                  modphot[nir[ib]].sb0fit gt 0 and modphot[nir[ib]].sb0fit_ivar gt 0,nmodgood)
;                amin_kpc = min(modphot[nir[ib]].majora[modgood])*pixscale*arcsec2kpc
;
;; the equivalent radius needs to be sorted (note: the ellipse
;; parameters in BCGMSTAR_ELLIPSE are monotonic in *semi-major axis*,
;; not equivalent radius!)  also add a minimum error floor to the
;; surface brightnesses
;                sb = modphot[nir[ib]].sb0fit[modgood]*1D
;                sb_var_floor = (sb*errfloor)^2.0
;                sb_ivar = 1D/(1D/modphot[nir[ib]].sb0fit_ivar[modgood]+sb_var_floor)
;                radius_kpc = modphot[nir[ib]].radius_kpc[modgood] ; [kpc]
;                
;                rmin_kpc = min(radius_kpc)
;                rmax_kpc = max(radius_kpc)
;                
;                out[nir[ib]].amin_pixels = amin_kpc/pixscale/arcsec2kpc
;                out[nir[ib]].amax_pixels = amax_kpc/pixscale/arcsec2kpc
;                out[nir[ib]].amin_kpc = amin_kpc
;                out[nir[ib]].amax_kpc = amax_kpc
;                out[nir[ib]].rmin_kpc = rmin_kpc
;                out[nir[ib]].rmax_kpc = rmax_kpc
;                
;                srt = sort(radius_kpc)
;                radius_kpc = radius_kpc[srt]
;                sb = sb[srt]
;                sb_ivar = sb_ivar[srt]
;             
;; fit with a single-Sersic and then a double-Sersic
;                bcgmstar_sersic, radius_kpc, sb, sersic, sb_ivar=sb_ivar, verbose=verbose
;                if dosersic2 then bcgmstar_sersic2, radius_kpc, sb, sersic2, $
;                  sb_ivar=sb_ivar, fixdevac=fixdevac, verbose=verbose
;
;; expand the output structure
;                if ib eq 0 then begin
;                   out = struct_addtags(out,im_empty_structure(sersic[0],ncopies=nfilt))
;                   if dosersic2 then out = struct_addtags(out,im_empty_structure(sersic2[0],ncopies=nfilt))
;                endif 
;             
;                out[nir[ib]] = im_struct_assign(sersic,out[nir[ib]],/nozero)
;                if dosersic2 then out[nir[ib]] = im_struct_assign(sersic2,out[nir[ib]],/nozero)
;             endfor
;          endif
;
;; near-IR, fitted together
;          if nir_multiband then begin
;             for ib = 0, nnir-1 do begin
;                amax_kpc = max(pp[nir[ib]].sma,mxindx) ; [kpc]
;                if pp[nir[ib]].mu[mxindx] gt modphot[nir[ib]].sblimit then $
;                  amax_kpc = pp[nir[ib]].sma[mxindx-1]
;
;                modgood = where(modphot[nir[ib]].majora*pixscale*arcsec2kpc le amax_kpc and $
;                  modphot[nir[ib]].sb0fit gt 0 and modphot[nir[ib]].sb0fit_ivar gt 0,nmodgood)
;                amin_kpc = min(modphot[nir[ib]].majora[modgood])*pixscale*arcsec2kpc
;
;                sb = modphot[nir[ib]].sb0fit[modgood]*1D
;                sb_var_floor = (sb*errfloor)^2.0
;                sb_ivar = 1D/(1D/modphot[nir[ib]].sb0fit_ivar[modgood]+sb_var_floor)
;                radius_kpc = modphot[nir[ib]].radius_kpc[modgood] ; [kpc]
;
;                rmin_kpc = min(radius_kpc)
;                rmax_kpc = max(radius_kpc)
;                
;                out[nir[ib]].amin_pixels = amin_kpc/pixscale/arcsec2kpc
;                out[nir[ib]].amax_pixels = amax_kpc/pixscale/arcsec2kpc
;                out[nir[ib]].amin_kpc = amin_kpc
;                out[nir[ib]].amax_kpc = amax_kpc
;                out[nir[ib]].rmin_kpc = rmin_kpc
;                out[nir[ib]].rmax_kpc = rmax_kpc
;             
;                srt = sort(radius_kpc)
;                radius_kpc = radius_kpc[srt]
;                sb = sb[srt]
;                sb_ivar = sb_ivar[srt]
;             
;                if ib eq 0 then begin
;                   fit_sb = sb
;                   fit_sb_ivar = sb_ivar
;                   fit_radius_kpc = radius_kpc
;                   fit_wave = replicate(modphot[nir[ib]].weff,nmodgood)
;                endif else begin
;                   fit_sb = [fit_sb,sb]
;                   fit_sb_ivar = [fit_sb_ivar,sb_ivar]
;                   fit_radius_kpc = [fit_radius_kpc,radius_kpc]
;                   fit_wave = [fit_wave,replicate(modphot[nir[ib]].weff,nmodgood)]
;                endelse
;             endfor
;             
;             bcgmstar_sersic_multiband, fit_radius_kpc, fit_sb, fit_wave, $
;               multisersic, sb_ivar=fit_sb_ivar, verbose=verbose
;             out = struct_addtags(out,im_empty_structure(multisersic[0],ncopies=nfilt))
;             
;             for ib = 0, nnir-1 do begin
;                out[nir[ib]] = im_struct_assign(multisersic[ib],out[nir[ib]],/nozero)
;                out[nir[ib]].sersic_covar = 0
;                out[nir[ib]].sersic_covar_nir = multisersic[ib].sersic_covar
;             endfor
;             
;             if dosersic2 then begin
;                bcgmstar_sersic2_multiband, fit_radius_kpc, fit_sb, fit_wave, $
;                  multisersic2, sb_ivar=fit_sb_ivar, sersicfit=sersicfit, $
;                  fixdevac=fixdevac, verbose=verbose
;                out = struct_addtags(out,im_empty_structure(multisersic2[0],ncopies=nfilt))
;                for ib = 0, nnir-1 do begin
;                   out[nir[ib]] = im_struct_assign(multisersic2[ib],out[nir[ib]],/nozero)
;                   out[nir[ib]].sersic2_covar = 0
;                   out[nir[ib]].sersic2_covar_nir = multisersic2[ib].sersic2_covar
;                endfor
;             endif
;          endif
;
;; optical
;          if nopt gt 0 then begin
;             for ib = 0, nopt-1 do begin
;                amax_kpc = max(pp[opt[ib]].sma,mxindx) ; [kpc]
;                if pp[opt[ib]].mu[mxindx] gt modphot[opt[ib]].sblimit then $
;                  amax_kpc = pp[opt[ib]].sma[mxindx-1]
;                
;                modgood = where(modphot[opt[ib]].majora*pixscale*arcsec2kpc le amax_kpc and $
;                  modphot[opt[ib]].sb0fit gt 0 and modphot[opt[ib]].sb0fit_ivar gt 0,nmodgood)
;                amin_kpc = min(modphot[opt[ib]].majora[modgood])*pixscale*arcsec2kpc
;                
;                sb = modphot[opt[ib]].sb0fit[modgood]*1D
;                sb_var_floor = (sb*errfloor)^2.0
;                sb_ivar = 1D/(1D/modphot[opt[ib]].sb0fit_ivar[modgood]+sb_var_floor)
;                radius_kpc = modphot[opt[ib]].radius_kpc[modgood] ; [kpc]
;
;                rmin_kpc = min(radius_kpc)
;                rmax_kpc = max(radius_kpc)
;
;                out[opt[ib]].amin_pixels = amin_kpc/pixscale/arcsec2kpc
;                out[opt[ib]].amax_pixels = amax_kpc/pixscale/arcsec2kpc
;                out[opt[ib]].amin_kpc = amin_kpc
;                out[opt[ib]].amax_kpc = amax_kpc
;                out[opt[ib]].rmin_kpc = rmin_kpc
;                out[opt[ib]].rmax_kpc = rmax_kpc
;             
;                srt = sort(radius_kpc)
;                radius_kpc = radius_kpc[srt]
;                sb = sb[srt]
;                sb_ivar = sb_ivar[srt]
;             
;                if ib eq 0 then begin
;                   fit_sb = sb
;                   fit_sb_ivar = sb_ivar
;                   fit_radius_kpc = radius_kpc
;                   fit_wave = replicate(modphot[opt[ib]].weff,nmodgood)
;                endif else begin
;                   fit_sb = [fit_sb,sb]
;                   fit_sb_ivar = [fit_sb_ivar,sb_ivar]
;                   fit_radius_kpc = [fit_radius_kpc,radius_kpc]
;                   fit_wave = [fit_wave,replicate(modphot[opt[ib]].weff,nmodgood)]
;                endelse
;             endfor
;             
;             bcgmstar_sersic_multiband, fit_radius_kpc, fit_sb, fit_wave, $
;               multisersic, sb_ivar=fit_sb_ivar, verbose=verbose
;             for ib = 0, nopt-1 do begin
;                out[opt[ib]] = im_struct_assign(multisersic[ib],out[opt[ib]],/nozero)
;                out[opt[ib]].sersic_covar = 0
;                out[opt[ib]].sersic_covar_opt = multisersic[ib].sersic_covar
;             endfor
;          
;             if dosersic2 then begin
;                bcgmstar_sersic2_multiband, fit_radius_kpc, fit_sb, fit_wave, $
;                  multisersic2, sb_ivar=fit_sb_ivar, sersicfit=sersicfit, $
;                  fixdevac=fixdevac, verbose=verbose
;                for ib = 0, nopt-1 do begin
;                   out[opt[ib]] = im_struct_assign(multisersic2[ib],out[opt[ib]],/nozero)
;                   out[opt[ib]].sersic2_covar = 0
;                   out[opt[ib]].sersic2_covar_opt = multisersic2[ib].sersic2_covar
;                endfor
;             endif
;          endif
;
;; blue
;          if nblue ne 0 then begin
;             for ib = 0, nblue-1 do begin
;                amax_kpc = max(pp[blue[ib]].sma,mxindx) ; [kpc]
;                if pp[blue[ib]].mu[mxindx] gt modphot[blue[ib]].sblimit then $
;                  amax_kpc = pp[blue[ib]].sma[mxindx-1]
;                
;                modgood = where(modphot[blue[ib]].majora*pixscale*arcsec2kpc le amax_kpc and $
;                  modphot[blue[ib]].sb0fit gt 0 and modphot[blue[ib]].sb0fit_ivar gt 0,nmodgood)
;                amin_kpc = min(modphot[blue[ib]].majora[modgood])*pixscale*arcsec2kpc
;                
;                sb = modphot[blue[ib]].sb0fit[modgood]*1D
;                sb_var_floor = (sb*errfloor)^2.0
;                sb_ivar = 1D/(1D/modphot[blue[ib]].sb0fit_ivar[modgood]+sb_var_floor)
;                radius_kpc = modphot[blue[ib]].radius_kpc[modgood] ; [kpc]
;
;                rmin_kpc = min(radius_kpc)
;                rmax_kpc = max(radius_kpc)
;
;                out[blue[ib]].amin_pixels = amin_kpc/pixscale/arcsec2kpc
;                out[blue[ib]].amax_pixels = amax_kpc/pixscale/arcsec2kpc
;                out[blue[ib]].amin_kpc = amin_kpc
;                out[blue[ib]].amax_kpc = amax_kpc
;                out[blue[ib]].rmin_kpc = rmin_kpc
;                out[blue[ib]].rmax_kpc = rmax_kpc
;             
;                srt = sort(radius_kpc)
;                radius_kpc = radius_kpc[srt]
;                sb = sb[srt]
;                sb_ivar = sb_ivar[srt]
;             
;                if ib eq 0 then begin
;                   fit_sb = sb
;                   fit_sb_ivar = sb_ivar
;                   fit_radius_kpc = radius_kpc
;                   fit_wave = replicate(modphot[blue[ib]].weff,nmodgood)
;                endif else begin
;                   fit_sb = [fit_sb,sb]
;                   fit_sb_ivar = [fit_sb_ivar,sb_ivar]
;                   fit_radius_kpc = [fit_radius_kpc,radius_kpc]
;                   fit_wave = [fit_wave,replicate(modphot[blue[ib]].weff,nmodgood)]
;                endelse
;             endfor
;             
;             bcgmstar_sersic_multiband, fit_radius_kpc, fit_sb, fit_wave, $
;               multisersic, sb_ivar=fit_sb_ivar, verbose=verbose
;             for ib = 0, nblue-1 do begin
;                out[blue[ib]] = im_struct_assign(multisersic[ib],out[blue[ib]],/nozero)
;                out[blue[ib]].sersic_covar = 0
;                out[blue[ib]].sersic_covar_blue = multisersic[ib].sersic_covar
;             endfor
;             
;             if dosersic2 then begin
;                bcgmstar_sersic2_multiband, fit_radius_kpc, fit_sb, fit_wave, $
;                  multisersic2, sb_ivar=fit_sb_ivar, sersicfit=sersicfit, $
;                  fixdevac=fixdevac, verbose=verbose
;                for ib = 0, nblue-1 do begin
;                   out[blue[ib]] = im_struct_assign(multisersic2[ib],out[blue[ib]],/nozero)
;                   out[blue[ib]].sersic2_covar = 0
;                   out[blue[ib]].sersic2_covar_blue = multisersic2[ib].sersic2_covar
;                endfor
;             endif
;          endif
;
;; write out
;          im_mwrfits, out, sersicpath+cluster+'-sersic.fits', clobber=clobber
;       endfor
;    endif
;    

;pro bcgmstar_sersic2_multiband, rr, sb, wave, scoeff, sb_ivar=sb_ivar, $
;  init_params=init_params, fixed=fixed, sersicfit=sersicfit, $
;  fixdevac=fixdevac, verbose=verbose
;; OBSOLETE!
;    
;; fit a double-Sersic function to multiple bands simultaneously; solve
;; for the best-fit half-light radius and Sersic n parameter, but allow
;; the surface brightness at re to vary 
;
;; parse the data
;    good = where(sb_ivar gt 0.0 and finite(sb),nn)
;    xx = rr[good]
;    xxsb = sb[good]
;    xxsb_ivar = sb_ivar[good] 
;
;; figure out how many bands we need to fit and then set up the
;; parameters 
;    nband = n_elements(uniq(wave,sort(wave)))
;    nparam = 2*nband+4
;
;    parinfo = replicate({value: 0D, fixed: 0, $
;      limited: [0,0], limits: [0D,0D]},nparam)
;
;; re1>0
;    parinfo[0].limited = [1,1]
;    parinfo[0].limits = [0.01D,500D]
;; 1<n1<10
;    parinfo[1].limited = [1,1]
;    parinfo[1].limits = [0.1D,10D]
;; re2>0
;    parinfo[2].limited = [1,1]
;    parinfo[2].limits = [0.01D,500D]
;; 1<n2<10
;    parinfo[3].limited = [1,1]
;    parinfo[3].limits = [0.1D,10D]
;
;; sbe1 - surface brightness at re
;    for ib = 0, nband-1 do begin
;       parinfo[4+ib].limited = [1,0]
;;      parinfo[4+ib].limits = [0D,0D]
;       parinfo[4+ib].limits = [1D-12,0D]
;    endfor
;
;; sbe2 - surface brightness at re
;    for ib = 0, nband-1 do begin
;       parinfo[4+nband+ib].limited = [1,0]
;;      parinfo[4+nband+ib].limits = [0D,0D]
;       parinfo[4+nband+ib].limits = [1D-12,0D]
;    endfor
;
;    if n_elements(fixed) eq nparam then parinfo.fixed = fixed
;    if n_elements(init_params) eq nparam then $
;      parinfo.value = init_params else $
;        parinfo.value = [1D,1D,30D,4D,replicate(median(xxsb),nband),$
;      0.1*replicate(median(xxsb),nband)]
;
;    if keyword_set(fixdevac) then begin
;       parinfo[3].value = 4D
;       parinfo[3].fixed = 1
;    endif
;
;;   struct_print, parinfo
;    params = mpfitfun('bcgmstar_sersic2_multiband_func',xx,xxsb,parinfo=parinfo,$
;      yfit=sersicfit,perror=perror,covar=covar,weights=xxsb_ivar,dof=dof,$
;      bestnorm=chi2,status=status,quiet=keyword_set(verbose) eq 0,$
;      functargs={parinfo: parinfo, wave: wave})
;
;    factor = sqrt(chi2/dof)
;    for ib = 0, nband-1 do begin
;       params1 = [params[4+ib],params[0],params[1]]
;       params2 = [params[4+nband+ib],params[2],params[3]]
;
;       perror1 = [perror[4+ib],perror[0],perror[1]]
;       perror2 = [perror[4+nband+ib],perror[2],perror[3]]
;       
;       scoeff1 = {$
;         sersic2_status: status,$
;         sersic2_chi2:     chi2,$
;         sersic2_dof:       dof,$
;         sersic2_covar:   covar,$
;
;         sersic2_sbe1: params1[0],$
;         sersic2_re1:  params1[1],$
;         sersic2_n1:   params1[2],$
;         sersic2_sbe1_err: perror1[0],$ ; *factor,$
;         sersic2_re1_err:  perror1[1],$ ; *factor,$
;         sersic2_n1_err:   perror1[2],$ ; *factor,$
;
;         sersic2_sbe2: params2[0],$
;         sersic2_re2:  params2[1],$
;         sersic2_n2:   params2[2],$
;         sersic2_sbe2_err: perror2[0],$ ; *factor,$
;         sersic2_re2_err:  perror2[1],$ ; *factor,$
;         sersic2_n2_err:   perror2[2]}  ; *factor,$
;       if ib eq 0 then scoeff = scoeff1 else scoeff = [scoeff,scoeff1]
;    endfor
;
;return
;end
;
;pro bcgmstar_sersic_multiband, rr, sb, wave, scoeff, sb_ivar=sb_ivar, $
;  init_params=init_params, fixed=fixed, sersicfit=sersicfit, verbose=verbose
;; OBSOLETE!
;    
;; fit multiple bands simultaneously; solve for the best-fit half-light
;; radius and Sersic n parameter, but allow the surface brightness at
;; re to vary 
;
;; parse the data; take the log
;    good = where(sb_ivar gt 0.0 and finite(sb),nn)
;    xx = rr[good]
;    xxsb = -2.5*alog10(sb[good])
;    xxsb_ivar = sb_ivar[good]*(alog(10)*sb[good]/2.5)^2.0
;    
;; figure out how many bands we need to fit and then set up the
;; parameters 
;    nband = n_elements(uniq(wave,sort(wave)))
;    nparam = nband+2
;
;    parinfo = replicate({value: 0D, fixed: 0, $
;      limited: [0,0], limits: [0D,0D]},nparam)
;
;; re>0
;    parinfo[0].limited = [1,1]
;    parinfo[0].limits = [0.1D,500D]
;; 0<n<10
;    parinfo[1].limited = [1,1]
;    parinfo[1].limits = [0.1D,10D]
;; sbe - surface brightness at re
;    for ib = 0, nband-1 do begin
;       parinfo[2+ib].limited = [1,1]
;       parinfo[2+ib].limits = [10D,35D]
;    endfor
;    
;    if n_elements(init_params) eq nparam then $
;      parinfo.value = init_params else $
;        parinfo.value = [10D, 4D, replicate(median(xxsb),nband)]
;    if n_elements(fixed) eq nparam then parinfo.fixed = fixed
;    
;    params = mpfitfun('bcgmstar_sersic_multiband_func',xx,xxsb,parinfo=parinfo,$
;      perror=perror,covar=covar,weights=xxsb_ivar,bestnorm=chi2,dof=dof,$
;      status=status,quiet=keyword_set(verbose) eq 0,yfit=sersicfit,$
;      functargs={parinfo: parinfo, wave: wave})
;
;    factor = sqrt(chi2/dof)
;    for ib = 0, nband-1 do begin
;       params1 = [params[2+ib],params[0],params[1]]
;;      covar1 = fltarr(3,3)
;;      for ii = 0, 2 do covar1[*,ii] = covar[[0,1,2+ib]]
;       
;       scoeff1 = {$
;         sersic_status: status,$
;         sersic_chi2:     chi2,$
;         sersic_dof:       dof,$
;         sersic_covar:   covar,$
;
;         sersic_sbe: params[2+ib],$
;         sersic_re:  params[0],$
;         sersic_n:   params[1],$
;         sersic_sbe_err: perror[2+ib],$ ; *factor,$
;         sersic_re_err:  perror[0],$ ; *factor,$
;         sersic_n_err:   perror[1]}  ; *factor,$
;
;;        sersic_total: cumsersic_total(params1)}
;       if ib eq 0 then scoeff = scoeff1 else scoeff = [scoeff,scoeff1]
;    endfor
;
;return
;end
;


;; -------------------------    
;; QAplot: color-radius plots, relative to F160W
;    if keyword_set(qaplot_colorradius) then begin
;       for ic = 0, ncl-1 do begin
;          cluster = strtrim(sample[ic].shortname,2)
;          splog, cluster
;          psfile = qapath+'qa_color_sersic_'+cluster+'.ps'
;          im_plotconfig, 0, pos, psfile=psfile, charsize=1.1
;
;          arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
;          
;          sersic = mrdfits(sersicpath+cluster+'-sersic.fits.gz',1,/silent)
;          modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
;          nfilt = n_elements(modphot)-1 ; relative to F160W
;          
;          nrow = ceil(nfilt/float(ncol))
;          pos = im_getposition(nx=ncol,ny=nrow,yspace=0.1,xspace=0.8*[1,1],$
;            xmargin=[1.0,0.2],width=1.9)
;          for ib = 1, nfilt-1 do begin
;             band = strtrim(strupcase(modphot[ib].band),2)
;
;             modgood = where($
;               modphot[ib].majora*pixscale*arcsec2kpc le sersic[ib].amax_kpc and $
;               modphot[ib].sb0fit gt 0 and modphot[ib].sb0fit_ivar gt 0 and $
;               modphot[0].majora*pixscale*arcsec2kpc le sersic[0].amax_kpc and $
;               modphot[0].sb0fit gt 0 and modphot[0].sb0fit_ivar gt 0)
;             modbad = where($
;               (modphot[ib].majora*pixscale*arcsec2kpc gt sersic[ib].amax_kpc or $
;               modphot[0].majora*pixscale*arcsec2kpc gt sersic[0].amax_kpc) and $
;               modphot[ib].sb0fit gt 0 and modphot[ib].sb0fit_ivar gt 0 and $
;               modphot[0].sb0fit gt 0 and modphot[0].sb0fit_ivar gt 0)
;             
;             radius_kpc = modphot[ib].radius_kpc[modgood]  ; [kpc]
;             c1 = -2.5*alog10(modphot[ib].sb0fit[modgood])
;             c2 = -2.5*alog10(modphot[0].sb0fit[modgood])
;             c1err = 2.5/(alog(10)*modphot[ib].sb0fit[modgood]*sqrt(modphot[ib].sb0fit_ivar[modgood]))
;             c2err = 2.5/(alog(10)*modphot[0].sb0fit[modgood]*sqrt(modphot[0].sb0fit_ivar[modgood]))
;;            c1err = modphot[ib].sb0fit_ivar[modgood]*(modphot[ib].sb0fit[modgood]*alog(10)/2.5)^2
;             
;             color = c1 - c2
;             colorerr = sqrt(c1err^2+c2err^2)
;
;             colorbad = -2.5*alog10(modphot[ib].sb0fit[modbad]/modphot[0].sb0fit[modbad])
;
;             yrange = [min(color-colorerr),max(color+colorerr)]
;;            splog, band, yrange
;             
;             sersiccolor = bcgmstar_sersic_func(rr,params=sersic[ib])-$
;               bcgmstar_sersic_func(rr,params=sersic[0])
;             if dosersic2 then sersic2color = -2.5*alog10(bcgmstar_sersic2_func(rr,params=sersic[ib])/$
;               bcgmstar_sersic2_func(rr,params=sersic[0]))
;             
;             if ib-1 eq 1 then title = strupcase(cluster) else delvarx, title
;             if ib-1 ge nfilt-4 then begin
;                xtitle = 'Equivalent Radius (kpc)'
;                delvarx, xtickname
;             endif else begin
;                xtitle = ''
;                xtickname = replicate(' ',10)
;             endelse
;;            if (ib-1 mod 3) eq 0 then begin
;;               delvarx, ytickname
;;            endif else begin
;;               ytickname = replicate(' ',10)
;;            endelse
;             
;             djs_plot, [0], [0], /nodata, /xlog, noerase=ib-1 gt 0, $
;               xrange=[0.3,200], xsty=1, yrange=yrange, position=pos[*,ib-1], $
;               xtickname=xtickname, ytickname=ytickname, title=title, $
;               symsize=0.5, ysty=1, ytitle=band+'-F160W', xtitle=xtitle
;             oploterror, radius_kpc, color, colorerr, psym=3, /nohat, $
;               errcolor=cgcolor('light grey')
;             djs_oplot, radius_kpc, color, psym=symcat(16), symsize=0.5
;
;;            djs_oplot, modphot[ib].radius_kpc[modbad], colorbad, $
;;              psym=symcat(6), symsize=0.3; color=cgcolor('medium grey')
;             
;;            im_legend, band+'-F160W', /right, /top, box=0, margin=0, charsize=1.0
;             djs_oplot, rr, sersiccolor, color=cgcolor('firebrick')
;             if dosersic2 then djs_oplot, rr, sersic2color, color=cgcolor('dodger blue')
;          endfor
;          
;;         xyouts, min(pos[0,*])-0.06, (max(pos[3,*])-min(pos[1,*]))/2.0+min(pos[1,*]), $
;;           textoidl(''), orientation=90, align=0.5, charsize=1.4, /norm
;;         xyouts, (max(pos[2,*])-min(pos[0,*]))/2.0+min(pos[0,*]), min(pos[1,*])-0.06, $
;;           textoidl('Equivalent Radius (kpc)'), align=0.5, charsize=1.1, /norm
;          im_plotconfig, psfile=psfile, /psclose, /pdf
;       endfor
;    endif
;
