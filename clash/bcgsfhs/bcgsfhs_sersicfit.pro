pro bcgsfhs_sersic2, rr, sb, scoeff, sb_err=sb_err, $
  sb_ivar=sb_ivar, init_params=init_params, fixed=fixed, $
  sersicfit=sersicfit
; fit a double Sersic function

    parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
      limits:[0D,0D]},6)

; sb01 > 0 
    parinfo[0].limited = [1,0]
    parinfo[0].limits = [0D,0D]
; k1>0
    parinfo[1].limited = [1,0]
    parinfo[1].limits = [0D,0D]
; 1<n1<10
    parinfo[2].limited = [1,1]
    parinfo[2].limits = [1D,10D]
; sb02 > 0
    parinfo[3].limited = [1,0]
    parinfo[3].limits = [0D,0D]
; k2>0
    parinfo[4].limited = [1,0]
    parinfo[4].limits = [0D,0D]
; 1<n2<10
    parinfo[5].limited = [1,1]
    parinfo[5].limits = [1D,10D]

    if n_elements(init_params) eq 0 then begin
       parinfo.value = [max(sb), 4.0, 4.0, 0.1*max(sb), 4.0, 4.0]
    endif else parinfo.value = init_params

    if n_elements(fixed) ne 0 then parinfo.fixed = fixed

    good = where(sb_ivar gt 0.0 and finite(sb) eq 1, nn)
    xx = rr[good]
    xxsb = sb[good]
    xxsb_ivar = sb_ivar[good]
    
    params = mpfitfun('bcgsfhs_sersic2_func',xx,xxsb,parinfo=parinfo,$
      yfit=sersicfit,perror=perror,covar=covar,weights=xxsb_ivar,$
      status=status,quiet=1)
    
    scoeff = {$
      sersic2_sb01: max(sb), sersic2_k1: 0D, sersic2_n1: 0D, $
      sersic2_sb01_err: 0.0, sersic2_k1_err: 0.0, sersic2_n1_err: 0.0, $
      sersic2_sb02: max(sb), sersic2_k2: 0D, sersic2_n2: 0D, $
      sersic2_sb02_err: 0.0, sersic2_k2_err: 0.0, sersic2_n2_err: 0.0, $
      sersic2_covar: fltarr(6,6), sersic2_status: status, $
      sersic2_total1: 0.0, sersic2_total2: 0.0}
    
    scoeff.sersic2_sb01 = params[0]
    scoeff.sersic2_k1 = params[1]
    scoeff.sersic2_n1 = params[2]
    scoeff.sersic2_sb02 = params[3]
    scoeff.sersic2_k2 = params[4]
    scoeff.sersic2_n2 = params[5]
    
    scoeff.sersic2_sb01_err = perror[0]
    scoeff.sersic2_k1_err = perror[1]
    scoeff.sersic2_n1_err = perror[2]
    scoeff.sersic2_sb02_err = perror[3]
    scoeff.sersic2_k2_err = perror[4]
    scoeff.sersic2_n2_err = perror[5]
    
    scoeff.sersic2_covar = covar
    scoeff.sersic2_total1 = cumsersic_total(params[0:2])
    scoeff.sersic2_total2 = cumsersic_total(params[3:5])
    
return
end

pro bcgsfhs_sersic, rr, sb, scoeff, sb_err=sb_err, $
  sb_ivar=sb_ivar, init_params=init_params, fixed=fixed, $
  sersicfit=sersicfit
; fit a single Sersic function
    
    parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
      limits:[0.D,0.D]}, 3)

; sb0 > max(sb0)
    parinfo[0].limited = [1,0]
    parinfo[0].limits = [max(sb),0D]
; k>0
    parinfo[1].limited = [1,0]
    parinfo[1].limits = [0D,0D]
; n>0, <10
    parinfo[2].limited = [1,1]
    parinfo[2].limited = [1D,10D]
   
    if n_elements(init_params) eq 0 then begin
       parinfo.value = [max(sb), 4.0, 4.0]
    endif else parinfo.value = init_params 

    if n_elements(fixed) ne 0 then parinfo.fixed = fixed

    good = where(sb_ivar gt 0. and finite(sb) eq 1, nn)
    xx = rr[good]
    xxsb = sb[good]
    xxsb_ivar = sb_ivar[good] 
 
    params = mpfitfun('bcgsfhs_sersic_func',xx,xxsb,parinfo=parinfo,$
      perror=perror,covar=covar,weights=xxsb_ivar,$
      status=status,quiet=1,yfit=sersicfit)
    scoeff = {sersic_sb0: max(sb), sersic_k:0.d, sersic_n: 0D, $
      sersic_sb0_err:0., sersic_k_err:0., sersic_n_err: 0.0, $
      sersic_covar:fltarr(3,3), sersic_status: status, $
      sersic_total: 0.}

    scoeff.sersic_sb0 = params[0]
    scoeff.sersic_k = params[1]
    scoeff.sersic_n = params[2]
    scoeff.sersic_sb0_err = perror[0]
    scoeff.sersic_k_err = perror[1]
    scoeff.sersic_n_err = perror[2]
    scoeff.sersic_covar = covar
    scoeff.sersic_total = cumsersic_total(params)

return
end

pro bcgsfhs_sersicfit, debug=debug
; jm13oct22siena - fit various Sersic models to the output of
; BCGSFHS_ELLIPSE 

; read the sample
    sample = read_bcgsfhs_sample(/noa2261)
;   struct_print, sample
    ncl = n_elements(sample)

    pixscale = 0.065D                      ; [arcsec/pixel]

; wrap on each cluster    
    for ic = 3, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       splog, 'Working on cluster '+cluster
       datapath = bcgsfhs_path(/bcg)+cluster+'/'

       modphot = mrdfits(datapath+cluster+'-ellipse-model.fits.gz',1,/silent)
       nfilt = n_elements(modphot)

       for ii = 0, nfilt-1 do begin ; just fit F160W
;      for ii = 0, 0 do begin ; just fit F160W
          band = strtrim(strupcase(modphot[ii].band),2)

          modgood = where(modphot[ii].sb0fit gt 10^(-0.4*modphot[ii].sblimit))
          radius = modphot[ii].radius[modgood]*pixscale ; [arcsec]
          radius_kpc = modphot[ii].radius_kpc[modgood]  ; [kpc]
          sb = modphot[ii].sb0fit[modgood]
          sb_ivar = modphot[ii].sb0fit_ivar[modgood]

; fit with a single-Sersic and then a double-Sersic
          bcgsfhs_sersic, radius_kpc, sb, sersic, sb_ivar=sb_ivar
;; this (working) code is to use the F160W to constrain the free
;; parameters of the Sersic model
;          if ii gt 0 then begin
;             bcgsfhs_sersic, radius_kpc, sb, sersic, sb_ivar=sb_ivar, $
;               init_params=[out[0].sersic_sb0,out[0].sersic_k,$
;               out[0].sersic_n], fixed=[0,1,1]
;          endif else begin
;             bcgsfhs_sersic, radius_kpc, sb, sersic, sb_ivar=sb_ivar
;          endelse
          bcgsfhs_sersic2, radius_kpc, sb, sersic2, sb_ivar=sb_ivar

; do photometry in radial apertures, using the Sersic model to
; extrapolate inward
          int_radius_kpc = [0,range(min(radius_kpc)*1E-3,min(radius_kpc)*0.95,30,/log)]
          int_sb = [bcgsfhs_sersic_func(int_radius_kpc,params=sersic),sb]
          tot = -2.5*alog10(2.0*!pi*im_integral(int_radius_kpc,int_radius_kpc*int_sb))
          print, band, tot
          
; pack into a structure          
          if ii eq 0 then out = struct_addtags(sersic,sersic2) else $
            out = [out,struct_addtags(sersic,sersic2)]
          
          if keyword_set(debug) then begin
;            help, sersic, sersic2, /str
;            splog, band, sersic.sersic_n, sersic2.sersic2_n1, sersic2.sersic2_n2

             rr = [0,range(0.01,200,500,/log)]
             djs_plot, radius_kpc, -2.5*alog10(sb), psym=8, /xlog, $
               xrange=[0.1,200], xsty=1, yrange=[30,16]
             djs_oplot, rr, -2.5*alog10(bcgsfhs_sersic_func(rr,params=out[ii])), $
               color=cgcolor('yellow')
             if ii gt 0 then djs_oplot, rr, -2.5*alog10(bcgsfhs_sersic_func(rr,$
               params=out[0])), color=cgcolor('forest green')
             
             djs_oplot, rr, -2.5*alog10(bcgsfhs_sersic2_func(rr,params=out[ii])), $
               color=cgcolor('red')
             djs_oplot, rr, -2.5*alog10(bcgsfhs_sersic_func(rr,[out[ii].sersic2_sb01,$
               out[ii].sersic2_k1,out[ii].sersic2_n1])), color=cgcolor('orange'), line=5
             djs_oplot, rr, -2.5*alog10(bcgsfhs_sersic_func(rr,[out[ii].sersic2_sb02,$
               out[ii].sersic2_k2,out[ii].sersic2_n2])), color=cgcolor('orange'), line=5
             cc = get_kbrd(1)
          endif
       endfor
; write out
       im_mwrfits, out, datapath+cluster+'-sersic.fits', clobber=clobber
    endfor

;;         bgt_ellipse_sersicradius, ellipse, outradius=outrad
;;         bgt_ellipse_radius() ; get the half-light radius

return
end
