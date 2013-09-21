;+
; NAME:
;   bgt_sersic_fit1d
; CALLING SEQUENCE:
; INPUT:
; OUTPUT:
; COMMENTS:
;   -- not in ln scale, R vs I
; REVISION HISTORY:
;   ??-NOV-2008, Started, Guangtun Zhu, NYU
;   15-Mar-2009, sort of documented, Guangtun Zhu, NYU
;-

; fit a Sersic + disk

pro streams_sersic_fitbd, rr, sb, scoeff, sb_err=sb_err, $
  sb_ivar=sb_ivar, ini_value=ini_value, fixed=fixed, $
  sersicfit=sersicfit

   parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                        limits:[0.D,0.D]}, 5)

;; sb0 ~ max(sb0)
   parinfo[0].limited = [1,0]
   parinfo[0].limits = [0.D,0.D]
;; k>0
   parinfo[1].limited = [1,0]
   parinfo[1].limits = [0.D,0.D]
;; 2<n<6
   parinfo[2].limited = [1,1]
   parinfo[2].limits = [1.D,10.D]
;; sb0 for exp
   parinfo[3].limited = [1,0]
   parinfo[3].limits = [0.0D,0.D]
;; kk(1/rd) for exp
   parinfo[4].limited = [1,0]
   parinfo[4].limits = [0.D,0.0D]

   if (not keyword_set(ini_value)) then begin
      parinfo[*].value = [max(sb), 4., 4., 0.1*max(sb), 1./30.]
   endif else begin
      parinfo[*].value = ini_value
   endelse
   if (not keyword_set(fixed)) then begin
      parinfo[*].fixed = [0,0,0,0,0]
   endif else begin
      parinfo[*].fixed = fixed
   endelse
   ii = where(sb_ivar gt 0. and finite(sb) eq 1, nn)
   if nn gt 2 then begin
      xx = rr[ii]
      xxsb = sb[ii]
      xxsb_ivar = sb_ivar[ii]*xx^1.5

      params= MPFITFUN('bd_func', xx, xxsb, xxsb_err, parinfo=parinfo, yfit=sersicfit, $
                      perror=perror, covar=covar, weights=xxsb_ivar, status=status, /quiet)
   endif else status = -1

   if(n_tags(scoeff) eq 0) then $
   scoeff={bd_sersic_sb0:max(sb), bd_sersic_k:0.d, bd_sersic_n:0.d, $
            bd_sersic_sb0_err:0., bd_sersic_k_err:0., bd_sersic_n_err:0., $
            bd_exp_sb0:0.1*max(sb), bd_exp_k:0.d, $
            bd_exp_sb0_err:0.1*max(sb), bd_exp_k_err:0.d, $
            bd_covar:fltarr(5,5), bd_status:status, $
            bd_sersic_total: 0., bd_exp_total:0.}

   if (status le 0) then begin
       params = parinfo[*].value
       perror = fltarr(5)
       covar = fltarr(5,5)
       status = -1
   endif
   scoeff.bd_sersic_sb0 = params[0]
   scoeff.bd_sersic_k = params[1]
   scoeff.bd_sersic_n = params[2]
   scoeff.bd_exp_sb0 = params[3]
   scoeff.bd_exp_k = params[4]
   scoeff.bd_sersic_sb0_err = perror[0]
   scoeff.bd_sersic_k_err = perror[1]
   scoeff.bd_sersic_n_err = perror[2]
   scoeff.bd_exp_sb0_err = perror[3]
   scoeff.bd_exp_k_err = perror[4]

   scoeff.bd_covar = covar
   cum_params = params[0:2]
   cum_params1 = [params[3:4],1]
;  cum_params[0] = exp(params[0])
   scoeff.bd_sersic_total = cumsersic_total(cum_params)
   scoeff.bd_exp_total = cumsersic_total(cum_params1)
end
