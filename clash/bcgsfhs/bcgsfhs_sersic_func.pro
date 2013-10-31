function bcgsfhs_sersic_func, xe, pe, params=params
; k > 0, n > 0
; sb = sb0*exp(-k*r^(1/n))
; pe = [sb0,k,n]
    if n_elements(params) eq 0 then begin
       use_sb0 = pe[0]
       use_k = pe[1]
       use_n = pe[2]
    endif else begin
       use_sb0 = params.sersic_sb0
       use_k = params.sersic_k
       use_n = params.sersic_n
    endelse
;   use_k = 1.0/Re (see Graham & Driver 2005)

;   model = alog(use_sb0)-use_k*xe^(1D/use_n)
;   model = use_sb0*exp(-(use_k*xe)^(1D/use_n))
    model = use_sb0*exp(-get_sersicb(use_n)*((xe/use_k)^(1D/use_n)-1))

return, model
end

