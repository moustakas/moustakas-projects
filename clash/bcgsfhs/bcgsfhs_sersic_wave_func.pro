function bcgsfhs_sersic_wave_func, xe, pe, params=params

    

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
    model = use_sb0*exp(-use_k*xe^(1.0/use_n))

    
    
return, model
end

