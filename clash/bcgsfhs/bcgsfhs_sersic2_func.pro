function bcgsfhs_sersic2_func, xe, pe, params=params
; sb = sb01*exp(-k1*r^(1/n1)) + sb02*exp(-k2*r^(1/n2))
; pe = [sb01,k1,n1,sb02,k2,n2]
    if n_elements(params) eq 0 then begin
       use_sb01 = pe[0]
       use_k1 = pe[1]
       use_n1 = pe[2]
       use_sb02 = pe[3]
       use_k2 = pe[4]
       use_n2 = pe[5]
    endif else begin
       use_sb01 = params.sersic2_sb01
       use_k1 = params.sersic2_k1
       use_n1 = params.sersic2_n1
       use_sb02 = params.sersic2_sb02
       use_k2 = params.sersic2_k2
       use_n2 = params.sersic2_n2
    endelse
    model = use_sb01*exp(-use_k1*xe^(1.0/use_n1)) + $
      use_sb02*exp(-use_k2*xe^(1.0/use_n2))
return, model
end

