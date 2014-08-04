function bcgmstar_sersic2_func, rr, pp, params=params, parinfo=parinfo, allbands=allbands
; pp = [sbe1,re1,n1,sbe2,re2,n2]

    if keyword_set(allbands) then begin
; support input from the bcgmstar_sersic2_allbands_func function; see
; bcgmstar_sersicfit, /qaplot_sbprofiles for proper usage
       if n_elements(params) ne 0 then pp = $
         [params.sersic2_all_sbe1,params.sersic2_all_re1,params.sersic2_all_n1,$
         params.sersic2_all_sbe2,params.sersic2_all_re2,params.sersic2_all_n2]
    endif else begin
       if n_elements(params) ne 0 then pp = $
         [params.sersic2_sbe1,params.sersic2_re1,params.sersic2_n1,$
         params.sersic2_sbe2,params.sersic2_re2,params.sersic2_n2]
    endelse

; make sure the parameters don't go outside the boundaries,
; since mpfit does not 
    use_pp = pp
    if n_elements(parinfo) ne 0 then begin
       for ii = 0, n_elements(pp)-1 do if parinfo[ii].limited[0] then $
         use_pp[ii] = use_pp[ii] > parinfo[ii].limits[0]
       for ii = 0, n_elements(pp)-1 do if parinfo[ii].limited[1] then $
         use_pp[ii] = use_pp[ii] < parinfo[ii].limits[1]
    endif

    use_sbe1 = pp[0]
    use_re1 = pp[1]
    use_n1 = pp[2]
    use_sbe2 = pp[3]
    use_re2 = pp[4]
    use_n2 = pp[5]

; fit the sum of two Sersic models     
    model = use_sbe1*exp(-get_sersicb(use_n1)*((rr/use_re1)^(1D/use_n1)-1D)) + $
      use_sbe2*exp(-get_sersicb(use_n2)*((rr/use_re2)^(1D/use_n2)-1D))

return, model
end

