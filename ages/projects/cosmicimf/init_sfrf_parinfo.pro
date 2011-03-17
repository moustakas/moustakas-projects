function init_sfrf_parinfo, fixslope=fixslope, schechter=schechter
; jm10mar18ucsd - initialize the SFRF fitting parameters

    if keyword_set(schechter) then begin
       parinfo = {value: 0.0D, limited: [0,0], $
         limits: [0.0D,0.0D], fixed: 0}
       parinfo = replicate(parinfo,3)
       
       parinfo[0].value = 3D-3 ; ; phi
       parinfo[0].limited[0] = 1
       parinfo[0].limits[0] = 1D-8
       
       parinfo[1].value = 5D
       parinfo[1].limited = 1
       parinfo[1].limits = [1D,12D]

       parinfo[2].value = -1.4D
       parinfo[2].fixed = keyword_set(fixslope)
    endif else begin ; double power-law
       parinfo = {value: 0.0D, limited: [0,0], $
         limits: [0.0D,0.0D], fixed: 0}
       parinfo = replicate(parinfo,4)
       
       parinfo[0].value = 3D-3 ; phi*
       parinfo[0].limited[0] = 1
       parinfo[0].limits[0] = 1D-8
       
       parinfo[1].value = 5D ; SFR*
       parinfo[1].limited = 1
       parinfo[1].limits = [1D,12D]
   
       parinfo[2].value = -0.4D ; alpha
       parinfo[2].fixed = keyword_set(fixslope)
   
       parinfo[3].value = -1.8D ; beta
       parinfo[3].fixed = 1 ; keyword_set(fixslope)
    endelse

return, parinfo
end
