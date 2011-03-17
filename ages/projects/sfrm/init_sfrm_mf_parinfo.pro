function init_sfrm_mf_parinfo, fixslope=fixslope, quiescent=quiescent,$
  active=active, double_schechter=double_schechter
; jm10feb15ucsd - initialize the MF fitting parameters with the
;   results from BGD08

; ToDo: give the option of starting with the SDSS or AGES results 

    if keyword_set(double_schechter) then begin
       parinfo = {value: 0.0D, limited: [0,0], $
         limits: [0.0D,0.0D], fixed: 0}
       parinfo = replicate(parinfo,5)
       
       parinfo[0].value = 4.04D-3 ; 4.26D-3 ; phi_1
       parinfo[0].limited[0] = 1
       parinfo[0].limits[0] = 1D-8
       
       parinfo[3].value = 3.69D-4 ; 0.58D-3 ; phi_2
       parinfo[3].limited[0] = 1
       parinfo[3].limits[0] = 1D-8
       parinfo[3].fixed = keyword_set(fixslope)
   
       parinfo[1].value = alog10(4.8D10) ; 10.0D^10.648
       parinfo[1].limited = 1
       parinfo[1].limits = [9,12]
   
       parinfo[2].value = -0.4155D ; -0.46D ; alpha
       parinfo[2].fixed = keyword_set(fixslope)
   
       parinfo[4].value = -1.6416D ; -1.58D ; alpha+
       parinfo[4].fixed = keyword_set(fixslope)
    endif else begin
       parinfo = {value: 0.0D, limited: [0,0], $
         limits: [0.0D,0.0D], fixed: 0}
       parinfo = replicate(parinfo,3)
       
       parinfo[0].value = 4.04D-3 ; 4.26D-3 ; phi_1
       parinfo[0].limited[0] = 1
       parinfo[0].limits[0] = 1D-8
       
       parinfo[1].value = alog10(4.8D10) ; 10.0D^10.648
       parinfo[1].limited = 1
       parinfo[1].limits = [9,12]

       parinfo[2].value = -1.07D ; -1.1D
       if keyword_set(quiescent) then parinfo[2].value = -0.1D ; -0.72D
       if keyword_set(active) then parinfo[2].value = -1.29D ; -1.2D
       parinfo[2].fixed = keyword_set(fixslope)
    endelse

return, parinfo
end
