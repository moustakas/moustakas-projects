function get_sfrd_evolution, inzaxis, zaxis=zaxis, zcut=zcut, $
  zform=zform, dz=dz, sfrdmin=sfrdmin, sfrdmax=sfrdmax
; jm10mar30ucsd - parameterize the Madau plot (SFRD vs redshift) in a
;   way that fits the latest 24-micron results from Rujopakarn; the
;   default is to return a "best fit", but we also compute a "maximum"
;   and a "minimum" set of curves that fit the data (1-sigma)

; include a break at z=ZCUT and adopt a formation redshift ZFORM=5 
    if (n_elements(zcut) eq 0) then zcut = 1.0
    if (n_elements(zform) eq 0) then zform = 5.0
    if (n_elements(dz) eq 0) then dz = 0.01
    zaxis = im_array(0.0,zform,dz) ; default redshift grid
    loz = where(zaxis le zcut,nloz,comp=hiz,ncomp=nhiz)

; maximum curve    
    maxsfrd = zaxis*0.0
    maxsfrd[loz] = poly(alog10(1+zaxis[loz]),[-1.95+0.1-alog10(0.66),3.4+0.25])
;   maxsfrd[loz] = poly(alog10(1+zaxis[loz]),[-1.95+0.08-alog10(0.66),3.4+0.2])
;   if (nhiz ne 0) then maxsfrd[hiz] = (-0.4)*alog10((1+zaxis[hiz])/(1+zaxis[nloz-1])) + maxsfrd[nloz-1]
    if (nhiz ne 0) then maxsfrd[hiz] = (+0.0)*alog10((1+zaxis[hiz])/(1+zaxis[nloz-1])) + maxsfrd[nloz-1]

; minimum curve    
    minsfrd = zaxis*0.0
    minsfrd[loz] = poly(alog10(1+zaxis[loz]),[-1.95-0.1-alog10(0.66),3.4-0.25])
;   minsfrd[loz] = poly(alog10(1+zaxis[loz]),[-1.95-0.08-alog10(0.66),3.4-0.2])
    if (nhiz ne 0) then minsfrd[hiz] = (-1.2)*alog10((1+zaxis[hiz])/(1+zaxis[nloz-1])) + minsfrd[nloz-1]
;   if (nhiz ne 0) then minsfrd[hiz] = (-0.5)*alog10((1+zaxis[hiz])/(1+zaxis[nloz-1])) + minsfrd[nloz-1]

; mean curve; use the coefficients from Wiphu's paper,
; scaled to Salpeter     
    sfrd = zaxis*0.0
    sfrd[loz] = poly(alog10(1+zaxis[loz]),[-1.95-alog10(0.66),3.4])
;   if (nhiz ne 0) then sfrd[hiz] = (+0.0)*alog10((1+zaxis[hiz])/(1+zaxis[nloz-1])) + sfrd[nloz-1]
    if (nhiz ne 0) then sfrd[hiz] = (-0.6)*alog10((1+zaxis[hiz])/(1+zaxis[nloz-1])) + sfrd[nloz-1]

    if keyword_set(sfrdmin) then sfrd = minsfrd
    if keyword_set(sfrdmax) then sfrd = maxsfrd
    if (n_elements(inzaxis) ne 0) then sfrd = interpol(sfrd,zaxis,inzaxis)
    
return, sfrd
end
    
    
