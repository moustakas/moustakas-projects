function desi_get_hizelg, cat, mostek=mostek, johan=johan, magcut=magcut, $
  sigma_kms=sigma_kms, zbandlimit=zbandlimit
; jm14jun28siena - select high-redshift galaxies using grz color-cuts;
; note: this code expects a catalog that has been processed through
; DEEP2_GET_UGRIZ() 

    rr = cat.cfhtls_r
    zz = cat.cfhtls_z
    gr = cat.cfhtls_g-cat.cfhtls_r
    rz = cat.cfhtls_r-cat.cfhtls_z
    if n_elements(magcut) eq 0 then magcut = 99.0
    if keyword_set(zbandlimit) then mag = zz else mag = rr

; fiducial cut    
    int1 = -0.1 & slope1 = 1.0
    int2 = 1.7 & slope2 = -1.0
    hiz = where(mag lt magcut and rz gt 0.2 and rz lt 1.4 and $
      ((rz le 0.9 and gr lt poly(rz,[int1,slope1])) or $
      (rz gt 0.9 and gr lt poly(rz,[int2,slope2]))))

; optionally also cut on line-width    
    if n_elements(sigma_kms) ne 0L then begin
       hiz = where(sigma_kms lt 250.0 and rr lt magcut and rz gt 0.2 and rz lt 1.4 and $
         ((rz le 0.9 and gr lt poly(rz,[int1,slope1])) or $
         (rz gt 0.9 and gr lt poly(rz,[int2,slope2]))))
    endif
    
;    hiz = where(rr lt magcut and rz gt 0.2 and rz lt 1.4 and $
;;     ((rz lt 0.2 and gr lt poly(0.2,[int,slope])) or $
;;     (rz gt 0.2 and gr lt poly(rz,[int,slope])<poly(0.9,[int,slope]))))
;      gr lt poly(rz,[int,slope])<poly(0.9,[int,slope]))

; Mostek's proposed cut    
    if keyword_set(mostek) then hiz = where(gr lt (0.68*rz-0.08) and $
      (rz gt 0.2) and (rz lt 1.3) and rr lt magcut)
        int = 0.1 & slope = 0.55/1.25
; Johan's proposed cut    
    if keyword_set(johan) then hiz = where(gr lt (0.55/1.25*rz+0.1) and $
      (rz gt 0.2) and (rz lt 1.25) and rr lt magcut)
return, hiz
end
