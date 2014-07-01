function desi_get_hizelg, ugriz, mostek=mostek, johan=johan, magcut=magcut, $
  sigma_kms=sigma_kms
; select high-redshift galaxies using a color-cut
    rr = ugriz[2,*]
    gr = ugriz[1,*]-ugriz[2,*]
    rz = ugriz[2,*]-ugriz[4,*]
    if n_elements(magcut) eq 0 then magcut = 99.0

; fiducial cut    
    int1 = -0.1 & slope1 = 1.0
    int2 = 1.7 & slope2 = -1.0
    hiz = where(rr lt magcut and rz gt 0.2 and rz lt 1.4 and $
      ((rz le 0.9 and gr lt poly(rz,[int1,slope1])) or $
      (rz gt 0.9 and gr lt poly(rz,[int2,slope2]))))

    if n_elements(sigma_kms) ne 0L then begin
       hiz = where(sigma_kms lt 250.0 and rr lt magcut and rz gt 0.2 and rz lt 1.4 and $
         ((rz le 0.9 and gr lt poly(rz,[int1,slope1])) or $
         (rz gt 0.9 and gr lt poly(rz,[int2,slope2]))))
    endif
    
;     (rz gt 0.2 and gr lt poly(rz,[int,slope])<poly(0.9,[int,slope]))))

;    hiz = where(rr lt magcut and rz gt 0.2 and rz lt 1.4 and $
;;     ((rz lt 0.2 and gr lt poly(0.2,[int,slope])) or $
;;     (rz gt 0.2 and gr lt poly(rz,[int,slope])<poly(0.9,[int,slope]))))
;      gr lt poly(rz,[int,slope])<poly(0.9,[int,slope]))

;   hiz = where(rr lt magcut and rz gt 0.1 and rz lt 1.8 and $
;     gr lt poly(rz,[-0.08,1.0])<poly(1.1,[-0.08,1.0]))

;   rzcut = 1.2
;   blue = where(rz lt rzcut,comp=red,nblue,ncomp=nred)
;   if nblue ne 0 then hiz_blue = where(gr[blue] lt (0.7*rz[blue]+0.05))
;   if nred ne 0 then hiz_red = where(gr[red] lt (1.4*rz[red]-0.79))
;   if nblue ne 0 and nred ne 0 then hiz = [blue[hiz_blue],red[hiz_red]]
;   if nblue ne 0 and nred eq 0 then hiz = blue[hiz_blue]
;   if nblue eq 0 and nred ne 0 then hiz = red[hiz_red]
; Mostek's cut    
    if keyword_set(mostek) then hiz = where(gr lt (0.68*rz-0.08) and $
      (rz gt 0.2) and (rz lt 1.3) and rr lt magcut)
        int = 0.1 & slope = 0.55/1.25
; Johan's cut    
    if keyword_set(johan) then hiz = where(gr lt (0.55/1.25*rz+0.1) and $
      (rz gt 0.2) and (rz lt 1.25) and rr lt magcut)
return, hiz
end
