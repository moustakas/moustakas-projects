function desi_get_hizelg, ugriz, mostek=mostek, magcut=magcut
; select high-redshift galaxies using a color-cut
    rr = ugriz[2,*]
    gr = ugriz[1,*]-ugriz[2,*]
    rz = ugriz[2,*]-ugriz[4,*]
    if n_elements(magcut) eq 0 then magcut = 99.0
; fiducial cut    
    hiz = where(rr lt magcut and rz gt 0.1 and rz lt 1.8 and $
      gr lt poly(rz,[-0.08,1.0])<poly(1.1,[-0.08,1.0]))
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
return, hiz
end
