pro benitez
; jm01mar27uofa

    z = findgen(101)/20. ; redshift
    m0 = findgen(8)+21.0

; from Table 1 of Benitez 2000
    
    ell = 0.35*exp(-0.147*(m0-20))
    sbc = 0.50*exp(-0.450*(m0-20))
    irr = 1.0-(ell+sbc)

    plot, m0, ell, line=0, xsty=3, ysty=3, yr=[0,1.0], $
      xtit='Apparent magnitude', ytit='Galaxy fraction'
    oplot, m0, sbc, line=1
    oplot, m0, irr, line=2
    legend, ['E/S0','Sbc,Scd','Irr'], line=[0,1,2], box=0, /left, /top

    zmtell = 0.431+0.0910*(m0-20)
    zmtsbc = 0.390+0.0636*(m0-20)
    zmtirr = 0.063+0.1230*(m0-20)
    
    pztm0ell = fltarr(n_elements(m0),n_elements(z))
    pztm0sbc = fltarr(n_elements(m0),n_elements(z))
    pztm0irr = fltarr(n_elements(m0),n_elements(z))

    for j = 0L, n_elements(m0)-1L do begin
       pztm0ell[j,*] = z^(2.46)*exp(-(z/zmtell[j])^(2.46))
       pztm0sbc[j,*] = z^(1.81)*exp(-(z/zmtsbc[j])^(1.81))
       pztm0irr[j,*] = z^(0.91)*exp(-(z/zmtirr[j])^(0.91))
    endfor
    
    pzm0 = fltarr(n_elements(m0),n_elements(z))

    for j = 0L, n_elements(m0)-1L do pzm0[j,*] = ell[j]*pztm0ell[j,*] + sbc[j]*pztm0sbc[j,*] + irr[j]*pztm0irr[j,*]

; attempting to remake figure 4 from benitez 2000 but didn't work
    
    plot, z, pzm0[0,*], line=0, xsty=3, ysty=3, yr=[0,0.15]
    for j = 1L, 7L do oplot, z, pzm0[j,*]
    
stop   

return
end
