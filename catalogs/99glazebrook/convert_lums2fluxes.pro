pro convert_lums2fluxes
; jm04dec31uofa
; convert the luminosities in Table 2 (lums.dat) of Glazebrook et
; al. (1999) to fluxes, which can be tabulated in 99glazebrook.dat,
; and subsequently parsed by write_99glazebrook

; cosmology: H0=50 km/s/Mpc    

    red, h100=0.5, omega_lambda=0.0, omega_matter=1.0
    
    readcol, 'lums.dat', cfrs, ha, eha, oii, eoii, z, $
      format='A,D,D,D,D,D', /silent, comment='#'

    ha = ha*1D41 & eha = eha*1D41
    oii = oii*1D41 & eoii = eoii*1D41

    d = dluminosity(z,/cm)
    area = 4.0 * !dpi * d * d

    fha = ha / area   & efha = eha / area
    foii = oii / area & efoii = eoii / area

    niceprint, cfrs, ha, 1.6D17*fha;, efha

return
end
    
