function read_irclusters, maggies=maggies, ivarmaggies=ivarmaggies
; jm10nov09ucsd - read the redshift and photometric catalog 
    irpath = irclusters_path()
    zcat = rsex(irpath+'FullSample_v1.cat')
    phot = rsex(irpath+'FullSample.flux_v1')
    cat = struct_addtags(zcat,phot)

    splog, 'Enforcing z>0.01 minimum redshift!'
    cat.z = cat.z>0.01D

; fluxes are in microJy, so to convert to maggies multiply by this factor:    
    factor = 10D^(-0.4*23.9)
    maggies = transpose([[cat.bw],[cat.r],[cat.i],[cat.j],[cat.h],$
      [cat.ks],[cat.ch1],[cat.ch2],[cat.ch3],[cat.ch4]])*factor
    dflux = transpose([[cat.dbw],[cat.dr],[cat.di],[cat.dj],[cat.dh],[cat.dks],$
      [cat.dch1],[cat.dch2],[cat.dch3],[cat.dch4]])
    ivarmaggies = 1.0/((dflux*factor)^2.0+(dflux eq 0))*(dflux ne 0)
    
return, cat
end
