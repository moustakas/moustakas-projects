pro test

    common sirtf_simulations, sirtf

;   sirtfz_initialize, filters=['U Johnson','B Kitt Peak','V Kitt Peak'], templates='devriendt'

    obands = filter_match(*sirtf.filters,sirtf.bandcube.bandnames)
    nbands = size(obands,/n_elements)

;   nseds = size(sirtf.sedcube,/n_elements)

    sedindx = [5,8,9,12]
    nseds = n_elements(sedindx)
    
    z = 0.0

    flux = fltarr(nbands,nseds)
    for i = 0L, nbands-1L do begin

       wband = *sirtf.bandcube[obands[i]].wband 
       rband = *sirtf.bandcube[obands[i]].rband

       for j = 0L, nseds-1L do begin
       
          wave = reform(*sirtf.sedcube[sedindx[j]].lambda)
          mlum = reform(*sirtf.sedcube[sedindx[j]].mlum)
    
          flux[i,j] = band_flux(wave,mlum,wband,rband,z=z)

       endfor

    endfor

    ub = -2.5*alog10(flux[0,*]/flux[1,*]) & print, 'U-B = ', ub
    bv = -2.5*alog10(flux[1,*]/flux[2,*]) & print, 'B-V = ', bv

;   plot, alog10(wave), alog10(mlum), /xsty, /ysty
;   oplot, alog10(*sirtf.bandcube[obands[1]].wband), *sirtf.bandcube[obands[1]].rband, line=2
;   oplot, alog10(*sirtf.bandcube[obands[2]].wband), *sirtf.bandcube[obands[2]].rband, line=2

;   jm = flux[*,0]
;   print, -2.5*alog10(jm[0]) - 14.08 + 2.5*alog10(jm[1])+13.00 ; should be 1.26
    
stop

return
end
