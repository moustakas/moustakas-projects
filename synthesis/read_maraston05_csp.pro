function read_csp, file
; jm11oct05ucsd 
    
    massfile = repstr(file,'.sed_agb','.mass_agb')
    readfast, massfile, allmass, skip=3
    
    readfast, file, data, skip=3
    allage = reform(data[0,*])
    age = allage[uniq(allage,sort(allage))]
    nage = n_elements(age)

    for ii = 0L, nage-1 do begin
       these = where(age[ii] eq allage,npix)
       if (ii eq 0) then begin
          mara = {$
            age:            age*1D9,$ ; [yr]
            mstar:     fltarr(nage),$ ; [Msun]
            wave: float(reform(data[1,these])),$
            flux:      fltarr(npix,nage)} ; [erg/s/AA/Msun]
       endif
       mara.flux[*,ii] = reform(data[2,these])
;      if ii eq 60 then stop
; get the total stellar mass (living stars plus remnants)  
       mara.mstar = interpol(allmass[1,*],allmass[0,*]*1D9,mara.age)
    endfor
    
return, mara
end

    
