pro smf_wrapper2, iter=iter

	if not keyword_set(iter) then iter=5000

; ----------------------------------------------------------------------
; KECK
; ----------------------------------------------------------------------

; 'NGC2366', 'holmbergii', ,'ic2574'

	galkeck = ['SextansB','NGC3109']
        ngkeck = n_elements(galkeck)
        minkeck = [20.5,20.5]
        
        for j = 0L, ngkeck-1L do begin

            halo = 1L
            if galkeck[j] eq 'NGC2976' then ccut = 0L else ccut = 1L
            edgesmf, galkeck[j], iter=iter, minmag=minkeck[j], halo=halo, ccut=ccut

        endfor
        
stop
return
end
