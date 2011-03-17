; jm00nov12uofa
; wrapper to bootstrap resample the SMF luminosity functions of all
; the galaxies with reasonable detections

pro smf_wrapper1, iter=iter

	if not keyword_set(iter) then iter=5000

; ----------------------------------------------------------------------
; HST
; ----------------------------------------------------------------------

; wait list:  'ugc07577'
        
        galhst = ['ngc1313','ugc03755','ugc06456']
        nghst = n_elements(galhst)
        minhst = [22.0,22.0,22.0]
        
        for k = 1L, nghst-1L do begin

            if (galhst[k] eq 'ugc06456') or (galhst[k] eq 'ugc03755') then ccut = 1L else ccut = 0L
            edgesmf, galhst[k], /hst, iter=iter, minmag=minhst[k], ccut=ccut

        endfor

; ----------------------------------------------------------------------
; KECK
; ----------------------------------------------------------------------

; 'NGC2366', 'holmbergii', ,'ic2574'

	galkeck = ['SextansB','NGC3109','NGC1560','HolmbergIX','NGC2976']
        ngkeck = n_elements(galkeck)
        minkeck = [20.5,20.5,20.5,20.5,20.5]
        
        for j = 0L, ngkeck-1L do begin

            halo = 1L
            if galkeck[j] eq 'NGC2976' then ccut = 0L else ccut = 1L
            edgesmf, galkeck[j], iter=iter, minmag=minkeck[j], halo=halo, ccut=ccut

        endfor
        
stop
return
end
