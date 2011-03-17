; jm00nov13uofa
; wrapper to bootstrap resample the maximum likelihood method for all
; the galaxies

pro mlike_wrapper, iter=iter

	if not keyword_set(iter) then iter=5000

; ----------------------------------------------------------------------
; KECK
; ----------------------------------------------------------------------

	galkeck = ['SextansB','NGC3109','NGC1560','NGC2366','holmbergii',$
                   'HolmbergIX','NGC2976','ic2574']
        ngkeck = n_elements(galkeck)
        minkeck = [20.5,20.5,20.0,20.0,19.0,20.0,20.0,20.0]
        
        for j = 0L, ngkeck-1L do begin

            halo = 1L
            if galkeck[j] eq 'NGC2976' then ccut = 0L else ccut = 1L
            trgb_mlike, galkeck[j], iter=iter, minmag=minkeck[j], halo=halo, ccut=ccut

        endfor
        
; ----------------------------------------------------------------------
; HST
; ----------------------------------------------------------------------

        galhst = ['ngc1313','ugc03755','ugc07577','ugc06456']
        nghst = n_elements(galhst)
        minhst = [22.0,22.0,22.0,22.0]
        
        for k = 0L, nghst-1L do begin

            if (galhst[k] eq 'ugc06456') or (galhst[k] eq 'ugc03755') then ccut = 1L else ccut = 0L
            trgb_mlike, galhst[k], /hst, iter=iter, minmag=minhst[k], ccut=ccut

        endfor

stop
return
end
