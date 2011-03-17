; jm00nov12uofa
; wrapper to bootstrap resample the SMF luminosity functions of all
; the galaxies with reasonable detections

pro smf_wrapper1, iter=iter

	if not keyword_set(iter) then iter=5000

; ----------------------------------------------------------------------
; HST
; ----------------------------------------------------------------------

; wait list:  'ugc07577'
; done: 'ugc03755'
        
        galhst = ['ugc06456','ngc1313','ugc07577']
        nghst = n_elements(galhst)
        minhst = [23.0,22.0,22.0]
        
        for k = 2L, nghst-1L do begin

            if (galhst[k] eq 'ugc06456') or (galhst[k] eq 'ugc03755') then ccut = 1L else ccut = 0L
            edgesmf, galhst[k], /hst, iter=iter, minmag=minhst[k], ccut=ccut, maxmag=24.5

        endfor

stop
return
end
