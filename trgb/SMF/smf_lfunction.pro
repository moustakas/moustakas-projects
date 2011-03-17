pro smf_lfunction, lf, error, phi, resplin, resplog, noise

        phi = lf.hist-lf.hist
        resplin = fltarr(lf.nhist)
        resplog = fltarr(lf.nhist)

        for j = 0L, lf.nhist-1L do begin
            phi[j] = total(1./(sqrt(2.*!pi)*lf.merr)* $
                           exp(-0.5*((lf.mags-lf.mag[j])/lf.merr)^2))
        endfor
            
;       doit = call_external('/deep1/ioannis/trgb/idl/SMF/smf_lf.so', 'smf_idl', $
;                            float(lf.mags), float(lf.merr), float(lf.mag), $
;                            long(lf.nhist), long(nstars), phi)

        for j = 1L, lf.nhist-2L do begin
            getelement_vector, lf.mag, (lf.mag[j]+error[j]), x1
            getelement_vector, lf.mag, (lf.mag[j]-error[j]), x2
            resplin[j] = phi[x1]-phi[x2]			 ; linear
            resplog[j] = alog10(phi[x1]>1.) - alog10(phi[x2]>1.) ; log
        endfor

	smth = smooth2(phi,0.2/lf.binsize) ; smooth the continous luminosity function

        noise = sqrt(2*smth>1.) ; Poisson noise per magnitude bin

;       resplin = resplin/noise
;       resplog = resplog*noise

return
end
