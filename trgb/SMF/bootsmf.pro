pro bootsmf, iter, lf, error, trgblinboot, trgblogboot

	print, 'Generating the bootstrap star lists . . . '
	boot_arrays, iter, lf, bootmags, booterrs

; bootstrap resample

        trgblinboot = fltarr(iter)
        trgblogboot = fltarr(iter)

        for k = 0L, iter-1L do begin

            trgb_lfunction, bootmags[*,k], booterrs[*,k], lfboot, binsize=lf.binsize, $
              minmag=lf.minmag, maxmag=lf.maxmag
            if k mod 100 eq 0L then print, 'Iteration '+strn(k)+'.'

            smf_lfunction, lfboot, error, phiboot, resplin, resplog, noise
            smf_edge, lfboot, phiboot, resplin, resplog, noise, dum1, dum2

            trgblinboot[k] = dum1
            trgblogboot[k] = dum2
            
        endfor

return
end
