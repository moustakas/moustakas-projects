function mpfunc, m, p

        model = alog10(p[3]*(10.^(p[0]*m)*(m ge 0.)+10.^(p[1]*m-p[2])*(m lt 0.)))
	return, model

end

pro chisqtest, objname

	colortable1

; good objects to run this test on
; ----------------------------------------------------------------------
        objname = 'SextansB'
        trgb = 21.74
        binsize = 0.01
        minmag = 20.
        maxmag = 23.4
; ----------------------------------------------------------------------
        
; read in the data

	trgb_readata, objname, datapath, data, infobase, /halo, /ccut

; generate a luminosity function and smooth it

        trgb_lfunction, data[3,*], data[4,*], lf, binsize=binsize, minmag=minmag, maxmag=maxmag

        linsmooth = smooth2(float(lf.hist),4)
        noise = sqrt(2.*(linsmooth>1.))		; Poisson noise in each bin
        weights = 1D/noise^2

        init = dblarr(4) ; best initial guesses for the fit

        init[0] = 0.3  ; alpha
        init[1] = 0.1  ; beta
        init[2] = 0.16 ; cutfactor
        init[3] = 10.  ; normalization
        
        params = mpfitfun('mpfunc', lf.mag-trgb, alog10(float(lf.hist)>1.), noise, $
                           init, bestnorm=bestnorm, covar=cov, $
                           maxiter=200, niter=ni, yfit=yfit, $
                           perror=perror, weights=weights)

        dof = lf.nhist-n_elements(params) ; degrees of freedom
        chi_dof = bestnorm/dof
        sigs = perror*sqrt(chi_dof)

        model = mpfunc(lf.mag-trgb,params)

        window, 0, xs=450, ys=450
        plot, lf.mag-trgb, alog10(float(lf.hist)>1.), col=7
        oplot, lf.mag-trgb, model, colo=5
  
        
stop
return
end
