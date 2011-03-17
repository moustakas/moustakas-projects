pro trgb_edge_detect, lf, cutarray, lnlike, alpha=alpha, beta=beta, gamma=gamma

        nstars = long(total(lf.hist))  ; number of stars
        ncut = n_elements(cutarray)    ; number of steps in the TRGB discontinuity
        total = fltarr(lf.nhist)   
        logtot = fltarr(lf.nhist)
        lnlike = fltarr(lf.nhist,ncut) ; 2D likelihood array

; bin the stellar magnitudes

        indx = fix((lf.mags-lf.minmag)*(1./lf.binsize))
        indx = indx[where((indx ge 0) and (indx lt lf.nhist),nmbins)]

; shift the origin of the TRGB in magnitude

        marray = fltarr(lf.nhist,lf.nhist)
        for k = 0L, lf.nhist-1L do marray[*,k] = lf.mag-lf.mag[k]

; iterate on the width of the discontinuity

        for j = 0L, ncut-1L do begin

            cutfactor = cutarray[j]

            g = form_g(marray,cutfactor,alpha=alpha,beta=beta,gamma=gamma)

            lnlike[*,j] = - nstars * alog10(total(g,1)) + total(alog10(g[indx,*]),1)

;           doit = call_external('/deep1/ioannis/trgb/idl/LIKELIHOOD/mlike.so', $
;                                'mlike_idl', float(g), long(indx), long(lf.nhist), $
;                                long(nmbins), total, logtot)
;           lnlike[*,j] = - nstars * alog10(total) + logtot

        endfor

        lnlike = lnlike - max(lnlike)	; normalize to zero

return
end

