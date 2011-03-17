pro trgb_edge_detect, lf, cutarray, betarray, error_data, lnlike, alpha=alpha

        nstars = long(total(lf.hist))  ; number of stars
        ncut = n_elements(cutarray)    ; number of steps in the TRGB discontinuity
        nbeta = n_elements(betarray)   ; number of steps in the bright end LF parameter

        lnlike = fltarr(lf.nhist,ncut,nbeta) ; 3D likelihood array

; bin the stellar magnitudes

        indx = round((lf.mags-lf.minmag)/lf.binsize)
        indx = indx[where((indx ge 0) and (indx lt lf.nhist),nmbins)]

; shift the origin of the TRGB in magnitude

        marray = fltarr(lf.nhist,lf.nhist)
        for k = 0L, lf.nhist-1L do marray[*,k] = lf.mag-lf.mag[k]

; iterate on the bright end LF parameter and on the width of the discontinuity

        for i = 0L, nbeta-1L do begin

            beta = betarray[i]

            for j = 0L, ncut-1L do begin

                cutfactor = cutarray[j]
                g = trgb_lfmodel(error_data,marray,cutfactor,beta,alpha=alpha)
                lnlike[*,j,i] = - nstars * alog10(total(g,1)) + total(alog10(g[indx,*]),1)

            endfor

        endfor

        lnlike = lnlike - max(lnlike)	; normalize to zero

return
end

