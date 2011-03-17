pro dum

        mag = findgen(500)*0.01+20.
        flux = 10^(-0.4*(30-mag))

        restore, '/deepscr1/ioannis/trgb/ic342/ic342_error_fit.dat'      

        sig = fltarr(n_elements(flux))
        si = fltarr(n_elements(flux),n_elements(a))

        for i = 0L, n_elements(flux)-1L do begin
            for j = 0L, n_elements(a)-1L do begin
                si[i,j] = a[j]*((flux[i])^float(j))
            endfor
            sig[i] = (2.5/alog(10.))*(sqrt(total(si[i,*]))/flux[i])
        endfor

        window, 0, xs=400, ys=400
        plot, mag, sig

stop
return
end
