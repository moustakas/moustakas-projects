function trgb_error_func, objname, datapath, mags, minmag=minmag, maxmag=maxmag, $
                          hst=hst, halo=halo, core=core, plot=plot

; generate a smoothed magnitude error function 

        flux = 10.^(-0.4*(mags-30.)) ; convert to data numbers

; read in the polynomial fit values

        if keyword_set(halo) then filename = datapath+'/'+objname+'_halo_error_fit.dat' else $
          if keyword_set(core) then filename = datapath+'/'+objname+'_core_error_fit.dat' else $
          filename = datapath+'/'+objname+'_error_fit.dat'

        restore, filename

; form the error function in terms of magnitude and magnitude errors

        sig = fltarr(n_elements(flux))
        si = fltarr(n_elements(flux),n_elements(a))
        for i = 0L, n_elements(flux)-1L do begin
            for j = 0L, n_elements(a)-1L do $
              si[i,j] = a[j]*((flux[i])^float(j))
            sig[i] = (2.5/alog(10.))*(sqrt(total(si[i,*]))/flux[i])
        endfor
        
        if keyword_set(plot) then begin
            window, 6, xs=400, ys=400
            plot, mags, sig
        endif
    
return, sig
end

