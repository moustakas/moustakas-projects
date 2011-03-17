pro sirtf_features
; jm01may23uofa
; at what redshift do various SED features move into each of the SIRTF
; bands?

    lambda0 = [0.658, 2.2, 3.6, 4.5, 5.8, 8.0, 24.0, 70.0, 160.0] ; SIRTf bands (mu)
    nbands = n_elements(lambda0)
    
    keen = [0.0912, 0.4000, $         ; lyman-limit and 4000-Angstrom break
            3.3, 6.2, 7.7, 8.6, 11.3] ; PAH features
    nkeen = n_elements(keen)
    
    zshow = fltarr(nbands,nkeen)

    for j = 0L, nkeen-1L do begin

       zshow[*,j] = lambda0/keen[j] - 1.0
       neg = where(zshow[*,j] lt 0.0,nneg)
       if nneg ne 0L then zshow[neg,j] = 0.0
       
    endfor

    print, zshow

stop

return
end
