function stellarlocus_project_kurucz, filterlist
; jm09jul24ucsd - project specified bandpasses onto the Kurucz model
; spectra  
    
    kuruczfile = getenv('IDLSPEC2D_DIR')+'/etc/kurucz_stds_v5.fit'
    flux = mrdfits(kuruczfile,0,hdr,/silent)
    info = mrdfits(kuruczfile,1,/silent)
    keep = where(info.feh eq -1.0 and info.g eq 4.5)
    flux = flux[*,keep]

    nx = n_elements(flux[*,0])
    wave = 10d^(sxpar(hdr, 'CRVAL1') + findgen(nx)*sxpar(hdr, 'CD1_1'))
    nband = n_elements(filterlist)
    nstar = n_elements(flux[0,*])

    lambda = k_lambda_to_edges(wave)
    maggies = fltarr(nband,nstar)
    for i = 0L, nstar-1L do maggies[*,i] = k_project_filters(lambda,$
      flux[*,i],filterlist=filterlist,/silent)
    
return, maggies
end
