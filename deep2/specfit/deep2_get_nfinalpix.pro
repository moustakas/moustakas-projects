pro deep2_get_nfinalpix

    spec1dpath = deep2_path(/dr4)

    if n_elements(zcat) eq 0L then zcat = read_deep2_zcat() ; Q34 sample
    nobj = n_elements(zcat)

    npix = lonarr(nobj)
    
    for iobj = 0, nobj-1 do begin
       spec = fill_gap(strtrim(spec1dpath+zcat[iobj].file,2),/silent)
       obswave = spec.lambda
       obsflux = spec.spec
       obsinvvar = spec.ivar
          
       srt = sort(obswave)
       obswave = obswave[srt]
       obsflux = obsflux[srt]

       good = where((obswave ge zcat[iobj].minwave) and $
         (obswave le zcat[iobj].maxwave),ngood)
       obswave = obswave[good]
       obsflux = obsflux[good]

       log_rebin, minmax(obswave), obsflux, flux, wave, velscale=velscale
       npix[iobj] = n_elements(flux)

       print, 'Object ', iobj, npix[iobj]
    endfor
       
stop

return
end
    
