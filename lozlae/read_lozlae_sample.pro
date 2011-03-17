function read_lozlae_sample
; jm10dec06ucsd - read the sample and pack into a useful format

    file = lozlae_path()+'hecto.fits'
    data = mrdfits(file,1)
    data.id[10] = 23 ; hack!

    ngal = n_elements(data.id)
    npix = n_elements(data.lambda)
    
    spec1d = replicate({id: 0, z: 0.0, wave: data.lambda, $
      flux: fltarr(npix), ferr: fltarr(npix)},ngal)
    spec1d.id = data.id
    spec1d.z = data.z

    for ii = 0, ngal-1 do begin
       spec1d[ii].flux = data.(tag_indx(data,'SPEC_EGS'+strtrim(spec1d[ii].id,2)))
       spec1d[ii].ferr = data.(tag_indx(data,'DSPEC_EGS'+strtrim(spec1d[ii].id,2)))
    endfor

return, spec1d
end
