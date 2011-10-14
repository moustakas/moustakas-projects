function read_postman, file
; jm11oct13ucsd - read a M. Postman style photometric catalog and put
; it into the standard CLASH catalog format

    if file_test(file) eq 0 then message, 'File '+file+' not found'
    readcol, file, xpos, ypos, rad, zpt, flux, fluxerr, mag, $
      magerr, wave, filename, format='F,F,F,F,F,F,F,F,F,A', $
      comment='#', /silent

    filt = clash_filterlist(short=short)
    nfilt = n_elements(filt)

    phot = {id: 0L}
    for ii = 0, nfilt-1 do phot = create_struct(phot,$
      short[ii]+'_mag',-99.0,short[ii]+'_magerr',-99.0)

; match the filters    
    weff = k_lambda_eff(filterlist=filt)
    get_element, weff, wave, windx
    for ii = 0, n_elements(windx)-1 do begin
       magindx = tag_indx(phot,short[windx[ii]]+'_mag')
       magerrindx = tag_indx(phot,short[windx[ii]]+'_magerr')
       phot.(magindx) = mag[ii]
       phot.(magerrindx) = magerr[ii]
    endfor

return, phot
end
