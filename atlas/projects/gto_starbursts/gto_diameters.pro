pro gto_diameters, diameters, write=write
; jm06dec12nyu - written
    
    analysis_path = gto_path(/ancillary)
    gto = mrdfits(analysis_path+'gto_ned.fits.gz',1,/silent)

    galaxy = strtrim(gto.nedgalaxy,2)
    ngalaxy = n_elements(galaxy)

    outfile = 'gto_diameters.fits'
    
    ned_webget_diameters, galaxy, diameters, outfile=outfile, $
      outpath=analysis_path, write=write

return
end
    
