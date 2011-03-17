pro gto_ancillary_data, gto, write=write
; jm06dec12nyu - written

    outpath = gto_path(/ancillary)

    basicname = 'gto_ned.fits.gz'
    photoname = 'gto_ned_photo.fits.gz'
    distname = 'gto_distances.fits.gz'
    diamname = 'gto_diameters.fits.gz'

;   leda = gto_read_leda()
    
    write_ancillary_data, gto, datapath=outpath, outpath=outpath, $
      basicname=basicname, photoname=photoname, distname=distname, $
      diamname=diamname;, leda=leda

; nice galaxy names

    readcol, outpath+'nice_galaxy_names.txt', galaxy, nice, format='A,A', delimiter='|', /silent

    match, strtrim(strupcase(gto.galaxy),2), strtrim(strupcase(galaxy),2), indx1, indx2
    gto = im_struct_trimtags(gto,select=tag_names(gto),newtags=repstr(tag_names(gto),'LEDA_GALAXY','NICE_GALAXY'))
    gto[indx1].nice_galaxy = nice[indx2]

; write out

    outname = 'gto_ancillary_data.fits'
    
    if keyword_set(write) then begin
       splog, 'Writing '+outpath+outname+'.'
       mwrfits, gto, outpath+outname, /create
       spawn, ['gzip -f '+outpath+outname], /sh
    endif

return
end
