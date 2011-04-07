pro build_sings_ppxf_templates, debug=debug
; jm10feb25ucsd - write out the basic set of BC03 templates for the
;   PPXF template fitting

    outpath = sings_path(/ppxf)
    version = sings_version(/ppxf_templates)

    light = 299792.458D
    lsun = 3.826D33         ; [erg/s]
    dist = 10.0*3.085678D18 ; 10 pc [cm]

    agegrid = [5.0,25.0,50.0,125.0,300.0,650.0,$
      1500.0,3000.0,7000.0,13000.0]
    metal = [0.004,0.02,0.05] ; metallicity grid
    nmodel = n_elements(agegrid)
    nmetal = n_elements(metal)

    bc03_Z004 = im_read_bc03(age=agegrid*1E6/1E9,$ ; chabrier IMF
      bc03_extras=extras,/silent,metallicity=2)
    bc03_Z02 = im_read_bc03(age=agegrid*1E6/1E9,$
      bc03_extras=extras,/silent,metallicity=4)
    bc03_Z05 = im_read_bc03(age=agegrid*1E6/1E9,$
      bc03_extras=extras,/silent,metallicity=5)
    bc03 = [bc03_Z004,bc03_Z02,bc03_Z05]

; resample to be constant in linear wavelength    
    minwave = 912.0D
    maxwave = 1D5 ; 3D4
;   minwave = 3500.0D ; v2.0 template wavelength range
;   maxwave = 7000.0D
    dwave = 1.0D
    wave = im_array(minwave,maxwave,dwave,/double)
    npix = n_elements(wave)

    outfile = outpath+'bc03_'+string(metal,format='(F5.3)')+$
      '_'+version+'.fits'
    for jj = 0, nmetal-1 do begin
       flux = fltarr(npix,nmodel)
       for ii = 0, nmodel-1 do begin
          linterp, bc03[jj].wave, bc03[jj].flux[*,ii], wave, newflux
          flux[*,ii] = lsun*newflux/extras[ii].m_/(4.0*!dpi*dist^2)
          if keyword_set(debug) then begin
             djs_plot, wave, newflux, xr=[1000,4E4], /xlog, /ylog
             djs_oplot, bc03[jj].wave, bc03[jj].flux[*,ii], color='orange'
             djs_oplot, wave, newflux1, color='blue'
             cc = get_kbrd(1)
          endif
       endfor
       
       info = {metallicity: metal[jj], age: fltarr(nmodel)}
       info.age = bc03[jj].age/1E9 ; [Gyr]
       
; instrumental velocity dispersion of the models    
       inst_vdisp = 60.0        ; [km/s]
    
; write out; generate an appropriate FITS header
       mkhdr, hdr, float(flux)
       sxdelpar, hdr, 'COMMENT'
       sxaddpar, hdr, 'CRVAL1', min(wave), ' central wavelength of first pixel [Angstrom]'
       sxaddpar, hdr, 'CRPIX1', 1, ' starting pixel (1-indexed)'
       sxaddpar, hdr, 'CD1_1', dwave, ' pixel size [Angstrom/pixel]'
       sxaddpar, hdr, 'VDISP', inst_vdisp, ' instrumental velocity dispersion [km/s]'
       sxaddpar, hdr, 'DATE', hogg_iso_date(), ' file creation date'
       
; write out
       splog, 'Writing '+outfile[jj]
       mwrfits, float(flux), outfile[jj], hdr, /create
       mwrfits, info, outfile[jj]
       spawn, 'gzip -f '+outfile[jj], /sh
    endfor
       
return
end
    
