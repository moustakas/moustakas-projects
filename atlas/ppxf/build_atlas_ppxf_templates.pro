pro build_atlas_ppxf_templates, debug=debug
; jm09dec17ucsd - write out the basic set of BC03 templates for the
;   PPXF template fitting

    outpath = getenv('RESEARCHPATH')+'/data/atlas/ppxf/'
    version = atlas_version(/ppxf_templates)
    outfile = outpath+'atlas_ppxf_bc03_templates_'+version+'.fits'

    light = 299792.458D
    lsun = 3.826D33         ; [erg/s]
    dist = 10.0*3.085678D18 ; 10 pc [cm]
    
    agegrid = [5.0,25.0,50.0,125.0,300.0,650.0,$
      1500.0,3000.0,7000.0,13000.0]
    bc03 = im_read_bc03(age=agegrid*1E6/1E9,$ ; chabrier IMF
      bc03_extras=extras,/silent)
    nmodel = n_elements(agegrid)

; resample to be constant in linear wavelength    
    minwave = 3500.0D
    maxwave = 7000.0D
    dwave = 1.0D
    wave = im_array(minwave,maxwave,dwave,/double)
    npix = n_elements(wave)

    flux = fltarr(npix,nmodel)
    for ii = 0, nmodel-1 do begin
       linterp, bc03.wave, bc03.flux[*,ii], wave, newflux
       flux[*,ii] = lsun*newflux/extras[ii].m_/(4.0*!dpi*dist^2)
       if keyword_set(debug) then begin
          djs_plot, wave, newflux, /ylog, xsty=3, ysty=3
          djs_oplot, bc03.wave, bc03.flux[*,ii], color='orange'
          djs_oplot, wave, newflux1, color='blue'
          cc = get_kbrd(1)
       endif
    endfor

    info = {age: fltarr(nmodel)}
    info.age = bc03.age/1E9

; instrumental velocity dispersion of the models    
    inst_vdisp = 60.0       ; [km/s]
    
; write out; generate an appropriate FITS header
    mkhdr, hdr, float(flux)
    sxdelpar, hdr, 'COMMENT'
    sxaddpar, hdr, 'CRVAL1', min(wave), ' central wavelength of first pixel [Angstrom]'
    sxaddpar, hdr, 'CRPIX1', 1, ' starting pixel (1-indexed)'
    sxaddpar, hdr, 'CD1_1', dwave, ' pixel size [Angstrom/pixel]'
    sxaddpar, hdr, 'VDISP', inst_vdisp, ' instrumental velocity dispersion [km/s]'
    sxaddpar, hdr, 'DATE', hogg_iso_date(), ' file creation date'

; write out
    splog, 'Writing '+outfile
    mwrfits, float(flux), outfile, hdr, /create
    mwrfits, info, outfile
    spawn, 'gzip -f '+outfile, /sh

return
end
    
