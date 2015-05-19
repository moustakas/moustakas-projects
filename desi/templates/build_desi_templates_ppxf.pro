pro build_desi_templates_ppxf, debug=debug
; jm10feb25ucsd - write out the basic set of FSPS templates for the
;   PPXF template fitting

    outpath = getenv('DESI_ROOT')+'/ppxf/'
    version = 'v1.0' ; atlas_version(/ppxf_templates)

    agegrid = [5.0,25.0,50.0,125.0,300.0,650.0,$
      1500.0,3000.0,7000.0,13000.0]
    metal = 'Z'+['0.0031','0.0190','0.0300'] ; metallicity grid
    nmodel = n_elements(agegrid)
    nmetal = n_elements(metal)

    fsps_Z003 = im_read_fsps(/miles,/chabrier,/flam,metallicity=metal[0])
    fsps_Z019 = im_read_fsps(/miles,/chabrier,/flam,metallicity=metal[1])
    fsps_Z030 = im_read_fsps(/miles,/chabrier,/flam,metallicity=metal[2])
    fsps = [fsps_Z003,fsps_Z019,fsps_Z030]

; resample to be constant in linear wavelength    
    minwave = 912.0D
    maxwave = 2D4 ; 3D4
    dwave = 0.5D
    npix = (maxwave-minwave)/dwave+1
    wave = dindgen(npix)*dwave+minwave

    outfile = outpath+'fsps_'+metal+'_'+version+'.fits'
    for jj = 0, nmetal-1 do begin
       ageindx = findex(fsps[jj].age,agegrid*1D6)
       flux = interpolate(fsps[jj].flux,findex(fsps[jj].wave,wave),$
         ageindx,/grid)/rebin(reform(interpolate(fsps[jj].mstar,ageindx),$
         1,nmodel),npix,nmodel)
       for ii = 0, nmodel-1 do begin
          if keyword_set(debug) then begin
             djs_plot, wave, flux[*,ii], xr=[3000,8000];, /xlog, /ylog
             djs_oplot, fsps[jj].wave, interpolate(fsps[jj].flux,ageindx[ii])/$
               interpolate(fsps[jj].mstar,ageindx[ii]), color='orange'
             cc = get_kbrd(1)
          endif
       endfor 
       
       info = {metallicity: metal[jj], age: fltarr(nmodel)}
       info.age = agegrid/1D3 ; [Gyr]
       
; instrumental velocity dispersion of the models    
       inst_vdisp = 70.0        ; [km/s]
    
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
    
