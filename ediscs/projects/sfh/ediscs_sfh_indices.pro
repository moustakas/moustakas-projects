pro ediscs_sfh_indices, debug=debug, clobber=clobber
; jm10may07ucsd - remeasure the Lick indices from spectra that have
; been smoothed to a common velocity dispersion

    sfhpath = ediscs_path(/projects)+'sfh/'
    outfile = sfhpath+'ediscs_sfh_smooth_indices.fits'
    if file_test(outfile) and (keyword_set(clobber) eq 0) then begin
       splog, 'Output file '+outfile+' exists; use /CLOBBER'
       return
    endif

    ppxf = read_ediscs(/ppxf)
    ngal = n_elements(ppxf)

; read the indices we want to measure
    version = ediscs_version(/ppxf_specfit)
    indexfile = ediscs_path(/ppxf)+'lick_indexlist_'+version+'.dat'
    
; group by cluster for computational ease    
    allcl = strtrim(ppxf.cluster,2)
    cl = allcl[uniq(allcl,sort(allcl))]
    nc = n_elements(cl)

    final_vdisp = 350.0 ; [km/s]

    t0 = systime(1)
    for ii = 0, nc-1 do begin
       splog, 'Cluster '+cl[ii]
       these = where(cl[ii] eq allcl,ngal)
;      these = these[0:10] & ngal = n_elements(these)
       spec = read_ediscs_gandalf_specfit(ppxf[these])

;      for jj = 4, ngal-1 do begin
       for jj = 0, ngal-1 do begin
          print, "Object ", jj, " of ", ngal, string(13b), $
            format='(A0,I0,A0,I0,A1,$)'

; get the pixel size in km/s and instrumental dispersion for this spectrum
          velscale = ediscs_ppxf_velscale(run34=ppxf[these[jj]].run eq 'run34')
          inst_vdisp = ediscs_ppxf_instvdisp(run34=ppxf[these[jj]].run eq 'run34')

          notzero = where(spec[jj].wave gt 0.0,npix)
          wave = spec[jj].wave[notzero]
          flux = spec[jj].flux[notzero]
          ferr = spec[jj].ferr[notzero]
          flux_nolines = spec[jj].flux[notzero] - spec[jj].linefit[notzero]
          continuum = spec[jj].continuum[notzero]

; observed spectrum
          ediscs_sfh_smooth, wave, flux, sflux, ferr=ferr, $
            smoothed_ferr=sferr, velscale=velscale, $
            vdisp=ppxf[these[jj]].vdisp, inst_vdisp=inst_vdisp, $
            final_vdisp=final_vdisp

          sivar = 1.0/ferr^2 ; note, unsmoothed errors!
          badpixels = where((ferr gt 1E5) or (abs(flux) gt 100.0),nbad)
          if (nbad ne 0) then sivar[badpixels] = 0.0
          
          indices_raw = spectral_indices(exp(wave),sflux,$
            ivar=sivar,indexfile=indexfile,debug=debug,/silent)

; emission-line subtracted spectrum
          ediscs_sfh_smooth, wave, flux_nolines, sflux_nolines, $
            velscale=velscale, vdisp=ppxf[these[jj]].vdisp, $
            inst_vdisp=inst_vdisp, final_vdisp=final_vdisp
          indices_cor = spectral_indices(exp(wave),sflux_nolines,$
            ivar=sivar,indexfile=indexfile,debug=debug,/silent)

; best-fitting model; if the continuum coefficients are all zero then
; just return the empty structure
          if (total(ppxf[these[jj]].continuum_coeff,/double) eq 0) then $
            empty_structure = 1 else empty_structure = 0
          ediscs_sfh_smooth, wave, continuum, scontinuum, $
            velscale=velscale, vdisp=ppxf[these[jj]].vdisp, $
            inst_vdisp=inst_vdisp, final_vdisp=final_vdisp

          modelsivar = sferr*0.0+1.0/median(sferr)^2 ; the model errors are illustrative
          indices_model = spectral_indices(exp(wave),scontinuum,$
            ivar=modelsivar,indexfile=indexfile,debug=debug,/silent,$
            empty_structure=empty_structure)

; pack everything together          
          indices_raw = struct_trimtags(temporary(indices_raw),except='INDICES')
          indices_raw = im_struct_trimtags(indices_raw,select=tag_names(indices_raw),$
            newtags=tag_names(indices_raw)+'_raw')

          indices_cor = struct_trimtags(temporary(indices_cor),except='INDICES')
          indices_cor = im_struct_trimtags(indices_cor,select=tag_names(indices_cor),$
            newtags=tag_names(indices_cor)+'_cor')
          
          indices_model = struct_trimtags(temporary(indices_model),except='INDICES')
          indices_model = im_struct_trimtags(indices_model,select=tag_names(indices_model),$
            newtags=tag_names(indices_model)+'_model')

          indices = struct_addtags(struct_addtags(temporary(indices_raw),$
            temporary(indices_cor)),temporary(indices_model))
;         newindices = [indkeep+'_raw',indkeep+'_cor',indkeep+'_model']
;         indices = struct_addtags(struct_addtags(struct_addtags({indices: newindices},$
;           temporary(indices_raw)),temporary(indices_cor)),temporary(indices_model))

          if (ii eq 0) and (jj eq 0) then out_indices = indices else $
            out_indices = [temporary(out_indices),indices]
       endfor
;      diff = out_indices.lick_hd_a_model[0]-ppxf[these].lick_hd_a_model[0]
;      plot, out_indices.d4000_narrow_model[0], diff, $
;        ps=6, yrange=[-1,1], xsty=3, ysty=3, xrange=[1.0,2.2]
    endfor
    splog, 'Total time = ', (systime(1)-t0)/60

; write out
    out = struct_addtags(struct_trimtags(ppxf,select=['EDISCS_ID']),out_indices)
    im_mwrfits, out, outfile, clobber=clobber

stop    
    
return
end
