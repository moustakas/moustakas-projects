function ediscs_sfh_measure_indices, lnwave1, lnflux1, $
  linear=linear, velscale=velscale, vdisp=vdisp, $
  inst_vdisp=inst_vdisp, final_vdisp=final_vdisp, $
  debug=debug
; jm10may07ucsd - measure spectral indices from an optionally smoothed
; spectrum; LNWAVE and LNFLUX should be in the rest frame and
; ln-binned, unless /LINEAR; KERNEL should be the velocity kernel in
; km/s 

    npix = n_elements(lnwave)
    version = ediscs_version(/ppxf_specfit)
    indexfile = ediscs_path(/ppxf)+'lick_indexlist_'+version+'.dat'

; rebin the spectra to be constant in ln-wavelength
    if keyword_set(linear) then begin
       log_rebin, minmax(lnwave1), lnflux1, lnflux, $
         lnwave, velscale=velscale
    endif else begin
       lnflux = lnflux1
       lnwave = lnwave1
    endelse

    if (n_elements(vdisp) eq 0) then vdisp = 0.0 ; [km/s]
    if (n_elements(inst_vdisp) eq 0) then inst_vdisp = 0.0 ; [km/s]
    if (n_elements(final_vdisp) eq 0) then final_vdisp = 350.0 ; [km/s]

    vkernel = sqrt((final_vdisp^2-vdisp^2-inst_vdisp^2)>0) ; [km/s]

    if (vkernel gt 0.0) then begin
       smoothing = vkernel/velscale ; pixels
       npix = long(4.0*ceil(smoothing))*2L+3
       klam = findgen(npix)-float(npix-1.0)/2.0
       kernel = exp(-0.5*(klam/smoothing)^2)/sqrt(2.0*!dpi)/smoothing
       kernel = kernel/total(kernel)
       lnflux_smooth = convol(lnflux,kernel,/edge_truncate)
    endif else lnflux_smooth = lnflux
    
    indices = spectral_indices(exp(lnwave),lnflux_smooth,$
      indexfile=indexfile,debug=debug,/silent) ; no weighting

return, indices
end
