function read_ediscs_ppxf_templates, out_tempwave, velscale=velscale, $
  inst_vdisp=inst_vdisp, ntemp=ntemp, nmetal=nmetal, npix=ntemppix, $
  linear_tempwave=linear_tempwave, linear_tempflux=linear_tempflux, $
  tempinfo=out_tempinfo, solar=solar
; jm10apr28ucsd - read the output from BUILD_EDISCS_PPXF_TEMPLATES and
;   optionally resample to a constant ln-lambda; also optionally
;   convolve to the fiducial EDISCS instrumental velocity dispersion;
;   note that INST_VDISP *requires* VELSCALE; also, if VELSCALE is
;   passed then OUT_TEMPWAVE is in natural log

;   common ppxf_templates, tempflux, tempwave, temphdr, tempinfo

    if (n_elements(tempflux) eq 0) then begin
       if keyword_set(solar) then metal = '0.020' else $
         metal = ['0.004','0.020','0.050']
       templatefile = ediscs_path(/ppxf)+'bc03_'+metal+$
         '_'+ediscs_version(/ppxf_templates)+'.fits.gz'
       for ii = 0, n_elements(templatefile)-1 do begin
          splog, 'Reading '+templatefile[ii]
          tempflux1 = mrdfits(templatefile[ii],0,temphdr,/silent)
          tempinfo1 = mrdfits(templatefile[ii],1,/silent)
          tempwave = float(make_wave(temphdr))
          if (ii eq 0) then begin
             tempflux = tempflux1
             tempinfo = tempinfo1
          endif else begin
             tempflux = [[[tempflux]],[[tempflux1]]]
             tempinfo = [tempinfo,tempinfo1]
          endelse
       endfor
    endif

    ndim = size(tempflux,/n_dim)
    tempsize = size(tempflux,/dim)
    ntemppix = tempsize[0]
    ntemp = tempsize[1]
    if keyword_set(solar) then nmetal = 1 else nmetal = tempsize[2]

; these variables are to get around the common-block restrictions    
    linear_tempflux = tempflux
    linear_tempwave = tempwave
    out_tempflux = tempflux
    out_tempwave = tempwave
    out_tempinfo = tempinfo
    
; resample logarithmically (natural log) for use with PPXF
    if (n_elements(velscale) ne 0) then begin
       log_rebin, minmax(linear_tempwave), tempflux[*,0,0], $
         lnflux1, out_tempwave, velscale=velscale
       ntemppix = n_elements(lnflux1)
       out_tempflux = fltarr(ntemppix,ntemp,nmetal)
       for jj = 0, nmetal-1 do begin
          for ii = 0, ntemp-1 do begin
             log_rebin, minmax(linear_tempwave), tempflux[*,ii,jj], $
               lnflux1, velscale=velscale
             out_tempflux[*,ii,jj] = lnflux1
          endfor
       endfor
    endif

; convolve to the EDISCS instrumental resolution
    if (n_elements(inst_vdisp) ne 0) and (n_elements(velscale) ne 0) then begin
       bc03_inst_vdisp = sxpar(temphdr,'VDISP')
       vdisp = sqrt(inst_vdisp^2.0-bc03_inst_vdisp^2)>0.0
       if (vdisp gt 1.0) then begin
          smoothing = vdisp/velscale ; [pixel]
;         kernel = psf_gaussian(npix=10.0*smoothing,$
;           fwhm=2.35*smoothing,/norm,ndim=1)
          nkpix = long(4.0*ceil(smoothing))*2L+3
          klam = findgen(nkpix)-float(nkpix-1.0)/2.0
          kernel = exp(-0.5*(klam/smoothing)^2)/sqrt(2.0*!dpi)/smoothing
          kernel = kernel/total(kernel)
          for jj = 0, nmetal-1 do for ii = 0, ntemp-1 do $
            out_tempflux[*,ii,jj] = convol(out_tempflux[*,ii,jj],$
            kernel,/edge_truncate)
       endif
    endif

return, out_tempflux
end

    
