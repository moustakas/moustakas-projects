function restore_ages_ppxf_bestfit, weights1, ebv=ebv, $
  bestwave=bestwave, velscale=velscale, inst_vdisp=inst_vdisp, $
  bestage=bestage, solar=solar
; jm09nov13ucsd - restore the best-fitting PPXF fit (in the rest
;   frame), ignoring the polynomials; deals with a 2D vector of
;   WEIGHTS and EBV values; note that the LOSVD is not handled here!
;   (and there is some funny-business in PPXF when the best-fit VDISP
;   is very small - see BUILD_AGES_PRIMUS_PARENTSAMPLE)

    ndim = size(weights1,/n_dim)
    if (ndim eq 0) then begin
       doc_library, 'restore_ages_ppxf_bestfit'
       return, -1
    endif

    if (ndim eq 2) then begin
       dims = size(weights1,/dim)
       if (n_elements(ebv) eq 0) then ebv1 = fltarr(dims[1]) else $
         ebv1 = ebv
       bestage = fltarr(dims[1])
       for jj = 0, dims[1]-1 do begin
          modelflux1 = restore_ages_ppxf_bestfit(weights1[*,jj],$
            ebv=ebv1[jj],bestwave=bestwave,velscale=velscale,$
            inst_vdisp=inst_vdisp,bestage=bestage1,solar=solar)
          if (jj eq 0) then modelflux = cmreplicate(modelflux1,dims[1])
          modelflux[*,jj] = modelflux1
          bestage[jj] = bestage1
       endfor
       return, modelflux
    endif
    
; read the templates and redden    
    tempflux = read_ages_ppxf_templates(bestwave,ntemp=ntemp,$
      velscale=velscale,inst_vdisp=inst_vdisp,tempinfo=tempinfo,$
      solar=solar)

    if (n_elements(ebv) ne 0) then begin
       if (n_elements(velscale) eq 0) then $
         kl = k_lambda(bestwave,/calzetti,/silent) else $
           kl = k_lambda(exp(bestwave),/calzetti,/silent)
       for ii = 0, ntemp-1 do tempflux[*,ii] = $
         tempflux[*,ii]*10.0^(-0.4*kl*ebv)
    endif

    weights = weights1
    bestfit = tempflux # weights

; compute the luminosity-weighted age    
    bestage = total(weights*tempinfo.age)/total(weights) ; [Gyr]

; if linear wavelength spacing is retained, then interpolate/rebin
; onto the native AGES pixel size
    if (n_elements(inst_vdisp) eq 0) and $
      (n_elements(velscale) eq 0) then begin
       dwave = ages_pixel_size()
       bestfit1 = bestfit
       bestwave1 = bestwave
       bestfit = im_linear_rebin(bestwave1,bestfit1,$
         dwave=dwave,outwave=bestwave)
    endif
       
return, bestfit
end
    
