function ages_vmax, imag, sample_zmin, sample_zmax, $
  ibright, ifaint, coeffs, area, actual_z, $
  q0=q0, q1=q1, qz0=qz0
; jm10feb01ucsd - AGES-specific wrapper on IM_CALC_VMAX() for
;   computing Vmax

;   if (n_params() ne 5) then begin
;      doc_library, 'ages_vmax'
;      return, -1
;   endif
    
;   filterlist = strtrim(ages_filterlist(),2)
;   ifilt = (where(strmatch(filterlist,'*ndwfs_I*')))[0]
;   v2ab = k_vega2ab(filterlist=filterlist[ifilt],/kurucz,/silent)
;   zpoffset = (ages_zpoffset())[ifilt]
;      
;   weff = k_lambda_eff(filterlist=filterlist[ifilt])
;   kl = k_lambda(weff,/odonnell,/silent)
;   glactc, anc.i_alpha_j2000, anc.i_delta_j2000, $
;     2000.0, gl, gb, 1, /deg
;   ebv = dust_getval(gl,gb,/interp,/noloop)
;
;   imag = -2.5*alog10(anc.maggies[ifilt]) - zpoffset + kl*ebv - v2ab ; Vega

; test code    
    aa = read_ages(/ancillary)
    vv = mrdfits(ages_path(/catalogs)+'catalog.Vmax.v3.fits.gz',1,rows=aa.ages_id)

    cut = where((aa.i_obs gt 15.0) and (aa.i_obs lt 20.0) and $
      (aa.z gt 0.05) and (aa.z lt 0.75) and (aa.mass gt 0.0) and $
      (vv.kcorr_flag eq 1))
    aa = aa[cut]
    vv = vv[cut]
    
    rr = im_calc_vmax(aa.i_obs,aa.coeffs,'ndwfs_I.par',$
      1.0,15.0,20.0,0.05,0.75,actual_z=aa.z,$
      q0=q0,q1=q1,qz0=qz0);,/silent)

return, result
end
    
