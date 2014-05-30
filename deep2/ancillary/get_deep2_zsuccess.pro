function get_zsuccess, mag, color1, color2, zmap_grid=zmap_grid, $
  zmap_err_grid=zmap_err_grid, attempt_grid=attempt_grid, got_grid=got_grid, $
  attempt1=attempt1, got1=got1, magaxis=magaxis, color1axis=color1axis, $
  color2axis=color2axis, zs_err=zs_err, justmag=justmag, debug=debug, $
  noextrap=noextrap
; compute the success rate

    minsnr = 1.0

    if keyword_set(justmag) then begin
       mindx = primus_findex(magaxis,mag)

       zs = interpolate(zmap_grid,mindx)
       zs_err = interpolate(zmap_err_grid,mindx)
       attempt1 = interpolate(attempt_grid,mindx)
       got1 = interpolate(got_grid,mindx)

;      zs = interpolate(zmap_grid,mindx,missing=-1)
;      zs_err = interpolate(zmap_err_grid,mindx,missing=-1)
;      attempt1 = interpolate(attempt_grid,mindx,missing=-1)
;      got1 = interpolate(got_grid,mindx,missing=-1)

       if (keyword_set(noextrap) eq 0) then begin
          doit = where((zs le 0.0) or (zs/zs_err lt minsnr),ndoit,comp=good)
          for ii = 0L, ndoit-1 do begin
             mindist = min(mag[doit[ii]]-mag[good],dindx)
             zs[doit[ii]] = zs[good[dindx]]
             zs_err[doit[ii]] = zs_err[good[dindx]]
          endfor
       endif
    endif else begin
       mindx = primus_findex(magaxis,mag)
       cindx1 = primus_findex(color1axis,color1)

; interpolate in 3D       
       if n_elements(color2) ne 0L then begin
          cindx2 = primus_findex(color2axis,color2)
          zs = interpolate(zmap_grid,mindx,cindx1,cindx2)
          zs_err = interpolate(zmap_err_grid,mindx,cindx1,cindx2)
          attempt1 = interpolate(attempt_grid,mindx,cindx1,cindx2)
          got1 = interpolate(got_grid,mindx,cindx1,cindx2)
       endif else begin
          zs = interpolate(zmap_grid,mindx,cindx1)
          zs_err = interpolate(zmap_err_grid,mindx,cindx1)
          attempt1 = interpolate(attempt_grid,mindx,cindx1)
          got1 = interpolate(got_grid,mindx,cindx1)
       endelse

;      zs = interpolate(zmap_grid,mindx,cindx,missing=-1)
;      zs_err = interpolate(zmap_err_grid,mindx,cindx,missing=-1)
;      attempt1 = interpolate(attempt_grid,mindx,cindx,missing=-1)
;      got1 = interpolate(got_grid,mindx,cindx,missing=-1)

       if keyword_set(debug) then djs_plot, mag, color1, ps=6, xsty=3, ysty=3

; optionally deal with poorly measured ZS values due to low S/N or
; because of out-of-bounds colors and/or magnitudes
       if (keyword_set(noextrap) eq 0) then begin
          doit = where((zs le 0.0) or (zs/zs_err lt minsnr),ndoit,comp=good)
          for ii = 0L, ndoit-1 do begin
             if n_elements(color2) ne 0L then begin
                mindist = min(sqrt((mag[doit[ii]]-mag[good])^2 + (color1[doit[ii]]-color1[good])^2 + $
                  (color2[doit[ii]]-color2[good])^2),dindx)
                zs[doit[ii]] = zs[good[dindx]]
                zs_err[doit[ii]] = zs_err[good[dindx]]
                if keyword_set(debug) then splog, 'Debugging plot must be coded!'
             endif else begin
                mindist = min(sqrt((mag[doit[ii]]-mag[good])^2 + (color1[doit[ii]]-color1[good])^2),dindx)
                zs[doit[ii]] = zs[good[dindx]]
                zs_err[doit[ii]] = zs_err[good[dindx]]
                if keyword_set(debug) then djs_oplot, [mag[doit[ii]],mag[good[dindx]]], $
                  [color1[doit[ii]],color1[good[dindx]]], color=djs_icolor('red')
             endelse
          endfor
;         plot, zs, zs/zs_err, ps=6, yr=[0,10]
       endif
    endelse
    
return, zs
end

function get_deep2_zsuccess, field, mag=mag, color1=color1, color2=color2, $
  flag=flag, attempt=attempt, got=got, notarg_weight=notarg_weight, $
  noextrap=noextrap, nosmooth=nosmooth
; jm10sep14ucsd - compute the statistical weight to correct for
;   redshift failures (i.e., the zconf cut)
; jm11jul19ucsd - major update to reflect smoothed completeness map

; notarg_weight - use the completeness map uncorrected for the targeting
;   weight (only useful for testing!)

;   common com_zsuccess, allcomp

; read the completeness map (see build_deep2_zsuccess)
    catpath = deep2_path(/catalogs)
    if (n_elements(allcomp) eq 0) then begin
       compfile = catpath+'deep2_zsuccess.fits.gz'
;      splog, 'Reading '+compfile
       allcomp = mrdfits(compfile,1)
    endif

    ngal = n_elements(mag)
    if (n_elements(field) ne 1) or (ngal eq 0L) then begin
       splog, '(Scalar) FIELD and MAG inputs required'
       return, -1
    endif

    if (n_elements(color1) eq 0L) then color1 = mag*0.0-999.0 else begin
       if (n_elements(color1) ne ngal) then begin
          splog, 'Dimensions of COLOR1 and MAG must match!'
          return, -1
       endif
    endelse
    if (n_elements(color2) eq 0L) then color2 = mag*0.0-999.0 else begin
       if (n_elements(color2) ne ngal) then begin
          splog, 'Dimensions of COLOR2 and MAG must match!'
          return, -1
       endif
    endelse

; pull out the right pieces of the completeness map and read the
; parameters for this field    
    params = get_deep2_completeness_params(field)
    comp = allcomp[where(strtrim(allcomp.field,2) eq field)]

    magaxis = deep2_grid_axis(params,0,/midbin)
    color1axis = deep2_grid_axis(params,1,/midbin)
    color2axis = deep2_grid_axis(params,2,/midbin)

;; renormalize the magnitudes and colors
;    magnorm = (mag-min(mag))/(max(mag)-min(mag))
;
;    color1norm = color1*0.0-999.0
;    gd = where(color1 gt -900.0,ngd)
;    if (ngd ne 0L) then color1norm[gd] = (color1[gd]-min(color1[gd]))/(max(color1[gd])-min(color1[gd]))
;    
;    color2norm = color2*0.0-999.0
;    gd = where(color2 gt -900.0,ngd)
;    if (ngd ne 0L) then color2norm[gd] = (color2[gd]-min(color2[gd]))/(max(color2[gd])-min(color2[gd]))
    
; if the redshift success is formally negative then it means the
; galaxy is in a weird part of color-magnitude space where we had no
; observed targets; assign the weight of the closest object in
; color-color-mag space (except for the offending objects!)
    zsuccess = mag*0.0-1.0
    zsuccess_err = mag*0.0-1.0
    attempt = abs(mag*0.0)
    got = abs(mag*0.0)
    flag = intarr(ngal)-1

    minsnr = 1.0
    
;; both colors; smoothed maps not supported! something is wrong with
;; this, so skip it for now
;    case0 = where(((zsuccess le 0.0) or (zsuccess_err le 0.0) or $
;      (zsuccess lt minsnr*zsuccess_err)) and (color1 gt -900.0) and $
;      (color2 gt -900.0),ncase0)
;    if (ncase0 ne 0L) then begin
;       flag[case0] = 0
;       if keyword_set(notarg_weight) then begin
;          zmap = comp.zsuccess_color1_color2
;          zmap_err = comp.zsuccess_color1_color2_err
;       endif else begin
;          zmap = comp.zsuccess_weighted_color1_color2
;          zmap_err = comp.zsuccess_weighted_color1_color2_err
;       endelse
;
;       zsuccess[case0] = get_zsuccess(mag[case0],color1[case0],color2[case0],$
;         magaxis=magaxis,color1axis=color1axis,color2axis=color2axis,zmap_grid=zmap,$
;         zmap_err_grid=zmap_err,attempt_grid=comp.attempt_color1_color2,$
;         got_grid=comp.got_color1_color2,attempt1=attempt1,got1=got1,$
;         zs_err=zs_err,debug=debug,noextrap=noextrap)
;       zsuccess_err[case0] = zs_err
;       attempt[case0] = attempt1
;       got[case0] = got1
;    endif

; first color    
    case1 = where(((zsuccess le 0.0) or (zsuccess_err le 0.0) or $
      (zsuccess lt minsnr*zsuccess_err)) and (color1 gt -900.0),ncase1)
    if (ncase1 ne 0L) then begin
       flag[case1] = 1
       if keyword_set(notarg_weight) then begin
          if keyword_set(nosmooth) then begin
             zmap = comp.zsuccess_color1
             zmap_err = comp.zsuccess_color1_err
          endif else begin
             zmap = comp.zsuccess_color1_smooth
             zmap_err = comp.zsuccess_color1_smooth_err
          endelse
       endif else begin
          if keyword_set(nosmooth) then begin
             zmap = comp.zsuccess_weighted_color1
             zmap_err = comp.zsuccess_weighted_color1_err
          endif else begin
             zmap = comp.zsuccess_weighted_color1_smooth
             zmap_err = comp.zsuccess_weighted_color1_smooth_err
          endelse
       endelse
       
       zsuccess[case1] = get_zsuccess(mag[case1],color1[case1],$
         magaxis=magaxis,color1axis=color1axis,zmap_grid=zmap,$
         zmap_err_grid=zmap_err,attempt_grid=comp.attempt_color1,$
         got_grid=comp.got_color1,attempt1=attempt1,got1=got1,$
         zs_err=zs_err,debug=debug,noextrap=noextrap)
       zsuccess_err[case1] = zs_err
       attempt[case1] = attempt1
       got[case1] = got1
    endif

; second color    
    case2 = where(((zsuccess le 0.0) or (zsuccess_err le 0.0) or $
      (zsuccess lt minsnr*zsuccess_err)) and (color2 gt -900.0),ncase2)
    if (ncase2 ne 0L) then begin
       flag[case2] = 2
       if keyword_set(notarg_weight) then begin
          if keyword_set(nosmooth) then begin
             zmap = comp.zsuccess_color2
             zmap_err = comp.zsuccess_color2_err
          endif else begin
             zmap = comp.zsuccess_color2_smooth
             zmap_err = comp.zsuccess_color2_smooth_err
          endelse
       endif else begin
          if keyword_set(nosmooth) then begin
             zmap = comp.zsuccess_weighted_color2
             zmap_err = comp.zsuccess_weighted_color2_err
          endif else begin
             zmap = comp.zsuccess_weighted_color2_smooth
             zmap_err = comp.zsuccess_weighted_color2_smooth_err
          endelse
       endelse

       zsuccess[case2] = get_zsuccess(mag[case2],color2[case2],$
         magaxis=magaxis,color1axis=color2axis,zmap_grid=zmap,$
         zmap_err_grid=zmap_err,attempt_grid=comp.attempt_color2,$
         got_grid=comp.got_color2,attempt1=attempt1,got1=got1,$
         zs_err=zs_err,debug=debug,noextrap=noextrap)
       zsuccess_err[case2] = zs_err
       attempt[case2] = attempt1
       got[case2] = got1
    endif

; magnitude    
    case3 = where(((zsuccess le 0.0) or (zsuccess_err le 0.0) or $
      (zsuccess lt minsnr*zsuccess_err)),ncase3)
    if (ncase3 ne 0L) then begin
       flag[case3] = 4
       if keyword_set(notarg_weight) then begin
          zmap = comp.zsuccess_mag
          zmap_err = comp.zsuccess_mag_err
       endif else begin
          zmap = comp.zsuccess_weighted_mag
          zmap_err = comp.zsuccess_weighted_mag_err
       endelse
       
       zsuccess[case3] = get_zsuccess(mag[case3],magaxis=magaxis,zmap_grid=zmap,$
         zmap_err_grid=zmap_err,attempt_grid=comp.attempt_mag,got_grid=comp.got_mag,$
         attempt1=attempt1,got1=got1,zs_err=zs_err,/justmag,noextrap=noextrap)
       zsuccess_err[case3] = zs_err
       attempt[case3] = attempt1
       got[case3] = got1
    endif

return, zsuccess
end
