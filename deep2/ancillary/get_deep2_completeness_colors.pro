function get_deep2_completeness_colors, maggies, targ_mag=targ_mag, $
  params=params, filterlist=filterlist
; jm10sep13ucsd - pull out the select magnitude and the observed-frame
; color(s) of interest, for the completeness corrections

    sz = size(maggies,/dim)
    ngal = sz[1]
    if (n_elements(targ_mag) ne ngal) then $
      message, 'Targeting magnitude required!'
    
    colors = replicate({mag: -999.0, color1: -999.0, $
      color2: -999.0, good1: 0, good2: 0},ngal)
    iselect = where(filterlist eq strtrim(params.select_filter,2))
    icolor0 = where(filterlist eq strtrim(params.color_filters[0,0],2)) ; e.g., r
    icolor1 = where(filterlist eq strtrim(params.color_filters[1,0],2)) ; e.g., i
    icolor2 = where(filterlist eq strtrim(params.color_filters[0,1],2)) ; e.g., g
    icolor3 = where(filterlist eq strtrim(params.color_filters[1,1],2)) ; e.g., i
    
; selection magnitude
    good = where((targ_mag gt 0.0) and (targ_mag lt 90.0),ngood)
    if (ngood ne 0L) then colors[good].mag = targ_mag[good]
;   good = where(photo.maggies[iselect] gt 0.0,ngood)
;   if (ngood ne 0L) then colors[good].mag = $
;     reform(-2.5*alog10(photo[good].maggies[iselect]))

; first color (e.g., r-i)
    m0 = reform(maggies[icolor0,*])
    m1 = reform(maggies[icolor1,*])
    good = where((m0 gt 0.0) and (m1 gt 0.0),ngood)
    if (ngood ne 0L) then colors[good].color1 = $
      reform(-2.5*alog10(m0[good]/m1[good]))

; second color (e.g., g-i)
    m2 = reform(maggies[icolor2,*])
    m3 = reform(maggies[icolor3,*])
    good = where((m2 gt 0.0) and (m3 gt 0.0),ngood)
    if (ngood ne 0L) then colors[good].color2 = $
      reform(-2.5*alog10(m2[good]/m3[good]))

; mask bit    
    colors.good1 = (colors.mag gt -900) and (colors.color1 gt -900)
    colors.good2 = (colors.mag gt -900) and (colors.color2 gt -900)

return, colors
end

