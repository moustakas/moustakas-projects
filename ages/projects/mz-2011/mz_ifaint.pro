function mz_ifaint, vega=vega, ibright=ibright, select_filter=select_filter, $
  select_vega2ab=select_vega2ab
; I-band limits of the MZ sample
    select_filter = 'ndwfs_I.par'
    select_vega2ab = k_vega2ab(filterlist=select_filter,/kurucz,/silent)
    if keyword_set(vega) then begin
       ibright = 15.0
       ifaint = 19.95
    endif else begin
       ibright = 15.0+select_vega2ab
       ifaint = 19.95+select_vega2ab
    endelse
return, ifaint
end
    
