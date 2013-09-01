function mzr_filterlist, sdss=sdss
; jm13aug28siena

    if keyword_set(sdss) then begin ; needed by mz_isedfit
       filt = [galex_filterlist(),sdss_filterlist(),(wise_filterlist())[0:1]]
    endif else begin
       filt = ages_filterlist()
       filt = filt[where(strmatch(filt,'*ch3*') eq 0 and strmatch(filt,'*ch4*') eq 0)]
    endelse

return, filt
end
