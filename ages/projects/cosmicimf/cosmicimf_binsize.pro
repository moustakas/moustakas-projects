function cosmicimf_binsize, histmin=histmin, histmax=histmax, lf24=lf24
; jm10mar18ucsd - Vmax binsize for the MFs, by default, or for the
; L(24) LFs if /lf24
    if keyword_set(lf24) then begin
       binsize = 0.1
       histmin = 9.0
       histmax = 12.0
    endif else begin
       binsize = 0.1
       histmin = 7.85
       histmax = 12.0
    endelse
return, binsize
end
    
