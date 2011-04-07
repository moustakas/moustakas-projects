function bcgs_filterlist, sdss=sdss
; jm11apr06ucsd - 
    if keyword_set(sdss) then begin
       filterlist = [sdss_filterlist(),twomass_filterlist()]
    endif else begin
       filterlist = bootes_filterlist()
       keep = where((strmatch(filterlist,'*ufilter*') eq 0) and $
         (strmatch(filterlist,'*bok*') eq 0))
       filterlist = filterlist[keep]
    endelse
return, filterlist
end
    
    
