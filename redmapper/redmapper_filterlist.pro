function redmapper_filterlist
; jm13mar28siena
    filterlist = [galex_filterlist(),sdss_filterlist(),$
      (wise_filterlist())[0:1]]
return, filterlist
end
