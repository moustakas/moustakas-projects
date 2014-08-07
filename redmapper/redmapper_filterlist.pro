function redmapper_filterlist
; jm13mar28siena
    filterlist = [sdss_filterlist(),wise_filterlist(/short)]
;   filterlist = [galex_filterlist(),sdss_filterlist(),wise_filterlist(/short)]
return, filterlist
end
