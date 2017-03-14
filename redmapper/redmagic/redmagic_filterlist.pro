function redmagic_filterlist, short=short
; jm17feb19siena
    filterlist = (sdss_filterlist())[1:4]
    short = ['G','R','I','Z']
return, filterlist
end
