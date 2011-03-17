pro nedscript_distance_catalog
; jm05feb23uofa

    parse_ned_literature, '04karachentsev.ned', outfile='04karachentsev_ned.fits'
    parse_ned_byname, '03karachentsev.ned', outfile='03karachentsev_ned.fits'
    parse_ned_byname, '01tonry.ned', outfile='01tonry_ned.fits'
    parse_ned_byname, '03whiting.ned', outfile='03whiting_ned.fits'
    parse_ned_byname, '01freedman.ned', outfile='01freedman_ned.fits'
    parse_ned_byname, '00ferrarese.ned', outfile='00ferrarese_ned.fits'
    parse_ned_byname, '01shapley.ned', outfile='01shapley_ned.fits'
    
return
end
    
