pro jd_pah_loiii, ho, ho_nodust, sings, nuclear
; jm06mar29uofa - retrieve the nuclear L(OIII) of the SINGS galaxies  

    outpath = sings_path(/projects)+'jd_pah/'
    typepath = sings_path(/projects)+'log12oh/'

    if (n_elements(ho) eq 0L) then ho = read_97ho(datanodust=ho_nodust)

    snrcut = 1.0
    dist_cut = 17.0 ; 17.0

    if (n_elements(ho) eq 0L) then ho = read_97ho(datanodust=ho_nodust)
    if (n_elements(sings) eq 0L) then sings = sings_read_info()
    if (n_elements(nuclear) eq 0L) then nuclear = read_sings(/nuclear)
    ngalaxy = n_elements(sings)

    result = {galaxy: '', type: '', loiii_ho: -1.0, phot_ho: '-', loiii: -1.0, phot: '-'}
    result = replicate(result,ngalaxy)
    result.galaxy = sings.galaxy

    types = rsex(typepath+'sings_types.txt')
;   niceprint, result.galaxy, types.galaxy
    result.type = types.type
    
; Ho et al. 1997
    
    doit = match_string(sings.galaxy,ho.galaxy,/exact,findex=homatch,/silent)
    good = where(homatch ne -1L,nmatch)
    niceprint, sings[good].galaxy, ho[homatch[good]].galaxy

    result[good].loiii_ho = ho[homatch[good]].oiii_5007_lum[0]
    result[good].phot_ho = ho[homatch[good]].photflag

; SINGS spectroscopy
    
    doit = match_string(sings.galaxy,nuclear.galaxy,/exact,findex=nucmatch,/silent)
    good = where(nucmatch ne -1L,nmatch)
    niceprint, sings[good].galaxy, nuclear[nucmatch[good]].galaxy

    result[good].loiii = nuclear[nucmatch[good]].oiii_5007_lum[0]
    result[good].phot = nuclear[nucmatch[good]].nuclear_photflag

    agn = where((strmatch(result.type,'*HII*',/fold) eq 0B),nagn)
    struct_print, result[agn]
    
stop
    
return
end
