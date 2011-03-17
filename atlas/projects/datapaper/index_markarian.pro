function index_markarian, data, markarian, verbose=verbose, searchrad=searchrad
; jm05jul25uofa
; given a data structure with ra and dec coordinates, return the
; indices of galaxies in the Markarian et al. (1989) subsample

    if n_elements(data) eq 0L then begin
       print, 'Syntax - index = index_markarian(data)'
       return, -1L
    endif

    if (n_elements(searchrad) eq 0L) then searchrad = 30.0 ; [arcsec]

    markarian = read_89markarian()

    ra = 15.0*im_hms2dec(markarian.ra)
    de = im_hms2dec(markarian.dec)

    raref = 15.0*im_hms2dec(data.ra)
    deref = im_hms2dec(data.dec)
    
; match by position

    ntot = im_djs_angle_match(raref,deref,ra,de,dtheta=searchrad/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
    index = where(mindx ne -1,ngood)
    if (ngood eq 0L) then return, -1L
    s = sort(mdist[index])

    if keyword_set(verbose) then niceprint, data[index[s]].galaxy, $
      data[index[s]].alt_galaxy, markarian[mindx[index[s]]].galaxy, $
      string(mdist[index[s]]*3600.0,format='(F12.3)')

    markarian = markarian[mindx[index]]
    
return, index
end
