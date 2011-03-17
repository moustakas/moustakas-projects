function index_gto_starbursts, data, gto, index_gto=index_gto, $
  verbose=verbose, searchrad=searchrad
; jm05jul25uofa
; given a data structure with ra and dec coordinates, return the
; indices of galaxies in the Gto_Starbursts et al. (1989) subsample

    if n_elements(data) eq 0L then begin
       print, 'Syntax - index = index_gto_starbursts(data)'
       return, -1L
    endif

    if (n_elements(searchrad) eq 0L) then searchrad = 15.0 ; [arcsec]
    
    gto = read_gto_starbursts()

    ra = 15.0*im_hms2dec(gto.ra)
    de = im_hms2dec(gto.dec)

    raref = 15.0*im_hms2dec(data.ra)
    deref = im_hms2dec(data.dec)
    
; match by position

    ntot = im_djs_angle_match(raref,deref,ra,de,dtheta=searchrad/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
    index = where(mindx ne -1,ngood)
    if (ngood eq 0L) then return, -1L
    s = sort(mdist[index])

    index_gto = mindx[index]
;   gto = gto[index_gto]
    
    if keyword_set(verbose) then niceprint, data[index[s]].galaxy, $
      gto[mindx[index[s]]].galaxy, data[index[s]].ned_galaxy, $
      gto[mindx[index[s]]].nedgalaxy, $
      string(mdist[index[s]]*3600.0,format='(F12.3)')
    
return, index
end
