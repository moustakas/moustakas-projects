function index_11mpc, data, verbose=verbose, searchrad=searchrad
; jm03feb11uofa
; given a data structure with ra and dec coordinates, return the
; indices of galaxies in the 11 Mpc sample

    if n_elements(data) eq 0L then begin
       print, 'Syntax - index = index_11mpc(data)'
       return, -1L
    endif

    if (n_elements(searchrad) eq 0L) then searchrad = 30.0 ; [arcsec]

;   mpc11 = rsex(atlas_path(/projects)+'11mpc/11mpc_sample.2003')
    mpc11 = rsex(atlas_path(/projects)+'11mpc/11mpc_sample.dat')
    ngalaxy = n_elements(mpc11)
    
    raref = 15.0*im_hms2dec(data.ra)
    deref = im_hms2dec(data.dec)
    
    ra = 15.0*im_hms2dec(mpc11.ra)
    de = im_hms2dec(mpc11.dec)

; match by position

    ntot = im_djs_angle_match(raref,deref,ra,de,dtheta=searchrad/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
    index = where(mindx ne -1,ngood)
    if (ngood eq 0L) then return, -1L
    s = lindgen(ngood)
;   s = sort(mdist[index])

    if keyword_set(verbose) then niceprint, data[index[s]].galaxy, $
      mpc11[mindx[index[s]]].galaxy, string(mdist[index[s]]*3600.0,format='(F12.3)')

return, index
end
