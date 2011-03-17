pro write_gto_hiiregions, hii, write=write
; jm06dec09nyu - 

    version = gto_log12oh_version()
    outpath = gto_path()

; read and cross-match the list of GTO galaxies against the
; HII-region database

    bighii = read_hii_regions()
    biggto = gto_read_ancillary()

; ---------------------------------------------------------------------------    
    ra = 15.0*im_hms2dec(biggto.ra)
    dec = im_hms2dec(biggto.dec)
    raref = 15.0*im_hms2dec(bighii.galaxy_ra)
    decref = im_hms2dec(bighii.galaxy_dec)
    
    searchrad = 10.0
    ntot = djs_angle_match(raref,decref,ra,dec,dtheta=searchrad/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
    index = where(mindx ne -1,ngood)
    index_gto = mindx[index]
 
    hii = bighii[index]
    gto = biggto[index_gto]
 
    splog, 'Replacing HII_GALAXY with GTO_GALAXY!'
    hii = struct_addtags(replicate({hii_galaxy_original: ''},ngood),hii)
    hii.hii_galaxy_original = hii.hii_galaxy
    hii.hii_galaxy = gto.galaxy

    s = sort(gto.galaxy)
    niceprint, gto[s].galaxy, hii[s].ned_galaxy, hii[s].hii_galaxy, $
      hii[s].hii_galaxy_original, hii[s].reference, $
      string(mdist[index[s]]*3600.0,format='(F12.3)')

    uindx = uniq(strtrim(hii.ned_galaxy,2),sort(strtrim(hii.ned_galaxy)))
    nmatch = n_elements(uindx)
    splog, 'Identified '+string(nmatch,format='(I0)')+' galaxies with fluxes for '+$
      string(ngood,format='(I0)')+' HII regions.'

; ---------------------------------------------------------------------------    

;   doit = match_string(biggto.ned_galaxy,bighii.ned_galaxy,index=index,/exact,/silent)
;   good = where(index ne -1L,ngood)
;   hii = bighii[index[good]]
;   uindx = uniq(strtrim(hii.ned_galaxy,2),sort(strtrim(hii.ned_galaxy)))
;   nmatch = n_elements(uindx)
;   splog, 'Identified '+string(nmatch,format='(I0)')+' galaxies with fluxes for '+$
;     string(ngood,format='(I0)')+' HII regions.'
;
;   for i = 0L, nmatch-1L do print, hii[uindx[i]].hii_galaxy, $
;     strtrim(hii[uindx[i]].ned_galaxy,2), hii[uindx[i]].texref, format='(A35,3A17)'

; write out    
    
    if keyword_set(write) then begin

;      hii = struct_addtags(im_struct_trimtags(gto,select=['GALAXY'],newtags=['GTO_GALAXY']),hii)
;      hii = struct_addtags(struct_trimtags(gto,select=['GTO_GALAXY']),hii)
       
       hiifile = 'gto_hiiregions_'+version+'.fits'
       splog, 'Writing '+outpath+hiifile+'.'
       mwrfits, hii, outpath+hiifile, /create
       spawn, ['gzip -f '+outpath+hiifile], /sh
       
    endif

return
end
    
