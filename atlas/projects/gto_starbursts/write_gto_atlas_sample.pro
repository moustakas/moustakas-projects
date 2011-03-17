pro write_gto_atlas_sample, atlas1, atlas_nodust1, write=write
; jm06dec09nyu

    version = gto_log12oh_version()
    outpath = gto_path()

    if (n_elements(atlas1) eq 0L) then atlas1 = read_integrated(atlasnodust=atlas_nodust1)
    
; ---------------------------------------------------------------------------    
; cross-match the GTO STARBURSTS and the spectral atlas
; ---------------------------------------------------------------------------    
       
    gto1 = gto_read_ancillary()

;   match, strtrim(gto1.ned_galaxy,2), strtrim(atlas1.ned_galaxy,2), indx_gto, indx_atlas

    ra = 15.0*im_hms2dec(gto1.ra)
    dec = im_hms2dec(gto1.dec)
    raref = 15.0*im_hms2dec(atlas1.ra)
    decref = im_hms2dec(atlas1.dec)
    
    searchrad = 10.0
    ntot = djs_angle_match(raref,decref,ra,dec,dtheta=searchrad/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
    index_atlas = where(mindx ne -1,ngood)
    index_gto = mindx[index_atlas]

    gto = gto1[index_gto]
    atlas = atlas1[index_atlas]
    atlas_nodust = atlas_nodust1[index_atlas]

    s = sort(gto.galaxy)
    niceprint, gto[s].galaxy, gto[s].ned_galaxy, atlas[s].ned_galaxy, $
      string(mdist[index_atlas[s]]*3600.0,format='(F12.3)')
    splog, 'Identified '+string(n_elements(index_gto),format='(I0)')+' GTO galaxies in the MK06 spectral atlas.'
;   niceprint, gto.galaxy, atlas.galaxy, gto.ned_galaxy, atlas.ned_galaxy

; replace the ATLAS NED_GALAXY with the GTO NED_GALAXY; these are
; different because my NED query for the ATLAS was so long ago that
; the NED_GALAXY name for at least one object (MRK0331) changed in the
; interim; these need to match because in GTO_LOG12OH I match on
; NED_GALAXY name

    atlas.ned_galaxy = gto.ned_galaxy
    
; ---------------------------------------------------------------------------    
; classifications
; ---------------------------------------------------------------------------    

;   gto_class, atlas, write=write

; ---------------------------------------------------------------------------    
; write out
; ---------------------------------------------------------------------------    

    if keyword_set(write) then begin
       splog, 'Writing '+outpath+'gto_atlas_speclinefit_'+version+'.fits.gz'
       mwrfits, atlas, outpath+'gto_atlas_speclinefit_'+version+'.fits', /create
       spawn, ['gzip -f '+outpath+'gto_atlas_speclinefit_'+version+'.fits'], /sh

       splog, 'Writing '+outpath+'gto_atlas_speclinefit_nodust_'+version+'.fits.gz'
       mwrfits, atlas_nodust, outpath+'gto_atlas_speclinefit_nodust_'+version+'.fits', /create
       spawn, ['gzip -f '+outpath+'gto_atlas_speclinefit_nodust_'+version+'.fits'], /sh
    endif

return
end
    
