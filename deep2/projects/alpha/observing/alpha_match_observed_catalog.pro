function alpha_match_observed_catalog, mike, setup=setup, $
  side=side, mikeobj=mikeobj, verbose=verbose
; jm08apr16nyu - given a MIKE structure (e.g., for a specific night of
;                data), match against the catalog of observed objects
;                (see ALPHA_PARSE_OBSERVED_CATALOG)

    common alpha_observed_catalog, observed
    
    if (n_elements(mike) eq 0L) then return, -1L
    if (n_elements(setup) eq 0L) then setup = 1L
    if (n_elements(side) eq 0L) then side = 2L ; red

    latest = 'sep08'
    path = deep2_path(/projects)+'alpha/'

    if (n_elements(observed) eq 0L) then observed = mrdfits(path+$
      'alpha_observed_catalog_'+latest+'.fits',1)

    gd = where((mike.setup eq setup) AND (mike.type eq 'OBJ') and $
      (mike.side eq side),ngd)
    objindx = gd[uniq(mike[gd].obj_id,sort(mike[gd].obj_id))]
    obj = mike[objindx].obj
    objid = mike[objindx].obj_id
    nobj = n_elements(objid)

    mikeobj = mike[objindx]
    
;   spherematch, mike[objindx].ra, mike[objindx].dec, $
;     observed.ra, observed.dec, 5.0/3600.0, m1, m2, $
;     maxmatch=0

    baseobj = repstr(repstr(repstr(repstr(repstr(repstr(obj,$
      '_a',''),'_b'),'_c',''),'_d',''),'_e',''),'obj_','')

    match = lonarr(nobj)
    for ii = 0L, nobj-1L do match[ii] = where($
      string(baseobj[ii],format='(I3.3)') eq $
      string(observed.parent_id,format='(I3.3)'))
    nomatch = where((match eq -1),nnomatch)
    if (nnomatch ne 0L) then message, 'Matching problem'

    info = struct_addtags(observed[match],struct_trimtags(mike[objindx],$
      select=['obj_id','obj','exp','am','img_root']))
    if keyword_set(verbose) then struct_print, info

return, info
end
    
