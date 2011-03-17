pro ediscs_parse_ircatalogs, clobber=clobber
; jm08may22nyu - parse Rose's latest IR catalogs
; jm10may04ucsd - major update    

    irpath = ediscs_path(/catalogs)+'ircatalogs/'
    mycatpath = ediscs_path(/mycatalogs)

    out_template = {$
      galaxy:       '',$
      cluster:      '',$
      photmembflag:  0,$
      specmembflag:  0,$
      matchflag24:   0,$
      log_lir:     0.0,$
      f24:         0.0,$
      f24_err:     0.0,$
      zspec:       0.0}
    
    allcat = file_search(irpath+'*.dat',count=ncat)
    for icat = 0, ncat-1 do begin
       splog, 'Parsing '+allcat[icat]
       readcol, allcat[icat], galaxy, photmembflag, specmembflag, $
         matchflag24, log_lir, f24, f24_err, zspec, comment='#', $
         format='A,I,I,I,F,F,F,F', /silent
       nobj = n_elements(galaxy)
       out1 = replicate(out_template,nobj)
       out1.cluster = strmid(file_basename(allcat[icat]),0,$
         strpos(file_basename(allcat[icat]),'.'))
       out1.galaxy = galaxy
       out1.photmembflag = photmembflag
       out1.specmembflag = specmembflag
       out1.matchflag24 = matchflag24
       out1.log_lir = log_lir
       out1.f24 = f24
       out1.f24_err = f24_err
       out1.zspec = zspec
       if (icat eq 0) then out = out1 else out = [out,out1]
    endfor

; build a line-matched spectroscopic IR catalog
    spec1d = read_ediscs(/spec1d)
    out_spec = im_empty_structure(out,empty_value=-999,$
      ncopies=n_elements(spec1d))
    match, strtrim(out.galaxy,2), strtrim(spec1d.galaxy,2), m1, m2
    out_spec[m2] = out[m1]

;   niceprint, spec1d.galaxy, out_spec.galaxy, spec1d.z, out_spec.zspec
    ww = where(spec1d.z-out_spec.zspec gt 0.1 and out_spec.zspec ne -999)              
;   niceprint, spec1d[ww].cluster, spec1d[ww].galaxy, out_spec[ww].galaxy, spec1d[ww].z, out_spec[ww].zspec
    struct_print, out_spec[ww]

; write out both catalogs
    im_mwrfits, out, mycatpath+'ediscs_ircatalogs.fits', clobber=clobber
    im_mwrfits, out_spec, mycatpath+'ediscs_ircatalogs_spec.fits', clobber=clobber
    
stop    
    
return
end
    
