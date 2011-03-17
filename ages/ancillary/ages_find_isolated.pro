pro ages_find_isolated
; jm10feb11ucsd - identify isolated galaxies in the AGES spectroscopic
;   sample to test the GALEX deblending issues

    common ages_isolated, allparent, allbootes, ingroup, mult, first, next

; read the parent bootes catalog
    if (n_elements(allbootes) eq 0L) then begin
       bootespath = getenv('RESEARCHPATH')+'/data/bootes/'
       bootesfile = bootespath+'bootes_I.fits.gz'
       allbootes = mrdfits(bootesfile,1)
    endif

; read the full AGES sample
    if (n_elements(ages) eq 0L) then begin
       photfile = ages_path(/mycatalogs)+'ages_photometry_'+$
         ages_version(/photo)+'.fits.gz'
       splog, 'Reading '+photfile
       allparent = mrdfits(photfile,1,/silent)
    endif

; just consider the objects with spectroscopic redshifts    
    allout = struct_addtags(struct_trimtags(allparent,$
      select=['ages_id','ra','dec','z']),$
      replicate({nfriends: -1},n_elements(allparent)))
    zspec = where(allparent.z gt 0.0,ngal)
    out = allout[zspec]
    parent = allparent[zspec]

; restrict the full BOOTES sample to bright objects; NUV-I~3, so
; I=24.5 corresponds to NUV=27.5, which is *way* lower than the NUV
; detection limit (~25.5 mag)
    bright = where(allbootes.mag_auto lt 24.5,nallgal)
    splog, 'Ngal (I<24.5) = ', nallgal
    bootes = allbootes[bright]

; now group the BOOTES catalog using a radius of 6"
    rad = 6.0D/3600.0D
    if (n_elements(ingroup) eq 0L) then begin
       ingroup = spheregroup(bootes.alpha_j2000,bootes.delta_j2000,$
         rad,multgroup=mult,firstgroup=first,nextgroup=next)
    endif

; for each group, find an object which is specprimary (if none exists,
; use the first object in the group) and set primary_indx to the
; position of the chosen object

; loop through each group
    for ii = 0L, nallgal-1L do begin
       if ((ii mod 1000) eq 0) then splog,' galaxy '+string(ii)
       jj = first[ii]
; loop through each member of the group; if any member matches one of
; the objects in PARENT then store the multiplicity of that group and
; move onto the next group
       for kk = 0L, mult[ii]-1L do begin
          match = where(bootes[jj].object_position eq parent.i_object_position,nmatch)
          if (nmatch gt 1) then stop
          if (match[0] ne -1) then begin
             splog, parent[match].i_alpha_j2000, parent[match].i_delta_j2000, $
               bootes[jj].alpha_j2000, bootes[jj].delta_j2000
;            if (friends[match[0]]) ne -1 then stop
             out[match].nfriends = mult[ii] ; number of neighbors
          endif
          jj = next[jj]
       endfor
    endfor

    allout[zspec] = out
    im_mwrfits, allout, ages_path(/mycatalogs)+'ages_nfriends.fits', /clobber
    
stop    
    
return
end
    
