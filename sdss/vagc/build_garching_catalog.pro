;+
; NAME:
;   build_garching_catalog
; PURPOSE:
;   puts together the garching catalog; among duplicates, it picks the
;     observation with the highest median S/N
;   
;-
;------------------------------------------------------------------------------
pro build_garching_catalog

    garching_dr7_in= hogg_mrdfits(sdss_path(/mpa_dr7)+'gal_info_dr7_v5_2.fit.gz',1,nrow=28800)

    garching_dr71= create_struct(garching_dr7_in[0], $
      'garching_dr7_tag', 0L, $
      'garching_dr7_tag_primary', 0L, $
      'garching_dr7_multgroup', 0L, $
      'garching_dr7_nextgroup', 0L, $
      'garching_dr7_firstgroup', 0L, $
      'garching_dr7_ingroup', 0L)
    garching_dr7=replicate(garching_dr71,n_elements(garching_dr7_in))
    struct_assign,garching_dr7_in, garching_dr7
    garching_dr7_in=0
    garching_dr7.garching_dr7_tag=lindgen(n_elements(garching_dr7))
    garching_dr7.garching_dr7_tag_primary=lindgen(n_elements(garching_dr7))
    garching_dr7.garching_dr7_multgroup=1
    garching_dr7.garching_dr7_nextgroup=-1
    garching_dr7.garching_dr7_firstgroup=lindgen(n_elements(garching_dr7))
    garching_dr7.garching_dr7_ingroup=lindgen(n_elements(garching_dr7))

    ingroup=spheregroup(garching_dr7.ra,garching_dr7.dec, $
      2.D/3600.D,chunksize=0.2,multgroup=multgroup, $
      firstgroup=firstgroup,nextgroup=nextgroup)

    primary_indx=lonarr(n_elements(garching_dr7))-1L
;;  for i=0L, n_elements(multgroup)-1L do begin
;;      j=firstgroup[i]
;;      tmp_primary_indx=j
;;      j=firstgroup[i]
;;      for k=0L, multgroup[i]-1L do begin
;;          primary_indx[j]=tmp_primary_indx
;;          j=nextgroup[j]
;;      endfor
;;  endfor

; jm08mar08nyu
    ugroup = ingroup[uniq(ingroup,sort(ingroup))] ; unique groups
    for ii = 0L, n_elements(ugroup)-1L do begin
       if ((ii mod 1000) eq 0L) then splog,' group '+string(ii)
       match = where((ingroup eq ugroup[ii]),nmatch)
       maxsnr = max(garching_dr7[match].sn_median,maxindx)
;      if (nmatch gt 1L) then stop
       primary_indx[match] = match[maxindx[0]]
    endfor

    garching_dr7_tag=lindgen(n_elements(garching_dr7))
    garching_dr7_tag_primary=garching_dr7_tag[primary_indx]

    garching_dr7.garching_dr7_tag= garching_dr7_tag
    garching_dr7.garching_dr7_tag_primary= garching_dr7_tag_primary
    garching_dr7.garching_dr7_multgroup= multgroup
    garching_dr7.garching_dr7_nextgroup= nextgroup
    garching_dr7.garching_dr7_firstgroup= firstgroup
    garching_dr7.garching_dr7_ingroup= ingroup

    spawn,'mkdir -p '+getenv('VAGC_REDUX')+'/garching_dr7'
    outfile=getenv('VAGC_REDUX')+'/garching_dr7/garching_dr7_catalog.fits'
    hdr=['']
    sxaddpar,hdr,'VAGCTIME',systime(), $
      'Time of creation of garching catalog file'
    sxaddpar,hdr,'VAGCVERS',vagc_version(), $
      'Version of vagc used'
    mwrfits,garching_dr7,outfile,hdr, /create

;   null=(get_catalog('garching_dr7/garching_dr7',range=[0,0],/unsigned))[0]
;   struct_assign,{junk:0},null
;   nullfile=getenv('VAGC_REDUX')+'/garching_dr7/null_garching_dr7_catalog.fits'
;   null.garching_dr7_tag=-1
;   null.garching_dr7_tag_primary=-1
;   null.targettype= ' '
;   null.spectrotype= ' '
;   null.subclass= ' '
;   null.release= ' '
;   mwrfits,null,nullfile,/create

return
end

