pro build_deep2_photo_catalog
; jm13jul14siena - match the Matthews+13 extended photometric catalog
;   to my revised redshift catalog of objects with good spectra
; jm14may15siena - add unWISE photometry

    catpath = deep2_path(/catalogs)

; read the unwise catalog    
    unwise = mrdfits(catpath+'deep2.dr4.unwise.fits.gz',1)

; read the ACS-GC catalog from Griffith+13
    acs = mrdfits(getenv('IM_DATA_DIR')+'/acs-gc/egs_v_i_public_catalog_V1.0.fits.gz',1)

; good objects
    phot = mrdfits(catpath+'zcat_ext.uniq.fits.gz',1)
    zcat = read_deep2_zcat() ; Q>=3
    match, zcat.objno, phot.objno, m1, m2
;   spherematch, zcat.ra, zcat.dec, $
;     phot.ra_deep, phot.dec_deep, 3D/3600.0, m1, m2
    srt = sort(m1)
    m1 = m1[srt]
    m2 = m2[srt]
    if n_elements(m1) ne n_elements(zcat) then message, 'Problem here!'
    zcat = zcat[m1]
    phot = phot[m2]
    ngal = n_elements(phot)
    
; a small number of objects do not have matching photometry in the
; Matthews catalog but do in the original ZCAT catalog, so fix those
; here
    deep2_to_maggies, phot, mm, ii
    fix = where(total(mm gt 0,1) eq 0)
    
    phot[fix].bestb = zcat[fix].magb
    phot[fix].bestr = zcat[fix].magr
    phot[fix].besti = zcat[fix].magi

    phot[fix].bestberr = zcat[fix].magberr
    phot[fix].bestrerr = zcat[fix].magrerr
    phot[fix].bestierr = zcat[fix].magierr

; add unwise
    match, phot.objno, unwise.objno, m1, m2
    phot = struct_addtags(temporary(phot),im_empty_structure($
      struct_trimtags(unwise[0],select=['w1_*','w2_*']),$
      ncopies=ngal))
    phot[m1] = im_struct_assign(unwise[m2],phot[m1],/nozero)

; add the ACS-GC catalog; check that the missing objects are just
; outside the footprint
;   egs = where(strmid(strtrim(phot.objno,2),0,1) eq '1',negs)
;   match, phot[egs].objno, acs.survey_id, m1, m2
;   djs_plot, phot[egs].ra_deep, phot[egs].dec_deep, psym=3, xstsy=3, ysty=3
;   djs_oplot, acs.ra, acs.dec, psym=3, color='orange'
;   djs_oplot, phot[egs[m1]].ra_deep, phot[egs[m1]].dec_deep, psym=3, color='green'
    
    match, phot.objno, acs.survey_id, m1, m2
    acs1 = struct_trimtags(acs,select=['*radius*','*galfit*',$
      'vis_morph'],except=['x_*','y_*'])

    phot = struct_addtags(temporary(phot),im_empty_structure($
      struct_trimtags(acs1[0]),ncopies=ngal,empty_value=-999))
    phot[m1] = im_struct_assign(acs1[m2],phot[m1],/nozero)

    im_mwrfits, phot, catpath+'photo.dr4.goodspec1d.Q34.fits', /clobber

; ---------------------------------------------------------------------------
; full catalog    
    phot = mrdfits(catpath+'zcat_ext.uniq.fits.gz',1)
    zcat = read_deep2_zcat(/all)
    match, zcat.objno, phot.objno, m1, m2
;   spherematch, zcat.ra, zcat.dec, $
;     phot.ra_deep, phot.dec_deep, 3D/3600.0, m1, m2
    srt = sort(m1)
    m1 = m1[srt]
    m2 = m2[srt]
    if n_elements(m1) ne n_elements(zcat) then message, 'Problem here!'
    zcat = zcat[m1]
    phot = phot[m2]
    ngal = n_elements(phot)

; a small number of objects do not have matching photometry in the
; Matthews catalog but do in the original ZCAT catalog, so fix those
; here
    deep2_to_maggies, phot, mm, ii
    fix = where(total(mm gt 0,1) eq 0)
    phot[fix].bestb = zcat[fix].magb
    phot[fix].bestr = zcat[fix].magr
    phot[fix].besti = zcat[fix].magi

    phot[fix].bestberr = zcat[fix].magberr
    phot[fix].bestrerr = zcat[fix].magrerr
    phot[fix].bestierr = zcat[fix].magierr

; add unwise
    match, phot.objno, unwise.objno, m1, m2
    phot = struct_addtags(temporary(phot),im_empty_structure($
      struct_trimtags(unwise[0],select=['w1_*','w2_*']),$
      ncopies=ngal))
    phot[m1] = im_struct_assign(unwise[m2],phot[m1],/nozero)

; add the ACS-GC catalog
    match, phot.objno, acs.survey_id, m1, m2
    acs1 = struct_trimtags(acs,select=['*radius*','*galfit*',$
      'vis_morph'],except=['x_*','y_*'])

    phot = struct_addtags(temporary(phot),im_empty_structure($
      struct_trimtags(acs1[0]),ncopies=ngal,empty_value=-999))
    phot[m1] = im_struct_assign(acs1[m2],phot[m1],/nozero)

    im_mwrfits, phot, catpath+'photo.dr4.goodspec1d.fits', /clobber

return
end
    
