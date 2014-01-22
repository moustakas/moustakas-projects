function read_a2744, photoz=photoz, bpz_dz=bpz_dz, bpz_redshift=bpz_redshift, $
  index=index
; jm14jan10siena - read the photometric catalog

    date = 'jan03'
    
    path = getenv('CLASH_PROJECTS')+'/hff/a2744/'
    cat = rsex(path+'flx_iso.'+date)
    bpz = rsex(path+'bpz/'+date+'.v2.bpz')
    if total(cat.id-bpz.id) ne 0 then message, 'Problem here!'

    pbpz = read_bpz_probs(path+'bpz/'+date+'.probs',$
      redshift=bpz_redshift,dz=bpz_dz)
    
; trim    
    keep = where(cat.id ne 1215 and cat.id ne 3117,ngal)
    cat = cat[keep]
    bpz = bpz[keep]
    pbpz = pbpz[keep]

    cat = struct_addtags(cat,struct_trimtags(bpz,except='id'))
    cat = struct_addtags(cat,replicate({group: -1, galaxy: '', $
      mu: -1.0, pofz: fltarr(n_elements(bpz_redshift))},ngal))
    cat.galaxy = 'ID'+strtrim(cat.id,2)
    cat.pofz = pbpz.pofz

    mag = rsex(path+'model_jan03.txt') ; magnifications
    match, cat.id, mag.id, m1, m2
    cat[m1].mu = mag[m2].nfw

    group1 = where($
      cat.id eq 383 or $
      cat.id eq 2509 or $
      cat.id eq 1880 or $
      cat.id eq 406 or $
      cat.id eq 400 or $
      cat.id eq 407 or $
      cat.id eq 652 or $
      cat.id eq 270 or $
      cat.id eq 575 or $
      cat.id eq 339 or $
      cat.id eq 292 or $
      cat.id eq 527 or $
      cat.id eq 120 or $
      cat.id eq 3903)
    
    group2 = where($
      cat.id eq 534 or $
      cat.id eq 690 or $
      cat.id eq 429 or $
      cat.id eq 523)

    group3 = where($
      cat.id eq 2521 or $
      cat.id eq 286 or $
      cat.id eq 2129 or $
      cat.id eq 2338)

    cat[group1].group = 1
    cat[group2].group = 2
    cat[group3].group = 3

; get rid of the tentative candidates    
    if keyword_set(photoz) eq 0 then begin
       index = where(cat.group ne 3)
       cat = cat[index]
    endif

return, cat
end    
