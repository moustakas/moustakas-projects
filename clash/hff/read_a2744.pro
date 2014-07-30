function read_a2744, photoz=photoz, bpz_dz=bpz_dz, bpz_redshift=bpz_redshift, $
  index=index
; jm14jan10siena - read the photometric catalog

    date = 'jul24'
;   date = 'feb24'
;   date2 = 'jan03'
    
    path = getenv('CLASH_PROJECTS')+'/hff/a2744/'
    cat = rsex(path+'flx_iso.'+date)
    bpz = rsex(path+'bpz/'+date+'.bpz')
;   bpz = rsex(path+'bpz/'+date+'.v2.bpz')
    if total(cat.id eq bpz.other) ne n_elements(bpz) then message, 'Problem here!'
;   if total(cat.id-bpz.id) ne 0 then message, 'Problem here!'

    pbpz = read_bpz_probs(path+'bpz/'+date+'.probs',$
      redshift=bpz_redshift,dz=bpz_dz)
    ngal = n_elements(bpz)

;; trim (jan03)
;   keep = where(cat.id ne 1215 and cat.id ne 3117,ngal)
;   cat = cat[keep]
;   bpz = bpz[keep]
;   pbpz = pbpz[keep]

    cat = struct_addtags(cat,struct_trimtags(bpz,except='id'))
    cat = struct_addtags(cat,replicate({group: 1, galaxy: '', $
      pofz: fltarr(n_elements(bpz_redshift))},ngal))
    cat.galaxy = strtrim(cat.id,2)
;   cat.galaxy = 'ID'+strtrim(cat.id,2)
    cat.pofz = pbpz.pofz

;    mag = rsex(path+'model_jan03.txt') ; magnifications
;    match, oldid, mag.id, m1, m2
;
;;   match, cat.id, mag.id, m1, m2
;    cat[m1].mu = mag[m2].nfw

;; old code, some of which might still be useful
;    newid = strarr(ngal)
;    oldid = strarr(ngal)
;    for ii = 0, ngal-1 do oldid[ii] = strmid(cat[ii].id,0,strpos(cat[ii].id,'_'))
;    for ii = 0, ngal-1 do newid[ii] = strmid(cat[ii].id,strpos(cat[ii].id,'_')+1)
;    
;; feb17
;    group1 = where($
;      newid eq 'YD1' or $
;      newid eq 'YD2' or $
;      newid eq 'YD4' or $
;      newid eq 'YD6' or $
;      newid eq 'YD7' or $
;      newid eq 'YD8' or $
;      newid eq 'ZD2' or $
;      newid eq 'ZD3' or $
;      newid eq 'ZD4' or $
;      newid eq 'ZD5' or $
;      newid eq 'ZD7' or $
;      newid eq 'ZD8' or $
;      newid eq 'ZD9')
;    
;    group2 = where($
;      newid eq 'YD3' or $
;      newid eq 'YD5' or $
;      newid eq 'YD9' or $
;      newid eq 'ZD1' or $
;      newid eq 'ZD6')
;
;    group3 = where($
;      newid eq 'PD1' or $
;      newid eq 'PD2' or $
;      newid eq 'PD3' or $
;      newid eq 'PD4')
;
;; change ID numbers - again!!!! and finalize galaxy names
;    cat[where(cat.id eq '407_YD7')].id = '407_YD6'
;    cat[where(cat.id eq '400_YD6')].id = '400_YD7'
;    cat[where(cat.id eq '527_ZD5')].id = '527_ZD4'
;    cat[where(cat.id eq '120_ZD8')].id = '120_ZD5'
;    cat[where(cat.id eq '292_ZD4')].id = '292_ZD6'
;    cat[where(cat.id eq '523_ZD6')].id = '523_ZD8'
;    
;    for ii = 0, ngal-1 do cat[ii].galaxy = strmid(cat[ii].id,strpos(cat[ii].id,'_')+1)

;; jan03    
;    group1 = where($
;      cat.id eq 383 or $
;      cat.id eq 2509 or $
;      cat.id eq 1880 or $
;      cat.id eq 406 or $
;      cat.id eq 400 or $
;      cat.id eq 407 or $
;      cat.id eq 652 or $
;      cat.id eq 270 or $
;      cat.id eq 575 or $
;      cat.id eq 339 or $
;      cat.id eq 292 or $
;      cat.id eq 527 or $
;      cat.id eq 120 or $
;      cat.id eq 3903)
;    
;    group2 = where($
;      cat.id eq 534 or $
;      cat.id eq 690 or $
;      cat.id eq 429 or $
;      cat.id eq 523)
;
;    group3 = where($
;      cat.id eq 2521 or $
;      cat.id eq 286 or $
;      cat.id eq 2129 or $
;      cat.id eq 2338)
    
;    cat[group1].group = 1
;    cat[group2].group = 2
;    cat[group3].group = 3

; get rid of the tentative candidates    
;   if keyword_set(photoz) eq 0 then begin
;      index = where(cat.group ne 3)
;      cat = cat[index]
;   endif

; get rid of the tentative candidates    
    keep = where(strtrim(cat.id,2) ne 'JD1' and $
      strtrim(cat.id,2) ne '820_YD3_534' and $
      strtrim(cat.id,2) ne '971_YD11_527')
    cat = cat[keep]
    
return, cat
end    
