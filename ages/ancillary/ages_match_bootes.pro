pro ages_match_bootes
; jm10jan06ucsd - match the AGES and BOOTES catalogs; all the BOOTES
;   catalogs are line-matched to the I-band catalog

    common match_bootes, ages, alliband

    suffix = '2011a' ; '2010b' ; '2009b'
    
    agespath = ages_path(/catalogs)
    catpath = ages_path(/mycatalogs)
    bootespath = getenv('IM_ARCHIVE_DIR')+'/data/bootes/'+suffix+'/'

    if (n_elements(ages) eq 0) or (n_elements(alliband) eq 0) then begin
       agesfile = agespath+'catalog.cat.noguidestars.fits.gz'
       splog, 'Reading '+agesfile
       ages = mrdfits(agesfile,1)
       ages.ra = 15.0D*ages.ra  ; convert to decimal degrees

       bootesfile = bootespath+'bootes_I.fits.gz'
       splog, 'Reading '+bootesfile
       alliband = mrdfits(bootesfile,1)
    endif
    ngal = n_elements(ages)
    
    splog, 'Matching...'
    m1 = im_spherematch(alliband,ages,match2=m2,radius=1.0,$
      ratagname1='alpha_j2000',dectagname1='delta_j2000')

; now check if any of the matches are duplicates and, if so, choose
; the one that has FLAG_DUPLICATE=0
    dup = where(alliband[m1].flag_duplicate eq 1,ndup)
    for ii = 0, ndup-1 do begin
       spherematch, alliband.alpha_j2000, alliband.delta_j2000, $
         alliband[m1[dup[ii]]].alpha_j2000, alliband[m1[dup[ii]]].delta_j2000, $
         1.0/3600.0, d1, d2, max=0
       splog, ii, m2[dup[ii]], n_elements(d1)
       this = where(alliband[d1].flag_duplicate eq 0,nthis)
       if (nthis ne 0) then m1[dup[ii]] = d1[this]
    endfor
       
    bootes = im_empty_structure(alliband[0],$
      empty_value=-999.0,ncopies=ngal)
    bootes[m2] = alliband[m1]

;   tosstags1 = ['*aper_05','*aper_06','*aper_07','*aper_08',$
;     '*aper_09','*aper_10','*aper_15','*aper_20']
    tosstags2 = ['*object_position*','*field*','alpha_j2000',$
      'delta_j2000','flag_duplicate','flux_*','fluxerr_*',$
      'imaflags_*','segflags_*']
;   tosstags2 = ['*aper_01','*aper_02','*aper_03','*aper_05',$
;     '*aper_07','*aper_08','*aper_09','*aper_10',$
;     '*aper_15','*aper_20','field','*object_position*',$
;     'alpha_j2000','delta_j2000','flag_*']
    bootes = struct_trimtags(temporary(bootes),except=tosstags1)
    bootes = im_struct_trimtags(bootes,select=tag_names(bootes),$
      newtags='I_'+tag_names(bootes))

; store the index number in the original catalog of 2.5 million
; objects (used in SFRM_FIND_ISOLATED)
;   moretags = replicate({bootes_object_position: -1L},ngal)
;   bootes = struct_addtags(moretags,temporary(bootes))
;   bootes[m2].bootes_object_position = m1 

;; are the missing objects near bright stars?
;    usnofile = agespath+'catalog.usno.fits.gz'
;    splog, 'Reading '+usnofile
;    usno = mrdfits(usnofile,1)
;
;    miss = where((bootes.i_mag_auto lt 0.0))
;    struct_print, usno[miss]

; now read the rest of the catalogs, but just the objects we care
; about 
    bands = ['u','Bw','R','z','y','J','H',$
      'Ks','ch1','ch2','ch3','ch4']
;   bands = ['FUV','NUV','Bw','R','J','H',$
;     'Ks','ch1','ch2','ch3','ch4']
    for ii = 0, n_elements(bands)-1 do begin
       bootesfile = bootespath+'bootes_'+bands[ii]+'.fits.gz'
       splog, 'Reading '+bootesfile
       allcat = mrdfits(bootesfile,1,rows=m1) 

       cat = im_empty_structure(allcat[0],$
         empty_value=-999.0,ncopies=ngal)
       cat[m2] = temporary(allcat)

       cat = struct_trimtags(temporary(cat),except=tosstags2)
;      cat = struct_trimtags(temporary(cat));,except=tosstags2)
       cat = im_struct_trimtags(cat,select=tag_names(cat),$
         newtags=bands[ii]+'_'+tag_names(cat))
       
       bootes = struct_addtags(temporary(bootes),temporary(cat))
    endfor

; write out    
    im_mwrfits, bootes, catpath+'ages_bootes_'+suffix+'.fits', /clobber
    
return
end
    
