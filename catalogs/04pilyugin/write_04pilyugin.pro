pro write_04pilyugin, data
; jm05mar11uofa

    Bmsun = 5.42
    Vmsun = 4.77
    
    root = '04pilyugin'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    ned = mrdfits(path+root+'_ned.fits.gz',1,/silent)
    nedphoto = mrdfits(path+root+'_ned_photo.fits.gz',1,/silent)

;   data = struct_addtags(ned,struct_trimtags(nedphoto,$
;      except=['*GALAXY*','IRAS*','*_FLAG','UV*']))
    data = ned

    table1 = rsex(path+root+'_table1.dat')
    table3 = rsex(path+root+'_table3.dat')
    table5 = rsex(path+root+'_table5.dat')
    table7 = rsex(path+root+'_table7.dat')
    tableA1 = rsex(path+root+'_tableA1.dat')

;   niceprint, data.galaxy, table1.galaxy, table3.galaxy

    data = struct_addtags(data,struct_trimtags(table1,except=['GALAXY']))
    data = struct_addtags(data,struct_trimtags(table3,except=['GALAXY','NHII_NH']))
    data = struct_addtags(data,struct_trimtags(table5,select=['OH_CHAR',$
      'GAS_FRACTION','MHI_LB','MH2_LB','VROT']))
    data = struct_addtags(data,struct_trimtags(tableA1,except=['GALAXY','LOGLB']))
    ngalaxy = n_elements(data)

; add a tag that gives the gradient in dex/kpc; also add a flag that
; identifies objects with "good" gradients: (1) not strongly barred;
; (2) more than 5 HII regions; (3) good radial coverage

    newtags = replicate({ohgradient_kpc: -999.0, ohgradient_good: 0L},ngalaxy)
    newtags.ohgradient_kpc = data.ohgradient_r25/(data.r25*!dtor*data.distance*1E3/60.0)
    data = struct_addtags(data,newtags)

    flag = where($
      (data.nhii_oh ge 5.0) and $
      (where(strmatch(data.morphology,'*SB*') eq 0B)) and $
      (where(strmatch(data.galaxy,'*4736*') eq 0B)) and $
      (where(strmatch(data.galaxy,'*4725*') eq 0B)))
    data[flag].ohgradient_good = 1L
    
; correct the RC3 photometry for foreground extinction

    glactc, 15.0D*im_hms2dec(data.ra), im_hms2dec(data.dec), 2000.0, gl, gb, 1, /degree
    ebv = dust_getval(gl,gb,/interp)

    photo = {$
      b:         -999.0, $
      b_err:     -999.0, $
      M_b:       -999.0, $
      M_b_err:   -999.0, $
      b_lum:     -999.0, $
      b_lum_err: -999.0, $
      v:         -999.0, $
      v_err:     -999.0, $
      M_v:       -999.0, $
      M_v_err:   -999.0, $
      v_lum:     -999.0, $
      v_lum_err: -999.0}
    photo = replicate(photo,n_elements(data))

    photo.b         = nedphoto.b_rc3 - ebv*k_lambda(4370.9,/odonnell)
    photo.b_err     = nedphoto.b_rc3_err
    photo.M_b       = photo.b - 5*alog10(data.distance) - 25.0
    photo.M_b_err   = photo.b_err
    b_lum = 10^(0.4*(Bmsun-photo.M_b))
    b_lum_err = alog(10.0)*0.4*b_lum*photo.M_b_err
    photo.b_lum     = alog10(b_lum)
    photo.b_lum_err = b_lum_err/b_lum/alog(10.0)

    good = where(nedphoto.v_rc3 gt -900)
    
    photo[good].v         = nedphoto[good].v_rc3 - ebv*k_lambda(5477.6,/odonnell)
    photo[good].v_err     = nedphoto[good].v_rc3_err
    photo[good].M_v       = photo[good].v - 5*alog10(data[good].distance) - 25.0
    photo[good].M_v_err   = photo[good].v_err
    v_lum = 10^(0.4*(Vmsun-photo[good].M_v))
    v_lum_err = alog(10.0)*0.4*v_lum*photo[good].M_v_err
    photo[good].v_lum     = alog10(v_lum)
    photo[good].v_lum_err = v_lum_err/v_lum/alog(10.0)
    
; compute stellar masses

    mass = compute_stellar_mass(photo)
    data = struct_addtags(struct_addtags(data,photo),mass)

;; which of these objects are either in the atlas or 11HUGS?    
;    
;    moredata = replicate({pa: 0.0, atlas: 0L, hugs: 0L},ngalaxy)
;    data = struct_addtags(data,moredata)
;
;    p04ra = 15.0D*im_hms2dec(data.ra)
;    p04dec = im_hms2dec(data.dec)
;
;    atlas = read_integrated()
;    ra = 15.0*im_hms2dec(atlas.ra)
;    dec = im_hms2dec(atlas.dec)
;    
;    ntot = im_djs_angle_match(ra,dec,p04ra,p04dec,dtheta=30.0/3600.0,$
;      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
;    match = where(mindx ne -1L,nmatch,comp=nomatch,ncomp=nnomatch)
;    if (nmatch ne 0L) then begin
;       niceprint, atlas[match].galaxy, data[mindx[match]].galaxy, data[mindx[match]].nedgalaxy
;       data[mindx[match]].atlas = 1L
;    endif
;
;; 11HUGS
;    
;    hugspath = '/home/ioannis/research/projects/11hugs/'
;    hugs = mrdfits(hugspath+'hugs_ned.fits.gz',1,/silent)
;    ra = 15.0*im_hms2dec(hugs.ra)
;    dec = im_hms2dec(hugs.dec)
;    
;    ntot = im_djs_angle_match(ra,dec,p04ra,p04dec,dtheta=30.0/3600.0,$
;      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
;    match = where(mindx ne -1L,nmatch,comp=nomatch,ncomp=nnomatch)
;    if (nmatch ne 0L) then begin
;       niceprint, hugs[match].galaxy, data[mindx[match]].galaxy, data[mindx[match]].nedgalaxy
;       data[mindx[match]].hugs = 1L
;    endif
;
;; grab the RC3 position angle
;
;    rc3 = read_rc3()
;    ra = 15.0D*im_hms2dec(rc3.ra)
;    dec = im_hms2dec(rc3.dec)
;
;    searchrad = 125.0 ; search radius
;
;    splog, 'RC3 search radius = '+string(searchrad,format='(G0.0)')+'".'
;    ntot = im_djs_angle_match(ra,dec,p04ra,p04dec,dtheta=searchrad/3600.0,$
;      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
;    match = where(mindx ne -1L,nmatch,comp=nomatch,ncomp=nnomatch)
;    if (nmatch ne 0L) then begin
;       niceprint, rc3[match].name, data[mindx[match]].galaxy, data[mindx[match]].nedgalaxy
;       data[mindx[match]].pa = rc3[match].pa
;    endif
    
; write out    

    splog, 'Writing '+path+root+'.fits.'
    mwrfits, data, path+root+'.fits', /create
    spawn, ['gzip -f ']+path+root+'.fits', /sh
    
return
end
    
