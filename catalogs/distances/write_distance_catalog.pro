;+
; NAME:
;       WRITE_DISTANCE_CATALOG
;
; PURPOSE:
;       Generate a single catalog of direct distance estimates.
;
; CALLING SEQUENCE:
;       write_distance_catalog, cat, /write
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       write - write the catalog to disk
;
; OUTPUTS:
;       cat - catalog; if WRITE=1 then write a binary FITS table
;             called 'distance_catalog.fits'. 
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       IM_READ_FMR(), RSEX()
;
; COMMENTS:
;       Distances are compiled in the following order of preference:
;       Cepheid, TRGB, Surface-Brightness Fluctuations, Tully-Fisher,
;       Group Membership (Membership), Fundamental Plane (FP), type II SN
;       (SNII), type I SN (SNI), Geometric, Average (an average of
;       many techniques), PNLF, and Brightest Stars.
; 
;       See associated routine READ_DISTANCE_CATALOG().
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Feb 23, U of A - written
;       jm05aug02uofa - added a BIBTEX reference field 
;-

pro write_distance_catalog, cat, write=write

    path = filepath('',root=getenv('CATALOGS_DIR'),subdirectory='distances')

; initialize the output catalog

    cat_template = {$
      galaxy:        '', $
      alt_galaxy:    '', $   ; alternate galaxy name
      ra:            '', $   ; J2000    
      dec:           '', $   ; J2000
      z:            -999.0,$ ; NED redshift
      distance:     -999.0,$ ; [Mpc]
      distance_err: -999.0,$ ; [Mpc]
      method:        '', $
      reference:     '', $
      texref:        ''}

; ###########################################################################
; Tonry et al. 2001
; ###########################################################################

    tonry = rsex(path+'01tonry.dat')
    ntonry = n_elements(tonry)
    cat1 = replicate(cat_template,ntonry)
    
    cat1.galaxy = strtrim(tonry.galaxy,2)
    cat1.alt_galaxy = cat1.galaxy

; use NED coordinates

;   openw, lun, '01tonry_galaxylist.dat', /get_lun
;   struct_print, struct_trimtags(cat1,sele=['ALT_GALAXY']), /no_head, lun=lun
;   free_lun, lun

;   cat1.ra = repstr(strtrim(im_dec2hms(tonry.ra/15.0),2),' ',':')
;   cat1.dec = repstr(strtrim(im_dec2hms(tonry.dec),2),' ',':')
    
    nedtonry = mrdfits(path+'01tonry_ned.fits.gz',1,/silent)

    cat1.alt_galaxy = strtrim(nedtonry.nedgalaxy,2)
    cat1.ra  = nedtonry.ra
    cat1.dec = nedtonry.dec
    cat1.z   = nedtonry.z

    dist = im_dmod2d(tonry.dmod,sigdmod=tonry.dmod_err,sigdist=dist_err) ; [Mpc]
    cat1.distance = dist
    cat1.distance_err = dist_err

    cat1.method = 'SBF'
    cat1.reference = 'Tonry et al. 2001, ApJ, 546, 681'
    cat1.texref = 'tonry01'

    cat = cat1 ; NOTE!
    
; ###########################################################################
; Freedman et al. 2001
; ###########################################################################

    freedman = rsex(path+'01freedman.dat')
    nfreedman = n_elements(freedman)
    cat1 = replicate(cat_template,nfreedman)
    
    cat1.galaxy = strtrim(freedman.galaxy,2)
    cat1.alt_galaxy = cat1.galaxy

;   cat1.ra = freedman.ra
;   cat1.dec = freedman.dec
    
    nedfreedman = mrdfits(path+'01freedman_ned.fits.gz',1,/silent)

    cat1.alt_galaxy = strtrim(nedfreedman.nedgalaxy,2)
    cat1.ra  = nedfreedman.ra
    cat1.dec = nedfreedman.dec
    cat1.z   = nedfreedman.z

; use the error in the uncorrected distance modulus    
    
    dist = im_dmod2d(freedman.dmod_cor,sigdmod=freedman.dmod_uncor_err,sigdist=dist_err) ; [Mpc]
;   dist = im_dmod2d(freedman.dmod_uncor,sigdmod=freedman.dmod_uncor_err,sigdist=dist_err)
    cat1.distance = dist
    cat1.distance_err = dist_err

    cat1.method = 'Cepheid'
    cat1.reference = 'Freedman et al 2001, ApJ, 553, 47'
    cat1.texref = 'freedman01'

    cat = [cat,cat1]
    
; ###########################################################################
; Ferrarese et al. 2000
; ###########################################################################

    ferrarese = rsex(path+'00ferrarese.dat')
    nferrarese = n_elements(ferrarese)
    cat1 = replicate(cat_template,nferrarese)
    
    cat1.galaxy = strtrim(ferrarese.galaxy,2)
    cat1.alt_galaxy = cat1.galaxy

;   openw, lun, '00ferrarese_galaxylist.dat', /get_lun
;   struct_print, struct_trimtags(cat1,sele=['ALT_GALAXY']), /no_head, lun=lun
;   free_lun, lun

;   cat1.ra = ferrarese.ra
;   cat1.dec = ferrarese.dec

    nedferrarese = mrdfits(path+'00ferrarese_ned.fits.gz',1,/silent)
 
    cat1.alt_galaxy = strtrim(nedferrarese.nedgalaxy,2)
    cat1.ra  = nedferrarese.ra
    cat1.dec = nedferrarese.dec
    cat1.z   = nedferrarese.z

    good = where(ferrarese.dmod_cepheid gt -900,ngood)
    if (ngood ne 0L) then begin

       dist = im_dmod2d(ferrarese[good].dmod_cepheid,sigdmod=ferrarese[good].dmod_cepheid_err,sigdist=dist_err) ; [Mpc]
       cat1[good].distance = dist
       cat1[good].distance_err = dist_err
       cat1[good].method = 'Cepheid'

    endif
    
    good = where(ferrarese.I_trgb gt -900,ngood)
    if (ngood ne 0L) then begin

       I_trgb = ferrarese[good].I_trgb 
       I_trgb_err = ferrarese[good].I_trgb_err
       
; correct for foreground Galactic extinction using the O'Donnell
; (1994) Galactic extinction law

       glactc, 15.0D*im_hms2dec(cat1[good].ra), im_hms2dec(cat1[good].dec), 2000.0, gl, gb, 1, /degree
       ebv_mw = dust_getval(gl,gb,/interp)
       
       I_trgb = I_trgb - ebv_mw*k_lambda(8020.0,/odonnell) ; lambda_eff from Bessell (1990)
       
       MI_trgb = -4.06    ; from Ferrarese et al. 2000, ApJ, 529, 745
       MI_trgb_err = 0.15 ; includes random and systematic uncertainty

       dmod = I_trgb - MI_trgb
       dmod_err = MI_trgb_err
       
       dist = im_dmod2d(dmod,sigdmod=dmod_err,sigdist=dist_err) ; [Mpc]
       cat1[good].distance = dist
       cat1[good].distance_err = dist_err
       cat1[good].method = 'TRGB'

    endif
    
    cat1.reference = 'Ferrarese et al 2000, ApJS, 128, 431'
    cat1.texref = 'ferrarese00'

    cat = [cat,cat1]
    
;;; ###########################################################################
;;; Karachentsev et al. 2003
;;; ###########################################################################
;;
;;    karachentsev = rsex(path+'03karachentsev.dat')
;;    nkarachentsev = n_elements(karachentsev)
;;    cat1 = replicate(cat_template,nkarachentsev)
;;    
;;; alter some of the galaxy names to be compatible with NED; use the
;;; NED coordinates 
;;
;;    cat1.galaxy = strtrim(karachentsev.galaxy,2)
;;    cat1.alt_galaxy = cat1.galaxy
;;
;;    replace = where((strmatch(cat1.alt_galaxy,'*KK*',/fold) eq 1B) and $
;;      (strmatch(cat1.alt_galaxy,'*KKH*',/fold) eq 0B) and $
;;      (strmatch(cat1.alt_galaxy,'*KKR*',/fold) eq 0B),nreplace)
;;    if (nreplace ne 0L) then cat1[replace].alt_galaxy = repstr(cat1[replace].alt_galaxy,'KK','[KK98]')
;;
;;    cat1[0].alt_galaxy  = 'ESO349-G031'
;;    cat1[3].alt_galaxy  = 'ESO294-G010'
;;    cat1[12].alt_galaxy = 'ESO245-G005'
;;    cat1[21].alt_galaxy = 'ESO115-G021'
;;    cat1[39].alt_galaxy = '[KM96]055451+072840'
;;    cat1[40].alt_galaxy = 'ESO489-G?056'
;;    cat1[41].alt_galaxy = 'ESO490-G017'
;;    cat1[42].alt_galaxy = 'ARGODWARF'
;;    cat1[60].alt_galaxy = 'ArpsLoop'
;;    cat1[67].alt_galaxy = 'AntliaDwarf'
;;    cat1[70].alt_galaxy = 'HIJASSJ1021+6842'
;;    cat1[71].alt_galaxy = '[HS98]117'
;;    cat1[82].alt_galaxy = 'ESO379-G007'
;;    cat1[85].alt_galaxy = 'ESO321-G014'
;;    cat1[109].alt_galaxy = 'ESO269-G058'
;;    cat1[111].alt_galaxy = 'ESO269-G?066'
;;    cat1[120].alt_galaxy = 'ESO324-G024'
;;    cat1[126].alt_galaxy = 'UGCA365'
;;    cat1[128].alt_galaxy = 'HIPASSJ1337-39'
;;    cat1[129].alt_galaxy = 'ESO444-G084'
;;    cat1[135].alt_galaxy = 'ESO325-G?011'
;;    cat1[136].alt_galaxy = 'HIPASSJ1348-37'
;;    cat1[137].alt_galaxy = 'HIPASSJ1351-47'
;;    cat1[141].alt_galaxy = 'ESO384-G016'
;;    cat1[146].alt_galaxy = 'UKS1424-460'
;;    
;;;   openw, lun, '03karachentsev_galaxylist.dat', /get_lun
;;;   struct_print, struct_trimtags(cat1,sele=['ALT_GALAXY']), /no_head, lun=lun
;;;   free_lun, lun
;;
;;; old code: precess the coordinates
;;
;;;   kra = 15.0*im_hms2dec(karachentsev.ra)
;;;   kdec = im_hms2dec(karachentsev.dec)
;;;   precess, kra, kdec, 1950.0, 2000.0
;;;
;;;   cat1.ra = repstr(strtrim(im_dec2hms(kra/15.0),2),' ',':')
;;;   cat1.dec = repstr(strtrim(im_dec2hms(kdec),2),' ',':')
;;
;;    nedkarachentsev = mrdfits(path+'03karachentsev_ned.fits.gz',1,/silent)
;;
;;    cat1.alt_galaxy  = strtrim(nedkarachentsev.nedgalaxy,2)
;;    cat1.ra  = nedkarachentsev.ra
;;    cat1.dec = nedkarachentsev.dec
;;    cat1.z   = nedkarachentsev.z
;;
;;; assume a 25% error in distances based on the brightest stars and a
;;; 15% error in all other techniques
;;
;;    dist = karachentsev.distance
;;    dist_err = dist*0.15
;;    
;;    stars = where(strmatch(karachentsev.method,'*stars*',/fold) eq 1B,nstars)
;;    if (nstars ne 0L) then dist_err[stars] = dist[stars]*0.20
;;
;;    cat1.distance = dist
;;    cat1.distance_err = dist_err
;;
;;    cat1.method = karachentsev.method
;;    cat1.reference = 'Karachentsev et al. 2003, A&A, 398, 479'
;;
;;    cat = [cat,cat1]

; ###########################################################################
; Karachentsev et al. 2004
; ###########################################################################

    karachentsev = im_read_fmr(path+'04karachentsev.dat')
    nkarachentsev = n_elements(karachentsev)
    cat1 = replicate(cat_template,nkarachentsev)

; parse the data table

    ra = string(karachentsev.rah,format='(I2.2)')+':'+string(karachentsev.ram,format='(I2.2)')+':'+$
      repstr(string(karachentsev.ras,format='(F4.1)'),' ','0')
    dec = karachentsev.de_+string(karachentsev.ded,format='(I2.2)')+':'+string(karachentsev.dem,format='(I2.2)')+':'+$
      repstr(string(karachentsev.des,format='(F4.1)'),' ','0')
    
    name = strtrim(karachentsev.name,2) & galaxy = name
    for i = 0L, nkarachentsev-1L do if strmatch(name[i],'*,*') then galaxy[i] = strmid(name[i],0,strpos(name[i],','))
    
    cat1.galaxy = galaxy
    cat1.alt_galaxy = cat1.galaxy

; repair the coordinates for NGC5194

    ra[364] = '13:29:52.7'
    dec[364] = '+47:11:43'
    
; cross-match the Karachentsev catalog with the NED literature search;
; a search radius of 325" gets everything, including the Saggitarius
; Dwarf; note: only 450/451 objects match, because for some reason
; Karachentsev includes the Milky Way as a galaxy in his sample; the
; distance? 0.01 Mpc!

    nedkarachentsev = mrdfits(path+'04karachentsev_ned.fits.gz',1,/silent)

    raref = im_hms2dec(ra)*15.0
    decref = im_hms2dec(dec)

    raned = im_hms2dec(nedkarachentsev.ra)*15.0
    decned = im_hms2dec(nedkarachentsev.dec)

    ntot = im_djs_angle_match(raref,decref,raned,decned,dtheta=325.0/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist,mmax=1)
    good = where(mindx ne -1L,ngood)

    srt = sort(mdist[good])
    niceprint, nedkarachentsev[mindx[good]].galaxy, name[good], mdist[good]*3600.0
;   niceprint, nedkarachentsev[mindx[good[srt]]].galaxy, name[good[srt]], mdist[good[srt]]*3600.0

    nedkarachentsev = nedkarachentsev[mindx[good]]
    karachentsev = karachentsev[good]
    cat1 = cat1[good]
    
;   bigindx = lindgen(nkarachentsev)
;   remove, good, bigindx
    
    cat1.galaxy  = strtrim(nedkarachentsev.galaxy,2)
    cat1.ra  = nedkarachentsev.ra
    cat1.dec = nedkarachentsev.dec
    cat1.z   = nedkarachentsev.z
    
; distance method; Note!  replace "fj" (="Faber-Jackson") with "FP"
; (="Fundamental Plane")

    method = karachentsev.n_dis
    cat1.method = strtrim(repstr(repstr(repstr(repstr(repstr(repstr(repstr(method,'rgb','TRGB'),$
      'cep','Cepheid'),'bs','Stars'),'mem','Membership'),'tf','TF'),'sbf','SBF'),'fj','FP'),2)

; only retain objects with direct distance estimates

    keep = where(strmatch(cat1.method,'h') eq 0B,nkeep)
    cat1 = cat1[keep]
    karachentsev = karachentsev[keep]
    
; assume a 25% error in distances based on the brightest stars and a
; 15% error in all other techniques

    dist = karachentsev.dis
    dist_err = dist*0.15
    
    stars = where(strmatch(cat1.method,'*stars*',/fold) eq 1B,nstars)
    if (nstars ne 0L) then dist_err[stars] = dist[stars]*0.20

    cat1.distance = dist
    cat1.distance_err = dist_err
    cat1.reference = 'Karachentsev et al. 2004, AJ, 127, 2031'
    cat1.texref = 'karachentsev04'

    cat = [cat,cat1]

;; ###########################################################################
;; Whiting 2003 - basically everything in Whiting is in Karachentsev 
;; ###########################################################################
;
;    whiting = im_read_fmr(path+'03whiting.dat')
;    nwhiting = n_elements(whiting)
;    cat1 = replicate(cat_template,nwhiting)
;
;; alter some of the galaxy names to be compatible with NED; use the
;; NED coordinates 
;
;    cat1.galaxy = strtrim(whiting.name,2)
;    cat1.alt_galaxy = cat1.galaxy
;
;    cat1[0].alt_galaxy = '[KM96]055451+072840'
;    cat1[[1,67]].alt_galaxy = strtrim(cat1[[1,67]].galaxy,2)+'Dwarf'
;    cat1[23:24].alt_galaxy = repstr(cat1[23:24].galaxy,'KK','[KK98]')
;    
;;   openw, lun, '03whiting_galaxylist.dat', /get_lun
;;   struct_print, struct_trimtags(cat1,sele=['ALT_GALAXY']), /no_head, lun=lun
;;   free_lun, lun
;
;; old code: convert from galactic to equatorial coordinates
;    
;;   glactc, ra, dec, 2000.0, whiting.glon, whiting.glat, 2, /degree, /supergalactic
;;   cat1.ra = repstr(strtrim(im_dec2hms(ra/15.0),2),' ',':')
;;   cat1.dec = repstr(strtrim(im_dec2hms(dec),2),' ',':')
;
;    nedwhiting = mrdfits(path+'03whiting_ned.fits.gz',1,/silent)
;
;    cat1.alt_galaxy  = strtrim(nedwhiting.nedgalaxy,2)
;    cat1.ra  = nedwhiting.ra
;    cat1.dec = nedwhiting.dec
;    cat1.z   = nedwhiting.z
;
;; assume a 25% error in distances based on the brightest stars and a
;; 15% error in all other techniques
;
;    dist = whiting.dist
;    dist_err = dist*0.15
;
;    stars = where(strmatch(whiting.meth,'*stars*',/fold) eq 1B,nstars)
;    if (nstars ne 0L) then dist_err[stars] = dist[stars]*0.20
;
;    cat1.distance = dist
;    cat1.distance_err = dist_err
;
;    cat1.method = repstr(repstr(whiting.meth,'stars','Stars'),'ceph','Cepheid')
;    cat1.reference = 'Whiting 2003, ApJ, 587, 186'
;
;; Note!  NGC3274 and UGC05721 are the same object and they appear
;; twice in Whitting's catalog - remove UGC05721
;
;    good = where(strmatch(cat1.galaxy,'*5721*') eq 0B)
;    cat1 = cat1[good]
;
;    cat = [cat,cat1]
    
; ###########################################################################
; Shapley et al. 2001 - two objects are in 01shapley.dat but *not* in
;                       01shapley_table1.dat for some unknown reason:
;                       NGC3660 and NGC5079
; ###########################################################################

    shapley = im_read_fmr(path+'01shapley.dat')
;   shapley2 = im_read_fmr(path+'01shapley_table1.dat')
    nshapley = n_elements(shapley)
    cat1 = replicate(cat_template,nshapley)

; alter some of the galaxy names to be compatible with NED; use the
; NED coordinates 

    cat1.galaxy = strtrim(shapley.name,2)
    cat1.alt_galaxy = cat1.galaxy

;   openw, lun, '01shapley_galaxylist.dat', /get_lun
;   struct_print, struct_trimtags(cat1,sele=['ALT_GALAXY']), /no_head, lun=lun
;   free_lun, lun

;   openw, lun, '01shapley_galaxylist2.dat', /get_lun
;   struct_print, struct_trimtags(shapley2,sele=['NAME']), /no_head, lun=lun
;   free_lun, lun

; old code: precess the coordinates

;   sra = 15.0*im_hms2dec(shapley.ra)
;   sdec = im_hms2dec(shapley.dec)
;   precess, sra, sdec, 1950.0, 2000.0
;
;   cat1.ra = repstr(strtrim(im_dec2hms(sra/15.0),2),' ',':')
;   cat1.dec = repstr(strtrim(im_dec2hms(sdec),2),' ',':')

    nedshapley = mrdfits(path+'01shapley_ned.fits.gz',1,/silent)

    cat1.alt_galaxy  = strtrim(nedshapley.nedgalaxy,2)
    cat1.ra  = nedshapley.ra
    cat1.dec = nedshapley.dec
    cat1.z   = nedshapley.z

; only retain objects with literature distances

    good = where(shapley.m_m gt -900.0,ngood)
    if (ngood ne 0L) then cat1 = cat1[good]

; read the methods used for each object    

;   openw, lun, '01shapley_methods.dat', /get_lun & niceprintf, lun, cat1.galaxy & free_lun, lun
    readcol, path+'01shapley_methods.dat', mgalaxy, method, format='A,A', /silent

    cat1.method = method
    cat1.reference = 'Shapley 2001, ApJS, 137, 139'
    cat1.texref = 'shapley01'

; assign the final distances and distance errors    
    
    cat1.distance = im_dmod2d(shapley[good].m_m)
    cat1.distance_err = cat1.distance*0.10

    stars = where(strmatch(cat1.method,'*stars*',/fold) eq 1B,nstars)
    if (nstars ne 0L) then cat1[stars].distance_err = cat1[stars].distance*0.15
    
    cat = [cat,cat1]

; write out
    
    if keyword_set(write) then begin

       srt = sort(im_hms2dec(cat.ra))
       cat = cat[srt]
       
       outfile = 'distance_catalog.fits'
       splog, 'Writing '+path+outfile+'.'
       mwrfits, cat, path+outfile, /create
       spawn, ['gzip -f '+path+outfile], /sh
       
    endif
    
return
end
    
