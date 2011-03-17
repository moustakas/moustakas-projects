pro parse_pilyugin04, result
; jm04aug02uofa

    light = 2.99792458D5 ; speed of light [km/s]

    litpath = atlas_path(/literature)
    datapath = litpath+'ABUNDANCES/'

    data = mrdfits(datapath+'2004_pilyugin_data.fits.gz',1,/silent)
    photo = mrdfits(datapath+'2004_pilyugin_photo.fits.gz',1,/silent)
    oh = rsex(datapath+'2004_pilyugin_table1.dat')
    ngalaxy = n_elements(photo)

    Umsun = +5.66
    Bmsun = +5.47
    Vmsun = +4.82

    result = {$
      galaxy:           ' ', $
      B:             -999.0, $
      B_err:         -999.0, $
      B_lum:         -999.0, $
      B_lum_err:     -999.0, $
      M_B:           -999.0, $
      M_B_err:       -999.0, $
      V:             -999.0, $
      V_err:         -999.0, $
      V_lum:         -999.0, $
      V_lum_err:     -999.0, $
      M_V:           -999.0, $
      M_V_err:       -999.0, $
      distance:      -999.0, $
      distance_err:     0.0}
    result = replicate(result,ngalaxy)

    result.galaxy = data.galaxy
    result.B = photo.B_rc3
    result.B_err = photo.B_rc3_err
    result.V = photo.V_rc3
    result.V_err = photo.V_rc3_err

; ---------------------------------------------------------------------------
; Freedman et al 2001, ApJ, 553, 47 (HST Key Project on the Hubble
; Constant); since there are no coordinates, match on galaxy name
; (mostly NGC galaxies) 
; ---------------------------------------------------------------------------

    print, format='("Reading the Freedman HST database . . . ",$)'
    readcol, litpath+'DISTANCES/freedman_distances.dat', /silent, comment='#', $
      name, fdmod, fdmoderr, format='A,X,X,X,X,X,X,F,F', delimiter=' '
    
    doit = match_string(data.galaxy,name,/exact,index=mindx,/silent)
    good = where(mindx ne -1L,ngood)
    niceprint, data[good].galaxy, name[mindx[good]], im_dmod2d(fdmod[mindx[good]])
    
    print, 'found '+string(ngood,format='(I3)')+'/'+$
      string(ngalaxy,format='(I3)')+' galaxies.'

    result[good].distance = im_dmod2d(fdmod[mindx[good]])
    
; ---------------------------------------------------------------------------
; Karachentsev et al. 2003 galaxy distance compilation (d < 5 Mpc).
; use a 45" search radius to get NGC2366.  assign a 12% error to every
; distance in this catalog as typical of the error in the TRGB
; distance technique (Karachentsev et al 2003)
; ---------------------------------------------------------------------------

    print, format='("Reading Karachentsev et al. 2003 distance database . . . ",$)'
    readcol, litpath+'DISTANCES/karachentsev03_distances.dat', /silent, comment='#', $
      name, kra1950, kde1950, kdist, format='A,A,A,X,X,X,X,F'

    raref = 15.0*im_hms2dec(data.ra)
    decref = im_hms2dec(data.dec)

    kra = 15.0*im_hms2dec(kra1950)
    kde = im_hms2dec(kde1950)
    precess, kra, kde, 1950.0, 2000.0
    
    need = where(result.distance lt -900.0,nneed)
    
    ntot = im_djs_angle_match(raref[need],decref[need],kra,kde,$
      dtheta=45.0/3600.0,units='degrees',mcount=mcount,mindx=mindx,$
      mdist=mdist)
    good = where(mindx ne -1L,ngood)
    niceprint, data[need[good]].galaxy, name[mindx[good]], kdist[mindx[good]]

    print, 'found '+string(ngood,format='(I3)')+'/'+$
      string(nneed,format='(I3)')+' galaxies.'

    result[need[good]].distance = kdist[mindx[good]]

; ---------------------------------------------------------------------------    
; Mould model flow distances
; ---------------------------------------------------------------------------    
    
    need = where(result.distance lt -900.0,nneed)

    dist = mould_distance(data[need].ra,data[need].dec,$
      data[need].z*light,object=data[need].galaxy)
    result[need].distance = dist[need].distance

    niceprint, result.galaxy, result.distance

; compute UBVJHKs absolute magnitudes and luminosities

    good = where((result.distance gt -900.0) and (result.B gt -900.0),ngood)
    if ngood ne 0L then begin

       result[good].M_B = result[good].B - 5.0*alog10(result[good].distance) - 25.0
       result[good].M_B_err = sqrt(result[good].B_err^2.0 + $
         ( (5.0/alog(10))*result[good].distance_err/result[good].distance )^2.0)

       B_lum = 10.0D^(0.4*Bmsun) * 10.0^(-0.4*result[good].M_B)     ; [L_sun]
       B_lum_err = alog(10.0) * 0.4D * B_lum * result[good].M_B_err ; [L_sun]
       
       result[good].B_lum_err = B_lum_err/B_lum/alog(10.0) ; log L_sun
       result[good].B_lum = alog10(B_lum)                  ; log L_sun

    endif

    good = where((result.distance gt -900.0) and (result.V gt -900.0),ngood)
    if ngood ne 0L then begin

       result[good].M_V = result[good].V - 5.0*alog10(result[good].distance) - 25.0
       result[good].M_V_err = sqrt(result[good].V_err^2.0 + $
         ( (5.0/alog(10))*result[good].distance_err/result[good].distance )^2.0)

       V_lum = 10.0D^(0.4*Vmsun) * 10.0^(-0.4*result[good].M_V)     ; [L_sun]
       V_lum_err = alog(10.0) * 0.4D * V_lum * result[good].M_V_err ; [L_sun]
       
       result[good].V_lum_err = V_lum_err/V_lum/alog(10.0) ; log L_sun
       result[good].V_lum = alog10(V_lum)                  ; log L_sun

    endif
    
; compute stellar masses and then append everything

    mass = compute_stellar_mass(result)

    result = struct_addtags(result,struct_trimtags(data,except='galaxy'))
    result = struct_addtags(struct_addtags(result,mass),struct_trimtags(oh,except='galaxy'))

    splog, 'Writing '+datapath+'2004_pilyugin_parsed.fits.gz'
    mwrfits, result, datapath+'2004_pilyugin_parsed.fits', /create
    spawn, ['gzip -f '+datapath+'2004_pilyugin_parsed.fits'], /sh
        
return
end
    
