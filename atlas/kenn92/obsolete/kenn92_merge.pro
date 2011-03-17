pro kenn92_merge, table, write=write
; jm02may14uofa
; merge the basic data output and the photometry output from ned and
; compute several more useful quantities
    
    forward_function dangular, redh100
    
    red
    h100 = redh100()

    mpc2m = 3.086D22            ; meter per Mpc
    ergs2watts = 1D7            ; erg/s per Watts
    lsun = 3.826D26             ; Watts
    msun = +5.47                ; the Sun's absolute B-band magnitude
    msun_bol = +4.74            ; the Sun's absolute bolometric magnitude
    light = 2.99792458D5        ; speed of light [km/s]
    
    datapath = atlas_path(/kenn92)+'analysis/'
    outpath = datapath
    nedpath = datapath

; retrieve the NED basic data and photometry on the atlas

    basicdata = mrdfits(nedpath+'kenn92_ned.fits.gz',1,/silent)
    photodata = mrdfits(nedpath+'kenn92_ned_photo.fits.gz',1,/silent)
    ngalaxy = n_elements(basicdata)

    if ngalaxy ne n_elements(photodata) then message, 'Incompatible dimensions!'

; we assume that BASICDATA and PHOTODATA are sorted correctly

;   niceprint, basicdata.galaxy, photodata.galaxy

; ---------------------------------------------------------------------------    
; initialize the output data table
; ---------------------------------------------------------------------------    

    table = {$
      galaxy: '',                                 $ ; galaxy
      ned_galaxy:  ' ',                             $ ; NED galaxy name
      alt_galaxy: ' ',                            $ ; alternate galaxy name
      nice_galaxy: ' ',                            $ ; nice galaxy name
      ra: '', dec: ' ',                            $ ; position
      driftdate:   '1985-1991',                   $
      driftexptime: 0.0,                          $
      driftcode:   -999L,                         $
      driftscan:    0.0,                          $
      driftap:      0.0,                          $
      driftpa:      0.0,                          $
      driftnote:    ' ',                          $
      driftphotflag: 'Y',                         $
      driftcomments: ' ',                         $
      cz: -999.0D, cz_err: -999.0D,               $ ; redshift and error
      xscale: -999.0,                             $ ; [pc/arcsec] at that distance
      ned_morphology: '',                         $ ; NED morphology and type
      rc3_type: '',                               $ ; RC3/LEDA morphological type
      rc3_t: -999.0,                              $ ; RC3/LEDA numerical type (continuous)
      ned_class: '',                              $ ; AGN according to NED? [U/Y/N] (U=Unknown)
      distance: -999.0,                           $ ; distance [Mpc]
      distance_err: -999.0,                       $ ; distance error [Mpc]
      distance_ref: '...',                        $ ; distance reference
      ebv_mw: -999.0,                             $ ; Galactic E(B-V) (SFD)
      U: -999.0, U_err: -999.0,                   $ ; RC3/LEDA total apparent U magnitude and error
      B: -999.0, B_err: -999.0,                   $ ; RC3/LEDA total apparent B magnitude and error
      V: -999.0, V_err: -999.0,                   $ ; RC3/LEDA total apparent V magnitude and error
      I: -999.0, I_err: -999.0,                   $ ; RC3/LEDA total apparent I magnitude and error
      J: -999.0,  J_err: -999.0,                  $ ; 2MASS J magnitude and error
      H: -999.0,  H_err: -999.0,                  $ ; 2MASS H magnitude and error
      Ks: -999.0,  Ks_err: -999.0,                $ ; 2MASS Ks magnitude and error
      M_B: -999.0, M_B_err: -999.0,               $ ; absolute B magnitude and error (normalized to B)
      B_lum: -999.0, B_lum_err: -999.0,           $ ; B luminosity and error (normalized to B)
      B_lum_bol: -999.0, B_lum_bol_err: -999.0,   $ ; B luminosity and error (bolometric)
      Ks_lum: -999.0, Ks_lum_err: -999.0,         $ ; Ks luminosity and error
      iras_12: -999.0, iras_12_err: -999.0,       $ ; IRAS 12 micron flux and error [Jy] (upper limits are negative)
      iras_25: -999.0, iras_25_err: -999.0,       $ ; IRAS 25 micron flux and error [Jy]
      iras_60: -999.0, iras_60_err: -999.0,       $ ; IRAS 60 micron flux and error [Jy]
      iras_100: -999.0, iras_100_err: -999.0,     $ ; IRAS 100 micron flux and error [Jy]
      fir_flux: -999.0,                           $ ; FIR flux (Helou et al 1988)
      fir_lum: -999.0, fir_lum_err: -999.0,       $ ; FIR luminosity
      ir_lum: -999.0, ir_lum_err: -999.0,         $ ; total IR luminosity (Calzetti et al 2001)
      ir_sfr: -999.0,                             $ ; IR star-formation rate (Kennicutt 1998)
      fir_opt_ratio: -999.0, fir_opt_ratio_err: -999.0,$ ; FIR-to-optical luminosity ratio and error
      ir_opt_ratio: -999.0, ir_opt_ratio_err:   -999.0,$ ; IR-to-optical luminosity ratio and error
      pa: -999.0,                                 $ ; position angle (N-->E)
      d25_maj: -999.0,                            $ ; D25 major axis [arcmin]
      d25_min: -999.0,                            $ ; d25 minor axis [arcmin]
      inclination: -999.0}                          ; inclination angle
    table = replicate(table,ngalaxy)
    
; fill the structure with preliminaries
    
    table.galaxy = strcompress(basicdata.galaxy,/remove)
    table.ra = basicdata.ra
    table.dec = basicdata.dec
    table.ned_galaxy = strcompress(basicdata.nedgalaxy,/remove)
    table.ned_morphology = strcompress(basicdata.morph,/remove)

    ra = 15.0D*im_hms2dec(table.ra)
    dec = im_hms2dec(table.dec)

; read some extra information and add it
    
    readcol, datapath+'kenn92_extra_info.txt', obj, nice_galaxy, alt_galaxy, $
      rc3type, aperture, scanlen, pa, code, comments, /silent, $
      comment='#', delimiter='|', format='A,A,A,A,F,F,F,L,A'

    table.nice_galaxy = strtrim(nice_galaxy,2)
    table.alt_galaxy = strtrim(alt_galaxy,2)
    table.driftscan = scanlen
    table.driftcode = code
    table.driftap = aperture
    table.driftpa = pa
    table.driftcomments = comments
    table.rc3_type = rc3type
    
; ---------------------------------------------------------------------------
; based on NED classification, determine which objects have nuclear AGN
; ---------------------------------------------------------------------------

    unknown = where(strcompress(basicdata.morph,/remove) eq '',nunknown,comp=known,ncomp=nknown)
    if nunknown ne 0L then table[unknown].ned_class = 'U'

    if nknown ne 0L then begin

;      ned_hii = where(stregex(basicdata[known].morph,'HII',/fold) gt -1L,nhii)
;      if nhii ne 0L then table[known[ned_hii]].ned_class = strjoin([table[known[ned_hii]].ned_class,'H'],', ')

       wagn = where(((stregex(basicdata[known].morph,'LINER',/fold) gt -1L) + $
         (stregex(basicdata[known].morph,'Sy[1-2]',/fold) gt -1L)) ge 1L,nagn)
       if nagn ne 0L then table[known[wagn]].ned_class = 'Y'
    endif

; ---------------------------------------------------------------------------
; determine the E(B-V) color excess in our Galaxy along the line of
; sight using the SFD dust maps.  convert RA and DEC to Galactic
; coordinates.
; ---------------------------------------------------------------------------

    glactc, 15.0D*im_hms2dec(table.ra), im_hms2dec(table.dec), 2000.0, gl, gb, 1, /degree
    table.ebv_mw = dust_getval(gl,gb,/interp)

; ---------------------------------------------------------------------------
; redshift information
; ---------------------------------------------------------------------------
    
    zinfo = where(basicdata.z gt -900.0,nz)
    if nz ne 0L then table[zinfo].cz = basicdata[zinfo].z*light ; km/s

    zerrinfo = where(basicdata.z_err gt -900.0,nzerr)
    if nzerr ne 0L then table[zerrinfo].cz_err = basicdata[zerrinfo].z_err*light ; km/s

; ---------------------------------------------------------------------------
; compute Mould model distances for all the galaxies; assign a model
; distance error that goes as D/(1+D) 
; ---------------------------------------------------------------------------

    czgood = where(table.cz gt -900.0,nczgood,comp=czbad,ncomp=nczbad)
    if (nczgood ne 0L) then begin

       mould = mould_distance(table[czgood].ra,table[czgood].dec,table[czgood].cz,$
         object=table[czgood].galaxy)
       table[czgood].distance = mould.distance
       table[czgood].distance_ref = 'Flow Model'

; assign a model distance error; choose an exponential model with an
; e-folding distance of 10 Mpc.  this corresponds to a fractional
; error of 20% at Virgo (D=15 Mpc)

       alpha = 1.0/10.0 ; [1/Mpc]

       table[czgood].distance_err = table[czgood].distance*exp(-alpha*table[czgood].distance)

    endif

; ---------------------------------------------------------------------------    
; compute the physical scale [pc/arcsec] at the distance of every
; galaxy 
; ---------------------------------------------------------------------------

    good = where((table.distance gt -900.0) and (table.cz gt -900.0),ngood)
    if ngood ne 0L then begin
       dz = h100*100.0*table[good].distance/light      ; cosmological redshift
       table[good].xscale = dangular(dz,/pc)/206265.0D ; [pc/arcsec]
;      niceprint, table[good].galaxy, dz, table[good].distance, table[good].xscale
    endif
    
; ---------------------------------------------------------------------------
; 2MASS near-infrared photometry.  correct for Galactic extinction.
; we use the Cardelli, Clayton, & Mathis (1989) extinction curve and
; the O'Donnell (1994) coefficients with the effective filter
; wavelengths tabulated in Allen 2000, Table 7.5
; ---------------------------------------------------------------------------

    good = where(photodata.J_2MASS gt -900.0,ngood)
    if ngood ne 0L then table[good].J  = photodata[good].J_2MASS - table[good].EBV_MW*k_lambda(12150.0,/odonnell)
    good = where(photodata.J_2MASS_ERR gt -900.0,ngood)
    if ngood ne 0L then table[good].J_ERR  = photodata[good].J_2MASS_ERR

    good = where(photodata.H_2MASS gt -900.0,ngood)
    if ngood ne 0L then table[good].H  = photodata[good].H_2MASS - table[good].EBV_MW*k_lambda(16540.0,/odonnell)
    good = where(photodata.H_2MASS_ERR gt -900.0,ngood)
    if ngood ne 0L then table[good].H_ERR  = photodata[good].H_2MASS_ERR

    good = where(photodata.K_2MASS gt -900.0,ngood)
    if ngood ne 0L then table[good].Ks  = photodata[good].K_2MASS - table[good].EBV_MW*k_lambda(21570.0,/odonnell)
    good = where(photodata.K_2MASS_ERR gt -900.0,ngood)
    if ngood ne 0L then table[good].Ks_ERR  = photodata[good].K_2MASS_ERR

; ---------------------------------------------------------------------------
; IRAS fluxes and errors
; ---------------------------------------------------------------------------
    
    table.IRAS_12 = photodata.IRAS_12
    table.IRAS_12_ERR = photodata.IRAS_12_ERR

    table.IRAS_25 = photodata.IRAS_25
    table.IRAS_25_ERR = photodata.IRAS_25_ERR

    table.IRAS_60 = photodata.IRAS_60
    table.IRAS_60_ERR = photodata.IRAS_60_ERR

    table.IRAS_100 = photodata.IRAS_100
    table.IRAS_100_ERR = photodata.IRAS_100_ERR

; ---------------------------------------------------------------------------
; radio data, errors, and associated quantities
; ---------------------------------------------------------------------------

    

; ---------------------------------------------------------------------------
; RC3 magnitudes from NED.  do not use "corrected" RC3 magnitudes 
; ---------------------------------------------------------------------------

    good = where(strmatch(strtrim(photodata.U_FLAG,2),'U_T') eq 1B,ngood)
    if ngood ne 0L then begin
       table[good].U = photodata[good].U_RC3
       table[good].U_err = photodata[good].U_RC3_err
    endif

    good = where((strmatch(strtrim(photodata.B_FLAG,2),'B_T') eq 1B) or $
      (strmatch(strtrim(photodata.B_FLAG,2),'m_B') eq 1B),ngood)
    if ngood ne 0L then begin
       table[good].B = photodata[good].B_RC3
       table[good].B_err = photodata[good].B_RC3_err
    endif

    good = where(strmatch(strtrim(photodata.V_FLAG,2),'V_T') eq 1B,ngood)
    if ngood ne 0L then begin
       table[good].V = photodata[good].V_RC3
       table[good].V_err = photodata[good].V_RC3_err
    endif

; ---------------------------------------------------------------------------
; read the RC3 catalog to retrieve the morphological type, position
; angle angle, inclination, R25, D25_MAJ, and D25_MIN of every galaxy.
; first we match on position to crop the RC3 catalog to a more
; manageable size.  then we hand-check the position matching and
; remove "contaminants";  all the galaxy matches were checked by hand
; (jm03sep15uofa)
; ---------------------------------------------------------------------------
    
    rc3 = read_rc3()
    rc3ra = 15.0D*im_hms2dec(rc3.ra)
    rc3dec = im_hms2dec(rc3.dec)

    searchrad = 125.0 ; search radius

    splog, 'RC3 search radius = '+string(searchrad,format='(G0.0)')+'".'
    ntot = im_djs_angle_match(ra,dec,rc3ra,rc3dec,dtheta=searchrad/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)

    match = where(mindx ne -1L,nmatch,comp=nomatch,ncomp=nnomatch)

;   s = sort(mdist[match])
;   for i = 0L, nmatch-1L do print, i, strn(table[match[s[i]]].galaxy), $
;     strcompress(rc3[mindx[match[s[i]]]].name,/remove), $
;     strcompress(rc3[mindx[match[s[i]]]].altname,/remove), $
;     mdist[match[s[i]]]*3600.0, mindx[match[s[i]]], format='(I5,A25,A15,A25,F10.3,I10)'

    big = where(mdist[match]*3600.0 gt 30.0,nbig)
    percent = 100.0*nbig/nmatch

    splog, 'The following objects ('+string(nbig,format='(I0)')+'/'+$
      string(nmatch,format='(I0)')+', '+strtrim(string(percent,format='(F4.1)'),2)+$
      '%) have position differences greater than 30".'

    s = sort(mdist[match[big]])
    for i = 0L, nbig-1L do print, i, strn(table[match[big[s[i]]]].galaxy), $
      strcompress(rc3[mindx[match[big[s[i]]]]].name,/remove), $
      strcompress(rc3[mindx[match[big[s[i]]]]].altname,/remove), $
      strcompress(rc3[mindx[match[big[s[i]]]]].pgc,/remove), $
      mdist[match[big[s[i]]]]*3600.0, mindx[match[big[s[i]]]], format='(I5,A25,A15,A25,A25,F10.3,I10)'

; now update TABLE with RC3 quantities.  note, some objects have
; sizes in the RC3 but no sizes in NED, for whatever ungodly reason.
    
;   niceprint, table[match].galaxy, rc3[mindx[match]].name, rc3[mindx[match]].altname
;   niceprint, table[match].galaxy, table[match].d25_maj, rc3[mindx[match]].d25_maj

    doit = where(rc3[mindx[match]].d25_maj gt -900.0,ndoit)
    if ndoit ne 0L then table[match[doit]].d25_maj = rc3[mindx[match[doit]]].d25_maj ; major axis diameter [arcmin]
    doit = where(rc3[mindx[match]].d25_min gt -900.0,ndoit)
    if ndoit ne 0L then table[match[doit]].d25_min = rc3[mindx[match[doit]]].d25_min ; minor axis diameter [arcmin]

    table[match].pa = rc3[mindx[match]].pa ; position angle [degrees]
    table[match].rc3_t = rc3[mindx[match]].T
;;  table[match].rc3_type = rc3[mindx[match]].type <-- NOTE!!
    table[match].inclination = rc3[mindx[match]].inclination

; ---------------------------------------------------------------------------
; check the NED results for additional size measurements
; ---------------------------------------------------------------------------

    nodmaj = where(table.D25_maj lt -900.0,nnodmaj)
    if (nnodmaj ne 0L) then $
      good = where((basicdata[nodmaj].dmaj gt -900.0) and (basicdata[nodmaj].dmin gt -900.0),ngood)
    if ngood ne 0L then begin
       table[nodmaj[good]].D25_maj = basicdata[nodmaj[good]].dmaj
       table[nodmaj[good]].D25_min = basicdata[nodmaj[good]].dmin
    endif
    
; ---------------------------------------------------------------------------
; correct the UV/optical photometric data for Galactic extinction.  we
; use the Cardelli, Clayton, & Mathis (1989) extinction curve and the
; O'Donnell (1994) coefficients with the UBV effective filter
; wavelengths tabulated in Table 9 of Fukugita et al (1995)
; ---------------------------------------------------------------------------

    good = where(table.U gt -900,ngood)
    if ngood ne 0L then table[good].U = table[good].U - table[good].EBV_MW*k_lambda(3652.0,/odonnell)
    good = where(table.B gt -900,ngood)
    if ngood ne 0L then table[good].B = table[good].B - table[good].EBV_MW*k_lambda(4448.0,/odonnell)
    good = where(table.V gt -900,ngood)
    if ngood ne 0L then table[good].V = table[good].V - table[good].EBV_MW*k_lambda(5505.0,/odonnell)

; ---------------------------------------------------------------------------
; include additional galaxy properties from the literature here
; ---------------------------------------------------------------------------

    

; ---------------------------------------------------------------------------
; for galaxies with photometric magnitudes but no errors, assign an
; error of 0.2 mag (there should not be any that satisfy the criteria
; below) 
; ---------------------------------------------------------------------------

    noerr = where((table.U gt -900) and (table.U_err lt -900),no)
    if no ne 0L then table[noerr].U_err = 0.2
    noerr = where((table.B gt -900) and (table.B_err lt -900),no)
    if no ne 0L then table[noerr].B_err = 0.2
    noerr = where((table.V gt -900) and (table.V_err lt -900),no)
    if no ne 0L then table[noerr].V_err = 0.2
    
; ---------------------------------------------------------------------------
; compute the total IR flux and log luminosity using Eric Bell's/Karl
; Gordon's CALC_IR() routine 
; ---------------------------------------------------------------------------

;    good = where((table.iras_12 gt -900.0) and $
;                 (table.iras_25 gt -900.0) and $
;                 (table.iras_60 gt -900.0) and $
;                 (table.iras_100 gt -900.0),ngood)
;    irasflux = transpose([[table[good].iras_12],[table[good].iras_25],$
;                          [table[good].iras_60],[table[good].iras_100]])
;    for j = 0L, ngood-1L do table[good[j]].ir_flux = calc_ir(reform(irasflux[*,j]))
;
;    mpc2cm = 3.086D24 ; cm per Mpc
;    good = where((table.distance gt -900.0) and (table.ir_flux gt -900.0),ngood)
;    if ngood ne 0L then table[good].ir_lum = $
;      alog10(table[good].ir_flux*4.0*!dpi*(table[good].distance*mpc2cm)^2.0/lsun) ; log L_sun

; ---------------------------------------------------------------------------
; compute the FIR flux and luminosity using the Helou et al. (1988)
; formula, which only needs 60 and 100 micron IRAS fluxes.  estimate
; the Calzetti et al. (2000) estimate: L(IR) ~ 1.75 L(FIR).  also
; estimate the infrared star-formation rate from Kewley et al. (2002)
; equation (5)
; ---------------------------------------------------------------------------

    good = where((table.iras_60 gt 0.0) and (table.iras_100 gt 0.0),ngood)
    table[good].fir_flux = 1.26D-14*(2.58*table[good].iras_60+table[good].iras_100)

    good = where((table.distance gt -900.0) and (table.fir_flux gt -900.0),ngood)

    if ngood ne 0L then begin

; compute the IR and FIR luminosity and error (only distance errors
; included) 

       fir_lum = table[good].fir_flux*4.0*!dpi*(table[good].distance*mpc2m)^2.0 ; [Watts]
       fir_lum_err = 4.0*!dpi*table[good].fir_flux*$
         (2*table[good].distance*mpc2m*table[good].distance_err*mpc2m) ; [Watts]

       ir_lum = 1.75*fir_lum         ; [Watts]
       ir_lum_err = 1.75*fir_lum_err ; [Watts]

; take the log of the quantities
       
       table[good].fir_lum_err = fir_lum_err/fir_lum/alog(10.0) ; log L_sun
       table[good].fir_lum = alog10(fir_lum/lsun)               ; log L_sun

       table[good].ir_lum_err = ir_lum_err/ir_lum/alog(10.0) ; log L_sun
       table[good].ir_lum = alog10(ir_lum/lsun)              ; log L_sun
       
; compute the SFR
       
       table[good].ir_sfr = float(4.5D-44*(10.0^table[good].ir_lum*lsun*ergs2watts)) ; [M_sun/yr]

    endif
    
; compute the absolute B-band magnitude and error

    good = where((table.distance gt -900.0) and (table.B gt -900.0),ngood)
    if ngood ne 0L then begin

       table[good].M_B = table[good].B - 5.0*alog10(table[good].distance) - 25.0
       table[good].M_B_err = sqrt(table[good].B_err^2.0 + $
         ( (5.0/alog(10))*table[good].distance_err/table[good].distance )^2.0)

    endif
       
; compute the optical luminosity and error (normalized to B)

    good = where(table.M_B gt -900.0,ngood)
    if ngood ne 0L then begin

       B_lum = 10.0D^(0.4*msun) * 10.0^(-0.4*table[good].M_B)      ; [L_sun]
       B_lum_err = alog(10.0) * 0.4D * B_lum * table[good].M_B_err ; [L_sun]
       
; take the log of the quantities       
       
       table[good].B_lum_err = B_lum_err/B_lum/alog(10.0) ; log L_sun
       table[good].B_lum = alog10(B_lum)                  ; log L_sun

; bolometric       
       
       B_lum_bol = 10.0D^(0.4*msun_bol) * 10.0^(-0.4*table[good].M_B)      ; [L_sun]
       B_lum_bol_err = alog(10.0) * 0.4D * B_lum_bol * table[good].M_B_err ; [L_sun]
       
       table[good].B_lum_bol_err = B_lum_bol_err/B_lum_bol/alog(10.0) ; log L_sun
       table[good].B_lum_bol = alog10(B_lum_bol)                      ; log L_sun

    endif

; compute the K-band luminosity and error

    good = where((table.distance gt -900.0) and (table.Ks gt -900.0) and $
      (table.Ks_err gt -900),ngood)
    if ngood ne 0L then begin

       M_K = table[good].Ks - 5.0*alog10(table[good].distance) - 25.0
       M_K_err = sqrt(table[good].Ks_err^2.0 + $
         ( (5.0/alog(10))*table[good].distance_err/table[good].distance )^2.0)

       Ks_lum = 10.0D^(0.4*msun) * 10.0^(-0.4*M_K)     ; [L_sun]
       Ks_lum_err = alog(10) * 0.4D * Ks_lum * M_K_err ; [L_sun]
       
; take the log of the quantities       

       table[good].Ks_lum_err = Ks_lum_err/Ks_lum/alog(10.0) ; log L_sun
       table[good].Ks_lum = alog10(Ks_lum)                   ; log L_sun

    endif
    
; optical/FIR and optical/IR ratios and errors

    good = where((table.B_lum_bol gt -900.0) and (table.fir_lum gt -900) and $
      (table.B_lum_bol_err gt -900.0) and (table.fir_lum_err gt -900.0),ngood)
    if ngood ne 0L then begin
       table[good].fir_opt_ratio = 10^(table[good].fir_lum-table[good].B_lum_bol)
       x1 = 10^table[good].fir_lum & x1err = table[good].fir_lum_err*alog(10.0)*x1
       y1 = 10^table[good].B_lum_bol   & y1err = table[good].B_lum_bol_err*alog(10.0)*y1
       table[good].fir_opt_ratio_err = im_compute_error(x1,x1err,y1,y1err,/quotient)
    endif

    good = where((table.B_lum_bol gt -900.0) and (table.ir_lum gt -900) and $
      (table.B_lum_bol_err gt -900.0) and (table.ir_lum_err gt -900.0),ngood)
    if ngood ne 0L then begin
       table[good].ir_opt_ratio = 10^(table[good].ir_lum-table[good].B_lum_bol)
       x1 = 10^table[good].ir_lum & x1err = table[good].ir_lum_err*alog(10.0)*x1
       y1 = 10^table[good].B_lum_bol   & y1err = table[good].B_lum_bol_err*alog(10.0)*y1
       table[good].ir_opt_ratio_err = im_compute_error(x1,x1err,y1,y1err,/quotient)
    endif

; ---------------------------------------------------------------------------
; sort by RA
; ---------------------------------------------------------------------------

    srtra = sort(im_hms2dec(table.ra))
    table = table[srtra]

; ---------------------------------------------------------------------------
; write out the full atlas table as a binary fits file if requested
; ---------------------------------------------------------------------------

    if keyword_set(write) then begin
    
       splog, 'Writing '+outpath+'kenn92_data.fits.gz'
       mwrfits, table, outpath+'kenn92_data.fits', /create
       spawn, ['gzip -f '+outpath+'kenn92_data.fits'], /sh

    endif

return
end
