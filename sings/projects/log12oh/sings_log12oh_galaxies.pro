pro sings_log12oh_galaxies, result, clobber=clobber
; jm10mar10ucsd - derive the mean abundances for the SINGS galaxies
; based on our nuclear, circumnuclear, and radial-strip spectra

; read the data    
    version = sings_log12oh_version()
    outpath = sings_path(/projects)+'log12oh/'
    singspath = sings_path(/analysis)

    outfile = 'sings_log12oh_'+version+'.fits'
    if file_test(outfile+'.gz') and (keyword_set(clobber) eq 0) then begin
       splog, 'Output file '+outfile+' exists; use /CLOBBER'
       return
    endif

    sings = sings_read_info()
    nuclear = read_sings_log12oh_samples(/nodust_nuclear)
    drift20 = read_sings_log12oh_samples(/nodust_drift20)
    drift56 = read_sings_log12oh_samples(/nodust_drift56)
    ngal = n_elements(sings)

    branchinfo = mrdfits(outpath+'sings_r23_branch_'+version+'.fits.gz',1,/silent)
    classinfo = mrdfits(outpath+'sings_class_'+version+'.fits',1,/silent)

; initialize the output data structure    
    result = {$
      sings_id:                                 0,$
      galaxy:                                  '',$
      ned_galaxy:                              '',$
      nice_galaxy:                             '',$
      ra:                                      '',$
      dec:                                     '',$
      ebv_mw:                              -999.0,$
      bvri:                         fltarr(4)-999,$
      bvri_err:                     fltarr(4)-999,$
      bvri_absmag:                  fltarr(4)-999,$
      bvri_absmag_err:              fltarr(4)-999,$
      bvri_ref:                                '',$
      mb:                                  -999.0,$
      mb_err:                              -999.0,$
      bv:                                  -999.0,$
      bv_err:                              -999.0,$
      
      distance:                            -999.0,$
      distance_err:                        -999.0,$
      distance_texref:                         '',$
      distance_method:                         '',$
      type:                                    '',$
      t:                                      0.0,$
      d25:                                 -999.0,$
      r25:                                 -999.0,$
      incl:                                -999.0,$
      pa:                                  -999.0,$

      class:                            '\nodata',$ ; final spectral class
      class_remarks:                    '\nodata',$
      ho_class:                         '\nodata',$ ; spectral class from Ho et al. 1997
      r23_branch:                             '?',$
      
      nuclear:                                  0, $ ; Boolean flag
      nuclear_aperture_kpc:                -999.0, $ ; physical aperture size
      nuclear_aperture_kpc_tex:         '\nodata', $ ; physical aperture size
      nuclear_rr25:                        -999.0, $ ; approximate characteristic radial position
      nuclear_fraction_b:         [-999.0,-999.0], $ ; B-band light fraction
      nuclear_hahb:               [-999.0,-999.0], $
      nuclear_ebv:                [-999.0,-999.0], $
      nuclear_class:                          '?', $ ; the default is insufficient data
      nuclear_broad:                            0, $ ; broad-line AGN?
      nuclear_n2flag:                           0, $ ; classification based on [NII]/Ha only?
      nuclear_bpt_d:                       -999.0, $
      nuclear_bpt_phi:                     -999.0, $
      nuclear_ewhb:                        -999.0, $
;     nuclear_ewhb_abs:                    -999.0, $
      nuclear_r23branch_kk04:                 '?', $
      nuclear_r23branch_pt05:                 '?', $
      nuclear_log12oh_kk04:       [-999.0,-999.0], $
      nuclear_log12oh_pt05:       [-999.0,-999.0], $
      nuclear_log12oh_o3n2:       [-999.0,-999.0], $
      nuclear_log12oh_n2o2:       [-999.0,-999.0], $
      nuclear_logu_kk04:          [-999.0,-999.0], $
      nuclear_r23:                [-999.0,-999.0], $
      nuclear_o32:                [-999.0,-999.0], $
      nuclear_p:                  [-999.0,-999.0], $
      nuclear_n2o2:               [-999.0,-999.0], $

      drift20:                                  0, $ ; Boolean flag
      drift20_aperture_kpc:                -999.0, $ ; physical aperture size
      drift20_aperture_kpc_tex:         '\nodata', $ ; physical aperture size
      drift20_rr25:                        -999.0, $ ; approximate characteristic radial position
      drift20_fraction_b:         [-999.0,-999.0], $ ; B-band light fraction
      drift20_hahb:               [-999.0,-999.0], $
      drift20_ebv:                [-999.0,-999.0], $
      drift20_class:                          '?', $
      drift20_broad:                            0, $ ; broad-line AGN?
      drift20_n2flag:                           0, $ ; classification based on [NII]/Ha only?
      drift20_bpt_d:                       -999.0, $
      drift20_bpt_phi:                     -999.0, $
      drift20_ewhb:                        -999.0, $
;     drift20_ewhb_abs:                    -999.0, $
      drift20_r23branch_kk04:                 '?', $
      drift20_r23branch_pt05:                 '?', $
      drift20_log12oh_kk04:       [-999.0,-999.0], $
      drift20_log12oh_pt05:       [-999.0,-999.0], $
      drift20_log12oh_o3n2:       [-999.0,-999.0], $
      drift20_log12oh_n2o2:       [-999.0,-999.0], $
      drift20_logu_kk04:          [-999.0,-999.0], $
      drift20_r23:                [-999.0,-999.0], $
      drift20_o32:                [-999.0,-999.0], $
      drift20_p:                  [-999.0,-999.0], $
      drift20_n2o2:               [-999.0,-999.0], $

      drift56:                                  0, $ ; Boolean flag
      drift56_aperture_kpc:                -999.0, $ ; physical aperture size
      drift56_aperture_kpc_tex:         '\nodata', $ ; physical aperture size
      drift56_rr25:                        -999.0, $ ; approximate characteristic radial position
      drift56_fraction_b:         [-999.0,-999.0], $ ; B-band light fraction
      drift56_hahb:               [-999.0,-999.0], $
      drift56_ebv:                [-999.0,-999.0], $
      drift56_class:                          '?', $
      drift56_broad:                            0, $ ; broad-line AGN?
      drift56_n2flag:                           0, $ ; classification based on [NII]/Ha only?
      drift56_bpt_d:                       -999.0, $
      drift56_bpt_phi:                     -999.0, $
      drift56_ewhb:                        -999.0, $
;     drift56_ewhb_abs:                    -999.0, $
      drift56_r23branch_kk04:                 '?', $
      drift56_r23branch_pt05:                 '?', $
      drift56_log12oh_kk04:       [-999.0,-999.0], $
      drift56_log12oh_pt05:       [-999.0,-999.0], $
      drift56_log12oh_o3n2:       [-999.0,-999.0], $
      drift56_log12oh_n2o2:       [-999.0,-999.0], $
      drift56_logu_kk04:          [-999.0,-999.0], $
      drift56_r23:                [-999.0,-999.0], $
      drift56_o32:                [-999.0,-999.0], $
      drift56_p:                  [-999.0,-999.0], $
      drift56_n2o2:               [-999.0,-999.0]}
    result = replicate(result,ngal)
    
    result.sings_id    = sings.sings_id
    result.galaxy      = sings.galaxy
    result.ned_galaxy  = sings.ned_galaxy
    result.nice_galaxy = sings.nice_galaxy
    result.ra          = sings.ra
    result.dec         = sings.dec
    result.ebv_mw      = sings.ebv_mw
    result.type        = strtrim(sings.lit_type,2)
    result.t           = sings.lit_t
    result.d25         = sings.d25_maj
    result.r25         = sings.d25_maj/2.0

; distances    
    result.distance        = sings.distance
    result.distance_err    = sings.distance_err
    result.distance_texref = sings.distance_texref
    result.distance_method = sings.distance_method

; physical aperture, latex style, and the characteristic radial
; position of each spectrum; note that I have to use D25 because the
; spectroscopic apertures are *diameters*, but the resulting quantity
; is an R/R25
    these = where(sings.nuclear)
    factor = sings[these].distance*1E3*!dtor/3600.0
    result[these].nuclear_aperture_kpc_tex = '$'+$
      string(factor*sings[these].nuclear_ap,format='(F5.2)')+'\times'+$
      string(factor*sings[these].nuclear_scan,format='(F5.2)')+'$'
    result[these].nuclear_aperture_kpc = factor^2.0*$
      sings[these].nuclear_ap*sings[these].nuclear_scan
    result[these].nuclear_rr25 = sqrt(sings[these].nuclear_ap*$
      sings[these].nuclear_scan)/result[these].d25/60.0 

    these = where(sings.drift20)
    factor = sings[these].distance*1E3*!dtor/3600.0
    result[these].drift20_aperture_kpc_tex = '$'+$
      string(factor*sings[these].drift20_ap,format='(F5.2)')+'\times'+$
      string(factor*sings[these].drift20_scan,format='(F5.2)')+'$'
    result[these].drift20_aperture_kpc = factor^2.0*$
      sings[these].drift20_ap*sings[these].drift20_scan
    result[these].drift20_rr25 = sqrt(sings[these].drift20_ap*$
      sings[these].drift20_scan)/result[these].d25/60.0 
    
    these = where(sings.drift56)
    factor = sings[these].distance*1E3*!dtor/3600.0
    result[these].drift56_aperture_kpc_tex = '$'+$
      string(factor*sings[these].drift56_ap,format='(F5.2)')+'\times'+$
      string(factor*sings[these].drift56_scan,format='(F5.2)')+'$'
    result[these].drift56_aperture_kpc = factor^2.0*$
      sings[these].drift56_ap*sings[these].drift56_scan
    result[these].drift56_rr25 = sqrt(sings[these].drift56_ap*$
      sings[these].drift56_scan)/result[these].d25/60.0 
    
; inclination and position angles: priority: RC3, then 2MASS; set the
; inclination angle of NGC5194 to be 20, the kinematic value (Tully
; 1974); also set the inclination and position angles according to Lee
; et al. (2006)
    result.incl = sings.inclination
    result.pa   = sings.posangle

    m51 = where(strmatch(sings.galaxy,'*5194*'),nm51)
    if (nm51 ne 0) then result[m51].incl = 20.0

    n6822 = where(strmatch(sings.galaxy,'*6822*'),nn6822)
    if (nn6822 ne 0) then begin
       result[n6822].incl = 50.1
       result[n6822].pa = 122.0
    endif
    
; photometry; start with Munoz-Mateos et al. 2009 (with some
; exceptions), and then use Dale et al. 2007; convert to Vega at the
; end 
;   result.bvri         = sings.bvri ; [Vega]
;   result.bvri_err     = sings.bvri_err

; Munoz-Mateos+09    
    table5 = rsex(getenv('CATALOGS_DIR')+'/09munoz/table5.sex')
    table6 = rsex(getenv('CATALOGS_DIR')+'/09munoz/table6.sex')

; first, use the recalibrated BVRI imaging    
    bvri = transpose([[table5.b],[table5.v],[table5.r],[table5.i]]) ; AB
    bvri_err = transpose([[table5.b_err],[table5.v_err],[table5.r_err],[table5.i_err]])
    match, strtrim(strlowcase(result.galaxy),2), strtrim(strlowcase(table5.galaxy),2), m1, m2
    
    result[m1].bvri = bvri[*,m2] ; AB
    result[m1].bvri_err = bvri_err[*,m2]
    result[m1].bvri_ref = 'munoz09a-ubvri'
;   niceprint, result[m1].galaxy, table5[m2].galaxy, result[m1].bvri[0], table5[m2].b

;   plot, result[m1].bvri[0]-im_d2dmod(result[m1].distance), $
;     result[m1].bvri[0]-result[m1].bvri[1], psym=6, $
;     xrange=[-12,-25], yrange=[-0.2,1.2], xsty=3, ysty=3
    
;; NGC2915 is special because it has V-band but no B-band photometry;
;; we want accurate colors, so zero out the Munoz-Mateos value and use
;; Dale+07
;    toss = where(strtrim(result.galaxy,2) eq 'NGC2915')
;    result[toss].bvri = 
    
; next, transform the Munoz-Mateos+09 SDSS/ugriz magnitudes to BVRI;
; the Lupton+05 transformations can be found here:
; http://www.sdss.org/dr5/algorithms/sdssUBVRITransform.html#Lupton2005
; but we will use the galaxy transformations in Blanton & Roweis+07

; for reference, here are the Chonis+08 equations for stars
;    b = table6.g + 0.327*(table6.g-table6.r) + 0.216
;    v = table6.g - 0.587*(table6.g-table6.r) - 0.011
;    r = table6.r - 0.272*(table6.g-table6.r) - 0.159
;    i = table6.i - 0.337*(table6.r-table6.i) - 0.370

; note the errors are negligible compared to the zeropoint error,
; which is 10%-15%; see Blanton & Roweis+07
    b = table6.g + 0.3915*(table6.g-table6.r-0.6102) + 0.2354 ; AB
    v = table6.g - 0.7585*(table6.g-table6.r-0.6102) - 0.3516 ; AB
    berr = sqrt(table6.g_err^2 + 0.3915*(table6.g_err^2+table6.r_err^2))
    verr = sqrt(table6.g_err^2 + 0.7585*(table6.g_err^2+table6.r_err^2))
;   b = table6.g + 0.3130*(table6.g-table6.r) + 0.2271 
;   v = table6.g - 0.5784*(table6.g-table6.r) - 0.0038
;   berr = sqrt(table6.g_err^2 + 0.3130*(table6.g_err^2+table6.r_err^2))
;   verr = sqrt(table6.g_err^2 + 0.5784*(table6.g_err^2+table6.r_err^2))
    
    need = where(result.bvri[0] lt -900.0)
    match, strtrim(strlowcase(result[need].galaxy),2), $
      strtrim(strlowcase(table6.galaxy),2), m1, m2
;   niceprint, result[need[m1]].galaxy, table6[m2].galaxy
    
    result[need[m1]].bvri[0:1] = transpose([[b[m2]],[v[m2]]]) ; AB
    result[need[m1]].bvri_err[0:1] = transpose([[berr[m2]],[verr[m2]]])
    result[need[m1]].bvri_ref = 'munoz09a-ugriz'

;   good = where(result.bvri[0] gt 0.0)
;   plot, result[good].bvri[0]-im_d2dmod(result[good].distance), $
;     result[good].bvri[0]-result[good].bvri[1], psym=6, $
;     xrange=[-12,-25], yrange=[-0.2,1.2], xsty=3, ysty=3
    
; finally, use Dale+07 for everything else
    flux = rsex(singspath+'sings_photometry_2006dale.sex') ; Dale+07
    ferr = rsex(singspath+'sings_photometry_2006dale.uncertainty.sex')
    dale_bvri = transpose([[flux.b],[flux.v],[flux.r],[flux.i]])
    dale_bvri_err = transpose([[ferr.b],[ferr.v],[ferr.r],[ferr.i]])
    dale_galaxy = repstr(repstr(repstr(repstr(repstr(flux.galaxy,'HoIX','HOLMBERGIX'),$
      'HoII','HOLMBERGII'),'HoI','HOLMBERGI'),'Mrk33','MRK0033'),'Tol89','TOLOLO89')
    for ii = 0, 3 do begin
       good = where(dale_bvri[ii,*] gt 0)
       dale_bvri[ii,good] = -2.5*alog10(dale_bvri[ii,good]*1D-23)-48.6 ; AB
       dale_bvri_err[ii,good] = 2.5*dale_bvri_err[ii,good]/dale_bvri[ii,good]/alog(10.0)
    endfor

    need = where(result.bvri[0] lt -900.0)
    match, strtrim(strlowcase(result[need].galaxy),2), strtrim(strlowcase(dale_galaxy),2), m1, m2
    result[need[m1]].bvri = dale_bvri[*,m2] ; AB
    result[need[m1]].bvri_err = dale_bvri_err[*,m2]
    result[need[m1]].bvri_ref = 'dale07a'

; convert to Vega    
    v2ab = k_vega2ab(filterlist='bessell_'+['B','V','R','I']+'.par',/kurucz,/silent)
    for ii = 0, 3 do begin
       good = where((result.bvri[ii] gt 0.0),ngood)
       if (ngood ne 0) then result[good].bvri[ii] = result[good].bvri[ii] - v2ab[ii] ; AB-->Vega
    endfor

;; check to make sure that the Dale BV photometry agrees with the RC3    
;    plot, result[need[m1]].bvri[0], sings[need[m1]].rc3_ubv[1], $
;      psym=6, xrange=[6,15], yrange=[6,15], xsty=3, ysty=3
;    plot, result[need[m1]].bvri[1], sings[need[m1]].rc3_ubv[2], $
;      psym=6, xrange=[6,15], yrange=[6,15], xsty=3, ysty=3
;    djs_oplot, !x.crange, !y.crange
;    gg = where(result[need[m1]].bvri[0] gt -900 and sings[need[m1]].rc3_ubv[1] gt -900)
;    ss = im_stats(result[need[m1[gg]]].bvri[0]-sings[need[m1[gg]]].rc3_ubv[1],/ver,sigrej=3.0)
;    gg = where(result[need[m1]].bvri[1] gt -900 and sings[need[m1]].rc3_ubv[2] gt -900)
;    ss = im_stats(result[need[m1[gg]]].bvri[1]-sings[need[m1[gg]]].rc3_ubv[2],/ver,sigrej=3.0)
    
; the Dale+07 calibration for NGC6822 is off (too faint) by ~0.37 mag
; compared to the RC3 or Karachentsev+04 (see the email thread with
; Danny); however, we still want the B-V color (there's no V-band in
; the RC3), so simply apply a zeropoint offset to Danny's photometry;
; the 9.32 magnitude comes from Lee+06 and Karachentsev+04
    n6822 = where(strtrim(result.galaxy,2) eq 'NGC6822')
    offset = 9.32-result[n6822].ebv_mw*k_lambda(4400.0,/odon)-result[n6822].bvri[0]
    splog, 'NGC6822 offset ', offset
    result[n6822].bvri = result[n6822].bvri + offset

; for convenience compute the B-V color; include a minimum 5% error
; (the statistical errors are tiny!) 
    minerr = 0.05
    result.bv = result.bvri[0]-result.bvri[1]
    result.bv_err = sqrt(result.bvri_err[0]^2+result.bvri_err[1]^2+minerr^2)

; compute absolute magnitudes; include a 10% minimum photometric floor
; on the observed photometry; the absolute magnitude error includes
; the floor plus the error on the distance modulus
    minerr = 0.1
    for ii = 0, 3 do begin
       good = where((result.distance gt -900.0) and (result.bvri[ii] gt -900.0),ngood)
       if (ngood ne 0) then begin
          dmod = im_d2dmod(result[good].distance,err_dmod=err_dmod,$
            err_dist=result[good].distance_err)
          result[good].bvri_err[ii] = sqrt(result[good].bvri_err[ii]^2 + minerr^2)
          result[good].bvri_absmag[ii] = result[good].bvri[ii] - dmod
          result[good].bvri_absmag_err[ii] = sqrt(result[good].bvri_err[ii]^2.0 + err_dmod^2.0)
       endif
    endfor
    result.mb = result.bvri_absmag[0]
    result.mb_err = result.bvri_absmag_err[0]

; classifications    
    result.class         = classinfo.class
    result.class_remarks = classinfo.class_remarks
    result.ho_class      = classinfo.ho_class

    result.nuclear_class   = classinfo.nuclear_class
    result.nuclear_broad   = classinfo.nuclear_broad
    result.nuclear_n2flag  = classinfo.nuclear_n2flag
    result.nuclear_bpt_d   = classinfo.nuclear_bpt_d
    result.nuclear_bpt_phi = classinfo.nuclear_bpt_phi

    result.drift20_class   = classinfo.drift20_class
    result.drift20_broad   = classinfo.drift20_broad
    result.drift20_n2flag  = classinfo.drift20_n2flag
    result.drift20_bpt_d   = classinfo.drift20_bpt_d
    result.drift20_bpt_phi = classinfo.drift20_bpt_phi

    result.drift56_class   = classinfo.drift56_class
    result.drift56_broad   = classinfo.drift56_broad
    result.drift56_n2flag  = classinfo.drift56_n2flag
    result.drift56_bpt_d   = classinfo.drift56_bpt_d
    result.drift56_bpt_phi = classinfo.drift56_bpt_phi

    result.r23_branch = branchinfo.r23_branch

; light-fractions
    result[where(sings.nuclear)].nuclear_fraction_b = $
      transpose([[sings[where(sings.nuclear)].nuclear_fraction_b],$
      [sings[where(sings.nuclear)].nuclear_fraction_b_err]])
    result[where(sings.drift20)].drift20_fraction_b = $
      transpose([[sings[where(sings.drift20)].drift20_fraction_b],$
      [sings[where(sings.drift20)].drift20_fraction_b_err]])
    result[where(sings.drift56)].drift56_fraction_b = $
      transpose([[sings[where(sings.drift56)].drift56_fraction_b],$
      [sings[where(sings.drift56)].drift56_fraction_b_err]])

; store the nuclear, circumnuclear, and radial-strip abundances

    psfile = repstr(outfile,'.fits','.ps') ; R23branch assignment QAplot
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.5
    
; nuclear    
    spherematch, 15.0*im_hms2dec(result.ra), im_hms2dec(result.dec), 15.0*im_hms2dec(nuclear.ra), $
      im_hms2dec(nuclear.dec), 1.0/3600.0, result_match, nuclear_match, maxmatch=1
    if (n_elements(nuclear_match) ne n_elements(nuclear)) then message, 'We have a problem, Dave.'
    result[result_match].nuclear = 1
    result[result_match].nuclear_ewhb     = nuclear[nuclear_match].h_beta_ew[0]
;   result[result_match].nuclear_ewhb_abs = nuclear[nuclear_match].babs_h_beta_ew[0]
    result[result_match].nuclear_hahb     = transpose([[nuclear[nuclear_match].hahb],[nuclear[nuclear_match].hahb_err]])
    result[result_match].nuclear_ebv      = transpose([[nuclear[nuclear_match].ebv_hahb],[nuclear[nuclear_match].ebv_hahb_err]])
    result[result_match].nuclear_r23      = transpose([[nuclear[nuclear_match].zstrong_r23],[nuclear[nuclear_match].zstrong_r23_err]])
    result[result_match].nuclear_o32      = transpose([[nuclear[nuclear_match].zstrong_o32],[nuclear[nuclear_match].zstrong_o32_err]])
    result[result_match].nuclear_p        = transpose([[nuclear[nuclear_match].zstrong_p],[nuclear[nuclear_match].zstrong_p_err]])
    result[result_match].nuclear_log12oh_o3n2 = transpose([[nuclear[nuclear_match].zstrong_12oh_oiiinii_pettini],$
      [nuclear[nuclear_match].zstrong_12oh_oiiinii_pettini_err]])

    result[result_match].nuclear_n2o2     = transpose([[nuclear[nuclear_match].zstrong_niioii],[nuclear[nuclear_match].zstrong_niioii_err]])
    result[result_match].nuclear_log12oh_n2o2 = transpose([[nuclear[nuclear_match].zstrong_12oh_kd02_niioii],$
      [nuclear[nuclear_match].zstrong_12oh_kd02_niioii_err]])

    splog, 'NUCLEAR abundances:'
    splog, '-------------------'
    splog, ' KK04:'
    kk04_nuclear_branch = sings_assign_r23branch(nuclear[nuclear_match],silent=0,$
      r23branch=result[result_match].r23_branch,/justflux,/kk04,/debug,$
      title='Nuclear/KK04')
    kk04_nuclear_good = where((strmatch(kk04_nuclear_branch.r23branch_kk04,'*rejected*',/fold) eq 0B) and $
      (kk04_nuclear_branch.zstrong_12oh_kk04 gt -900.0),nkk04_nuclear_good)
    if (nkk04_nuclear_good ne 0L) then begin
       result[result_match[kk04_nuclear_good]].nuclear_r23branch_kk04  = strtrim(kk04_nuclear_branch[kk04_nuclear_good].r23branch_kk04,2)
       result[result_match[kk04_nuclear_good]].nuclear_log12oh_kk04[0] = kk04_nuclear_branch[kk04_nuclear_good].zstrong_12oh_kk04
       result[result_match[kk04_nuclear_good]].nuclear_log12oh_kk04[1] = kk04_nuclear_branch[kk04_nuclear_good].zstrong_12oh_kk04_err
       result[result_match[kk04_nuclear_good]].nuclear_logu_kk04[0]    = kk04_nuclear_branch[kk04_nuclear_good].zstrong_logu_kk04
       result[result_match[kk04_nuclear_good]].nuclear_logu_kk04[1]    = kk04_nuclear_branch[kk04_nuclear_good].zstrong_logu_kk04_err
    endif
    splog, ' PT05:'
    pt05_nuclear_branch = sings_assign_r23branch(nuclear[nuclear_match],silent=0,$
      r23branch=result[result_match].r23_branch,/justflux,/pt05,/debug,$
      title='Nuclear/PT05')
    pt05_nuclear_good = where((strmatch(pt05_nuclear_branch.r23branch_pt05,'*rejected*',/fold) eq 0B) and $
      (pt05_nuclear_branch.zstrong_12oh_pt05 gt -900.0),npt05_nuclear_good)
    if (npt05_nuclear_good ne 0L) then begin
       result[result_match[pt05_nuclear_good]].nuclear_r23branch_pt05  = strtrim(pt05_nuclear_branch[pt05_nuclear_good].r23branch_pt05,2)
       result[result_match[pt05_nuclear_good]].nuclear_log12oh_pt05[0] = pt05_nuclear_branch[pt05_nuclear_good].zstrong_12oh_pt05
       result[result_match[pt05_nuclear_good]].nuclear_log12oh_pt05[1] = pt05_nuclear_branch[pt05_nuclear_good].zstrong_12oh_pt05_err
    endif

; drift20    
    
    spherematch, 15.0*im_hms2dec(result.ra), im_hms2dec(result.dec), 15.0*im_hms2dec(drift20.ra), $
      im_hms2dec(drift20.dec), 1.0/3600.0, result_match, drift20_match, maxmatch=1
    if (n_elements(drift20_match) ne n_elements(drift20)) then message, 'We have a problem, Dave.'
    result[result_match].drift20 = 1
    result[result_match].drift20_ewhb     = drift20[drift20_match].h_beta_ew[0]
;   result[result_match].drift20_ewhb_abs = drift20[drift20_match].babs_h_beta_ew[0]
    result[result_match].drift20_hahb     = transpose([[drift20[drift20_match].hahb],[drift20[drift20_match].hahb_err]])
    result[result_match].drift20_ebv      = transpose([[drift20[drift20_match].ebv_hahb],[drift20[drift20_match].ebv_hahb_err]])
    result[result_match].drift20_r23      = transpose([[drift20[drift20_match].zstrong_r23],[drift20[drift20_match].zstrong_r23_err]])
    result[result_match].drift20_o32      = transpose([[drift20[drift20_match].zstrong_o32],[drift20[drift20_match].zstrong_o32_err]])
    result[result_match].drift20_p        = transpose([[drift20[drift20_match].zstrong_p],[drift20[drift20_match].zstrong_p_err]])
    result[result_match].drift20_n2o2     = transpose([[drift20[drift20_match].zstrong_niioii],[drift20[drift20_match].zstrong_niioii_err]])
    result[result_match].drift20_log12oh_o3n2 = transpose([[drift20[drift20_match].zstrong_12oh_oiiinii_pettini],$
      [drift20[drift20_match].zstrong_12oh_oiiinii_pettini_err]])
    result[result_match].drift20_log12oh_n2o2 = transpose([[drift20[drift20_match].zstrong_12oh_kd02_niioii],$
      [drift20[drift20_match].zstrong_12oh_kd02_niioii_err]])

    splog, 'DRIFT20 abundances:'
    splog, '-------------------'
    splog, ' KK04:'
    kk04_drift20_branch = sings_assign_r23branch(drift20[drift20_match],silent=0,$
      r23branch=result[result_match].r23_branch,/justflux,/kk04,/debug,$
      title='Circumnuclear/KK04')
    kk04_drift20_good = where((strmatch(kk04_drift20_branch.r23branch_kk04,'*rejected*',/fold) eq 0B) and $
      (kk04_drift20_branch.zstrong_12oh_kk04 gt -900.0),nkk04_drift20_good)
    if (nkk04_drift20_good ne 0L) then begin
       result[result_match[kk04_drift20_good]].drift20_r23branch_kk04  = strtrim(kk04_drift20_branch[kk04_drift20_good].r23branch_kk04,2)
       result[result_match[kk04_drift20_good]].drift20_log12oh_kk04[0] = kk04_drift20_branch[kk04_drift20_good].zstrong_12oh_kk04
       result[result_match[kk04_drift20_good]].drift20_log12oh_kk04[1] = kk04_drift20_branch[kk04_drift20_good].zstrong_12oh_kk04_err
       result[result_match[kk04_drift20_good]].drift20_logu_kk04[0]    = kk04_drift20_branch[kk04_drift20_good].zstrong_logu_kk04
       result[result_match[kk04_drift20_good]].drift20_logu_kk04[1]    = kk04_drift20_branch[kk04_drift20_good].zstrong_logu_kk04_err
    endif
    splog, ' PT05:'
    pt05_drift20_branch = sings_assign_r23branch(drift20[drift20_match],silent=0,$
      r23branch=result[result_match].r23_branch,/justflux,/pt05,/debug,$
      title='Circumnuclear/PT05')
    pt05_drift20_good = where((strmatch(pt05_drift20_branch.r23branch_pt05,'*rejected*',/fold) eq 0B) and $
      (pt05_drift20_branch.zstrong_12oh_pt05 gt -900.0),npt05_drift20_good)

    if (npt05_drift20_good ne 0L) then begin
       result[result_match[pt05_drift20_good]].drift20_r23branch_pt05  = strtrim(pt05_drift20_branch[pt05_drift20_good].r23branch_pt05,2)
       result[result_match[pt05_drift20_good]].drift20_log12oh_pt05[0] = pt05_drift20_branch[pt05_drift20_good].zstrong_12oh_pt05
       result[result_match[pt05_drift20_good]].drift20_log12oh_pt05[1] = pt05_drift20_branch[pt05_drift20_good].zstrong_12oh_pt05_err
    endif
stop

; drift56    
    spherematch, 15.0*im_hms2dec(result.ra), im_hms2dec(result.dec), 15.0*im_hms2dec(drift56.ra), $
      im_hms2dec(drift56.dec), 1.0/3600.0, result_match, drift56_match, maxmatch=1
    if (n_elements(drift56_match) ne n_elements(drift56)) then message, 'We have a problem, Dave.'
    result[result_match].drift56 = 1
    result[result_match].drift56_ewhb     = drift56[drift56_match].h_beta_ew[0]
;   result[result_match].drift56_ewhb_abs = drift56[drift56_match].babs_h_beta_ew[0]
    result[result_match].drift56_hahb     = transpose([[drift56[drift56_match].hahb],[drift56[drift56_match].hahb_err]])
    result[result_match].drift56_ebv      = transpose([[drift56[drift56_match].ebv_hahb],[drift56[drift56_match].ebv_hahb_err]])
    result[result_match].drift56_r23      = transpose([[drift56[drift56_match].zstrong_r23],[drift56[drift56_match].zstrong_r23_err]])
    result[result_match].drift56_o32      = transpose([[drift56[drift56_match].zstrong_o32],[drift56[drift56_match].zstrong_o32_err]])
    result[result_match].drift56_p        = transpose([[drift56[drift56_match].zstrong_p],[drift56[drift56_match].zstrong_p_err]])
    result[result_match].drift56_n2o2     = transpose([[drift56[drift56_match].zstrong_niioii],[drift56[drift56_match].zstrong_niioii_err]])
    result[result_match].drift56_log12oh_o3n2 = transpose([[drift56[drift56_match].zstrong_12oh_oiiinii_pettini],$
      [drift56[drift56_match].zstrong_12oh_oiiinii_pettini_err]])
    result[result_match].drift56_log12oh_n2o2 = transpose([[drift56[drift56_match].zstrong_12oh_kd02_niioii],$
      [drift56[drift56_match].zstrong_12oh_kd02_niioii_err]])

    splog, 'DRIFT56 abundances:'
    splog, '-------------------'
    splog, ' KK04:'
    kk04_drift56_branch = sings_assign_r23branch(drift56[drift56_match],silent=0,$
      r23branch=result[result_match].r23_branch,/justflux,/kk04,/debug,$
      title='Radial-Strip/KK04')
    kk04_drift56_good = where((strmatch(kk04_drift56_branch.r23branch_kk04,'*rejected*',/fold) eq 0B) and $
      (kk04_drift56_branch.zstrong_12oh_kk04 gt -900.0),nkk04_drift56_good)
    if (nkk04_drift56_good ne 0L) then begin
       result[result_match[kk04_drift56_good]].drift56_r23branch_kk04  = strtrim(kk04_drift56_branch[kk04_drift56_good].r23branch_kk04,2)
       result[result_match[kk04_drift56_good]].drift56_log12oh_kk04[0] = kk04_drift56_branch[kk04_drift56_good].zstrong_12oh_kk04
       result[result_match[kk04_drift56_good]].drift56_log12oh_kk04[1] = kk04_drift56_branch[kk04_drift56_good].zstrong_12oh_kk04_err
       result[result_match[kk04_drift56_good]].drift56_logu_kk04[0]    = kk04_drift56_branch[kk04_drift56_good].zstrong_logu_kk04
       result[result_match[kk04_drift56_good]].drift56_logu_kk04[1]    = kk04_drift56_branch[kk04_drift56_good].zstrong_logu_kk04_err
    endif
    splog, ' PT05:'
    pt05_drift56_branch = sings_assign_r23branch(drift56[drift56_match],silent=0,$
      r23branch=result[result_match].r23_branch,/justflux,/pt05,/debug,$
      title='Radial-Strip/PT05')
    pt05_drift56_good = where((strmatch(pt05_drift56_branch.r23branch_pt05,'*rejected*',/fold) eq 0B) and $
      (pt05_drift56_branch.zstrong_12oh_pt05 gt -900.0),npt05_drift56_good)
    if (npt05_drift56_good ne 0L) then begin
       result[result_match[pt05_drift56_good]].drift56_r23branch_pt05  = strtrim(pt05_drift56_branch[pt05_drift56_good].r23branch_pt05,2)
       result[result_match[pt05_drift56_good]].drift56_log12oh_pt05[0] = pt05_drift56_branch[pt05_drift56_good].zstrong_12oh_pt05
       result[result_match[pt05_drift56_good]].drift56_log12oh_pt05[1] = pt05_drift56_branch[pt05_drift56_good].zstrong_12oh_pt05_err
    endif

    im_plotconfig, psfile=psfile, /psclose, /gzip

; write out
    im_mwrfits, result, outfile, /clobber

return
end
