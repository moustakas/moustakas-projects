pro build_sfrm_sample, anc, sdss=sdss, iauto=iauto
; jm09nov30ucsd - build the AGES/SFRM parent samples

    common sfrm_sample, phot, kcorr

;   jj = where((ww.spec_yesno ge 1) and (ww.z_yesno ge 1) and (kk.i_obs le 19.95) and $
;     (window eq 1) and (cc.gshort gt 0) and (cc.bgood eq 1) and (cc.rgood eq 1) and (ww.target_weight gt 0))
    
; select the parent sample    
    agespath = ages_path(/analysis)
    sfrmpath = ages_path(/projects)+'sfrm/'

    if (n_elements(phot) eq 0) then begin
       phot_version = ages_version(/photometry)
       phot = mrdfits(agespath+'ages_photometry_'+phot_version+'.fits.gz',1)
    endif
    if (n_elements(kcorr) eq 0) then begin
       kcorr_version = ages_version(/kcorrect)
       kcorr = mrdfits(agespath+'ages_photo_kcorr_'+kcorr_version+'.fits.gz',1)
    endif

;   window = ages_isin_window(phot.ra,phot.dec)
;   parent = where((anc1.i_obs le ifaint) and $
;     (anc1.sdss_star eq 0) and (anc1.target_weight gt 0.0) and $
;     (window eq 1) and (anc1.gshort gt 0) and $ ;(anc1.mass gt 0.0) and $
;     (anc1.bgood gt 0) and (anc1.rgood gt 0),ngal)
    
    
    ifaint = 19.95 ; [Vega]
    stars = mrdfits(agespath+'catalog.sdss_photometry.fits.gz',1)
    istar = (stars.type eq 6) and (stars.psfmag_r lt 19.0)

    window = ages_isin_window(phot.ra,phot.dec)
    kk = mrdfits(agespath+'catalog.kcorr.v3.fits.gz',1)
    parent = where($
      (kk.flag eq 1) and $
      (phot.i_obs le ifaint) and $
      (istar eq 0) and $
      (phot.gshort gt 0) and (phot.gshort ne 2048) and $
      (window eq 1),ngal)
    splog, 'Ngal = ', ngal
    
    cc = mrdfits(agespath+'catalog.codes.fits.gz',1)
    ww = mrdfits(agespath+'catalog.spectweight.fits.gz',1)
    kk = mrdfits(agespath+'catalog.kcorr.v3.fits.gz',1)
    dan = where((ww.main_weight ge 1) and (ww.spec_yesno ge 1) and (ww.z_yesno ge 1) and (kk.flag eq 1))
    m1 = cmset_op(parent,'and',dan)       ; objects in common
    m2 = cmset_op(parent,'and',dan,/not2) ; in my sample, not in Dan's
    m3 = cmset_op(parent,'and',dan,/not1) ; in Dan's sample, not in mine
    help, m1, m2, m3, dan, parent

    djs_plot, phot.ra, phot.dec, psym=3, xsty=3, ysty=3, xrange=[220,216], yr=[32,36], color='cyan'
    djs_oplot, phot[parent].ra, phot[parent].dec, psym=6, sym=0.1, color='orange' 
    djs_oplot, phot[dan].ra, phot[dan].dec, psym=3, color='red' 
    djs_oplot, phot[m3].ra, phot[m3].dec, psym=6, color='yellow'

    read_mangle_polygons, ages_path(/window)+'ages_bsmask.ply', stars
    flag = is_in_window(stars,ra=phot.ra,dec=phot.dec)
    pp = phot[where(flag eq 1 and phot.galaxy eq 1)]
    help, where(pp.galaxy eq 1)
    spherematch, pp.ra, pp.dec, usno.ra, usno.dec, 1.0/3600.0, m1, m2
    
    
    
    jj = is_in_window(stars,ra=phot[m3].ra,dec=phot[m3].dec)
    
stop    

    usno = mrdfits(ages_path(/analysis)+'ages_usno.fits.gz',1)
    keep = where(usno.rmag lt 15.0,nkeep)
    usno = usno[keep]

    djs_plot, phot.ra, phot.dec, psym=3, xsty=3, ysty=3, xrange=[219,219.2], yr=[33.8,34.0], color='cyan'
    djs_oplot, phot[parent].ra, phot[parent].dec, psym=6, sym=0.2, color='orange' 
    djs_oplot, phot[dan].ra, phot[dan].dec, psym=6, sym=0.2, color='red' 
    djs_oplot, phot[m3].ra, phot[m3].dec, psym=6, color='yellow'
    djs_oplot, usno.ra, usno.dec, ps=7



    read_mangle_polygons, ages_path(/window)+'ages_window_final.ply', balkans
    plot_poly, balkans, /over, color=djs_icolor('green')

    read_mangle_polygons, ages_path(/window)+'ages_fields.ply', balkans
    plot_poly, balkans, /over, color=djs_icolor('green')
    read_mangle_polygons, ages_path(/window)+'ages_bsmask.ply', stars
    plot_poly, stars, /over, color=djs_icolor('orange')
    
    plot, phot[dan].ra, phot[dan].dec, psym=4, xsty=3, ysty=3, xrange=[220,216], yr=[32,36] ;, xr=[216.3,216.6], yr=[33.3,34.1]
    djs_oplot, phot[m3].ra, phot[m3].dec, psym=4, color='red'
    
    
stop    
    
    
; basic parameters of the selection    
    sample_zmin = 0.05
    sample_zmax = 0.75

    vname = 'default.nolines'
    ifilt = 'ndwfs_I.par'
    ibright = 15.0 ; [Vega]
    ifaint = 19.95 ; [Vega]
    
; define the parent sample; the (soft) cut on MASS ensures that
; K-corrections have been measured (see the cuts in
; READ_AGES_PHOTOMETRY)
    window = ages_isin_window(anc1.ra,anc1.dec)

    parent = where((anc1.i_obs le ifaint) and $
      (anc1.sdss_star eq 0) and (anc1.target_weight gt 0.0) and $
      (window eq 1) and (anc1.gshort gt 0) and $ ;(anc1.mass gt 0.0) and $
      (anc1.bgood gt 0) and (anc1.rgood gt 0),ngal)

;   parent = where((anc1.i_tot gt ibright) and (anc1.i_tot lt ifaint) and $
;     (anc1.z gt sample_zmin) and (anc1.z lt sample_zmax) and $
;     (anc1.sdss_star eq 0) and (anc1.target_weight gt 0.0) and $
;     (window eq 1) and (anc1.mass gt 0.0) and $
;     (anc1.bgood gt 0) and (anc1.rgood gt 0),ngal)

    splog, 'Ngal = ', ngal

    anc = anc1[parent]
;   ispec = ispec1[parent]

    dan = where((anc1.main_weight ge 1) and (anc1.spec_yesno ge 1) and (anc1.z_yesno ge 1))
    m1 = cmset_op(parent,'and',dan)
    m2 = cmset_op(parent,'and',dan,/not2)
    m3 = cmset_op(parent,'and',dan,/not1)
    help, m1, m2, m3, dan, parent



    ww = mrdfits(ages_path(/analysis)+'catalog.spectweight.fits.gz',1)
    kk = mrdfits(ages_path(/analysis)+'catalog.kcorr.v3.fits.gz',1)
    cc = mrdfits(ages_path(/analysis)+'catalog.codes.fits.gz',1)
    all = mrdfits(ages_path(/analysis)+'catalog.cat.noguidestars.fits.gz',1)
    zz = mrdfits(ages_path(/analysis)+'catalog.zmerge.fits.gz',1)
    ddan = where((ww.main_weight ge 1) and (ww.spec_yesno ge 1) and (ww.z_yesno ge 1) and $
      (kk.flag eq 1))
    print, minmax(cc[ddan].gshort), minmax(kk[ddan].i_obs), minmax(kk[ddan].z)
    

stop    
    
    in = where(anc.main_weight eq 1,comp=out)
    plot, anc.ra, anc.dec, psym=3, xsty=3, ysty=3, xr=[217,218], yr=[34,34.5]
    djs_oplot, anc[in].ra, anc[in].dec, psym=4, color='red'
    djs_oplot, anc[out].ra, anc[out].dec, psym=4, color='blue'

    
    
    
stop    
    
    me = where((win eq 1) and (jj.sdss_star eq 0),comp=out)
    
    ww = where(anc.main_weight ge 1)
    jj = anc[ww]
    
    plot, jj.ra, jj.dec, psym=3, xsty=3, ysty=3, xr=[217,218], yr=[34,34.5]
    djs_oplot, jj[in].ra, jj[in].dec, psym=4, color='red'
    djs_oplot, jj[out].ra, jj[out].dec, psym=4, color='blue'
    
    
stop    
    
; compute ZMIN and ZMAX; be sure to use the targetting magnitude,
; I_obs, not the best magnitude, I_tot; also compute the final weight
; as the product of all the various selection terms
    zminzmax = im_zminzmax(anc.z,anc.i_obs,anc.coeffs,$
      bright=ibright,faint=ifaint,filter=ifilt,vname=vname)
    moretags = replicate({zmin: 0.0, zmax: 0.0, final_weight: 0.0},ngal)
    anc = struct_addtags(temporary(anc),moretags)
    anc.zmin = zminzmax.zmin
    anc.zmax = zminzmax.zmax
    anc.final_weight = anc.spec_weight*anc.target_weight*$
      anc.fiber_weight*anc.catalog_weight

; add the stellar masses
    isedpath = ages_path(/isedfit)

;   isedfile1 = isedpath+'BwRIJHKsirac_chab_calzetti_sfhgrid01.fits.gz'
;   splog, 'Reading '+isedfile1
;   ised1 = mrdfits(isedfile1,1,rows=anc.ages_id)

    isedfile2 = isedpath+'BwRIJHKsirac_chab_calzetti_sfhgrid02.fits.gz'
    splog, 'Reading '+isedfile2
    ised2 = mrdfits(isedfile2,1,rows=anc.ages_id)

    isedfile3 = isedpath+'BwRIJHKsirac_chab_calzetti_sfhgrid03.fits.gz'
    splog, 'Reading '+isedfile3
    ised3 = mrdfits(isedfile3,1,rows=anc.ages_id)

    tags = ['mass50','mass_err','age50','age_err','tau50','tau_err','chi2']
    newtags = 'isedfit_'+['mass','mass_err','age','age_err','tau','tau_err','chi2']
    isedout = im_struct_trimtags(ised2,select=tags,newtags=newtags)

; write out       
    anc = struct_addtags(temporary(anc),isedout)
;   final = where(anc.isedfit_chi2 lt 1E6,ngal)
;   anc = anc[final]
;   splog, 'Ngal = ', ngal
    im_mwrfits, anc, sfrmpath+'ages_parent_ancillary.fits', /clobber
    
return
end
    

;; compute Vmax
;    q0 = 0.0 & q1 = 0.0 & qz0 = 0.0 ; no evolution
;    noevol = im_calc_vmax(anc.(itag),anc.coeffs,ifilt,$
;      area,ibright,ifaint,sample_zmin,sample_zmax,$
;      actual_z=anc.z,q0=q0,q1=q1,qz0=qz0,/silent)
;
;    q0 = 1.6 & q1 = 0.0 & qz0 = 0.1 ; following Eisenstein+09
;    evol = im_calc_vmax(anc.(itag),anc.coeffs,ifilt,$
;      area,ibright,ifaint,sample_zmin,sample_zmax,$
;      actual_z=anc.z,q0=q0,q1=q1,qz0=qz0,/silent)
;
;    vmaxtags = {vol: 0.0, zmin_evol: 0.0, zmax_evol: 0.0, vmax_evol: 0.0, $
;      zmin_noevol: 0.0, zmax_noevol: 0.0, vmax_noevol: 0.0}
;    vmaxtags = replicate(vmaxtags,ngal)
;    anc = struct_addtags(temporary(anc),vmaxtags)
;    anc.vol         = noevol.vol
;    anc.zmin_evol   = evol.zmin
;    anc.zmax_evol   = evol.zmax
;    anc.vmax_evol   = evol.vmax
;    anc.zmin_noevol = noevol.zmin
;    anc.zmax_noevol = noevol.zmax
;    anc.vmax_noevol = noevol.vmax
