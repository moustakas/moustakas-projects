pro build_desi_deep2egs
; jm14mar01siena - build a sample of DEEP2 galaxies in the Extended
; Groth Strip (Field 1) in order to optimize DESI/ELG targeting

    templatepath = getenv('IM_PROJECTS_DIR')+'/desi/templates/v1.1/'
    isedfit_paramfile = templatepath+'desi_deep2_paramfile.par'
    targpath = getenv('IM_PROJECTS_DIR')+'/desi/targeting/'
    technotepath = getenv('IM_SVNREPOS')+'/desi/technotes/targeting/elg-deep2/trunk/'
    winpath = deep2_path(/window)
    catpath = deep2_path(/cat)
    
    area = 0.4342 ; deg^2
;   area = 0.6 ; [deg^2]

    brightcut = 18.5
    faintcut = 24.0
    oiicut1 = 8D-17 ; [erg/s/cm2]
    oiisnrcut = 3.0
    
    allphot = mrdfits(catpath+'deep2.pcat_ext.fits.gz',1)
    unwise = mrdfits(catpath+'deep2.dr4.unwise.fits.gz',1)
    match, allphot.objno, unwise.objno, m1, m2
    allphot = struct_addtags(allphot,im_empty_structure($
      struct_trimtags(unwise[0],select=['w1_*','w2_*']),$
      ncopies=n_elements(allphot)))
    allphot[m1] = im_struct_assign(unwise[m2],allphot[m1],/nozero)
       
    egs = where(strmid(strtrim(allphot.objno,2),0,1) eq '1',negs)
    phot_egs = deep2_get_ugriz(allphot[egs],/unwise)
    pointing_egs = strmid(strtrim(phot_egs.objno,2),0,2) ; test
    
; select objects in the spectroscopic footprint in a reasonable
; magnitude range with PGAL>0.5
    win = mrdfits(winpath+'windowf.egs.fits.gz',0,hdr)
    extast, hdr, astr
    area = total(win gt 0)*determ(astr.cd) 
    splog, 'Spectroscopic area (deg^2) = ', area
    
    ad2xy, phot_egs.ra, phot_egs.dec, astr, xx, yy
    maskweight = interpolate(win,xx,yy,missing=0)

    keep = where(maskweight gt 0.25 and phot_egs.pgal ge 0.5 and $
      phot_egs.badflag eq 0 and $
      (phot_egs.bestr + 2.5*alog10(!pi*(3*phot_egs.rg*0.207)^2)) le 26.5 and $
      phot_egs.ugriz[2,*] gt brightcut and phot_egs.ugriz[2,*] lt faintcut,nphot)
    phot = phot_egs[keep]

    nmissing = total(phot.ugriz[2,*] eq 0) ; missing r-band photometry --> none!
    splog, 'Parent photometric sample = ', nphot, negs
    splog, 'Missing r-band photometry = ', nmissing
    
    im_mwrfits, phot, targpath+'deep2egs-photparent.fits', /clobber

; also write out a sample of stars so we can investigate the stellar
; locus
    star = where(phot_egs.pgal lt 0.2 and phot_egs.badflag eq 0 and $
      phot_egs.ugriz[2,*] gt brightcut and phot_egs.ugriz[2,*] lt faintcut,nstar)
    splog, 'Sample of stars', nstar
    im_mwrfits, phot_egs[star], targpath+'deep2egs-photstars.fits', /clobber
       
; now build our parent spectroscopic sample of galaxies in the EGS
; field (apply the cut on Q>=3 below)
    allzcat = read_deep2_zcat(phot=allzcat_phot,weight=allzcat_weight,/all)
    allzcat = struct_addtags(struct_addtags(allzcat,struct_trimtags(allzcat_weight,select=$
      ['*_weight'])),replicate({cfhtls_u: 0.0, cfhtls_g: 0.0, cfhtls_r: 0.0, cfhtls_i: 0.0, $
      cfhtls_z: 0.0, mask_weight: -1.0},n_elements(allzcat)))
    egs = where(strmid(strtrim(allzcat.objno,2),0,1) eq '1',negs)
    
    zcat_egs = allzcat[egs]
    zcat_phot_egs = deep2_get_ugriz(allzcat_phot[egs],/unwise)

; this cut is to choose objects that are within the spectroscopic
; footprint of Field 1       
    pointing1 = strmid(strtrim(zcat_egs.objno,2),0,2)
    ad2xy, zcat_egs.ra, zcat_egs.dec, astr, xx, yy
    zcat_maskweight1 = interpolate(win,xx,yy,missing=0)
    zcat_egs.mask_weight = zcat_maskweight1
    
    keep = where(zcat_maskweight1 gt 0 and zcat_egs.targ_weight gt 0 and $
      zcat_phot_egs.ugriz[2,*] gt brightcut and zcat_phot_egs.ugriz[2,*] lt faintcut,nspec)
    zcat = zcat_egs[keep]
    zcat_phot = zcat_phot_egs[keep]
    zcat = struct_addtags(zcat,struct_trimtags(zcat_phot,$
      select=['ugriz','ugriz_err','bri','bri_err','wise','wise_err','*galfit*','*radius*']))
    zcat.cfhtls_u = zcat.ugriz[0]
    zcat.cfhtls_g = zcat.ugriz[1]
    zcat.cfhtls_r = zcat.ugriz[2]
    zcat.cfhtls_i = zcat.ugriz[3]
    zcat.cfhtls_z = zcat.ugriz[4]
    
    nspec_weighted = long(total(zcat.targ_weight))
    nmissing = total(zcat_phot.ugriz[2,*] eq 0) ; missing r-band photometry --> none!
    splog, 'Parent spectroscopic sample = ', negs, nspec, nspec_weighted
    splog, 'Missing r-band photometry = ', nmissing
    
    im_mwrfits, zcat, targpath+'deep2egs-zcatparent.fits', /clobber
    
; also write out a Q>=3 sample with the corresponding [OII] fluxes and
; upper limits 
    ppxf = read_deep2(/ppxf,/fixoii)
    match, zcat.objno, ppxf.objno, m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    zcat_phot_q34 = zcat_phot[m1]
    zcat_q34 = zcat[m1]
    ppxf_q34 = ppxf[m2]

; get the redshift success rate for this sample       
    deep2_to_maggies, zcat_phot_q34, maggies, filterlist=filters, /unwise
    params = get_deep2_completeness_params('EGS')
    colors = get_deep2_completeness_colors(maggies,params=params,$
      targ_mag=zcat_phot_q34.r,filterlist=filters)
    zweight = get_deep2_zsuccess('EGS',mag=colors.mag,$
      color1=colors.color1,color2=colors.color2)

    zcat_q34.zsuccess_weight = zweight
    zcat_q34.final_weight = 1/zweight

; build the [OII] flux
    kised = mrdfits(templatepath+'desi_deep2_fsps_v2.4_miles_'+$
      'chab_charlot_sfhgrid01_kcorr.z0.0.fits.gz',1,$
      rows=m2)
    oii = deep2_get_oiiflux(ppxf_q34,cflux_3727_rest=kised.cflux_3727)
    
; possibly come up with an algorithm to assign an [OII] flux to an
; object without one 
    zcat_q34 = struct_addtags(zcat_q34,oii)
    zcat_q34 = struct_addtags(zcat_q34,replicate({oii: 0.0},n_elements(zcat_q34)))
    zcat_q34.oii = oii.oii_3727[0] ; total flux so we can run ELGTUNE
    
    im_mwrfits, zcat_q34, targpath+'deep2egs-zcatparent.Q34.fits', /clobber
    
return
end
    
