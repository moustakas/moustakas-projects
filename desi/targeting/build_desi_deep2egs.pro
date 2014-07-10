pro build_desi_deep2egs, build_phot=build_phot, build_spec=build_spec
; jm14mar01siena - build a sample of DEEP2 galaxies in the Extended
; Groth Strip (Field 1) in order to optimize DESI/ELG target selection
; of ELGs 

    targpath = getenv('IM_PROJECTS_DIR')+'/desi/targeting/'
    winpath = deep2_path(/window)
    catpath = deep2_path(/cat)

    version = desi_deep2_template_version()
    templatepath = getenv('IM_PROJECTS_DIR')+'/desi/templates/'+version+'/'

    brightcut = 18.5
    faintcut = 24.0
    oiicut1 = 8D-17 ; [erg/s/cm2]
    oiisnrcut = 3.0

; read the DEEP2/EGS window function    
    win = mrdfits(winpath+'windowf.egs.fits.gz',0,hdr)
    extast, hdr, astr
    area = total(win gt 0)*determ(astr.cd) ; =0.4342 deg^2
    splog, 'Spectroscopic area (deg^2) = ', area
    
; ##################################################
; first build a parent photometric catalog that we will use to test
; our number counts
    if keyword_set(build_phot) then begin
    
; join the Matthews+13 photometric catalog with Dustin's unWISE catalog
       allphot = mrdfits(catpath+'deep2.pcat_ext.fits.gz',1)
       unwise = mrdfits(catpath+'deep2.dr4.unwise.fits.gz',1)
       match, allphot.objno, unwise.objno, m1, m2
       allphot = struct_addtags(allphot,im_empty_structure($
         struct_trimtags(unwise[0],select=['w1_*','w2_*']),$
         ncopies=n_elements(allphot)))
       allphot[m1] = im_struct_assign(unwise[m2],allphot[m1],/nozero)
       
; select objects in the wider EGS field and parse the photometry (also
; renaming some annoying structure tags); NB!  for objects detected at
; <1-sigma, compute the upper limit and set the uncertainty to zero;
; finally, degrade the grz photometry to the anticipated depth of our
; DECam imaging 
       egs = where(strmid(strtrim(allphot.objno,2),0,1) eq '1',negs)
       phot_egs = deep2_get_ugriz(allphot[egs],/unwise,/degrade)
       
; next, select objects in the spectroscopic footprint in a fiducial 
; magnitude range with PGAL>0.5
       ad2xy, phot_egs.ra, phot_egs.dec, astr, xx, yy
       maskweight = interpolate(win,xx,yy,missing=0)

       keep = where(maskweight gt 0.25 and phot_egs.pgal ge 0.5 and $
         phot_egs.badflag eq 0 and $
         (phot_egs.r + 2.5*alog10(!pi*(3*phot_egs.rg*0.207)^2)) lt 26.5 and $ ; SB cut
         phot_egs.cfhtls_r gt brightcut and phot_egs.cfhtls_r lt faintcut,nphot)
       phot = phot_egs[keep]

       splog, 'Parent photometric sample = ', nphot, negs
       splog, 'Surface density (#/deg^2) to r=22.8, 23, 24', $
         total(phot.cfhtls_r lt 22.8)/area, total(phot.cfhtls_r lt 23)/area, $
         total(phot.cfhtls_r lt 24)/area
       splog, 'Surface density (#/deg^2) in grz color box to r=22.8, 23, 24', $
         n_elements(desi_get_hizelg(phot,magcut=22.8))/area, $
         n_elements(desi_get_hizelg(phot,magcut=23.0))/area, $
         n_elements(desi_get_hizelg(phot,magcut=24.0))/area
       im_mwrfits, phot, targpath+'deep2egs-phot.fits', /clobber

; also write out a sample of stars so we can investigate the position
; of the stellar locus
       star = where(phot_egs.pgal lt 0.2 and phot_egs.badflag eq 0 and $
         phot_egs.cfhtls_r ge brightcut and phot_egs.cfhtls_r le faintcut,nstar)
       splog, 'Sample of stars', nstar
       im_mwrfits, phot_egs[star], targpath+'deep2egs-stars.fits', /clobber
    end

; ##################################################
; next, build a parent *spectroscopic* sample of galaxies in the
; DEEP2/EGS field
    if keyword_set(build_spec) then begin
       allzcat = read_deep2_zcat(phot=allzcat_phot,weight=allzcat_weight,/all)
       allzcat = struct_addtags(allzcat,struct_trimtags(allzcat_weight,$
         select='*_weight'))
       egs = where(strmid(strtrim(allzcat.objno,2),0,1) eq '1',negs)
       
       zcat_egs = allzcat[egs]
       zcat_phot_egs = deep2_get_ugriz(allzcat_phot[egs],/unwise,/degrade)

; this cut is to choose objects that are within the spectroscopic
; footprint of Field 1       
       pointing1 = strmid(strtrim(zcat_egs.objno,2),0,2)
       ad2xy, zcat_egs.ra, zcat_egs.dec, astr, xx, yy
       zcat_maskweight1 = interpolate(win,xx,yy,missing=0)
       
       keep = where(zcat_maskweight1 gt 0 and zcat_egs.targ_weight gt 0 and $
         zcat_phot_egs.cfhtls_r gt brightcut and zcat_phot_egs.cfhtls_r lt faintcut,nspec)

       zcat = zcat_egs[keep]
       zcat_phot = zcat_phot_egs[keep]
       zcat = struct_addtags(zcat,struct_trimtags(zcat_phot,select=$
         ['b','r','i','berr','rerr','ierr','cfhtls_*','w1*','w2*',$
         '*galfit*','*radius*']))

; cut to Q>=3 and then cross-match with the [OII] flux catalog
       zcat_q34 = zcat[where(zcat.zquality ge 3,ngal)]

       oiicat = mrdfits(templatepath+'deep2_oiicat_'+version+'.fits.gz',1)
       match, zcat_q34.objno, oiicat.objno, m1, m2
       srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
       if total(zcat_q34.objno-oiicat[m2].objno) ne 0 then message, 'Problem!'

       zcat_q34 = struct_addtags(zcat_q34,struct_trimtags(oiicat[m2],$
         select=['sigma_kms','oii_*']))

; assign [OII] fluxes to the missing objects by finding the nearest
; object in rest-frame color-magnitude space
       
       
       
stop       

       
;      nspec_weighted = long(total(zcat.targ_weight))
;      splog, 'Parent spectroscopic sample = ', negs, nspec, nspec_weighted
       
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
;   oii = deep2_get_oiiflux(ppxf_q34,cflux_3727_rest=kised.cflux_3727)
    
; possibly come up with an algorithm to assign an [OII] flux to an
; object without one 
       zcat_q34 = struct_addtags(zcat_q34,oii)
       zcat_q34 = struct_addtags(zcat_q34,replicate({oii: 0.0},n_elements(zcat_q34)))
       zcat_q34.oii = oii.oii_3727[0] ; total flux so we can run ELGTUNE
       
       im_mwrfits, zcat_q34, targpath+'deep2egs-zcatparent.Q34.fits', /clobber
    endif
    
return
end
    
