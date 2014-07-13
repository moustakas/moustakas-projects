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
       phot_egs = deep2_get_ugriz(allphot[egs],/unwise,degrade=0)
       
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
       zcat_phot_egs = deep2_get_ugriz(allzcat_phot[egs],/unwise,degrade=0)

; this cut is to choose objects that are within the spectroscopic
; footprint of Field 1 and which also have Q>=3
       pointing1 = strmid(strtrim(zcat_egs.objno,2),0,2)
       ad2xy, zcat_egs.ra, zcat_egs.dec, astr, xx, yy
       zcat_maskweight1 = interpolate(win,xx,yy,missing=0)
       
       keep = where(zcat_maskweight1 gt 0 and zcat_egs.zquality ge 3 and $
         zcat_egs.targ_weight gt 0 and zcat_phot_egs.cfhtls_r gt brightcut and $
         zcat_phot_egs.cfhtls_r lt faintcut,nspec)

       zcat = zcat_egs[keep]
       zcat_phot = zcat_phot_egs[keep]
       zcat = struct_addtags(zcat,struct_trimtags(zcat_phot,select=$
         ['b','r','i','berr','rerr','ierr','cfhtls_*','w1*','w2*',$
         '*galfit*','*radius*']))

; cross-match with the [OII] flux catalog
       kcorr = mrdfits(templatepath+'desi_deep2_fsps_v2.4_miles_'+$
         'chab_charlot_sfhgrid01_kcorr.z0.0.fits.gz',1)
       oiicat = mrdfits(templatepath+'deep2_oiicat_'+version+'.fits.gz',1)
       match, zcat.objno, oiicat.objno, m1, m2
       srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
       if total(zcat.objno-oiicat[m2].objno) ne 0 then message, 'Problem!'

       zcat = struct_addtags(zcat,struct_trimtags(oiicat[m2],$
         select=['sigma_kms','oii_*']))
       kcorr = kcorr[m2]

; assign [OII] fluxes to the missing objects by finding the nearest
; object in rest-frame color-magnitude space
       zmin = 0.6               ; 0.8
       zmax = 1.4

       refindx = where(zcat.zbest gt 0.8 and zcat.zbest lt 1.1 and zcat.oii_3727_err gt 0,nrefindx)
;      refindx = where(zcat.zbest gt zmin and zcat.zbest lt zmax and zcat.oii_3727_err gt 0,nrefindx)
;      refindx = where(zcat.zbest gt zmin and zcat.zbest lt zmax and $
;        zcat.oii_3727_err gt 0 and zcat.oii_3727_ew/zcat.oii_3727_ew_err gt 1.0,nrefindx)
       mrref = kcorr[refindx].absmag[2]
       umrref = kcorr[refindx].absmag[0]-mrref
       oiiref = zcat[refindx].oii_3727
       
       needindx = where(zcat.zbest gt zmin and zcat.zbest lt zmax and zcat.oii_3727_err lt 0.0,nneedindx)
;      needindx = where(zcat.zbest gt zmin and zcat.zbest lt zmax and $
;        (zcat.oii_3727_err eq 0.0 or zcat.oii_3727_ew/zcat.oii_3727_ew_err le 1.0),nneedindx)
       mrneed = kcorr[needindx].absmag[2]
       umrneed = kcorr[needindx].absmag[0]-mrneed
       oiineed = zcat[needindx].oii_3727

       keep = lonarr(nneedindx)
       for ii = 0, nneedindx-1 do begin
          dist = sqrt((mrref-mrneed[ii])^2+(umrref-umrneed[ii])^2)
          mindist = min(dist,thisref)
;         print, rref[thisref], rneed[ii], rzref[thisref], rzneed[ii], $
;           grref[thisref], grneed[ii], zref[thisref], zneed[ii]
          keep[ii] = thisref

          junk = struct_trimtags(zcat[refindx[thisref]],select='oii_*')
          zcat[needindx[ii]] = im_struct_assign(junk,zcat[needindx[ii]],/nozero)
       endfor

; here's some code to get the redshift success rate for this
; sample, but J. Newmann suggests that we don't apply this 
;       deep2_to_maggies, zcat_phot_q34, maggies, filterlist=filters, /unwise
;       params = get_deep2_completeness_params('EGS')
;       colors = get_deep2_completeness_colors(maggies,params=params,$
;         targ_mag=zcat_phot_q34.r,filterlist=filters)
;       zweight = get_deep2_zsuccess('EGS',mag=colors.mag,$
;         color1=colors.color1,color2=colors.color2)
;       zcat.zsuccess_weight = zweight

; assign the final weight as the targeting weight
       zcat.final_weight = zcat.targ_weight

; quick QAplots       
       djs_plot, mrref, umrref, psym=3, xsty=3, ysty=3, color='cyan', $
         xrange=[-18,-24], yrange=[-0.3,2.7]
       djs_oplot, mrneed, umrneed, psym=6, symsize=0.1, color='green'

;       ww = where(zcat.oii_3727 gt 8D-17)
;       im_plothist, zcat.z, bin=0.05
;       im_plothist, zcat[ww].z, bin=0.05, /over, /fill
       
       im_mwrfits, zcat, targpath+'deep2egs-oii.fits', /clobber
       
    endif
    
return
end
    
