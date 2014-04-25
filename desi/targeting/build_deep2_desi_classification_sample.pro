function get_ugriz, cat, bri=bri, errbri=brierr, ugrizerr=ugrizerr
; compute dust-corrected ugriz magnitudes
    deep2_to_maggies, cat, mm, ii
    mag = maggies2mag(mm,ivarmaggies=ii,magerr=magerr)
    bri = mag[0:2,*]
    brierr = magerr[0:2,*]
    ugriz = mag[3:7,*]
    ugrizerr = magerr[3:7,*]
return, ugriz
end

function get_ugriz, cat
; compute dust-corrected ugriz magnitudes
    deep2_to_maggies, cat, mm, ii
    offset = 3 ; offset from the BRI photometry
    ugriz = fltarr(5,n_elements(cat))-99.0
    for ii = 0, 4 do begin
       good = where(mm[offset+ii,*] gt 0)
       ugriz[ii,good] = -2.5*alog10(mm[offset+ii,good])
    endfor
return, ugriz
end

pro build_deep2_desi_classification_sample
; jm14mar30siena - build a sample of DEEP2/EGS galaxies so that Ross
; Fadely can work his Bayesian classification magic on it 

    templatepath = getenv('IM_PROJECTS_DIR')+'/desi/templates/'
    isedfit_paramfile = templatepath+'desi_deep2_paramfile.par'
    winpath = deep2_path(/window)
    outpath = getenv('IM_PROJECTS_DIR')+'/desi/targeting/classification/'

    brightcut = 18.0
    faintcut = 24.0
    oiicut1 = 8D-17 ; [erg/s/cm2]
    
    allphot = mrdfits(deep2_path(/cat)+'deep2.pcat_ext.fits.gz',1)    
    allugriz = get_ugriz(allphot)

    egs = where(strmid(strtrim(allphot.objno,2),0,1) eq '1',negs)
    phot1 = allphot[egs]
    ugriz1 = allugriz[*,egs]

; select objects in the spectroscopic footprint in a reasonable
; magnitude range with PGAL>0.5
    win = mrdfits(winpath+'windowf.egs.fits.gz',0,hdr)
    extast, hdr, astr
    area = total(win gt 0)*determ(astr.cd) 
    splog, 'Spectroscopic area (deg^2) = ', area
       
    ad2xy, phot1.ra_deep, phot1.dec_deep, astr, xx, yy
    weight = interpolate(win,xx,yy,missing=0)
       
    keep = where(weight gt 0 and ugriz1[2,*] lt faintcut,nphot) $ ; phot1.pgal ge 0.5 and $
    phot = phot1[keep]
    ugriz = ugriz1[*,keep]
    
    splog, 'Parent photometric sample = ', nphot, negs

stop    
    
; write out the parent photometric sample
    phot = struct_addtags(phot,replicate({ugriz: fltarr(5)},nphot))
    phot.ugriz = ugriz
    im_mwrfits, phot, outpath+'egs_parent_phot.fits', /clobber

; select and write out the parent spectroscopic sample       
    allzcat = read_deep2_zcat(phot=allzcat_phot,weight=allzcat_weight)
    parentindx = lindgen(n_elements(allzcat))
    allzcat_ugriz = get_ugriz(allzcat_phot)

    egs = where(strmid(strtrim(allzcat_phot.objno,2),0,1) eq '1',negs)
    parentindx = parentindx[egs]
    zcat1 = allzcat[egs]
    zcat_phot1 = allzcat_phot[egs]
    zcat_weight1 = allzcat_weight[egs]
    zcat_ugriz1 = allzcat_ugriz[*,egs]
       
    keep = where(zcat_weight1.weight gt 0 and zcat_ugriz1[2,*] gt brightcut and $
      zcat_ugriz1[2,*] lt faintcut,nspec)
    parentindx = parentindx[keep]
    zcat = zcat1[keep]
    zcat_phot = zcat_phot1[keep]
    zcat_weight = zcat_weight1[keep]
    zcat_ugriz = zcat_ugriz1[*,keep]

    nspec_weighted = long(total(1.0/zcat_weight.weight))
    nmissing = total(zcat_ugriz[2,*] eq 0) ; missing r-band photometry --> none!
    splog, 'Missing r-band photometry = ', nmissing
    splog, 'Parent spectroscopic sample = ', negs, nspec, nspec_weighted
    
    out = struct_addtags(zcat,replicate({ugriz: fltarr(5), $
      weight: 0.0, parentindx: 0L},nspec))
    out.parentindx = parentindx
    out.ugriz = zcat_ugriz
    out.weight = zcat_weight.weight
    im_mwrfits, out, outpath+'egs_parent_deep2.fits', /clobber
    
    
    oii = deep2_get_oiiflux(line,cflux_3727_obs=cflux_3727_obs,$
      cflux_3727_rest=cflux_3727_rest)
       
       
return
end
    
