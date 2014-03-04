pro desi_lris14mar
; jm14feb26siena - select candidate z>1.4 galaxies that DEEP2 failed
; to get redshifts for, for follow-up with Keck/LRIS, all for training
; DESI targeting

    outpath = getenv('IM_PROJECTS_DIR')+'/desi/lris14mar/'

    deep3 = mrdfits(getenv('IM_DATA_DIR')+'/deep3/deep3.egs.2012jun13.fits.gz',1)
    deep2 = mrdfits(deep2_path(/catalogs)+'zcat.deep2.dr4.uniq.fits.gz',1)
    deep23 = im_empty_structure(deep3,ncopies=n_elements(deep2))
    deep23 = im_struct_assign(deep2,deep23)
    all = [deep3,deep23]
    
;   all = mrdfits(deep2_path(/catalogs)+'zcat.deep2.dr4.uniq.fits.gz',1)
    br = all.magb-all.magr
    ri = all.magr-all.magi
    qq = all.zquality
    zz = all.z
    
; select likely z>1.4 galaxies using criteria recommended by J. Newman
    fail = where((qq gt 0 and qq lt 3) and (br lt 0.6) and (ri lt 0.5) and $
      strmid(strtrim(all.objno,2),0,1) eq '1',nfail)
    splog, nfail
    out = struct_trimtags(all[fail],select=['ra','dec','magb','magr','magi'])

    im_mwrfits, out, outpath+'lris14mar_deep2egs_candidates.fits', /clobber

; build a comparison catalog for KG with just magnitude and ra,dec
    train = where(qq ge 3 and zz lt 1.4 and $
      strmid(strtrim(all.objno,2),0,1) eq '1',ntrain)
    splog, ntrain
    out = struct_trimtags(all[train],select=['ra','dec','magb','magr','magi'])

    im_mwrfits, out, outpath+'lris14mar_deep2egs_lowztrain.fits', /clobber

return
end



