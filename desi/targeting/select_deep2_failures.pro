pro select_deep2_failures
; jm14feb26siena - select candidate z>1.3 galaxies that DEEP2 failed
; to get redshifts for, for follow-up with Keck/LRIS, all for training
; DESI targeting

    outpath = getenv('IM_PROJECTS_DIR')+'/desi/deep2/lris14mar/'
    
    all = mrdfits(deep2_path(/catalogs)+'zcat.deep2.dr4.uniq.fits.gz',1)
    br = all.magb-all.magr
    ri = all.magr-all.magi
    qq = all.zquality
    zz = all.zbest

    need = where((qq gt 0 and qq lt 3) and (br lt 0.6) and (ri lt 0.5) and $
      strmatch(all.objname,'11*'),nneed)
    
; build a comparison catalog for KG with just magnitude and ra,dec    
    good = where(qq ge 3 and zz lt 1.4,ngood)
    out = struct_trimtags(all[good],select=['ra','dec','magb','magr','magi'])

    im_mwrfits, out, outpath+
    

return
end



