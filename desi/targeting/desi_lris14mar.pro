pro desi_lris14mar
; jm14feb26siena - select candidate z>1.4 galaxies that DEEP2 failed
; to get redshifts for, for follow-up with Keck/LRIS, all for training
; DESI targeting

    outpath = getenv('IM_PROJECTS_DIR')+'/desi/lris14mar/'
    
    all = mrdfits(deep2_path(/catalogs)+'zcat.deep2.dr4.uniq.fits.gz',1)
    br = all.magb-all.magr
    ri = all.magr-all.magi
    qq = all.zquality
    zz = all.zbest
    
; select likely z>1.4 galaxies using criteria recommended by J. Newman
    need = where((qq gt 0 and qq lt 3) and (br lt 0.6) and (ri lt 0.5) and $
      strmatch(all.objname,'11*'),nneed)
    out = struct_trimtags(all[need],select=['ra','dec','magb','magr','magi'])
    
    im_mwrfits, out, outpath+'lris14mar_deep2_candidates.fits', /clobber

stop    
    
; build a comparison catalog for KG with just magnitude and ra,dec
    good = where(qq ge 3 and zz lt 1.4 and strmatch(all.objname,'11*'),ngood)
    out = struct_trimtags(all[good],select=['ra','dec','magb','magr','magi'])

    im_mwrfits, out, outpath+'lris14mar_deep2_lowz.fits', /clobber
    

return
end



