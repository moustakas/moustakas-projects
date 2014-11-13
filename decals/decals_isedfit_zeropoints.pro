pro decals_isedfit_zeropoints

    prefix = 'decals'
    isedfit_dir = getenv('DECALS_DIR')+'/isedfit/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

; fitted sample    
    cat = mrdfits(isedfit_dir+'decals_edr.fits.gz',1)
    these = where(cat.z ge 0.05 and cat.z le 0.7 and cat.zwarning eq 0 and $
      cat.brick_primary eq 'T' and strtrim(cat.class,2) eq 'GALAXY',ngal)
    cat = cat[these]

; isedfit results    
    ised = read_isedfit(isedfit_paramfile,isedfit_dir=isedfit_dir,outprefix='nodecals')

stop    

return
end
    
