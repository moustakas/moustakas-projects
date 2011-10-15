pro clash_bcg_isedfit, supergrid, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber
; jm11oct14ucsd - fit the CLASH BCGs

; echo "clash_bcg_isedfit, /model, /ised, /clob" | idl > & ~/clash.bcg.log & 
    
    isedpath = clash_path(/ised)
    catpath = clash_path(/cat)
    isedfit_sfhgrid_dir = clash_path(/monte)
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/clash_sfhgrid.par'

; read the supergrid parameter file
    super = get_clash_supergrid(supergrid,nsuper=nsuper,/bcg)
    struct_print, super

; gather the BCG photometry - will eventually need a wrapper for this 
    bcgfile = catpath+'a2261_bcg_apmag_'+['3arcsec','9arcsec','total']+'.txt'
    cat = [read_postman(bcgfile[0]),read_postman(bcgfile[1]),read_postman(bcgfile[2])]
    cat = struct_addtags(cat,replicate({z: 0.22331, galaxy: ''},3))
    cat.galaxy = 'A2261 BCG - '+['3 arcsec','9 arcsec','Total']
    
    igm = '0'
    nzz = '3'
    zlog = '0'
    prefix = 'bcgs'
    zminmax = [0.22,0.24]

; loop on each supergrid
    for gg = 0, nsuper-1 do begin
       splog, 'working on grid '+strtrim(super[gg].supergrid,2)

       paramfile = isedpath+prefix+'_supergrid'+string(super[gg].supergrid,$
         format='(i2.2)')+'_isedfit.par'
       clash_write_paramfile, paramfile, prefix=prefix, zminmax=zminmax, $
         nzz=nzz, zlog=zlog, igm=igm, super=super[gg]
       
; build the models
       if keyword_set(models) then begin
          isedfit_models, paramfile, iopath=isedpath, clobber=clobber, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir
       endif

; do the fitting!
       if keyword_set(isedfit) then begin
          clash_to_maggies, cat, maggies, ivarmaggies
          isedfit, paramfile, maggies, ivarmaggies, cat.z, iopath=isedpath, $
            clobber=clobber, sfhgrid_paramfile=sfhgrid_paramfile, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir
       endif       

; make some QAplots
       if keyword_set(qaplot) then begin
          isedfit_qaplot, paramfile, result, iopath=isedpath, galaxy=cat.galaxy, $
            index=index, clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
            outprefix=outprefix
       endif
    endfor

return
end
