pro macs0329_z6arcs_isedfit, supergrid, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber
; jm11nov08ucsd - fit the z~6 arcs in macs0329

    isedpath = clash_path(/ised)
    datapath = clash_path(/macs0329_z6arcs)

    isedfit_sfhgrid_dir = clash_path(/monte)
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/clash_sfhgrid.par'

; read the supergrid parameter file
    supergrid = 2 ; current supergrid for this project
    super = get_clash_supergrid(2,nsuper=nsuper)
    struct_print, super

; gather the photometry
    alladi = rsex(datapath+'macs0329_z6arcs.cat')
    allcat = read_clash_catalog('macs0329',/arcs)
    spherematch, allcat.ra, allcat.dec, 15D*hms2dec(alladi.ra), hms2dec(alladi.dec), 1D/3600, m1, m2
    cat = allcat[m1]
    adi = alladi[m2]
    struct_print, adi

    igm = '1'
    nzz = '3'
    zlog = '0'
    prefix = 'macs0329_z6arcs'
    zminmax = [6.0,6.2]

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
          isedfit, paramfile, maggies, ivarmaggies, adi.z, iopath=isedpath, $
            clobber=clobber, sfhgrid_paramfile=sfhgrid_paramfile, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir;, index=index
       endif       

; make some QAplots
       if keyword_set(qaplot) then begin
          isedfit_qaplot, paramfile, result, iopath=isedpath, galaxy=adi.id, $
            index=index, clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
            outprefix=outprefix
       endif
    endfor

return
end
