pro irclusters_isedfit, supergrid, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber, use_zcluster=use_zcluster
; jm12jun08ucsd - measure stellar masses for galaxies in Brodwin's sample of IR-selected clusters

; echo "irclustesr_isedfit, /model, /ised, /clob" | idl > & ~/irclusters.log & 

    isedpath = irclusters_path(/isedfit)
    isedfit_sfhgrid_dir = irclusters_path(/monte)
    sfhgrid_paramfile = getenv('IRCLUSTERS_DIR')+'/irclusters_supergrid.par'

; read the supergrid parameter file
    super = get_irclusters_supergrid(supergrid,nsuper=nsuper)
    struct_print, super

    igm = '1'
    nzz = '150'
    zlog = '0'
    zminmax = [0.01D,3.7D]
    prefix = 'irclusters'

; loop on each supergrid
    for gg = 0, nsuper-1 do begin
       splog, 'Working on grid '+strtrim(super[gg].supergrid,2)

       paramfile = isedpath+prefix+'_supergrid'+string(super[gg].supergrid,$
         format='(i2.2)')+'_isedfit.par'
       irclusters_write_paramfile, paramfile, prefix=prefix, zminmax=zminmax, $
         nzz=nzz, zlog=zlog, igm=igm, super=super[gg], filters=filters
       
; build the models
       if keyword_set(models) then begin
          isedfit_models, paramfile, iopath=isedpath, clobber=clobber, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir
       endif
       
; optionally fit at the cluster redshift  
       if keyword_set(isedfit) then begin
          cat = read_irclusters(maggies=maggies,ivarmaggies=ivarmaggies)

          toss = where(strmatch(filters,'*ch3*') or strmatch(filters,'*ch4*'))
          ivarmaggies[toss,*] = 0.0

          delvarx, index, outprefix
          if keyword_set(use_zcluster) then begin
             zobj = cat.zclust
             outprefix = prefix+'_zclust'
          endif else begin
; only fit objects at 0<z<2
;            index = where(total(maggies gt 0,1) lt 3) ; test!!
;            index = [index,index+1]
             zobj = cat.z
             index = where((cat.z ge zminmax[0]) and (cat.z le zminmax[1]),nindex)
          endelse
          isedfit, paramfile, maggies, ivarmaggies, zobj, iopath=isedpath, $
            clobber=clobber, sfhgrid_paramfile=sfhgrid_paramfile, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, index=index, outprefix=outprefix
       endif       

; make some QAplots of a random subset of objects
       if keyword_set(qaplot) then begin
          cat = read_irclusters(maggies=maggies,ivarmaggies=ivarmaggies)
          if keyword_set(use_zcluster) then begin
             outprefix = prefix+'_zclust'
             qaindx = shuffle_indx(n_elements(cat),num=30)
          endif else begin
             index = where((cat.z ge zminmax[0]) and (cat.z le zminmax[1]),nindex)
             qaindx = index[shuffle_indx(nindex,num=30)]
          endelse
          galaxy = 'Galaxy '+string(qaindx,format='(I5.5)')
          isedfit_qaplot, paramfile, result, iopath=isedpath, galaxy=galaxy, $
            index=qaindx, clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
            outprefix=outprefix, /xlog
       endif
    endfor

return
end
