pro z11_write_paramfile, paramfile, prefix=prefix, zminmax=zminmax, $
  nzz=nzz, zlog=zlog, igm=igm, super=super, filters=filters
; jm12aug14siena - support routine for the various z11_*_isedfit
; routines 

    isedpath = file_dirname(paramfile)
    if file_test(isedpath) eq 0 then begin
       splog, 'Making isedfit directory '+isedpath
       spawn, 'mkdir -p '+isedpath, /sh
    endif
    
    filters = z11_filterlist()

    sfhgridstring = strtrim(super.sfhgrid,2)
    redcurvestring = redcurve2string(super.redcurve)
    synthmodels = strtrim(super.synthmodels,2)
    imf = strtrim(super.imf,2)
    h100 = clash_h100(omega0=omega0,omegal=omegal)
    
    splog, 'Writing '+paramfile
    zrange = strtrim(string(zminmax[0],format='(F5.2)'),2)+','+strtrim(string(zminmax[1],$
      format='(F5.2)'),2)+','+strtrim(nzz,2)+','+strtrim(zlog,2)+' # [minz,maxz,dz,log?]'
    openw, lun, paramfile, /get_lun
    printf, lun, 'h100                 '+string(h100,format='(F4.2)')
    printf, lun, 'omega0               '+string(omega0,format='(F4.2)')
    printf, lun, 'omegal               '+string(omegal,format='(F4.2)')
    printf, lun, 'synthmodels          '+synthmodels
    printf, lun, 'imf                  '+imf
    printf, lun, 'sfhgrid              '+sfhgridstring
    printf, lun, 'redcurve             '+redcurvestring
    printf, lun, 'prefix               '+prefix
    printf, lun, 'redshift             '+zrange
    printf, lun, 'igm                  '+strtrim(igm,2)+' # [0=no, 1=yes]'
    printf, lun, 'maxold               0 # [0=no, 1=yes]'
    printf, lun, 'filterlist           '+strjoin(filters,',')
    free_lun, lun 

return
end

