pro clash_write_paramfile, paramfile, prefix=prefix, zminmax=zminmax, $
  nzz=nzz, zlog=zlog, igm=igm, super=super
; jm11oct14ucsd - support routine for the various clash_*_isedfit
; routines
    
    filters = clash_filterlist()

    sfhgridstring = strtrim(super.sfhgrid,2)
    redcurvestring = redcurve2string(super.redcurve)
    synthmodels = strtrim(super.synthmodels,2)
    imf = strtrim(super.imf,2)
    h100 = clash_h100(omega0=omega0,omegal=omegal)
    
    splog, 'Writing '+paramfile
    zrange = string(zminmax[0],format='(F4.2)')+','+string(zminmax[1],$
      format='(F4.2)')+','+nzz+','+zlog+' # [minz,maxz,dz,log?]'
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
    printf, lun, 'igm                  '+igm+' # [0=no, 1=yes]'
    printf, lun, 'maxold               0 # [0=no, 1=yes]'
    printf, lun, 'filterlist           '+strjoin(filters,',')
    free_lun, lun 

return
end

