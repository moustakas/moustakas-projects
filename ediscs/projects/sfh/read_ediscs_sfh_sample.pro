function read_ediscs_sfh_sample, field=field, cluster=cluster, all=all, silent=silent
; jm10may03ucsd - read the output from BUILD_EDISCS_SFH_SAMPLE

;   common sfh_sample, ediscs_field, ediscs_cluster

    if keyword_set(all) then suffix = '_all' else suffix = ''
    datapath = ediscs_path(/projects)+'sfh/'
    if keyword_set(field) then begin
       thisfile = datapath+'ediscs_sfh_field'+suffix+'.fits.gz'
       if (size(ediscs_field,/type) ne 8) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          ediscs_field = mrdfits(thisfile,1,silent=0)
       endif else if (keyword_set(silent) eq 0) then $
         splog, 'Restoring '+file_basename(thisfile)
       return, ediscs_field
    endif

    if keyword_set(cluster) then begin
       thisfile = datapath+'ediscs_sfh_cluster'+suffix+'.fits.gz'
       if (size(ediscs_cluster,/type) ne 8) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          ediscs_cluster = mrdfits(thisfile,1,silent=0)
       endif else if (keyword_set(silent) eq 0) then $
         splog, 'Restoring '+file_basename(thisfile)
       return, ediscs_cluster
    endif

end
