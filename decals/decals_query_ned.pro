;+
; NAME:
;   ned_query_atlas
; PURPOSE:
;   Get all the low redshift objects in NED
; CALLING SEQUENCE:
;   ned_query_atlas
; COMMENTS:
;   Emails some batch queries to NED.
; REVISION HISTORY:
;   31-Mar-2004  MRB, NYU
;-
;------------------------------------------------------------------------------
pro decals_query_ned, st=st, version=version

    if (NOT keyword_set(st)) then st=0L

    outdir = getenv('IM_RESEARCH_DIR')+'decals/catalogs/ned'
    file_mkdir, outdir

    chunk=20L
    for chunkst=st, 359L, chunk do begin
       for id=chunkst+0L, ((chunkst+chunk-1L)<359L) do begin
          rast=strtrim(string(id, f='(f40.2)'),2)+'d'
          rand=strtrim(string(id+1, f='(f40.2)'),2)+'d'
          
          filename='blanton_atlas_'+string(id,f='(i3.3)')+'.txt'
          if(NOT file_test(outdir+'/'+filename+'.gz')) then begin
             splog, filename
             queryfile='query_'+filename
             openw,unit, outdir+'/'+queryfile, /get_lun 
             printf, unit, 'OUTPUT_FILENAME   '+filename
             printf, unit, 'OUTPUT_OPTION    compact'
             printf, unit, 'COMPRESS_OPTION    gzip'
             printf, unit, 'OUTPUT_COORDINATE_SYSTEM      equatorial'
             printf, unit, 'OUTPUT_EQUINOX                 J2000.0'
             printf, unit, 'REDSHIFT_VELOCITY     10000000.0'
             printf, unit, 'FIND_OBJECTS_BY_PARAMETERS'
             printf, unit, 'RA_between '+rast+', '+rand
             printf, unit, 'REDSHIFT Between -0.01, 0.055'
             printf, unit, 'UNIT z'
             printf, unit, 'END_OF_DATA'
             printf, unit, 'END_OF_REQUESTS'
             free_lun, unit
             spawn, 'mail nedbatch@ned.ipac.caltech.edu < '+outdir+'/'+queryfile
          endif
       endfor
       wait, 900.
    endfor

return
end
