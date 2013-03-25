function read_clash_catalog, cluster, redshift=redshift, arcs=arcs, ir=ir, $
  path=path
; jm11oct14ucsd - read SExtractor (default) or redshift catalog for a
; given cluster

    path = clash_path(cluster,redshift=redshift,arcs=arcs,ir=ir,/catalogs)

; optionally read the redshift catalog...
    if keyword_set(redshift) then begin

; temporarily deal with MACS1206 - everything is different
       if strmatch(cluster,'*1206*') then begin
          file = path+'m1206_specz.cat'
          if file_test(file) eq 0 then file = file+'.gz' ; try gzipped
       endif else begin
          file = path+repstr(strlowcase(cluster),'_','')+'_redshifts.dat'
          if file_test(file) eq 0 then file = file+'.gz' ; try gzipped

          if file_test(file) eq 0 then begin
             splog, 'Redshift catalog '+file+' not found!'
             return, -1
          endif
       endelse
       if strmatch(file,'*.gz*') then compress = 1

       if file_test(file) eq 0 then begin
          splog, 'Redshift catalog '+file+' not found!'
          return, -1
       endif

       if strmatch(cluster,'*1206*') then begin
          readcol, file, id, ra, dec, z, flag, $
            format='A,D,D,F,F', comment='#', /silent, compress=compress
          zerr = z*0.0
          zspec_R = z*0.0
          alt_id = strarr(n_elements(id))
          refcode = strarr(n_elements(id))
       endif else begin
          readcol, file, ra, dec, z, zerr, zspec_R, id, alt_id, refcode, $
            format='D,D,F,F,F,A,A,A', comment='#', /silent, compress=compress
       endelse
       ngal = n_elements(ra)
       
       cat = replicate({ra: 0D, dec: 0D, z: 0.0, zerr: 0.0, $
         zspec_R: 0.0, id: '', alt_id: '', refcode: ''},ngal)
       cat.ra = ra
       cat.dec = dec
       cat.z = z
       cat.zerr = zerr
       cat.zspec_R = zspec_R
       cat.id = id
       cat.alt_id = alt_id
       cat.refcode = refcode
       
       return, cat
    endif

; optionally read the arcs photometric catalog...
    if keyword_set(arcs) then begin
       file = path+'arcs.cat'
       if file_test(file) eq 0 then begin
          splog, 'Arcs SE catalog '+file+' not found!'
          return, -1
       endif
       cat = rsex(file)
       return, cat
    endif
       
; ...otherwise read the photometric catalog
    thiscluster = repstr(strlowcase(cluster),'_','')
    if keyword_set(ir) then suffix = 'IR' else suffix = 'ACS_IR'
    file = path+thiscluster+'_'+suffix+'.cat'
    if file_test(file) eq 0 then begin
       file = path+thiscluster+'_'+suffix+'.cat.gz'
       if file_test(file) eq 0 then begin
          file = path+strlowcase(cluster)+'_'+suffix+'.cat.gz' ; check different name
          if file_test(file) eq 0 then begin
             splog, 'SE catalog '+file+' not found!'
             return, -1
          endif
       endif 
    endif else begin
       file = path+strlowcase(cluster)+'_'+suffix+'.cat' ; check different name
       if file_test(file) eq 0 then begin
          splog, 'SE catalog '+file+' not found!'
          return, -1
       endif
    endelse 

    if strmatch(file,'*.gz*') then gz = 1 else gz = 0
    if gz then spawn, 'gunzip '+file, /sh
    cat = rsex(repstr(file,'.gz',''))
    if gz then spawn, 'gzip '+repstr(file,'.gz',''), /sh
    
return, cat    
end
