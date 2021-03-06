function clash_path, cluster, catalogs=catalogs, redshift=redshift, mosaics=mosaics, $
  ir=ir, arcs=arcs, isedfit=isedfit, montegrids=montegrids, bcgimf=bcgimf, $
  macs0329_z6arcs=macs0329_z6arcs, santorini=santorini, megaspitzer=megaspitzer, $
  lensedvolumes=lensedvolumes, z11=z11, bcgmodels=bcgmodels, mas30=mas30
; jm11apr18ucsd - 

    clashpath = getenv('CLASH_ARCHIVE')+'/'
;   if keyword_set(catalogs) then clashpath = clashpath+'catalogs/'
    if keyword_set(isedfit) then clashpath = clashpath+'isedfit/'
    if keyword_set(montegrids) then clashpath = clashpath+'montegrids/'
    if keyword_set(bcgimf) then clashpath = clashpath+'projects/bcgimf/'
    if keyword_set(macs0329_z6arcs) then clashpath = clashpath+'projects/macs0329_z6arcs/'
    if keyword_set(santorini) then clashpath = clashpath+'projects/santorini/'
    if keyword_set(megaspitzer) then clashpath = clashpath+'projects/megaspitzer/'
    if keyword_set(lensedvolumes) then clashpath = clashpath+'projects/lensedvolumes/'
    if keyword_set(z11) then clashpath = clashpath+'projects/z11/'

; if CLUSTER has been given, then build the relevant ARCHIVE path;
; assume the user wants the 'catalogs' directory by default; deal with
; the possibility that the cluster name has an underscore (or not)
    if n_elements(cluster) ne 0 then begin
       archivepath = clashpath
;      archivepath = clashpath+'archive/'
       thiscluster = strlowcase(cluster)
       if file_test(archivepath+thiscluster,/dir) eq 0 then begin
          thiscluster = repstr(strlowcase(cluster),'_','')
          if file_test(archivepath+thiscluster,/dir) eq 0 then begin
             thiscluster = repstr(strlowcase(cluster),'abell','abell_')
          endif
       endif

       if keyword_set(catalogs) then begin
          clashpath = archivepath+thiscluster+'/HST/catalogs/mosaicdrizzle_image_pipeline/'
          if keyword_set(ir) then clashpath = clashpath+'IR_detection/' else $
            clashpath = clashpath+'ACS_IR_detection/'
;         clashpath = clashpath+'SExtractor/'
       endif

       if keyword_set(mosaics) then begin
          if keyword_set(mas30) then suff = 'scale_30mas/' else suff = 'scale_65mas/'
          clashpath = archivepath+thiscluster+'/HST/images/'+$
            'mosaicdrizzle_image_pipeline/'+suff
       endif
       if keyword_set(bcgmodels) then begin
          if keyword_set(mas30) then suff = '30mas/' else suff = ''
          clashpath = archivepath+thiscluster+'/HST/galaxy_subtracted_images/marc/'+suff
          if cluster eq 'macs0329' then clashpath = repstr(clashpath,'marc','marc_subt')
       endif
       
       if keyword_set(redshift) then clashpath = archivepath+thiscluster+'/redshifts/'
       if keyword_set(arcs) then clashpath = archivepath+thiscluster+$
         '/HST/PhotoZ/mosaicdrizzle_image_pipeline/IR_detection/html/'
          
       if file_test(clashpath,/dir) eq 0 then splog, 'Directory '+clashpath+' not found!'
    endif
    
return, clashpath
end
