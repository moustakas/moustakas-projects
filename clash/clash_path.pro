function clash_path, cluster, catalogs=catalogs, redshift=redshift, $
  ir=ir, arcs=arcs, isedfit=isedfit, montegrids=montegrids, bcgimf=bcgimf, $
  macs0329_z6arcs=macs0329_z6arcs, santorini=santorini
; jm11apr18ucsd - 

    clashpath = getenv('CLASH_DATA')+'/'
    if keyword_set(catalogs) then clashpath = clashpath+'catalogs/'
    if keyword_set(isedfit) then clashpath = clashpath+'isedfit/'
    if keyword_set(montegrids) then clashpath = clashpath+'montegrids/'
    if keyword_set(bcgimf) then clashpath = clashpath+'projects/bcgimf/'
    if keyword_set(macs0329_z6arcs) then clashpath = clashpath+'projects/macs0329_z6arcs/'
    if keyword_set(santorini) then clashpath = clashpath+'projects/santorini/'

; if CLUSTER has been given, then build the relevant ARCHIVE path;
; assume the user wants the 'catalogs' directory by default; deal with
; the possibility that the cluster name has an underscore (or not)
    if n_elements(cluster) ne 0 then begin
       archivepath = clashpath+'archive/'
       thiscluster = strlowcase(cluster)
       if file_test(archivepath+thiscluster,/dir) eq 0 then begin
          thiscluster = repstr(strlowcase(cluster),'_','')
          if file_test(archivepath+thiscluster,/dir) eq 0 then begin
             thiscluster = repstr(strlowcase(cluster),'abell','abell_')
          endif
       endif

       clashpath = archivepath+thiscluster+'/HST/catalogs/mosaicdrizzle_image_pipeline/'
       if keyword_set(ir) then clashpath = clashpath+'IR_detection/' else $
         clashpath = clashpath+'ACS_IR_detection/'

       if keyword_set(redshift) then clashpath = archivepath+thiscluster+'/redshifts/'
       if keyword_set(arcs) then clashpath = archivepath+thiscluster+'/HST/PhotoZ/mosaicdrizzle_image_pipeline/IR_detection/html/'
       if file_test(clashpath,/dir) eq 0 then message, 'Directory '+clashpath+' not found!'
    endif
    
return, clashpath
end
