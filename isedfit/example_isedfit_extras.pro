pro example_isedfit_extras, build_parent=build_parent, tar_ssps=tar_ssps, $
  sync_webpage=sync_webpage
; jm13aug13siena - do some extra stuff for the isedfit_example like
; building the parent sample and making some plots
    
    prefix = 'example'
    htmlpath = getenv('IM_PROJECTS_DIR')+'/isedfit/html/'
    isedfit_dir = htmlpath+'example/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

; ---------------------------------------------------------------------------
; copy the webpage over to SoS
    if keyword_set(sync_webpage) then begin
       profile = getenv('IDL_PROJECTS_DIR')+'/isedfit/example_isedfit.pro'
       
       pushd, htmlpath
       spawn, 'rsync -auvn --delete * sos:"public_html/isedfit/"', /sh
       spawn, 'rsync -auvn '+profile+' sos:"public_html/isedfit/example/"', /sh
       popd

       splog, 'Did that look OK [Y/N]?'
       cc = get_kbrd(1)
       if strupcase(cc) eq 'Y' then begin
          pushd, getenv('IM_PROJECTS_DIR')+'/isedfit/html/'
          spawn, 'rsync -auv --progress --delete * sos:"public_html/isedfit/"', /sh
          spawn, 'rsync -auv --progress '+profile+' sos:"public_html/isedfit/example/"', /sh
          popd
       endif
    endif

; ---------------------------------------------------------------------------
; make a tarball of the SSPs
    if keyword_set(tar_ssps) then begin
       ver = 'v1.0'
       tarfile = htmlpath+'isedfit_ssp_dir_'+ver+'.tar.gz'
       pushd, getenv('IM_DATA_DIR')
       spawn, 'tar czvf '+tarfile+' isedfit_ssp', /sh
       popd
    endif

; ---------------------------------------------------------------------------
; build the parent sample for the iSEDfit example    
    if keyword_set(build_parent) then begin
       sample = 'dr72' & letter = 'bsafe' & poststr = '0'
       post = read_vagc_garching(sample=sample,$
         letter=letter,poststr=poststr,/postlss)
       
       keep = where((post.z gt 0.05) and (post.z lt 0.2),ngal)
       keep = keep[0:9999]
       ngal = n_elements(keep)
       post = post[keep]
       
       vagcpath = getenv('VAGC_REDUX')+'/'
       sdssphot = mrdfits(vagcpath+'object_sdss_imaging.fits.gz',$
         row=post.object_position,1)
;      sdssspec = mrdfits(vagcpath+'object_sdss_spectro.fits.gz',$
;        row=post.object_position,1)
       galex = mrdfits(vagcpath+'object_galex_gr6.fits.gz',$
         row=post.object_position,1)
       wise = mrdfits(vagcpath+'object_wise.fits.gz',$
         row=post.object_position,1)
;      if (n_elements(sdsstwomass) eq 0L) then sdsstwomass = $
;        mrdfits(vagcpath+'object_twomass.fits.gz',$
;        row=phot.object_position,1)

       sdss_to_maggies, maggies, ivarmaggies, calib=sdssphot, flux='cmodel'
       im_galex_to_maggies, galex, gmaggies, givarmaggies
       wise_to_maggies, wise, wmaggies, wivarmaggies, /mpro

       parent = struct_addtags(sdssphot,struct_trimtags(post,$
         select=['z','absm','object_position']))
       parent = struct_addtags(temporary(parent),replicate($
         {maggies: fltarr(9), ivarmaggies: fltarr(9)},ngal))
       parent.maggies = [gmaggies,maggies,wmaggies[0:1,*]]
       parent.ivarmaggies = [givarmaggies,ivarmaggies,wivarmaggies[0:1,*]]
    
       im_mwrfits, parent, isedfit_dir+'exampledata.fits', /clobber
    endif

return
end
