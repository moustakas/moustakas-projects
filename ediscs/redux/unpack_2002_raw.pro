;+
; NAME:
;      UNPACK_2002_RAW
;
; PURPOSE:
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;      J. Moustakas, 2004 Jun 19-21, U of A
;-

pro unpack_2002_raw, wfits=wfits
; there is a lot of junky data in the RAW subdirectory.  we will not
; need chip2 images (no data) nor calibration data that are not part
; of our science program (see DELETE_CHIP2.PRO)

; distinguish between MOS and LSS observations    

    rootpath = '/d0/ioannis/ediscs/2002/'
    
    cwdpath = rootpath+'redux/'
    rawpath = rootpath+'raw/'
    outpath = rootpath+'data/'

    pushd, rawpath
    flist = file_search('*.fits',count=fcount)
    popd
    
    for i = 0L, fcount-1L do begin

       print, format='("Unpacking file ",I3,"/",I3,".",A1,$)', $
         i+1, fcount, string(13b)
       
       h = headfits(rawpath+flist[i])

; parse the FITS header into a data structure and then trim to a list
; of useful tag names

       hinfo = eso_header_forage(h)
       hinfo = create_struct({rawfile: flist[i]},hinfo)

       subhinfo = struct_trimtags(hinfo,select=['RAWFILE','DATE_OBS','EXPTIME','RA','DEC',$
         'EQUINOX','UTC','OBS_NAME','OBS_TARG_NAME','TPL_ID','TPL_NAME','DPR_CATG','DPR_TYPE',$
         'SEQ_SPEC_TARG','INS_MODE','INS_SLIT_NAME','INS_SLIT_RA','INS_SLIT_DEC','INS_SLIT_WID',$
         'INS_GRIS1_NAME','INS_GRIS1_DISP','INS_GRIS1_WLEN','ORIGFILE'])
       if (i eq 0L) then bighinfo = subhinfo else bighinfo = struct_append(bighinfo,subhinfo)

; ---------------------------------------------------------------------------       
; MOS observations
; ---------------------------------------------------------------------------       

       if (strtrim(subhinfo.ins_mode,2) eq 'MOS') then begin
          
; MOS: standard stars

          if (strtrim(subhinfo.dpr_type,2) eq 'STD') then begin

             if (n_elements(mos_objinfo) eq 0L) then mos_objinfo = subhinfo else $
               mos_objinfo = [ [mos_objinfo], [subhinfo] ]
             
          endif

; MOS: dome flats

          if (strtrim(subhinfo.dpr_type,2) eq 'FLAT,LAMP') then begin
             
             if (n_elements(mos_flatinfo) eq 0L) then mos_flatinfo = subhinfo else $
               mos_flatinfo = [ [mos_flatinfo], [subhinfo] ]
             
          endif

; MOS: arc lamps
          
          if (strtrim(subhinfo.dpr_type,2) eq 'WAVE,LAMP') then begin

             if (n_elements(mos_arcinfo) eq 0L) then mos_arcinfo = subhinfo else $
               mos_arcinfo = [ [mos_arcinfo], [subhinfo] ]
             
          endif
          
       endif 

; ---------------------------------------------------------------------------       
; LSS observations
; ---------------------------------------------------------------------------       

       if (strtrim(subhinfo.ins_mode,2) eq 'LSS') then begin
          
; LSS: hot stars

          if (strtrim(subhinfo.dpr_type,2) eq 'STD') then begin

             if (n_elements(lss_objinfo) eq 0L) then lss_objinfo = subhinfo else $
               lss_objinfo = [ [lss_objinfo], [subhinfo] ]
             
          endif

; LSS: dome flats

          if (strtrim(subhinfo.dpr_type,2) eq 'FLAT,LAMP') then begin
             
             if (n_elements(lss_flatinfo) eq 0L) then lss_flatinfo = subhinfo else $
               lss_flatinfo = [ [lss_flatinfo], [subhinfo] ]
             
          endif

; LSS: arc lamps
          
          if (strtrim(subhinfo.dpr_type,2) eq 'WAVE,LAMP') then begin

             if (n_elements(lss_arcinfo) eq 0L) then lss_arcinfo = subhinfo else $
               lss_arcinfo = [ [lss_arcinfo], [subhinfo] ]
             
          endif
          
       endif 

    endfor 
    print

; ---------------------------------------------------------------------------    
; assorted LSS variables
; ---------------------------------------------------------------------------    

    nlss_arc = n_elements(lss_arcinfo)
    nlss_flat = n_elements(lss_flatinfo)
    nlss_obj = n_elements(lss_objinfo)
    
    splog, 'NLSS_ARC  = '+string(nlss_arc,format='(I0)')
    splog, 'NLSS_FLAT = '+string(nlss_flat,format='(I0)')
    splog, 'NLSS_OBJ  = '+string(nlss_obj,format='(I0)')
    print

; ---------------------------------------------------------------------------    
; assorted calibration file variables
; ---------------------------------------------------------------------------    

    indx_mos_arc_left = where(strtrim(mos_arcinfo.seq_spec_targ,2) eq 'MOS_center-offset',nmos_arc_left)
    indx_mos_arc_center = where(strtrim(mos_arcinfo.seq_spec_targ,2) eq 'MOS_center',nmos_arc_center)
    indx_mos_arc_right = where(strtrim(mos_arcinfo.seq_spec_targ,2) eq 'MOS_center+offset',nmos_arc_right)

    indx_mos_flat_left = where(strtrim(mos_flatinfo.seq_spec_targ,2) eq 'MOS_center-offset',nmos_flat_left)
    indx_mos_flat_center = where(strtrim(mos_flatinfo.seq_spec_targ,2) eq 'MOS_center',nmos_flat_center)
    indx_mos_flat_right = where(strtrim(mos_flatinfo.seq_spec_targ,2) eq 'MOS_center+offset',nmos_flat_right)

    splog, 'NMOS_ARC_LEFT    = '+string(nmos_arc_left,format='(I0)')
    splog, 'NMOS_ARC_CENTER  = '+string(nmos_arc_center,format='(I0)')
    splog, 'NMOS_ARC_RIGHT   = '+string(nmos_arc_right,format='(I0)')
    print
    
    splog, 'NMOS_FLAT_LEFT   = '+string(nmos_flat_left,format='(I0)')
    splog, 'NMOS_FLAT_CENTER = '+string(nmos_flat_center,format='(I0)')
    splog, 'NMOS_FLAT_RIGHT  = '+string(nmos_flat_right,format='(I0)')
    print

    maximage = 10L
    splog, 'MAXIMAGE = '+string(maximage,format='(I0)')

    nlss_arc = nlss_arc   ; < maximage
    nlss_flat = nlss_flat   < maximage

    nmos_arc_left = nmos_arc_left     ; < maximage
    nmos_arc_center = nmos_arc_center ; < maximage
    nmos_arc_right = nmos_arc_right   ; < maximage
    
    nmos_flat_left = nmos_flat_left     < maximage
    nmos_flat_center = nmos_flat_center < maximage
    nmos_flat_right = nmos_flat_right   < maximage
    
; ---------------------------------------------------------------------------    
; assorted object file variables
; ---------------------------------------------------------------------------    

    indx_mos_obj_left = where(strtrim(mos_objinfo.seq_spec_targ,2) eq 'MOS_center-offset',nmos_obj_left)
    indx_mos_obj_center = where(strtrim(mos_objinfo.seq_spec_targ,2) eq 'MOS_center',nmos_obj_center)
    indx_mos_obj_right = where(strtrim(mos_objinfo.seq_spec_targ,2) eq 'MOS_center+offset',nmos_obj_right)

    splog, 'NMOS_OBJ_LEFT    = '+string(nmos_obj_left,format='(I0)')
    splog, 'NMOS_OBJ_CENTER  = '+string(nmos_obj_center,format='(I0)')
    splog, 'NMOS_OBJ_RIGHT   = '+string(nmos_obj_right,format='(I0)')
    print

; ---------------------------------------------------------------------------
; copy the bias frames for the hot stars and the standards to the
; output directory 
; ---------------------------------------------------------------------------

    splog, 'Updating bias file header [bias_lss.fits].'
    if (file_test(cwdpath+'bias_lss.fits',/regular) eq 0L) then begin
       splog, 'Bias file '+cwdpath+'bias_lss.fits not found.'
       return
    endif

    h = headfits(cwdpath+'bias_lss.fits')
    hinfo = eso_header_forage(h)
    
    newh = update_eso_2002_header(h,hinfo,/bias)
    
    if keyword_set(wfits) then begin

       splog, 'Writing '+outpath+'a.bias_lss.fits'+'.'
       spawn, ['/bin/cp -f '+cwdpath+'bias_lss.fits '+outpath+'a.bias_lss.fits'], /sh
       djs_modfits, outpath+'a.bias_lss.fits', 0, newh

    endif

; --------------------    
    
    splog, 'Updating bias file header [bias_mos.fits].'
    if (file_test(cwdpath+'bias_mos.fits',/regular) eq 0L) then begin
       splog, 'Bias file '+cwdpath+'bias_mos.fits not found.'
       return
    endif

    h = headfits(cwdpath+'bias_mos.fits')
    hinfo = eso_header_forage(h)
    
    newh = update_eso_2002_header(h,hinfo,/bias)
    
    if keyword_set(wfits) then begin

       splog, 'Writing '+outpath+'a.bias_mos.fits'+'.'
       spawn, ['/bin/cp -f '+cwdpath+'bias_mos.fits '+outpath+'a.bias_mos.fits'], /sh
       djs_modfits, outpath+'a.bias_mos.fits', 0, newh

    endif
    
; ---------------------------------------------------------------------------
; generate a co-added LSS_ARC image
; ---------------------------------------------------------------------------

    splog, 'Generating a mean LSS arc lamp from '+$
      string(nlss_arc,format='(i0)')+' images.'

    if keyword_set(wfits) then begin

       lss_arc = mrdfits(rawpath+lss_arcinfo[0].rawfile,0,h,/silent,/fscale)
       for i = 1L, nlss_arc-1L do lss_arc = [ [ [lss_arc] ], $
         [ [mrdfits(rawpath+lss_arcinfo[i].rawfile,0,h,/silent,/fscale)] ] ]

       lss_arc_avg = lss_arc ; <-- NOTE!
;      lss_arc_avg = total(lss_arc,3)
;      lss_arc_avg = djs_avsigclip(lss_arc,3,sigrej=3.0,outmask=lss_arcmask)
       sxaddpar, h, 'NCOMBINE', nlss_arc, before='HISTORY'

       newh = update_eso_2002_header(h,lss_arcinfo[0],/arclamp)
       outfile = 'a.lss_arclamp.fits'
       splog, 'Writing '+outpath+outfile+'.'
       
       mwrfits, lss_arc_avg, outpath+outfile, newh, /create

    endif

; ---------------------------------------------------------------------------
; generate a sigma-clipped average LSS_FLAT image
; ---------------------------------------------------------------------------
    
    splog, 'Generating a mean LSS dome flat from '+$
      string(nlss_flat,format='(i0)')+' images.'
    
    if keyword_set(wfits) then begin

       lss_flat = mrdfits(rawpath+lss_flatinfo[0].rawfile,0,h,/silent,/fscale)
       for i = 1L, nlss_flat-1L do lss_flat = [ [ [lss_flat] ], $
         [ [mrdfits(rawpath+lss_flatinfo[i].rawfile,0,h,/silent,/fscale)] ] ]

       lss_flat_avg = djs_avsigclip(lss_flat,3,sigrej=3.0,outmask=lss_flatmask)
       sxaddpar, h, 'NCOMBINE', nlss_flat, before='HISTORY'

       newh = update_eso_2002_header(h,lss_flatinfo[0],/domeflat)
       outfile = 'a.lss_domeflat.fits'
       splog, 'Writing '+outpath+outfile+'.'
       
       mwrfits, lss_flat_avg, outpath+outfile, newh, /create

    endif
    
; ---------------------------------------------------------------------------
; generate a sigma-clipped average MOS_ARC_LEFT image
; ---------------------------------------------------------------------------

    splog, 'Generating a mean MOS LEFT arc lamp from '+$
      string(nmos_arc_left,format='(i0)')+' images.'
    
    if keyword_set(wfits) then begin

       mos_arc_left = mrdfits(rawpath+mos_arcinfo[indx_mos_arc_left[0]].rawfile,0,h,/silent,/fscale)
       for i = 1L, nmos_arc_left-1L do mos_arc_left = [ [ [mos_arc_left] ], $
         [ [mrdfits(rawpath+mos_arcinfo[indx_mos_arc_left[i]].rawfile,0,h,/silent,/fscale)] ] ]

       mos_arc_left_avg = mos_arc_left
;      mos_arc_left_avg = total(mos_arc_left,3)
;      mos_arc_left_avg = djs_avsigclip(mos_arc_left,3,sigrej=3.0,outmask=mos_arc_leftmask)
       sxaddpar, h, 'NCOMBINE', nmos_arc_left, before='HISTORY'

       newh = update_eso_2002_header(h,mos_arcinfo[indx_mos_arc_left[0]],/arclamp,/left)
       outfile = 'a.mos_arclamp_left.fits'
       splog, 'Writing '+outpath+outfile+'.'
       
       mwrfits, mos_arc_left_avg, outpath+outfile, newh, /create

    endif
    
; ---------------------------------------------------------------------------
; generate a sigma-clipped average MOS_ARC_CENTER image
; ---------------------------------------------------------------------------

    splog, 'Generating a mean MOS CENTER arc lamp from '+$
      string(nmos_arc_center,format='(i0)')+' images.'
    
    if keyword_set(wfits) then begin

       mos_arc_center = mrdfits(rawpath+mos_arcinfo[indx_mos_arc_center[0]].rawfile,0,h,/silent,/fscale)
       for i = 1L, nmos_arc_center-1L do mos_arc_center = [ [ [mos_arc_center] ], $
         [ [mrdfits(rawpath+mos_arcinfo[indx_mos_arc_center[i]].rawfile,0,h,/silent,/fscale)] ] ]

       mos_arc_center_avg = mos_arc_center
;      mos_arc_center_avg = total(mos_arc_center,3)
;      mos_arc_center_avg = djs_avsigclip(mos_arc_center,3,sigrej=3.0,outmask=mos_arc_centermask)
       sxaddpar, h, 'NCOMBINE', nmos_arc_center, before='HISTORY'

       newh = update_eso_2002_header(h,mos_arcinfo[indx_mos_arc_center[0]],/arclamp,/center)
       outfile = 'a.mos_arclamp_center.fits'
       splog, 'Writing '+outpath+outfile+'.'
       
       mwrfits, mos_arc_center_avg, outpath+outfile, newh, /create

    endif

; ---------------------------------------------------------------------------
; generate a sigma-clipped average MOS_ARC_RIGHT image
; ---------------------------------------------------------------------------

    splog, 'Generating a mean MOS RIGHT arc lamp from '+$
      string(nmos_arc_right,format='(i0)')+' images.'
    
    if keyword_set(wfits) then begin

       mos_arc_right = mrdfits(rawpath+mos_arcinfo[indx_mos_arc_right[0]].rawfile,0,h,/silent,/fscale)
       for i = 1L, nmos_arc_right-1L do mos_arc_right = [ [ [mos_arc_right] ], $
         [ [mrdfits(rawpath+mos_arcinfo[indx_mos_arc_right[i]].rawfile,0,h,/silent,/fscale)] ] ]

       mos_arc_right_avg = mos_arc_right
;      mos_arc_right_avg = total(mos_arc_right,3)
;      mos_arc_right_avg = djs_avsigclip(mos_arc_right,3,sigrej=3.0,outmask=mos_arc_rightmask)
       sxaddpar, h, 'NCOMBINE', nmos_arc_right, before='HISTORY'

       newh = update_eso_2002_header(h,mos_arcinfo[indx_mos_arc_right[0]],/arclamp,/right)
       outfile = 'a.mos_arclamp_right.fits'
       splog, 'Writing '+outpath+outfile+'.'
       
       mwrfits, mos_arc_right_avg, outpath+outfile, newh, /create

    endif

; ---------------------------------------------------------------------------
; generate a sigma-clipped average MOS_FLAT_LEFT image
; ---------------------------------------------------------------------------

    splog, 'Generating a mean MOS LEFT flat lamp from '+$
      string(nmos_flat_left,format='(i0)')+' images.'
    
    if keyword_set(wfits) then begin

       mos_flat_left = mrdfits(rawpath+mos_flatinfo[indx_mos_flat_left[0]].rawfile,0,h,/silent,/fscale)
       for i = 1L, nmos_flat_left-1L do mos_flat_left = [ [ [mos_flat_left] ], $
         [ [mrdfits(rawpath+mos_flatinfo[indx_mos_flat_left[i]].rawfile,0,h,/silent,/fscale)] ] ]

       mos_flat_left_avg = djs_avsigclip(mos_flat_left,3,sigrej=3.0,outmask=mos_flat_leftmask)
       sxaddpar, h, 'NCOMBINE', nmos_flat_left, before='HISTORY'

       newh = update_eso_2002_header(h,mos_flatinfo[indx_mos_flat_left[0]],/domeflat,/left)
       outfile = 'a.mos_domeflat_left.fits'
       splog, 'Writing '+outpath+outfile+'.'
       
       mwrfits, mos_flat_left_avg, outpath+outfile, newh, /create

    endif
    
; ---------------------------------------------------------------------------
; generate a sigma-clipped average MOS_FLAT_CENTER image
; ---------------------------------------------------------------------------

    splog, 'Generating a mean MOS CENTER flat lamp from '+$
      string(nmos_flat_center,format='(i0)')+' images.'
    
    if keyword_set(wfits) then begin

       mos_flat_center = mrdfits(rawpath+mos_flatinfo[indx_mos_flat_center[0]].rawfile,0,h,/silent,/fscale)
       for i = 1L, nmos_flat_center-1L do mos_flat_center = [ [ [mos_flat_center] ], $
         [ [mrdfits(rawpath+mos_flatinfo[indx_mos_flat_center[i]].rawfile,0,h,/silent,/fscale)] ] ]

       mos_flat_center_avg = djs_avsigclip(mos_flat_center,3,sigrej=3.0,outmask=mos_flat_centermask)
       sxaddpar, h, 'NCOMBINE', nmos_flat_center, before='HISTORY'

       newh = update_eso_2002_header(h,mos_flatinfo[indx_mos_flat_center[0]],/domeflat,/center)
       outfile = 'a.mos_domeflat_center.fits'
       splog, 'Writing '+outpath+outfile+'.'
       
       mwrfits, mos_flat_center_avg, outpath+outfile, newh, /create

    endif
    
; ---------------------------------------------------------------------------
; generate a sigma-clipped average MOS_FLAT_RIGHT image
; ---------------------------------------------------------------------------

    splog, 'Generating a mean MOS RIGHT flat lamp from '+$
      string(nmos_flat_right,format='(i0)')+' images.'
    
    if keyword_set(wfits) then begin

       mos_flat_right = mrdfits(rawpath+mos_flatinfo[indx_mos_flat_right[0]].rawfile,0,h,/silent,/fscale)
       for i = 1L, nmos_flat_right-1L do mos_flat_right = [ [ [mos_flat_right] ], $
         [ [mrdfits(rawpath+mos_flatinfo[indx_mos_flat_right[i]].rawfile,0,h,/silent,/fscale)] ] ]

       mos_flat_right_avg = djs_avsigclip(mos_flat_right,3,sigrej=3.0,outmask=mos_flat_rightmask)
       sxaddpar, h, 'NCOMBINE', nmos_flat_right, before='HISTORY'

       newh = update_eso_2002_header(h,mos_flatinfo[indx_mos_flat_right[0]],/domeflat,/right)
       outfile = 'a.mos_domeflat_right.fits'
       splog, 'Writing '+outpath+outfile+'.'
       
       mwrfits, mos_flat_right_avg, outpath+outfile, newh, /create

    endif
    
; ---------------------------------------------------------------------------
; LSS hot stars
; ---------------------------------------------------------------------------

    if keyword_set(wfits) then begin

       splog, 'Writing LSS hot stars.'

       for i = 0L, nlss_obj-1L do begin
 
          image = mrdfits(rawpath+lss_objinfo[i].rawfile,0,h,/silent,/fscale)
          newh = update_eso_2002_header(h,lss_objinfo[i],/object)

          outname = 'a.lss_'+strlowcase(lss_objinfo[i].obs_targ_name)
          files = file_search(outpath+outname+'*.fits',count=fcount)
 
          if (fcount eq 1L) then begin
             spawn, ['\mv '+files[0]+' '+outpath+outname+'_1.fits'], /sh
             outname = outname+'_2'
          endif
          
          if (fcount gt 1L) then outname = outname+'_'+string(fcount+1,format='(I0)')
 
          splog, 'Writing '+outpath+outname+'.fits.'
          mwrfits, image, outpath+outname+'.fits', newh, /create
 
       endfor 
    
    endif

; ---------------------------------------------------------------------------
; MOS_LEFT standard stars
; ---------------------------------------------------------------------------

    if keyword_set(wfits) then begin

       splog, 'Writing MOS LEFT standard stars.'

       for i = 0L, nmos_obj_left-1L do begin
 
          image = mrdfits(rawpath+mos_objinfo[indx_mos_obj_left[i]].rawfile,0,h,/silent,/fscale)
          newh = update_eso_2002_header(h,mos_objinfo[indx_mos_obj_left[i]],/object,/left)

          outname = 'a.mos_'+strlowcase(strtrim(mos_objinfo[indx_mos_obj_left[i]].obs_targ_name,2))+'_left'
          files = file_search(outpath+outname+'*.fits',count=fcount)
 
          if (fcount eq 1L) then begin
             spawn, ['\mv '+files[0]+' '+outpath+outname+'_1.fits'], /sh
             outname = outname+'_2'
          endif
          
          if (fcount gt 1L) then outname = outname+'_'+string(fcount+1,format='(I0)')
 
          splog, 'Writing '+outpath+outname+'.fits.'
          mwrfits, image, outpath+outname+'.fits', newh, /create
 
       endfor 
    
    endif

; ---------------------------------------------------------------------------
; MOS_CENTER standard stars
; ---------------------------------------------------------------------------

    if keyword_set(wfits) then begin

       splog, 'Writing MOS CENTER standard stars.'

       for i = 0L, nmos_obj_center-1L do begin
 
          image = mrdfits(rawpath+mos_objinfo[indx_mos_obj_center[i]].rawfile,0,h,/silent,/fscale)
          newh = update_eso_2002_header(h,mos_objinfo[indx_mos_obj_center[i]],/object,/center)

          outname = 'a.mos_'+strlowcase(strtrim(mos_objinfo[indx_mos_obj_center[i]].obs_targ_name,2))+'_center'
          files = file_search(outpath+outname+'*.fits',count=fcount)
 
          if (fcount eq 1L) then begin
             spawn, ['\mv '+files[0]+' '+outpath+outname+'_1.fits'], /sh
             outname = outname+'_2'
          endif
          
          if (fcount gt 1L) then outname = outname+'_'+string(fcount+1,format='(I0)')
 
          splog, 'Writing '+outpath+outname+'.fits.'
          mwrfits, image, outpath+outname+'.fits', newh, /create
 
       endfor 
    
    endif

; ---------------------------------------------------------------------------
; MOS_RIGHT standard stars
; ---------------------------------------------------------------------------

    if keyword_set(wfits) then begin

       splog, 'Writing MOS RIGHT standard stars.'

       for i = 0L, nmos_obj_right-1L do begin
 
          image = mrdfits(rawpath+mos_objinfo[indx_mos_obj_right[i]].rawfile,0,h,/silent,/fscale)
          newh = update_eso_2002_header(h,mos_objinfo[indx_mos_obj_right[i]],/object,/right)

          outname = 'a.mos_'+strlowcase(strtrim(mos_objinfo[indx_mos_obj_right[i]].obs_targ_name,2))+'_right'
          files = file_search(outpath+outname+'*.fits',count=fcount)
 
          if (fcount eq 1L) then begin
             spawn, ['\mv '+files[0]+' '+outpath+outname+'_1.fits'], /sh
             outname = outname+'_2'
          endif
          
          if (fcount gt 1L) then outname = outname+'_'+string(fcount+1,format='(I0)')
 
          splog, 'Writing '+outpath+outname+'.fits.'
          mwrfits, image, outpath+outname+'.fits', newh, /create
 
       endfor 
    
    endif

stop        
    
return
end
