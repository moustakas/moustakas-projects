;+
; NAME:
;   dchildren_atlas
; PURPOSE:
;   deblend children in NASA-Sloan Atlas imagees
; CALLING SEQUENCE:
;   dchildren_atlas [, /noclobber ]
; COMMENTS:
;   Requires dparents_atlas, dpsf_atls, dstargal_atlas to be run
;   Uses file:
;     atlases/#/[base]-nimage-#.fits - star-subtracted images
;   which has just the galaxies in it. Gets list of galaxies from:
;     atlases/#/[base]-#-sgset.fits - locations of stars, galaxies
;   Writes files:
;     atlases/#/[base]-acat-#.fits - catalog of children
;     atlases/#/[base]-parent-#.fits - parent images
;     atlases/#/[base]-ivar-#.fits - inverse var images
;     atlases/#/[base]-#-atlas-#.fits - atlas images
;     atlases/#/[base]-#-templates-#.fits - template images
; BUGS:
;   Current version ONLY works if template is used from 
;   same band every time.
; REVISION HISTORY:
;   11-Jan-2006  Written by Blanton, NYU
;-
;------------------------------------------------------------------------------
pro streams_children, base, sersic=sersic, noclobber=noclobber, debug=debug

    subdir='atlases'

; default to use base name same as directory name
    if n_elements(base) eq 0 then begin
       spawn, 'pwd', cwd
       base=(file_basename(cwd))[0]
    endif

; read in pset
    pset= gz_mrdfits(base+'-pset.fits',1)
    tref = pset.ref
    imfiles=strtrim(pset.imfiles,2)
    tuse= replicate(tref, n_elements(imfiles))
    puse=pset.puse
    nim= n_elements(imfiles)
    nx=lonarr(nim)
    ny=lonarr(nim)
    nimages=ptrarr(nim)

; loop on each parent
    pcat=gz_mrdfits(base+'-pcat.fits',1)

    splog, 'HACK!!'
;   for iparent = 55, 55 do begin
    for iparent = 0L, n_elements(pcat)-1 do begin
       splog, 'Parent ', iparent
       
; read in star and galaxy locations
       sgsetfile=subdir+'/'+strtrim(string(iparent),2)+'/'+base+'-'+ $
         strtrim(string(iparent),2)+'-sgset.fits'
       sgset= gz_mrdfits(sgsetfile, 1)

       if (sgset.ngals eq 0) then continue
       
; file for input images
       nimfile=subdir+'/'+strtrim(string(iparent),2)+'/'+base+ $
         '-nimage-'+strtrim(string(iparent),2)+'.fits'
    
; file for input ivars
       pfile='parents'+'/'+base+'-parent-'+strtrim(string(iparent),2)+'.fits'
    
; if result exists, skip out
       acatfile=subdir+'/'+strtrim(string(iparent),2)+ $
         '/'+base+'-acat-'+strtrim(string(iparent),2)+ $
         '.fits'
       if(gz_file_test(acatfile) gt 0 AND $
         keyword_set(noclobber) gt 0) then continue

; create acat structure to store children in
       acat=replicate({pid:iparent, $
         tuse:tuse, $
         aid:-1L, $
         racen:0.D, $
         deccen:0.D, $
         type:0, $
         bgood:lonarr(nim), $
         good:0L}, sgset.ngals)
       acat.aid=lindgen(sgset.ngals)
       acat.racen=sgset.ra_gals[0:sgset.ngals-1]
       acat.deccen=sgset.dec_gals[0:sgset.ngals-1]

       if keyword_set(debug) then begin
          window, 0
          window, 1
       endif
       
; loop over filters
       templates= ptrarr(nim)
       hdrs= ptrarr(nim)
       tnx= lonarr(nim)
       tny= lonarr(nim)
       for k=0L, nim-1L do begin
          if(k eq 0) then first=1 else first=0
          kuse=tuse[k]
          
          if(keyword_set(templates[kuse]) eq 0) then begin
             
; read in image to use for templates
             if file_test(nimfile+'*') eq 0 then stop
             timage= gz_mrdfits(nimfile, kuse, thdr,/silent)
             adxy, thdr, acat.racen, acat.deccen, xgals, ygals

;            dobjects, timage, obj=tobj, plim=10.0, nlevel=1.0
;            cgimage, tobj, /keep;, stretch=2

; make galaxy templates
             splog, 'Making basic templates ...'
             dtemplates, timage, xgals, ygals, templates=curr_templates, $
               parallel=parallel, sersic=sersic, ikept=ikept
;              parallel=0.6, sersic=sersic, ikept=ikept

; debugging plot
             if keyword_set(debug) then begin
                wset, 0
                cgimage, timage, /keep_aspect, /save, stretch=2, /negative
                djs_oplot, xgals, ygals, psym=8, color=cgcolor('orange')
                djs_oplot, xgals[ikept], ygals[ikept], psym=symcat(9,thick=1), $
                  symsize=1.5, color=cgcolor('red')
                cc = get_kbrd(1)
             endif

             acat=acat[ikept]
             xgals= xgals[ikept]
             ygals= ygals[ikept]
             
;            splog, 'Smoothing templates ...'
             for i=0L, n_elements(acat)-1L do begin
                tmp_template= curr_templates[*,*,i]
                sig= dsigma(tmp_template, sp=10)
                
; smooth out low signficance areas
                dobjects, tmp_template, obj=tobjects, plim=10.
                stemplates= dsmooth(tmp_template, 5.)
                cobject= tobjects[long(xgals[i]), long(ygals[i])]
                iout=where(tobjects ne cobject, nout)
                if(nout gt 0) then $
                  tmp_template[iout]= stemplates[iout]
                
; clip VERY negative areas
                ilow=where(tmp_template lt -3.*sig, nlow)
                if(nlow gt 0) then $
                  tmp_template[ilow]=-3.*sig
                
                curr_templates[*,*,i]=tmp_template
             endfor 
             
             templates[kuse]=ptr_new(curr_templates)
             hdrs[kuse]=ptr_new(thdr)
             tnx[kuse]= (size(curr_templates, /dim))[0]
             tny[kuse]= (size(curr_templates, /dim))[1]
          endif

          if keyword_set(debug) then begin
             for bb = 0, n_elements(ikept)-1 do begin
                wset, 0
                plots, xgals[bb], ygals[bb], psym=symcat(9,thick=1), $
                  symsize=4.0, color=cgcolor('green')
                
                wset, 1
                cgimage, curr_templates[*,*,bb], /save, /keep_aspect, stretch=2, /negative
                plots, xgals[bb], ygals[bb], psym=7, $
                  symsize=2, color=cgcolor('red')
                cc = get_kbrd(1)
             endfor
          endif
          
 ; read in current image for galaxies
          cimage= gz_mrdfits(nimfile, k, hdr,/silent)
          civar= gz_mrdfits(pfile, k*2L+1,/silent)
          nx=(size(cimage, /dim))[0]
          ny=(size(cimage, /dim))[1]
          
          if(nx ne tnx[kuse] OR $
            ny ne tny[kuse]) then begin
             splog, 'Mapping templates ...'
             extast, hdr, k_ast
             extast, *hdrs[kuse], kuse_ast
             ctemplates=fltarr(nx, ny, n_elements(acat))
             for i=0L, n_elements(acat)-1L do begin
                tmp_template=fltarr(nx, ny)
                smosaic_remap, (*templates[kuse])[*,*,i], kuse_ast, k_ast, $
                  refimage=tmp_template
                ctemplates[*,*,i]=tmp_template
             endfor
          endif else begin
             ctemplates=(*templates[kuse])
          endelse

;         splog, 'Finding weights ...'
          dweights, cimage, civar, ctemplates, weights=weights, /nonneg
          
;         splog, 'Finding fluxes ...'
          dfluxes, cimage, ctemplates, weights, xgals, ygals, children=children

;         splog, 'Outputting results ...'
          for i=0L, n_elements(acat)-1L do begin
             aid=acat[i].aid
             if(total(children[*,*,i]) gt 0) then acat[i].bgood[k]= 1
             afile= subdir+'/'+ strtrim(string(iparent),2)+ $
               '/'+base+'-'+strtrim(string(iparent),2)+ $
               '-atlas-'+strtrim(string(aid),2)+'.fits'
             mwrfits, children[*,*,i], afile, hdr, create=first, /silent
             tfile= subdir+'/'+ strtrim(string(iparent),2)+ $
               '/'+base+'-'+strtrim(string(iparent),2)+ $
               '-templates-'+strtrim(string(aid),2)+'.fits'
             mwrfits, ctemplates[*,*,i], tfile, hdr, create=first, /silent
          endfor
          
          pbase=base+'-parent-'+strtrim(string(iparent),2)+'.fits'
          pfile= 'parents/'+pbase
          opfile= subdir+'/'+ strtrim(string(iparent),2)+ $
            '/'+pbase
          pim= gz_mrdfits(pfile, 2*k, phdr,/silent)
          mwrfits, pim, opfile, phdr, create=first, /silent
          oifile= subdir+'/'+ strtrim(string(iparent),2)+ $
            '/'+base+'-ivar-'+strtrim(string(iparent),2)+'.fits'
          mwrfits, civar, oifile, phdr, create=first, /silent
       endfor

       if(n_tags(acat) gt 0) then begin
          acat.good= total(acat.bgood, 1) gt 0
          dhdr= dimage_hdr()
          mwrfits, acat, acatfile, dhdr, /create
       endif
    endfor
    
; free memory
    heap_free, nimages
    heap_free, templates
    heap_free, hdrs
    
return 
end
;------------------------------------------------------------------------------
