function atlas_header_forage, allfits, sings=sings, atlas=atlas
; jm09dec28ucsd - grab info from the header
    
    info1 = {date: '', rawfile: '', galaxy: '', scanlen: 0.0, $
      posangle: 0.0, aperwid: 0.0, apercen: 0.0, pixscale: 0.0, $
      naxis: 0, atlas: keyword_set(atlas), sings: keyword_set(sings)}

    nfits = n_elements(allfits)
    for ii = 0, nfits-1 do begin
       hdr = headfits(allfits[ii])
       thisinfo = info1
       thisinfo.date = sxpar(hdr,'date-obs')
       if (strmatch(thisinfo.date,'*T*') eq 0) then $
         thisinfo.date = strtrim(thisinfo.date,2)+'T'+sxpar(hdr,'ut')
       thisinfo.galaxy = sxpar(hdr,'galaxy')
       thisinfo.scanlen = sxpar(hdr,'scanlen')
       thisinfo.posangle = sxpar(hdr,'posangle')
       thisinfo.apercen = sxpar(hdr,'apercen')
       thisinfo.aperwid = sxpar(hdr,'aperwid')
       thisinfo.pixscale = sxpar(hdr,'cd2_2')
       thisinfo.naxis = sxpar(hdr,'naxis1')
       ncomb = sxpar(hdr,'INCOMB')
       if (ncomb gt 0) then begin ; crsplit
          thisinfo = replicate(thisinfo,ncomb)
          for jj = 0, ncomb-1 do begin
             thisinfo[jj].rawfile = sxpar(hdr,'iraw'+string(jj+1,format='(I2.2)'))
          endfor
       endif else begin
          thisinfo.rawfile = sxpar(hdr,'iraw')
       endelse
; smash everything together
       if (n_elements(info) eq 0) then info = thisinfo else $
         info = [info,thisinfo]
    endfor

return, info
end

function atlas_match_info, plan, info, allinfo=allinfo, atlas=atlas, sings=sings
; now match the extraction info structure to the input plan file    

    nap = 2
    nobj = n_elements(plan)
    plan = struct_addtags(plan,replicate({match: -1},nobj))
    
    apinfo_template = {filename: '', galaxy: strarr(nap), hand_x: fltarr(nap), $
      hand_y: fltarr(nap), hand_fwhm: fltarr(nap), box_rad: fltarr(nap)}

    allfiles = strtrim(info.rawfile,2)
    allfiles = repstr(repstr(allfiles,'.fits',''),'.gz')
    for ii = 0, nobj-1 do begin
       thisfile = strtrim(plan[ii].filename,2)
       thisfile = repstr(repstr(thisfile,'.fits',''),'.gz')

       match = where(thisfile eq allfiles,nmatch)
       if (nmatch ne 0) then begin
          plan[ii].match = 1
;         struct_print, info[match], /no_head

; multiple matches can occur for two reasons: (1) the same spectrum
; was placed into both ATLAS *and* SINGS; if the apertures are the
; same (e.g., for a nuclear extraction, then just keep the first
; occurance); if the apertures are *different* (e.g., in SINGS we
; extracted 50", but in the ATLAS we extracted 90"), then keep the
; smaller extraction aperture of the two; or (2) one or more objects
; was was extracted from the same 2D spectrum, in which case we care
; about all occurances
          if (nmatch gt 1) then begin
             if (total(info[match].galaxy eq info[match[0]].galaxy) ne nmatch) then message, 'Hey!'

             daper = max(info[match].aperwid)-min(info[match].aperwid)
             if (daper eq 0.0) then begin ; keep the first one
                match = match[0]
                nmatch = 1
             endif else begin
                if (nmatch gt 1) then message, 'Multiple matches!'
; pick the largest aperture
                maxap = max(info[match].aperwid,keep)
;               minap = min(info[match].aperwid,keep)
;               keep = where((info[match].sings eq 1),nkeep)
                match = match[keep]
                nmatch = 1
                splog, 'Multiple matches with different apertures! Selecting:'
                struct_print, info[match], /no_head
                print
             endelse
          endif
          if (nmatch gt 1) then message, 'Multiple matches!'

          apinfo1 = apinfo_template
          apinfo1.filename = plan[ii].filename
          for jj = 0, nmatch-1 do begin
             apinfo1.galaxy[jj] = strtrim(info[match[jj]].galaxy,2)
             apinfo1.hand_x[jj] = info[match[jj]].apercen
             apinfo1.hand_y[jj] = info[match[jj]].naxis/2.0
             apinfo1.hand_fwhm[jj] = info[match[jj]].aperwid/info[match[jj]].pixscale
             apinfo1.box_rad[jj] = info[match[jj]].aperwid/info[match[jj]].pixscale/2.0
          endfor
          if (n_elements(apinfo) eq 0) then apinfo = apinfo1 else $
            apinfo = [apinfo,apinfo1]
          if (n_elements(allinfo) eq 0) then allinfo = info[match] else $
            allinfo = [allinfo,info[match]]
       endif 
    endfor
    
return, apinfo
end

function atlas_get_apertures, plan, sings=sings, atlas=atlas
; get the extraction apertures for each object

    common atlas_apertures, atlasinfo, singsinfo

; grab extraction info on all the objects
    if keyword_set(atlas) then begin
       if (n_elements(atlasinfo) eq 0) then begin
          atlasfits = [file_search(atlas_path(/atlas1d)+'*.fits'),$
            file_search(atlas_path(/atlas1d)+'repeaters/*.fits')]
          atlasinfo = atlas_header_forage(atlasfits,/atlas)
       endif
       info = atlasinfo
    endif
    if keyword_set(sings) then begin
       if (n_elements(singsinfo) eq 0) then begin
          singsfits = [file_search(sings_path(/spec1d)+'*.fits'),$
            file_search(sings_path(/spec1d)+'repeaters/*.fits')]
          singsinfo = atlas_header_forage(singsfits,/sings)
       endif
       info = singsinfo
    endif

; now match the PLAN against the master aperture list
    apinfo = atlas_match_info(plan,info,allinfo=allinfo)
    plan = plan[where((strtrim(plan.flavor,2) ne 'science') or (plan.match eq 1))]
    
;   newplan = struct_addtags(newplan,apinfo)
;   struct_print, newplan

return, apinfo
end
    
