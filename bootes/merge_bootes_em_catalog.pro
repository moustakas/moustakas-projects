;+
; NAME:
;   merge_bootes_em_catalog
; PURPOSE:
;   Merge the final, output EMphot photometric catalog.
; CALLING SEQUENCE:
;   merge_bootes_em_catalog
; INPUTS: 
; KEYWORDS:
;   unpack - unpack the original [N,F]UV catalogs, resolve duplicates,
;     and write out merged catalogs
;   combine - combine the NUV and FUV catalogs (where available) into
;     a single catalog (may also add additional info from the nominal
;     GALEX pipeline) 
; OUTPUTS: 
;   Output catalogs are written to:
;   ${IM_ARCHIVE_DIR}/bootes/galex/
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Aug 04, UCSD
;   jm10aug26ucsd - updated to the most recent EM catalogs
;-

pro merge_bootes_em_catalog, gr=gr, unpack=unpack, combine=combine, clobber=clobber

    if (n_elements(gr) eq 0) then gr = 'gr6'
    galexpath = getenv('IM_ARCHIVE_DIR')+'/bootes/galex/'

    field = 'NGPDWS'
    
; --------------------------------------------------
; unpack the raw catalogs for each tile
    if keyword_set(unpack) then begin
       xsize = 3840 ; GALEX image size
       ysize = 3840
       pixscale = 1.5 ; [arcsex/pixel] (see Morrissey+07)

       outfile_nuv = galexpath+strlowcase(field)+'_emphot_nuv.fits'
       outfile_fuv = repstr(outfile_nuv,'nuv','fuv')
; parse each tile individually
       alltiles = file_search(galexpath+field+'/*',count=ntile,/test_dir)
       splog, 'Found '+strtrim(ntile,2)+' tile(s)'
       niceprint, alltiles & print
       delvarx, nuv, fuv
       for jj = 0, ntile-1 do begin
          root = file_basename(alltiles[jj])
          splog, 'Parsing tile '+root
; NUV: each file is row-matched to the input prior catalog, so just
; keep observations with non-negative identification numbers; and add
; some tags we'll use below
          nuvfile = alltiles[jj]+'/EM/'+root+'-nd_v1-mix.fits'
          nuv1 = mrdfits(nuvfile,1)
          ngood = where(nuv1.ident ne -9999,nnuv1)
          nuv1 = nuv1[ngood]
          nuv1 = struct_addtags(replicate({tile: root, radius: 0.0, $
            ndup: 0, isbest_snr: 1, isbest_radius: 1},nnuv1),$
            temporary(nuv1))
          if (jj eq 0) then nuv = nuv1 else $
            nuv = [temporary(nuv),nuv1]
; FUV, if it exists
          fuvfile = repstr(nuvfile,'-nd','-fd')
          if (file_test(fuvfile) eq 0) then begin
             splog, 'No FUV photometry'
          endif else begin
             fuv1 = mrdfits(fuvfile,1)
             fgood = where(fuv1.ident ne -9999,nfuv1)
             fuv1 = fuv1[fgood]
             fuv1 = struct_addtags(replicate({tile: root, radius: 0.0, $
               ndup: 0, isbest_snr: 1, isbest_radius: 1},nfuv1),$
               temporary(fuv1))
             if (jj eq 0) then fuv = fuv1 else $
               fuv = [temporary(fuv),fuv1]
          endelse 
       endfor 
; compute the radius from the center of the FOV
       if (total((nuv.xc_em lt 0.0) or (nuv.xc_em gt xsize-1) or $
         (nuv.yc_em lt 0.0) or (nuv.yc_em gt ysize-1)) ne 0.0) then $
           message, 'This should not be possible!'
       nuv.radius = sqrt((nuv.xc_em-xsize/2)^2+(nuv.yc_em-ysize/2)^2)*1.5/3600.0 ; [deg]
       if (n_elements(fuv) ne 0L) then begin
          if (total((fuv.xc_em lt 0.0) or (fuv.xc_em gt xsize-1) or $
            (fuv.yc_em lt 0.0) or (fuv.yc_em gt ysize-1)) ne 0.0) then $
              message, 'This should not be possible!'
          fuv.radius = sqrt((fuv.xc_em-xsize/2)^2+(fuv.yc_em-ysize/2)^2)*1.5/3600.0 ; [deg]
       endif
; histogram the identification number to identify duplicates (due to
; overlapping tiles) and to pick the best observation; we do this
; instead of spheregroup because the galex positions are based on the
; input prior list
       hh = histogram(nuv.ident,bin=1,reverse=rr)
       dup = where(hh gt 1,ndup)
       for idup = 0L, ndup-1 do begin
          these = rr[rr[dup[idup]]:rr[dup[idup]+1]-1]
          nuv[these].ndup = n_elements(these)
; pick the best observation based on S/N *and* radius
          nuv[these].isbest_snr = 0
          maxsnr = max(nuv[these].flux_em/nuv[these].flux_err_em,maxindx)
          nuv[these[maxindx]].isbest_snr = 1

          nuv[these].isbest_radius = 0
          minradius = min(nuv[these].radius,minindx)
          nuv[these[minindx]].isbest_radius = 1
       endfor
       if (n_elements(fuv) ne 0L) then begin
          hh = histogram(fuv.ident,bin=1,reverse=rr)
          dup = where(hh gt 1,ndup)
          for idup = 0L, ndup-1 do begin
             these = rr[rr[dup[idup]]:rr[dup[idup]+1]-1]
             fuv[these].ndup = n_elements(these)

             fuv[these].isbest_snr = 0
             maxsnr = max(fuv[these].flux_em/fuv[these].flux_err_em,maxindx)
             fuv[these[maxindx]].isbest_snr = 1
             
             fuv[these].isbest_radius = 0
             minradius = min(fuv[these].radius,minindx)
             fuv[these[minindx]].isbest_radius = 1
          endfor
       endif
; write out
       im_mwrfits, nuv, outfile_nuv, clobber=clobber
       if (n_elements(fuv) ne 0L) then $
         im_mwrfits, fuv, outfile_fuv, clobber=clobber
    endif

; --------------------------------------------------
; combine into a single catalog for each field; note that every FUV
; measurement should have a corresponding NUV measurement (although
; there can be edge effects; we ignore them here)
    if keyword_set(combine) then begin
       delvarx, nuv, fuv
       nuvfile = galexpath+strlowcase(field)+'_emphot_nuv.fits.gz'
       fuvfile = repstr(nuvfile,'nuv','fuv')
       allnuv = mrdfits(nuvfile,1)
       nuv = allnuv[where(allnuv.isbest_snr,nnuv)] ; use S/N
;      allident = nuv.ident          
       if file_test(fuvfile) then begin
          allfuv = mrdfits(fuvfile,1)
          fuv = allfuv[where(allfuv.isbest_snr,nfuv)]
;         allident = [allident,fuv.ident] ; crazy?!?
       endif
; sort by (unique) identification number
;      ident = allident[uniq(allident,sort(allident))]
       ident = nuv.ident
       nobj = n_elements(ident)
       splog, field+': number of objects = '+strtrim(nobj,2)
; build the line-matched output catalog by juggling the tag names!
       nuvtags = tag_names(nuv)
       nuvtags = nuvtags[where(strmatch(nuvtags,'*fuv*',/fold) eq 0)]
       newnuvtags = nuvtags
       these = where(strmatch(newnuvtags,'*nuv*',/fold) eq 0)
       newnuvtags[these] = 'NUV_'+newnuvtags[these]
       nuv_out = zerod_empty_structure(nuv[0],ncopies=nobj)
       match, ident, nuv.ident, m1, m2
       srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
       nuv_out[m1] = nuv[m2]
       nuv_out = im_struct_trimtags(temporary(nuv_out),$
         select=nuvtags,newtags=newnuvtags)
       if file_test(fuvfile) then begin
          fuvtags = tag_names(fuv)
          fuvtags = fuvtags[where(strmatch(fuvtags,'*nuv*',/fold) eq 0)]
          newfuvtags = fuvtags
          these = where(strmatch(newfuvtags,'*fuv*',/fold) eq 0)
          newfuvtags[these] = 'FUV_'+newfuvtags[these]
          fuv_out = zerod_empty_structure(fuv[0],ncopies=nobj)
          match, ident, fuv.ident, m1, m2
          srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
          fuv_out[m1] = fuv[m2]
          fuv_out = im_struct_trimtags(temporary(fuv_out),$
            select=fuvtags,newtags=newfuvtags)
       endif else begin
          fuv_out = zerod_empty_structure(nuv_out,ncopies=nobj)
          tags = tag_names(nuv_out)
          newtags = repstr(tags,'NUV','FUV')
          fuv_out = im_struct_trimtags(temporary(fuv_out),$
            select=tags,newtags=newtags)
       endelse
; combine and write out!          
       outfile = galexpath+strlowcase(field)+'_emphot_final.fits'
       out = struct_addtags(temporary(nuv_out),temporary(fuv_out))
       im_mwrfits, out, outfile, clobber=clobber
    endif
    
stop    

return
end
