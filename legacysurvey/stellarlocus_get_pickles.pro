function stellarlocus_get_pickles, maggies, bestmaggies, filterlist=filterlist, info=info
; called by the various stellarlocus plotting routines
    
    ngal = (size(maggies,/dim))[1]
    nfilt = n_elements(filterlist)
    filters = strarr(nfilt)
    for ii = 0, nfilt-1 do filters[ii] = repstr(strmid(filterlist[ii],$
      strpos(filterlist[ii],'_',/reverse_search)+1),'.par','')

; Pickles     
    pickles = mrdfits(getenv('IM_PROJECTS_DIR')+'/decals/dr1/qaplots/pickles.fits',1,/silent)
    plambda = k_lambda_to_edges(pickles[0].wave)
    pmaggies = fltarr(nfilt,n_elements(pickles))
    npickles = n_elements(pickles)
    for ii = 0, npickles-1 do pmaggies[*,ii] = k_project_filters($
      plambda,pickles[ii].flux,filterlist=filterlist,/silent)

; Kurucz    
    kuruczfile = getenv('IDLSPEC2D_DIR')+'/etc/kurucz_stds_v5.fit'
    kflux = mrdfits(kuruczfile,0,hdr,/silent)
    kinfo = mrdfits(kuruczfile,1,/silent)
    keep = where(kinfo.feh eq -1.0 and kinfo.g eq 4.5,nkurucz)
    kinfo = kinfo[keep]
    kflux = kflux[*,keep]
    kwave = 10D^(sxpar(hdr,'CRVAL1')+findgen(n_elements(kflux[*,0]))*sxpar(hdr,'CD1_1'))

    klambda = k_lambda_to_edges(kwave)
    kmaggies = fltarr(nfilt,nkurucz)
    for ii = 0, nkurucz-1 do kmaggies[*,ii] = k_project_filters(klambda,$
      kflux[*,ii],filterlist=filterlist,/silent)

; get the observed colors for the data and for Pickles and Kurucz
    colors = {filterlist: strtrim(filterlist,2), filters: filters, $
      colors: strarr(nfilt-1), pickles_feh: pickles.feh, pickles_type: pickles.type, $
      kurucz_feh: kinfo.feh, kurucz_teff: kinfo.teff, kurucz_g: kinfo.g}
    for ii = 1, nfilt-1 do begin
       tag = filters[ii-1]+filters[ii]
       colors.colors[ii-1] = tag
       colors = create_struct(colors,tag,fltarr(ngal)-999.0,tag+'_best',fltarr(ngal)-999.0,$
         tag+'_pickles',fltarr(npickles)-999.0,tag+'_kurucz',fltarr(nkurucz)-999.0)
; data, both with and without the zeropoint offsets
       indx = tag_indx(colors,tag)
       good = where((maggies[ii-1,*] gt 0.0) and (maggies[ii,*] gt 0.0),ngood)
       if (ngood ne 0L) then begin
          junk = colors.(indx)
          junk[good] = reform(-2.5*alog10(maggies[ii-1,good]/maggies[ii,good]))
          colors.(indx) = junk
       endif

       indx = tag_indx(colors,tag+'_best')
       good = where((bestmaggies[ii-1,*] gt 0.0) and (bestmaggies[ii,*] gt 0.0),ngood)
       if (ngood ne 0L) then begin
          junk = colors.(indx)
          junk[good] = reform(-2.5*alog10(bestmaggies[ii-1,good]/bestmaggies[ii,good]))
          colors.(indx) = junk
       endif
; Pickles
       indx = tag_indx(colors,tag+'_pickles')
       good = where((pmaggies[ii-1,*] gt 0.0) and (pmaggies[ii,*] gt 0.0),ngood)
       if (ngood ne 0L) then begin
          junk = colors.(indx)
          junk[good] = reform(-2.5*alog10(pmaggies[ii-1,good]/pmaggies[ii,good]))
          colors.(indx) = junk
       endif
; Kurucz
       indx = tag_indx(colors,tag+'_kurucz')
       good = where((kmaggies[ii-1,*] gt 0.0) and (kmaggies[ii,*] gt 0.0),ngood)
       if (ngood ne 0L) then begin
          junk = colors.(indx)
          junk[good] = reform(-2.5*alog10(kmaggies[ii-1,good]/kmaggies[ii,good]))
          colors.(indx) = junk
       endif
    endfor

return, colors
end

