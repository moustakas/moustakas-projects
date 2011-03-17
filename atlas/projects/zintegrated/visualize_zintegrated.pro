pro visualize_zintegrated, int_dust, hii, postscript=postscript, encapsulated=encapsulated
; jm05aug23uofa
; jm05sep09uofa - excised into its own routine and streamlined     

    htmlbase = 'zintegrated'
    pspath = atlas_path(/web)+'analysis/'+htmlbase+'/'

    if keyword_set(postscript) then begin
       postthick = 8.0
       postthick2 = 3.0
    endif else begin
       postthick = 2.0
       postthick2 = postthick
       im_window, 0, xratio=0.5, /square
       im_window, 2, xratio=0.5, yratio=0.2
    endelse
    
; read the galaxy and HII-region data

    if (n_elements(int_dust) eq 0L) or (n_elements(hii) eq 0L) then $
      int_dust = read_zintegrated_sample(int_nodust=int_nodust,hii=hii)

; ##################################################    
; SPIRAL GALAXIES
; ##################################################    
    
    galaxy = strtrim(int_dust.galaxy,2)
    nicegalaxy = strtrim(int_dust.nice_galaxy,2)
    ngalaxy = n_elements(int_dust)

    for k = 0L, ngalaxy-1L do begin
       
; retrieve the HII region measurements for this object to get the
; radius of the last HII region

       match = where(galaxy[k] eq strtrim(hii.atlas_galaxy,2),nmatch)
       keep = where(hii[match].radius gt -9000.0,nkeep,comp=nodata,ncomp=nnodata)

       hiispiral = hii[match[keep]]

; the following HII regions are being excluded, but should be
; investigated 
       
       case galaxy[k] of
          'NGC1058': good = where($
            (strmatch(hiispiral.hii_region,'*FGW1058G*',/fold) eq 0B) and $
            (strmatch(hiispiral.hii_region,'*FGW1058H*',/fold) eq 0B),comp=donotfit,ncomp=ndonotfit)
          'NGC2903': good = where($
            (strmatch(hiispiral.hii_region,'*-062-085*',/fold) eq 0B) and $
            (strmatch(hiispiral.hii_region,'*-065-073*',/fold) eq 0B) and $
            (strmatch(hiispiral.hii_region,'*-067-061*',/fold) eq 0B) and $
            (strmatch(hiispiral.hii_region,'*+171+243*',/fold) eq 0B) and $
            (strmatch(hiispiral.hii_region,'*H16*',/fold) eq 0B),comp=donotfit,ncomp=ndonotfit)
          'NGC4651': good = where($
            (strmatch(hiispiral.hii_region,'*+131+021*',/fold) eq 0B),comp=donotfit,ncomp=ndonotfit)
          else: begin
             good = lindgen(n_elements(hiispiral))
             ndonotfit = 0L
          end
       endcase 

       hiispiral = hiispiral[good]       
       circle_radius = max(hiispiral.radius)
       
       splog, 'Found '+string(nkeep,format='(I0)')+' good HII regions'
       struct_print, struct_trimtags(hiispiral,select=['HII_GALAXY','HII_REGION',$
         'ZSTRONG_12OH_ZKH94','RADIUS','REFERENCE']), /no_head

; visualize the galaxy; overlay a circle at the maximum HII-region
; radius           

       if (not keyword_set(postscript)) then wset, 0
       atlas_display_image, int_dust[k], imagepath=imagepath, $
         lcharsize=1.8, pcharsize=pcharsize, /preserve_aspect, $
         pspath=pspath, labeltype=1L, _extra=extra, postscript=postscript, $
         circle_radius=max(hiispiral.radius), encapsulated=encapsulated
       
; visualize the galaxy image and the spectrum

       if (not keyword_set(postscript)) then wset, 2
       atlas_display_image, int_dust[k], imagepath=imagepath, $
         lcharsize=1.0, pcharsize=1.0, /preserve_aspect, $
         pspath=pspath, labeltype=-1L, _extra=extra, postscript=postscript, $
         /with_spectrum, /nobar, /nolabelbar, encapsulated=encapsulated

       if (not keyword_set(postscript)) then cc = get_kbrd(1)

    endfor

stop    
    
return
end
