pro ages_render_thumb, id, band=band, detections=detections, psfile=psfile
; jm10jan31ucsd - render a thumbnail image of an AGES galaxy; ID
;   should be the AGES_ID number

    common ages_thumb, bootes
    
    ngal = n_elements(id)
    if (ngal eq 0) then begin
       doc_library, 'ages_render_thumb'
       return
    endif
    
    if (n_elements(band) eq 0) then band = 'I'
    path = ages_path(/thumb)+strupcase(band)+'/'

; load the full BOOTES catalog    
    if keyword_set(detections) then begin
       if (n_elements(bootes) eq 0L) then begin
          bootespath = getenv('RESEARCHPATH')+'/data/bootes/'
          bootesfile = bootespath+'bootes_I.fits.gz'
          splog, 'Reading '+bootesfile
          bootes = mrdfits(bootesfile,1)
       endif
    endif
    
    im_plotconfig, 0, pos, psfile=psfile
    for ii = 0, ngal-1 do begin
       fits = path+'ages_'+string(id[ii],format='(I5.5)')+'.fits.gz'
       if (file_test(fits,/reg) eq 0) then begin
          splog, 'Thumbnail image '+fits+' not found!'
          continue
       endif
       image = mrdfits(fits,0,hdr,/silent)
; scale and display       
       img = asinhscl(image,/neg);,alpha=5.0,beta=16.0,omin=0,omax=250)
       plotimage, img, /normal, /noaxes, /preserve, $
         position=pos
       if keyword_set(detections) then begin
          extast, hdr, astr
          ad2xy, bootes.alpha_j2000, bootes.delta_j2000, astr, xx, yy
          these = where((xx-astr.crpix[0] gt -astr.naxis[0]) and $
            (xx-astr.crpix[0] lt astr.naxis[0]) and $
            (yy-astr.crpix[1] gt -astr.naxis[1]) and $
            (yy-astr.crpix[1] lt astr.naxis[1]),nthese)
          djs_oplot, xx[these]-astr.crpix[0], yy[these]-astr.crpix[1], $
            psym=7, color='red', symsize=1.5

;         pscale = abs(astr.cd[0,0])
;         these = where((bootes.alpha_j2000 gt astr.crval[0]-$
;           pscale*1.05*astr.naxis[0]) and (bootes.alpha_j2000 lt $
;           astr.crval[0]+pscale*astr.naxis[0]) and $
;           (bootes.delta_j2000 gt astr.crval[1]-$
;           pscale*astr.naxis[1]) and (bootes.delta_j2000 lt $
;           astr.crval[1]+1.05*pscale*astr.naxis[0]),nthese)
;         ad2xy, bootes[these].alpha_j2000, $
;           bootes[these].delta_j2000, astr, xx, yy
;         djs_oplot, xx-astr.crpix[0], yy-astr.crpix[1], psym=7, $
;           color='red', symsize=1.5
;         print, astr & print
       endif
          
       if (n_elements(psfile) eq 0) then cc = get_kbrd(1)
    endfor

    im_plotconfig, psfile=psfile, /gzip, /psclose

return
end
    
