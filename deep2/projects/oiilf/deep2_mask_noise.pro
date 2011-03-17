pro deep2_mask_noise, write=write
; jm08sep04nyu - build and write out the noise model for each mask  

    ii = read_deep2(/ispec)

    file = strtrim(ii.file,2) ; strtrim(ii.maskname,2)
    mask = strmid(file,0,4)
    umask = mask[uniq(mask,sort(mask))]
    nmask = n_elements(umask)

    outpath = deep2_path(/projects)+'oiilf/mask_noise/'
    
    if keyword_set(write) then begin
       im_plotfaves, /postscript
       psfile = outpath+'deep2_mask_noise.ps'
       dfpsplot, psfile, /landscape, /color
    endif

    npix = 4096*2
;   for imask = 0L, 1L do begin
    for imask = 0L, nmask-1L do begin

       print, format='("Reading mask ",I0,"/",I0,".",A1,$)', $
         imask+1, nmask, string(13b)
       
       these = where(umask[imask] eq mask,nthese)
       ss = read_deep2_specfit(ii[these].file,/obswave,/silent)
       allwave = reform(ss[*,0,*])
       allferr = reform(1.0/sqrt((ss[*,5,*] + (ss[*,5,*] le 0.0))) * (ss[*,5,*] ne 0.0))

       meanwave = 7500.0 ; mean(ss[*,0,*])
       scale = fltarr(nthese)

       for jj = 0L, nthese-1L do begin
          ww = where((allwave[*,jj] gt meanwave-25.0) and (allwave[*,jj] lt meanwave+25.0))
          scale[jj] = median(allferr[ww,jj])
       endfor
       meanscale = mean(scale)

; interpolate onto a common wavelength vector       
       
       minwave = min(allwave) & maxwave = max(allwave)
       dwave = median(allwave-shift(allwave,1))
       wave = findgen((maxwave-minwave)/dwave+1)*dwave+minwave
       allnoise = fltarr(n_elements(wave),nthese)

       for jj = 0L, nthese-1L do begin
          linterp, allwave[*,jj], allferr[*,jj]*meanscale/scale[jj], $
            wave, allnoise1, missing=0.0
          allnoise[*,jj] = allnoise1
       endfor
;      noise = djs_median(allnoise,2)
       noise = wave*0.0-1.0
       for kk = 0L, n_elements(wave)-1L do begin
          good = where(allnoise[kk,*] gt 0.0,ngood)
          if (ngood gt 10L) then noise[kk] = median(allnoise[kk,good]) ; require 10 good spectra
       endfor

       good = where(noise gt 0.0)
       wave = wave[good] & noise = noise[good]
       
; make a plot       

;      xrange = [7000,7050]
       xrange = minmax(wave)
       yrange = minmax(noise)
       
       djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, xsty=3, ysty=3, $
         xtitle='Observed Wavelength (\AA)', ytitle='Scaled Noise (ADU)'
       for jj = 0L, 10L<(nthese-1L) do djs_oplot, allwave[*,jj], $
;      for jj = 0L, nthese-1L do djs_oplot, allwave[*,jj], $
         allferr[*,jj]*meanscale/scale[jj], ps=3
       djs_oplot, wave, noise, color='red', thick=1.5
       im_legend, 'Mask/'+strtrim(umask[imask],2), /left, /top, box=0

;      help, wave, noise
       if keyword_set(write) then begin
          outfile = outpath+'mask_'+strtrim(umask[imask],2)+'.fits'
          splog, 'Writing '+outfile
          mwrfits, {mask: fix(umask[imask]), nobj: fix(nthese), $
            wave: wave, noise: noise}, outfile, /create
       endif
       
    endfor

    if keyword_set(write) then begin
       dfpsclose & spawn, 'gzip -f '+psfile
;      spawn, 'ps2pdf '+psfile+' '+repstr(psfile,'.ps','.pdf')
;      spawn, 'rm -f '+psfile
       im_plotfaves
    endif

stop    
    
return
end
    
