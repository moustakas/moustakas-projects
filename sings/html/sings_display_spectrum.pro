;+
; NAME:
;       SINGS_DISPLAY_SPECTRUM
;
; PURPOSE:
;       Display an sings spectrum and the best-fitting model.  
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Sep 29, U of A
;       jm07nov27nyu - cleaned up significantly
;-

pro sings_display_spectrum, sings, specfit=specfit, nsmooth=nsmooth, maxwave=maxwave, minwave=minwave, $
  psname=psname, position=position, labeltype=labeltype, no_specfit=no_specfit, only_continuum=only_continuum, $
  yrangetype=yrangetype, xrange=xrange, plotobswave=plotobswave, postscript=postscript, bigpostscript=bigpostscript, $
  silent=silent, nowindow=nowindow, nuclear=nuclear, drift20=drift20, drift56=drift56

    ngalaxy = n_elements(sings)
    if (ngalaxy eq 0L) then begin
       doc_library, 'sings_display_spectrum'
       return
    endif

    if (n_elements(nsmooth) eq 0L) then nsmooth = 0L
    if (n_elements(labeltype) eq 0L) then labeltype = 1L
    if (n_elements(yrangetype) eq 0L) then yrangetype = 1L

    if (not keyword_set(nuclear)) and (not keyword_set(drift20)) and $
      (not keyword_set(drift56)) then begin
       splog, 'Either NUCLEAR *or* DRIFT20 *or* DRIFT56 keyword must be set!'
       return
    endif

    if (keyword_set(nuclear) and keyword_set(drift20)) or $
      (keyword_set(nuclear) and keyword_set(drift56)) or $
      (keyword_set(drift20) and keyword_set(drift56)) then begin
       splog, 'Only one keyword (NUCLEAR, DRIFT20, or DRIFT56) can be set at the same time!'
       return
    endif

; call this routine recursively    
    
    if (ngalaxy gt 1L) then begin
       if (n_elements(psname) eq 0L) then psname = 'sings_spectrum.ps'
       if keyword_set(postscript) then begin
          bigpostscript = 1L
          dfpsplot, psname, /color, /landscape
       endif else im_window, 0, xratio=0.9, yratio=0.7
       for igal = 0L, ngalaxy-1L do begin
          if (n_elements(specfit) ne 0L) then specfit1 = reform(specfit[*,*,igal])
          sings_display_spectrum, sings[igal], specfit=specfit1, nsmooth=nsmooth, maxwave=maxwave, minwave=minwave, $
            psname=psname, position=position, labeltype=labeltype, no_specfit=no_specfit, only_continuum=only_continuum, $
            yrangetype=yrangetype, xrange=xrange, plotobswave=plotobswave, bigpostscript=bigpostscript, /silent, /nowindow, $
            nuclear=nuclear, drift20=drift20, drift56=drift56
          if (n_elements(specfit) eq 0L) then delvarx, specfit1
          if (igal le ngalaxy-2L) and (not keyword_set(postscript)) then begin
             prompt:
             print, 'Galaxy #'+string(igal,format='(I0)')+', '+strtrim(sings[igal].galaxy,2)+' [Options: b,g,q]'
             cc = strupcase(get_kbrd(1))
             case strlowcase(strcompress(cc,/remove)) of
                'b': igal = (igal-2L) ; back
                'g': begin      ; goto 
                   number = ''
                   read, number, prompt='Goto galaxy (0-'+string(ngalaxy-1L,format='(I0)')+'): '
                   number = 0 > long(number-1L) < (ngalaxy-2L)
                   igal = number
                end
                'r': igal = igal-1L ; redraw
                'q': begin
                   return       ; quit
                end
                else: 
             endcase
          endif 
       endfor
       if keyword_set(postscript) then begin
          dfpsclose
          spawn, 'gzip -f '+psname, /sh
       endif
       return
    endif

; plotting variables    

    if (n_elements(psname) eq 0L) then psname = 'sings_spectrum.ps'
    
    if keyword_set(postscript) or keyword_set(bigpostscript) then begin
       postthick1 = 2.0
       postthick2 = 3.0
       plotthick1 = 1.5
       plotthick2 = 2.0
       if keyword_set(postscript) then dfpsplot, psname, /color, /landscape
    endif else begin
       postthick1 = 3.0
       postthick2 = 2.0
       plotthick1 = 2.0
       plotthick2 = 3.0
       if (not keyword_set(nowindow)) then im_window, 0, xratio=0.9, yratio=0.7
    endelse
    
    scale = 1E15
    ytitle = 'Flux (10^{-15} '+flam_units()+')'
    obsxtitle = 'Observed Wavelength (\AA)'
    restxtitle = 'Rest Wavelength (\AA)'

    if (n_elements(position) eq 0L) then begin
       if keyword_set(plotobswave) then $
         pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, xmargin=[1.2,0.2], $
         ymargin=[0.7,1.4], position=pos, /normal else $
         pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, xmargin=[1.2,0.2], $
         ymargin=[0.3,1.4], position=pos, /normal
    endif else pos = position

; read the data

    galaxy = strtrim(sings.galaxy,2)
    nicegalaxy = strtrim(sings.nice_galaxy,2)
    zobj = sings.z
    zobj_str = strtrim(string(sings.z,format='(F12.3)'),2)
    if (n_elements(specfit) eq 0L) then specfit = read_sings_specfit(galaxy,$
      silent=1,nuclear=nuclear,drift20=drift20,drift56=drift56)

    dims = size(specfit,/dimension)
    ndims = size(specfit,/n_dimension)

    if (dims[1] ne 6L) and (dims[0] ne 4650L) then begin
       splog, 'Improper SPECFIT input.'
       return
    endif

    if (ndims eq 2L) then begin
       if (ngalaxy ne 1L) then begin
          splog, 'Dimensions of SPECFIT and GALAXY do not agree.'
          return
       endif
    endif else begin
       if (dims[2] ne ngalaxy) then begin
          splog, 'Dimensions of SPECFIT and GALAXY do not agree.'
          return
       endif 
    endelse 

; format the spectrum

    restwave = reform(specfit[*,0])
    good = where((restwave ne 0.0),ngood)
    if (ngood eq 0L) then message, 'Problem here!' else restwave = restwave[good]

    obswave = restwave*(1.0+zobj)
    if keyword_set(plotobswave) then begin
       plotwave = obswave
       zfactor = (1.0+zobj) 
    endif else begin
       plotwave = restwave
       zfactor = 1.0
    endelse

    flux = scale*reform(specfit[good,1])/zfactor
    cflux = scale*reform(specfit[good,2])/zfactor
    eflux = (scale*reform(specfit[good,3])/zfactor) > 0.0 ; NOTE!
    mcflux = scale*reform(specfit[good,4])/zfactor

    cnolinesflux = flux - eflux
    fluxnoc = flux - cflux - mcflux
    modelflux1 = cflux + eflux
    modelflux2 = cflux + eflux + mcflux
    modelflux3 = cflux + mcflux

    smoothflux = smooth(flux,nsmooth,/nan) ; for the xyranges
    smoothfluxnoc = smooth(fluxnoc,nsmooth,/nan)
    smoothmodelflux1 = smooth(modelflux1,nsmooth,/nan)
    smoothmodelflux2 = smooth(modelflux2,nsmooth,/nan)
    smoothmodelflux3 = smooth(modelflux3,nsmooth,/nan)
    smootheflux = smooth(eflux,nsmooth,/nan)
    smoothcflux = smooth(cflux,nsmooth,/nan)
    smoothmcflux = smooth(mcflux,nsmooth,/nan)
    smoothcnolinesflux = smooth(cnolinesflux,nsmooth,/nan)
    
    stats = im_stats(smoothflux,sigrej=5.0)
    if (n_elements(xrange) eq 0L) then xrange = minmax(plotwave)
    if (n_elements(maxwave) ne 0L) then xrange[1] = xrange[1]<maxwave
    if (n_elements(minwave) ne 0L) then xrange[0] = xrange[0]>minwave

    case yrangetype of
       1L: begin
          get_element, plotwave, xrange, xx
          yrange = minmax(modelflux2[xx[0]:xx[1]])
          yrange[0] = yrange[0]*0.95
          yrange[1] = yrange[1]*1.05
       end
       else: yrange = minmax(flux)
    endcase

; plot it     
    
    if keyword_set(plotobswave) then begin
       plot, [0], [0], /nodata, xthick=postthick1, xsty=11, ysty=3, $
         ythick=postthick1, xtitle=textoidl(obsxtitle), ytitle=textoidl(ytitle), $
         charsize=1.7, charthick=postthick2, yrange=yrange, xrange=xrange, $
         thick=postthick1, position=pos[*,0], xtickname=xtickname, ytickname=yticknamep
       axis, /xaxis, xtitle=textoidl(restxtitle), charsize=1.7, charthick=postthick2, $
         xthick=postthick1, xrange=interpol(restwave,obswave,!x.crange), xsty=3
    endif else begin
       plot, [0], [0], /nodata, xthick=postthick1, xsty=3, ysty=3, $
         ythick=postthick1, xtitle=textoidl(restxtitle), ytitle=textoidl(ytitle), $
         charsize=1.7, charthick=postthick2, yrange=yrange, xrange=xrange, $
         thick=postthick1, position=pos[*,0], xtickname=xtickname, ytickname=yticknamep
    endelse

    if keyword_set(no_specfit) then begin
       djs_oplot, plotwave, flux, ps=10, color='grey', thick=plotthick
    endif else begin
       if keyword_set(bigpostscript) or keyword_set(postscript) then $
         scolors = ['grey','red','dark green','blue'] else $
           scolors = ['grey','red','yellow','blue']
       djs_oplot, plotwave, smoothflux, ps=10, color=scolors[0], thick=plotthick1
       djs_oplot, plotwave, smoothmodelflux2, ps=10, color=scolors[1], thick=plotthick1
       djs_oplot, plotwave, smoothcflux, ps=10, color=scolors[2], thick=plotthick1
       djs_oplot, plotwave, smoothmcflux, ps=10, color=scolors[3], thick=plotthick1
    endelse

; legend    
    
    label = nicegalaxy

    if tag_exist(sings,'oii_3727_ew') then begin
       if (sings.oii_3727_ew[1] lt 0.0) then ewoii = 'EW([O II]) = ...' else $
         ewoii = 'EW([O II]) = '+strtrim(string(sings.oii_3727_ew[0],format='(F12.1)'),2)+'\pm'+$
         strtrim(string(sings.oii_3727_ew[1],format='(F12.1)'),2)+' \AA, S/N = '+$
         strtrim(string(sings.oii_3727[0]/sings.oii_3727[1],format='(F12.1)'),2)
    endif else ewoii = ''
    if tag_exist(sings,'oiii_5007_ew') then begin
       if (sings.oiii_5007_ew[1] lt 0.0) then ewoiii = 'EW([O III]) = ...' else $
         ewoiii = 'EW([O III]) = '+strtrim(string(sings.oiii_5007_ew[0],format='(F12.1)'),2)+'\pm'+$
         strtrim(string(sings.oiii_5007_ew[1],format='(F12.1)'),2)+' \AA, S/N = '+$
         strtrim(string(sings.oiii_5007[0]/sings.oiii_5007[1],format='(F12.1)'),2)
    endif else ewoiii = ''
    if tag_exist(sings,'h_beta_ew') then begin
       if (sings.h_beta_ew[1] lt 0.0) then ewhb = 'EW(H\beta) = ...' else $
         ewhb = 'EW(H\beta) = '+strtrim(string(sings.h_beta_ew[0],format='(F12.1)'),2)+'\pm'+$
         strtrim(string(sings.h_beta_ew[1],format='(F12.1)'),2)+' \AA, S/N = '+$
         strtrim(string(sings.h_beta[0]/sings.h_beta[1],format='(F12.1)'),2)
    endif else ewhb = ''
    if tag_exist(sings,'h_alpha_ew') then begin
       if (sings.h_alpha_ew[1] lt 0.0) then ewha = 'EW(H\alpha) = ...' else $
         ewha = 'EW(H\alpha) = '+strtrim(string(sings.h_alpha_ew[0],format='(F12.1)'),2)+'\pm'+$
         strtrim(string(sings.h_alpha_ew[1],format='(F12.1)'),2)+' \AA, S/N = '+$
         strtrim(string(sings.h_alpha[0]/sings.h_alpha[1],format='(F12.1)'),2)
    endif else ewha = ''

    case labeltype of
       1L: label = [nicegalaxy] ; ,'z = '+zobj_str]
       2L: label = [nicegalaxy,ewoii,ewhb]
       3L: label = [nicegalaxy,ewoii,ewoiii,ewhb,ewha]
       else: 
    endcase

    if keyword_set(bigpostscript) or keyword_set(postscript) then lcharsize = 1.3 else lcharsize = 1.5
    legend, textoidl(label), /top, box=0, charsize=lcharsize, charthick=postthick2, $
      spacing=1.6, /normal, margin=0, clear=keyword_set(postscript)
       
    if keyword_set(postscript) then dfpsclose

return
end
    
