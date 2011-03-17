;+
; NAME:
;       DEEP2_DISPLAY_SPECTRUM
;
; PURPOSE:
;       Display an deep2 spectrum and the best-fitting model.  
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Sep 29, U of A
;       jm07nov27nyu - cleaned up significantly
;-

pro deep2_display_spectrum, deep2, specfit=specfit, nsmooth=nsmooth, maxwave=maxwave, minwave=minwave, $
  psname=psname, position=position, labeltype=labeltype, yrangetype=yrangetype, plotobswave=plotobswave, $
  postscript=postscript, bigpostscript=bigpostscript, silent=silent, nowindow=nowindow

    ngalaxy = n_elements(deep2)
    if (ngalaxy eq 0L) then begin
       doc_library, 'deep2_display_spectrum'
       return
    endif

    if (n_elements(nsmooth) eq 0L) then nsmooth = 10L
    if (n_elements(labeltype) eq 0L) then labeltype = 1L
    if (n_elements(yrangetype) eq 0L) then yrangetype = 5L

; call this routine recursively    
    
    if (ngalaxy gt 1L) then begin
       if (n_elements(psname) eq 0L) then psname = 'deep2_spectrum.ps'
       if keyword_set(postscript) then begin
          bigpostscript = 1L
          dfpsplot, psname, /color, /landscape
       endif else im_window, 0, xratio=0.9, yratio=0.7
       for igal = 0L, ngalaxy-1L do begin
          if (n_elements(specfit) ne 0L) then specfit1 = reform(specfit[*,*,igal])
          deep2_display_spectrum, deep2[igal], specfit=specfit1, nsmooth=nsmooth, maxwave=maxwave, minwave=minwave, $
            psname=psname, position=position, labeltype=labeltype, yrangetype=yrangetype, plotobswave=plotobswave, $
            bigpostscript=bigpostscript, /silent, /nowindow
          if (n_elements(specfit) eq 0L) then delvarx, specfit1
          if (igal le ngalaxy-2L) and (not keyword_set(postscript)) then begin
             prompt:
             print, 'Galaxy #'+string(igal,format='(I0)')+', '+strtrim(deep2[igal].galaxy,2)+' [Options: b,g,q]'
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

    if (n_elements(psname) eq 0L) then psname = 'deep2_spectrum.ps'
    
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
    
    scale = 1.0
    ytitle = 'Counts (ADU)'
;   ytitle = 'Flux (10^{-17} '+flam_units()+')'
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

    galaxy = strtrim(deep2.galaxy,2)
    nicegalaxy = repstr(strtrim(strmid(galaxy,5),2),'_','/')
    zobj = deep2.z
    zobj_str = strtrim(string(deep2.z,format='(F12.3)'),2)
    if (n_elements(specfit) eq 0L) then begin
       datapath = deep2_path(/dr3)
;      specfit = read_deep2_specfit(galaxy,silent=silent)
       spec1 = mrdfits(datapath+strtrim(deep2.file,2),1,/silent)
       spec2 = mrdfits(datapath+strtrim(deep2.file,2),2,/silent)
       good = where(([spec1.ivar,spec2.ivar] gt 0.0),npix)
       specfit = fltarr(npix,3)
       specfit[*,0] = ([spec1.lambda,spec2.lambda])[good]/(1.0+zobj)
       specfit[*,1] = ([spec1.spec,spec2.spec])[good]*(1.0+zobj)
       ivar = ([spec1.ivar,spec2.ivar])[good]
       specfit[*,2] = (1.0+zobj)^2.0/sqrt(ivar)
    endif

;   dims = size(specfit,/dimension)
;   ndims = size(specfit,/n_dimension)
;   if (dims[1] ne 5L) and (dims[0] ne 4650L) then begin
;      splog, 'Improper SPECFIT input.'
;      return
;   endif
;   if (ndims eq 2L) then begin
;      if (ngalaxy ne 1L) then begin
;         splog, 'Dimensions of SPECFIT and GALAXY do not agree.'
;         return
;      endif
;   endif else begin
;      if (dims[2] ne ngalaxy) then begin
;         splog, 'Dimensions of SPECFIT and GALAXY do not agree.'
;         return
;      endif 
;   endelse 

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
    smoothflux = smooth(flux,nsmooth,/nan)

    stats = im_stats(smoothflux,sigrej=5.0)
    xrange = minmax(plotwave)
    if (n_elements(maxwave) ne 0L) then xrange[1] = xrange[1]<maxwave
    if (n_elements(minwave) ne 0L) then xrange[0] = xrange[0]>minwave

    case yrangetype of
       1L: begin ; for cropping the red leak
          get_element, obswave, [min(obswave),8500.0], xx
          minstats = im_stats(smoothflux[xx[0]:xx[1]],sigrej=3.0)
          maxstats = im_stats(smoothflux[xx[0]:xx[1]],sigrej=5.0)
          yrange = [minstats.minrej,1.1*maxstats.maxrej]
       end
       2L: begin
          stats = im_stats(cflux,sigrej=8.0)
          yrange = [stats.minrej,factorymax*stats.maxrej]
          yrange[0] = (yrange[0]<(min(cflux)*offsetscale))>0
       end
       3L: begin
          get_element, plotwave, xrange, xx
          stats = im_stats(smoothflux[xx[0]:xx[1]],sigrej=5.0)
          yrange = [stats.minrej,max(modelflux[xx[0]:xx[1]])]
;         yrange = minmax(modelflux[xx[0]:xx[1]])
          yrange[0] = (yrange[0]<(min(smoothcflux[xx[0]:xx[1]])*offsetscale))>0
       end
       4L: begin                
          stats = im_stats(smoothflux,sigrej=8.0)
          yrange = [stats.minrej,stats.maxrej]
          yrange[0] = (yrange[0]<(min(smoothcflux)*offsetscale))>0
       end
       5L: begin
          get_element, plotwave, xrange, xx
          yrange = minmax(smoothflux[xx[0]:xx[1]])
       end
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

    djs_oplot, plotwave, smoothflux, ps=10, color='grey', thick=plotthick

; for examining the [O II] fits    
    
    if (labeltype eq 2L) then begin

    endif

; legend    
    
    label = nicegalaxy
    case labeltype of
       1L: label = ['DEEP2 '+nicegalaxy,'z = '+zobj_str]
       2L: begin
          if tag_exist(deep2,'oii_3727_1_ew') then begin
             if (deep2.oii_3727_1_ew[1] lt 0.0) then ewoii_1 = 'EW([O II]) \lambda3726 = ...' else $
               ewoii_1 = 'EW([O II]) \lambda3726 = '+strtrim(string(deep2.oii_3727_1_ew[0],format='(F12.1)'),2)+'\pm'+$
               strtrim(string(deep2.oii_3727_1_ew[1],format='(F12.1)'),2)+' \AA' 
          endif else ewoii_1 = ''
          if tag_exist(deep2,'oii_3727_2_ew') then begin
             if (deep2.oii_3727_2_ew[1] lt 0.0) then ewoii_2 = 'EW([O II]) \lambda3729 = ...' else $
               ewoii_2 = 'EW([O II]) \lambda3729 = '+strtrim(string(deep2.oii_3727_2_ew[0],format='(F12.1)'),2)+'\pm'+$
               strtrim(string(deep2.oii_3727_2_ew[1],format='(F12.1)'),2)+' \AA' 
          endif else ewoii_2 = ''
          if tag_exist(deep2,'h_beta_ew') then begin
             if (deep2.h_beta_ew[1] lt 0.0) then ewhb = 'EW(H\beta) = ...' else $
               ewhb = 'EW(H\beta) = '+strtrim(string(deep2.h_beta_ew[0],format='(F12.1)'),2)+'\pm'+$
               strtrim(string(deep2.h_beta_ew[1],format='(F12.1)'),2)+' \AA' 
          endif else ewhb = ''
          label = ['DEEP2 '+nicegalaxy,'z = '+zobj_str,ewoii_1,ewoii_2,ewhb]
       end
    endcase

    if (labeltype le 2L) then begin
       legend, textoidl(label), /top, box=0, charsize=1.5, charthick=postthick2, $
         spacing=1.6, /normal, margin=0, clear=keyword_set(postscript)
    endif
       
    if keyword_set(postscript) then dfpsclose

return
end
    
