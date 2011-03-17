;+
; NAME:
;       AGES_DISPLAY_SPECTRUM
;
; PURPOSE:
;       Display an ages spectrum and the best-fitting model.  
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Sep 29, U of A
;       jm07nov27nyu - cleaned up significantly
;-

pro ages_display_spectrum, ages, specfit=specfit, nsmooth=nsmooth, xrange=xrange1, $
  psname=psname, position=position, labeltype=labeltype, no_specfit=no_specfit, $
  only_continuum=only_continuum, yrangetype=yrangetype, plotobswave=plotobswave, $
  postscript=postscript, bigpostscript=bigpostscript, silent=silent, nowindow=nowindow, $
  postthick1=postthick1, postthick2=postthick2, plotthick1=plotthick1, plotthick2=plotthick2, $
  speccolor=speccolor, bc03color=bc03color, bestfitcolor=bestfitcolor, especcolor=especcolor, $
  plottype=plottype, noupperaxis=noupperaxis, pdf=pdf, _extra=extra

    ngalaxy = n_elements(ages)
    if (ngalaxy eq 0L) then begin
       doc_library, 'ages_display_spectrum'
       return
    endif

    if (n_elements(nsmooth) eq 0L) then nsmooth = 1L
    if (n_elements(labeltype) eq 0L) then labeltype = 1L
    if (n_elements(plottype) eq 0L) then plottype = 1L
    if (n_elements(yrangetype) eq 0L) then yrangetype = 1L

; call this routine recursively    
    
    if (ngalaxy gt 1L) then begin
       if (n_elements(psname) eq 0L) then psname = 'ages_spectrum.ps'
       if keyword_set(postscript) then begin
          bigpostscript = 1L
          dfpsplot, psname, /color, /landscape
       endif else im_window, 0, xratio=0.9, yratio=0.7
       if (n_elements(specfit) eq 0L) then specfit = $
         read_ages_specfit(ages.galaxy,silent=silent,_extra=extra)
       for igal = 0L, ngalaxy-1L do begin
          if (n_elements(specfit) ne 0L) then specfit1 = reform(specfit[*,*,igal])
          ages_display_spectrum, ages[igal], specfit=specfit1, nsmooth=nsmooth, $
            xrange=xrange1, psname=psname, position=position, $
            labeltype=labeltype, no_specfit=no_specfit, only_continuum=only_continuum, $
            yrangetype=yrangetype, plotobswave=plotobswave, bigpostscript=bigpostscript, $
            /silent, /nowindow, postthick1=postthick1, postthick2=postthick2, $
            plotthick1=plotthick1, plotthick2=plotthick2, speccolor=speccolor, $
            bc03color=bc03color, bestfitcolor=bestfitcolor, especcolor=especcolor, $
            plottype=plottype, noupperaxis=noupperaxis, pdf=pdf, _extra=extra
          if (n_elements(specfit) eq 0L) then delvarx, specfit1
          if (igal le ngalaxy-2L) and (not keyword_set(postscript)) then begin
             prompt:
             print, 'Galaxy #'+string(igal,format='(I0)')+', '+strtrim(ages[igal].galaxy,2)+' [Options: b,g,q]'
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
          if keyword_set(pdf) then begin
             pdfname = repstr(psname,'.ps','.pdf')
             splog, 'Building PDF: '+file_basename(psname)+'-->'+file_basename(pdfname)
             spawn, 'ps2pdf '+psname+' '+pdfname
             spawn, 'pdftk '+pdfname+' cat 1-endW output /tmp/junk.pdf'
             spawn, '/bin/mv -f /tmp/junk.pdf '+pdfname
             spawn, '/bin/rm -f '+psname
          endif else spawn, 'gzip -f '+psname, /sh
       endif
       return
    endif

; plotting variables    

    if (n_elements(psname) eq 0L) then psname = 'ages_spectrum.ps'
    
    if keyword_set(postscript) or keyword_set(bigpostscript) then begin
       im_plotfaves, /postscript
       if (n_elements(speccolor) eq 0L) then speccolor = djs_icolor('grey')
       if (n_elements(bc03color) eq 0L) then bc03color = djs_icolor('red')
       if (n_elements(medianfitcolor) eq 0L) then medianfitcolor = djs_icolor('default')
       if (n_elements(bestfitcolor) eq 0L) then bestfitcolor = djs_icolor('blue')
       if (n_elements(especcolor) eq 0L) then especcolor = djs_icolor('dark green')
       if keyword_set(postscript) then dfpsplot, psname, /color, /landscape
    endif else begin
       im_plotfaves
       if (n_elements(speccolor) eq 0L) then speccolor = djs_icolor('grey')
       if (n_elements(bc03color) eq 0L) then bc03color = djs_icolor('red')
       if (n_elements(medianfitcolor) eq 0L) then medianfitcolor = djs_icolor('default')
       if (n_elements(bestfitcolor) eq 0L) then bestfitcolor = djs_icolor('yellow')
       if (n_elements(especcolor) eq 0L) then especcolor = djs_icolor('cyan')
       if (not keyword_set(nowindow)) then im_window, 0, xratio=0.9, yratio=0.7
    endelse
    
    if (n_elements(position) eq 0L) then begin
       if keyword_set(plotobswave) then $
         pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, xmargin=[1.0,0.2], $
         ymargin=[0.9,1.4], position=pos, /normal else $
         pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, xmargin=[1.0,0.2], $
         ymargin=[0.3,1.4], position=pos, /normal
    endif else pos = position

; read the data

    galaxy = strtrim(ages.galaxy,2)
    nicegalaxy = repstr(strtrim(strmid(galaxy,5),2),'_','/')
    agesid = string(ages.ages_id,format='(I0)')
    zobj = ages.z
    spzobj = ages.z
    zobj_str = strtrim(string(ages.z,format='(F12.3)'),2)
    spzobj_str = strtrim(string(ages.spzbest,format='(F12.3)'),2)
    if (n_elements(specfit) eq 0L) then specfit = $
      read_ages_specfit(galaxy,silent=silent,_extra=extra)
    if tag_exist(ages,'phot_i') then imag = $
      strtrim(string(ages.phot_i,format='(F12.2)'),2) else imag = '?'

    dims = size(specfit,/dimension)
    ndims = size(specfit,/n_dimension)

    if (dims[1] ne 5L) and (dims[0] ne 4650L) then begin
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

; x- and y-titles
    
    obsxtitle = 'Observed Wavelength (\AA)'
    restxtitle = 'Rest Wavelength (\AA)'
    power = ceil(abs(alog10(median(specfit[*,1,*]))))
    if (power lt 10.0) then begin
       scale = 10.0^(-power)
       ytitle = 'Flux (10^{-'+string(power,format='(I0)')+'} ADU)'
    endif else begin
       scale = 10.0^power
       ytitle = 'Flux (10^{-'+string(power,format='(I0)')+'} '+flam_units()+')'
    endelse
;   scale = 1E17
;   ytitle = 'Flux (10^{-17} '+flam_units()+')'

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
    if (n_elements(xrange1) eq 0L) then xrange = minmax(plotwave) else xrange = xrange1

    case yrangetype of
       1L: begin                ; optimized for the continuum+emission-line spectrum
          get_element, plotwave, xrange, xx
          maxstats = im_stats(smoothflux[xx[0]:xx[1]],sigrej=3.0)
          yrange = [-0.12*maxstats.maxrej,1.5*max(modelflux2[xx[0]:xx[1]])]
;         yrange = [-0.12,1.15]*maxstats.maxrej
       end
       2L: begin                ; full range
          get_element, plotwave, xrange, xx
          yrange = minmax(smoothflux[xx[0]:xx[1]])
       end
       3L: begin                ; like 1L, but use the data minimum
          get_element, obswave, [min(obswave),8500.0], xx
          maxstats = im_stats(smoothflux[xx[0]:xx[1]],sigrej=5.0)
          yrange = [min(smoothflux[xx[0]:xx[1]]),1.1*maxstats.maxrej]
       end
       4L: begin                ; for pure emission lines
          get_element, plotwave, xrange, xx
          yrange = [-3.5*djsig(smoothfluxnoc[xx[0]:xx[1]]),max(eflux[xx[0]:xx[1]])]
       end
       5L: begin                ; for cropping the red leak
          get_element, obswave, [min(obswave),8500.0], xx
          minstats = im_stats(smoothflux[xx[0]:xx[1]],sigrej=3.0)
          maxstats = im_stats(smoothflux[xx[0]:xx[1]],sigrej=5.0)
;         yrange = [minstats.minrej,1.1*maxstats.maxrej]
          yrange = [-0.12,1.15]*maxstats.maxrej
       end
       else: message, 'YRANGETYPE not recognized!'
    endcase
       
; plot it     
    
    if keyword_set(noupperaxis) then thisxstyle = 3 else thisxstyle = 9
    if keyword_set(plotobswave) then begin
       plot, [0], [0], /nodata, xsty=thisxstyle, ysty=1, xtitle=textoidl(obsxtitle), $
         ytitle=textoidl(ytitle), charsize=1.7, yrange=yrange, xrange=xrange, $
         position=pos[*,0], xtickname=xtickname, ytickname=ytickname, _extra=extra
       if (not keyword_set(noupperaxis)) then begin
          axis, /xaxis, xtitle=textoidl(restxtitle), charsize=1.7, $
            xrange=interpol(restwave,obswave,!x.crange), xsty=1, _extra=extra
       endif
    endif else begin
       plot, [0], [0], /nodata, xsty=1, ysty=1, xtitle=textoidl(restxtitle), $
         ytitle=textoidl(ytitle), charsize=1.7, yrange=yrange, xrange=xrange, $
         position=pos[*,0], xtickname=xtickname, ytickname=ytickname, _extra=extra
    endelse

    case plottype of
       0L: begin ; just the spectrum
          djs_oplot, plotwave, smoothflux, ps=10, color=speccolor
       end
       1L: begin ; spectrum + best fit
          djs_oplot, plotwave, smoothflux, ps=10, color=speccolor
          djs_oplot, plotwave, smoothmodelflux3, ps=10, color=djs_icolor('red')
          djs_oplot, plotwave, smoothcflux, ps=10, color=bestfitcolor
       end
       2L: begin ; spectrum + best fit
          djs_oplot, plotwave, smoothflux, ps=10, color=speccolor
;         djs_oplot, plotwave, smoothmodelflux2, ps=10, color=djs_icolor('red'), line=2
          djs_oplot, plotwave, smoothmodelflux2, ps=10, color=djs_icolor('red'), line=0
          djs_oplot, plotwave, smoothcflux, ps=10, color=bestfitcolor
       end
       3L: begin ; emission-line (residual) spectrum + best fit
          djs_oplot, plotwave, smoothfluxnoc, ps=10, color=speccolor
          djs_oplot, plotwave, smootheflux, ps=10, color=especcolor
       end
       4L: begin ; spectrum + BC03 fit + smooth residual model
          djs_oplot, plotwave, smoothflux, ps=10, color=speccolor
;         djs_oplot, !x.crange, [0,0], line=0, thick=1.0
          djs_oplot, plotwave, smoothcflux, ps=10, color=bc03color
          djs_oplot, plotwave, mcflux, ps=10, color=medianfitcolor
          djs_oplot, plotwave, smoothmodelflux3, ps=10, color=bestfitcolor
;         djs_oplot, plotwave, smootheflux, ps=10, color=especcolor
       end
       else: message, 'PLOTTYPE not recognized!'
    endcase
    
; legend    
    
    if tag_exist(ages,'mgii_2800_ew') then begin
       if (ages.mgii_2800_ew[1] lt 0.0) then ewmgii = 'EW([O II])=...' else $
         ewmgii = 'EW([Mg II])='+strtrim(string(ages.mgii_2800_ew[0],format='(F12.1)'),2)+'\pm'+$
         strtrim(string(ages.mgii_2800_ew[1],format='(F12.1)'),2)+' \AA; '+'S/N='+$
         strtrim(string(ages.mgii_2800[0]/ages.mgii_2800[1],format='(F12.2)'),2)+'; '+$
         '\sigma='+strtrim(string(ages.mgii_2800_sigma[0],format='(F12.1)'),2)+' km/s'
    endif else ewmgii = ''
    if tag_exist(ages,'oii_3727_ew') then begin
       if (ages.oii_3727_ew[1] lt 0.0) then ewoii = 'EW([O II])=...' else $
         ewoii = 'EW([O II])='+strtrim(string(ages.oii_3727_ew[0],format='(F12.1)'),2)+'\pm'+$
         strtrim(string(ages.oii_3727_ew[1],format='(F12.1)'),2)+' \AA; '+'S/N='+$
         strtrim(string(ages.oii_3727[0]/ages.oii_3727[1],format='(F12.2)'),2)+'; '+$
         '\sigma='+strtrim(string(ages.oii_3727_sigma[0],format='(F12.1)'),2)+' km/s'
    endif else ewoii = ''
    if tag_exist(ages,'oiii_5007_ew') then begin
       if (ages.oiii_5007_ew[1] lt 0.0) then ewoiii = 'EW([O III])=...' else $
         ewoiii = 'EW([O III])='+strtrim(string(ages.oiii_5007_ew[0],format='(F12.1)'),2)+'\pm'+$
         strtrim(string(ages.oiii_5007_ew[1],format='(F12.1)'),2)+' \AA; S/N='+$
         strtrim(string(ages.oiii_5007[0]/ages.oiii_5007[1],format='(F12.2)'),2)+'; '+$
         '\sigma='+strtrim(string(ages.oiii_5007_sigma[0],format='(F12.1)'),2)+' km/s'
    endif else ewoiii = ''
    if tag_exist(ages,'nii_6584_ew') then begin
       if (ages.nii_6584_ew[1] lt 0.0) then ewnii = 'EW([N II])=...' else $
         ewnii = 'EW([N II])='+strtrim(string(ages.nii_6584_ew[0],format='(F12.1)'),2)+'\pm'+$
         strtrim(string(ages.nii_6584_ew[1],format='(F12.1)'),2)+' \AA; S/N='+$
         strtrim(string(ages.nii_6584[0]/ages.nii_6584[1],format='(F12.2)'),2)+'; '+$
         '\sigma='+strtrim(string(ages.nii_6584_sigma[0],format='(F12.1)'),2)+' km/s'
    endif else ewnii = ''
    if tag_exist(ages,'h_beta_ew') then begin
       if (ages.h_beta_ew[1] lt 0.0) then ewhb = 'EW(H\beta)=...' else $
         ewhb = 'EW(H\beta)='+strtrim(string(ages.h_beta_ew[0],format='(F12.1)'),2)+'\pm'+$
         strtrim(string(ages.h_beta_ew[1],format='(F12.1)'),2)+' \AA; S/N='+$
         strtrim(string(ages.h_beta[0]/ages.h_beta[1],format='(F12.2)'),2)+'; '+$
         '\sigma='+strtrim(string(ages.h_beta_sigma[0],format='(F12.1)'),2)+' km/s'
    endif else ewhb = ''
    if tag_exist(ages,'h_alpha_ew') then begin
       if (ages.h_alpha_ew[1] lt 0.0) then ewha = 'EW(H\alpha)=...' else $
         ewha = 'EW(H\alpha)='+strtrim(string(ages.h_alpha_ew[0],format='(F12.1)'),2)+'\pm'+$
         strtrim(string(ages.h_alpha_ew[1],format='(F12.1)'),2)+' \AA; S/N='+$
         strtrim(string(ages.h_alpha[0]/ages.h_alpha[1],format='(F12.2)'),2)+'; '+$
         '\sigma='+strtrim(string(ages.h_alpha_sigma[0],format='(F12.1)'),2)+' km/s'
    endif else ewha = ''

    label = nicegalaxy
    case labeltype of
       0: label1 = ''
       1: label1 = ['AGES '+nicegalaxy,'z = '+zobj_str]
       2: begin
          label1 = ['AGES '+nicegalaxy,'z = '+zobj_str]
          label2 = [ewmgii,ewoii,ewhb,ewoiii,ewha,ewnii]
       end
       3: label1 = ['AGES '+nicegalaxy,'ID = '+agesid,'z = '+zobj_str]
       4: label1 = ['AGES '+nicegalaxy,'z = '+zobj_str,'spz = '+spzobj_str]
       5: label1 = ['AGES '+nicegalaxy,'ID = '+agesid,'z = '+zobj_str,'I = '+imag]
       else: label1 = ''
    endcase

    if keyword_set(bigpostscript) or keyword_set(postscript) then $
      lcharsize = 1.3 else lcharsize = 1.6
    legend, textoidl(label1), /left, /top, box=0, charsize=lcharsize, spacing=1.6, $
      /normal, margin=0, clear=keyword_set(postscript) or keyword_set(bigpostscript)
    if (n_elements(label2) ne 0L) then legend, textoidl(label2), /right, /top, box=0, $
      charsize=lcharsize, spacing=1.6, /normal, margin=0, charthick=1.3, $
      clear=keyword_set(postscript) or keyword_set(bigpostscript)
       
    if keyword_set(postscript) then begin
       im_plotfaves
       dfpsclose
    endif

return
end
    
