;+
; NAME:
;       MULTIPLOT_BC03_TEMPLATES
;
; PURPOSE:
;       Generate a multi-panel plot for the atlas paper.
;
; CALLING SEQUENCE:
;       multiplot_bc03_templates, /salpeter, /postscript
;
; INPUTS:
;       None
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       postscript - generate multi-page postscript output 
;
; OUTPUTS:
;
; PROCEDURES USED:
;       REPSTR(), MRDFITS(), IM_OPENCLOSE, DJS_PLOT, OPLOT,
;       MINMAX(), LEGEND
;
; COMMENTS:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Apr 12, U of A - written
;-

pro multiplot_bc03_templates, salpeter=salpeter, postscript=postscript, paper=paper
    
    Z = ['Z004','Z02','Z05']
    Zstr = 'Z = '+['0.004','0.02','0.05']
    nZ = n_elements(Z)

    if keyword_set(salpeter) then begin
       splog, 'Selecting SALPETER IMF.'
       imf = 'salpeter'
    endif else begin
       splog, 'Selecting CHABRIER IMF.'
       imf = 'chabrier'
    endelse

    eigendir = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='templates')
    eigenfile = 'BC03_'+Z+'_'+imf+'_templates.fits'

    psname = 'BC03_templates.ps'
    if keyword_set(paper) then begin
       postscript = 1
       pspath = atlas_path(/papers)+'atlas/FIG_ATLAS/' 
    endif else pspath = atlas_path(/plots)

    im_openclose, pspath+psname, postscript=postscript, /color

    pagemaker, nx=1, ny=nZ, position=pos, xmargin=[0.5,0.1], $
      ymargin=[0.1,1.0], yspace=0.05, /normal
    
    loadct, 13

    xrange = [3200,7300]
    yrange = [0.1,9]
    
    for i = 0L, nZ-1L do begin

       if file_test(eigendir+eigenfile[i]) eq 0B then begin

          splog, 'Template file '+eigendir+eigenfile[i]+' does not exist.'

       endif else begin
          
          info = mrdfits(eigendir+eigenfile[i],1,/silent)
          flux = mrdfits(eigendir+eigenfile[i],2,header,/silent)
          wave = make_wave(header)
          ntemplate = info.ntemplate
          age = info.template_age/1E6
          strage = string(round(info.template_age/1E6),format='(I5)')

          newflux = fltarr(5000,ntemplate)
          for it = 0L, ntemplate-1L do newflux[*,it] = frebin(im_normalize(flux[*,it],wave,normwave=5500.0),5000)

          wave = frebin(wave,5000)
          flux = newflux

          colors = fix(255*findgen(ntemplate)/ntemplate)

          case i of
             0L: begin
                djs_plot, wave, flux[*,0], xsty=3, ysty=3, charsize=1.5, charthick=2.0, $
                  xthick=2.0, ythick=2.0, yrange=yrange, xrange=xrange, xtickname=replicate(' ',10), $
                  position=pos[*,0], thick=4.0, color=colors[0], /ylog
                for k = 1L, ntemplate-1L do djs_oplot, wave, flux[*,k], ps=10, $
                  thick=4.0, color=colors[k]
                legend, strage+' Myr', /right, box=0, charsize=1.5, charthick=2.0, $
                  /top, textcolors=colors
             end
             1L: begin
                djs_plot, wave, flux[*,0], xsty=3, ysty=3, charsize=1.5, charthick=2.0, $
                  xthick=2.0, ythick=2.0, yrange=yrange, xrange=xrange, xtickname=replicate(' ',10), $
                  position=pos[*,1], thick=4.0, color=colors[0], /noerase, /ylog
                for k = 1L, ntemplate-1L do djs_oplot, wave, flux[*,k], ps=10, $
                  thick=4.0, color=colors[k]
             end
             2L: begin
                djs_plot, wave, flux[*,0], xsty=3, ysty=3, charsize=1.5, charthick=2.0, $
                  xthick=2.0, ythick=2.0, yrange=yrange, xrange=xrange, /noerase, /ylog, $
                  position=pos[*,2], thick=4.0, color=colors[0], xtitle='Wavelength (\AA)'
                for k = 1L, ntemplate-1L do djs_oplot, wave, flux[*,k], ps=10, $
                  thick=4.0, color=colors[k]
             end
          endcase

          legend, Zstr[i], /left, /top, box=0, charsize=1.5, charthick=2.0

;            ytitle='Normalized Flux Density'

       endelse 

    endfor 

    im_openclose, postscript=postscript, /close

return
end
