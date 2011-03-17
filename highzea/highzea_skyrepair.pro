;+
; NAME:
;       HIGHZEA_SKYREPAIR
;
; PURPOSE:
;       Repair sky subtraction residuals in the 1D spectra by
;       identifying residuals near known sky wavelengths between the
;       data and the best-fitting BC03 spectrum.  For galaxies that
;       have not been fitted (ie, AGN), repair sky-subtraction
;       residuals in the traditional way using ISKYMASK(). 
;
; CALLING SEQUENCE:
;       highzea_skyrepair, highzeainfo, /maskoi, /repair, /postscript
;
; INPUTS:
;
; OPTIONAL INPUTS:
;       highzeainfo - output from HIGHZEA_READ_INFO()
;
; KEYWORD PARAMETERS:
;       repair     - write out new FITS files
;       postscript - write postscript output
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       HIGHZEA_PATH(), HIGHZEA_READ_INFO(), DFPSPLOT, IM_WINDOW,
;       RD1DSPEC(), EMISSION_MASK(), READ_HIGHZEA_SPECFIT(), ISKYMASK(),
;       SXADDHIST, DJS_MODFITS, DJS_PLOT, DJS_OPLOT, LEGEND, ICLEANUP,
;       DFPSCLOSE 
;
; COMMENTS:
;       This routine, of course, assumes that the extracted spectra
;       have been fitted using ISPECLINEFIT() first.  Also, this
;       routine should only be run once.  See also ATLAS1D_SKYREPAIR
;       and SINGS_SKYREPAIR.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2006 July 7, U of A - based on SINGS_SKYREPAIR
;
;-

pro highzea_skyrepair, highzeainfo, repair=repair, postscript=postscript

    datapath = highzea_path(/spec1d)
    analysis_path = highzea_path(/analysis)
    backuppath = highzea_path(/spec1d)+'beforeskyrepair/'

    skywaves = [5577.345,6300.32]
;   skywaves = [5577.345,5889.950,6300.32,6363.81]

    specfile = file_basename(file_search(datapath+'*.ms.fits'))
    highzeainfo1 = iforage(datapath+specfile)
    ngalaxy = n_elements(highzeainfo1)

    if keyword_set(repair) then begin
       pushd, datapath
       spawn, ['tar czvf '+backuppath+'highzea_beforeskyrepair.tar.gz '+strjoin(specfile,' ')], /sh
       popd
    endif

; setup some plotting variables

    if keyword_set(repair) then postscript = 1L
    
    if keyword_set(postscript) then begin
       dfpsplot, analysis_path+'highzea_skyrepair.ps', /color, xsize=8.5, ysize=8.5;, /square
       postthick = 5.0
       postthick2 = 2.0
    endif else begin
       im_window, 0, xratio=0.5, /square
       postthick = 2.0
       postthick2 = 1.0
    endelse
    
    plotsym, 8, 0.5, /fill
    scale = 1E17

    pagemaker, nx=1, ny=2, xspace=0.0, yspace=0.0, xmargin=[1.1,0.4], ymargin=[0.4,1.1], $
      xpage=8.5, ypage=8.5, width=7.0, height=[3.5,3.5], position=pos, /normal
    
    xtitle = 'Wavelength [\AA]'
    ytitle='Relative Flux'

;   xrange = minmax(skywaves)+[-100,+100]
;   xrange = minmax(skywaves)+[-100,+100]
    
; loop on each galaxy in the highzea    
    
;   for i = 16, 20 do begin
    for i = 0L, ngalaxy-1L do begin

       skywaves1 = skywaves

       galaxy = strtrim(highzeainfo1[i].galaxy,2)

       scube = rd1dspec(specfile[i],datapath=datapath)
       wave = scube.wave
       flux = scube.spec
       ferr = scube.sigspec
       mask = bytarr(scube.npix)
       header = scube.header

;      zobj = highzeainfo1[i].z
       zobj = 0.0

; mask emission lines and the region around Na D
       
;      emask = emission_mask(wave,z=zobj,width=width,/noskymask)
;      pix = where(emask eq 1B,comp=epix)
;      inputmask = mask[pix]
       
       epix = -1L
       pix = lindgen(n_elements(wave))
       inputmask = mask[pix]
       outmask = mask

; check to see if this object was fitted (AGN were not); if so, then
; restore the best-fitting continuum and compute the residuals;
; otherwise use the observed spectrum

       agnflag = 1L
       
;      if keyword_set(nuclear) then $
;        agnflag = highzeainfo1[i].nuclear_agnflag else $
;        agnflag = highzeainfo1[i].drift20_agnflag

       if (agnflag eq 0L) then begin

;         specdata = read_highzea_specfit(galaxy,nuclear=nuclear,$
;           drift20=drift20,drift56=drift56,/silent)
          bestfit = specdata[*,2]/(1+zobj)
;         bestfit = (specdata[*,2]+specdata[*,3])/(1+zobj)
          inputflux = flux - bestfit

          newflux = flux
          tempflux = iskymask(inputflux[pix],ferr[pix],wave[pix],$
            mask=inputmask,skywaves=skywaves1,nsig=1.0,doplot=debug)
          outmask[pix] = inputmask
          thesepix = where(outmask ne 0B,nthesepix,comp=otherpix)

; add some noise to make it look real          

          if (nthesepix ne 0L) then newflux[thesepix] = bestfit[thesepix] + $
            randomn(seed,nthesepix)*djs_median(ferr[otherpix])
          
       endif else begin

          inputflux = flux
          newflux = flux

          tempflux = iskymask(inputflux[pix],ferr[pix],wave[pix],$
            mask=inputmask,skywaves=skywaves1,nsig=1.0,doplot=debug)
          thesepix1 = where(inputmask ne 0B,nthesepix1)
          if (nthesepix1 ne 0L) then newflux[pix[thesepix1]] = tempflux[thesepix1]

          outmask[pix] = inputmask
          thesepix = where(outmask ne 0B,nthesepix)
          
       endelse

; write out

       if keyword_set(repair) then begin
          sxaddhist, "'Sky-subtraction residuals repaired "+im_today()+"'", header
          djs_modfits, datapath+specfile[i], float(newflux), header, exten_no=0L
       endif
       
; make a diagnostic plot       
       
       djs_plot, wave, scale*flux, xsty=3, ysty=3, ps=10, xrange=xrange, $
         xthick=postthick, ythick=postthick, charsize=1.5, charthick=postthick, $
         thick=postthick2, xtickname=replicate(' ',10), xtitle='', ytitle=ytitle, $
         position=pos[*,0];, color='grey'
       if (nthesepix ne 0L) then djs_oplot, wave[thesepix], scale*flux[thesepix], ps=8, color='red'
       if (epix[0] ne -1L) then djs_oplot, wave[epix], scale*flux[epix], ps=8, color='dark blue'
       legend, [galaxy], /right, /top, box=0, charsize=1.6, charthick=postthick

       djs_plot, wave, scale*newflux, xsty=3, ysty=3, ps=10, /noerase, xrange=xrange, $
         xthick=postthick, ythick=postthick, charsize=1.5, charthick=postthick, $
         thick=postthick2, xtitle=xtitle, ytitle=ytitle, position=pos[*,1], color='grey'

       if (not keyword_set(postscript)) then cc = get_kbrd(1)
       icleanup, scube
       
    endfor

    if keyword_set(postscript) then begin
       dfpsclose
       spawn, ['gzip -f '+analysis_path+'highzea_skyrepair.ps'], /sh
    endif
    
return
end
