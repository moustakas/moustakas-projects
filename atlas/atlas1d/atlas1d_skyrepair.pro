;+
; NAME:
;       ATLAS1D_SKYREPAIR
;
; PURPOSE:
;       Repair sky subtraction residuals in the 1D spectra by
;       identifying residuals near known sky wavelengths between the
;       data and the best-fitting BC03 spectrum.  For galaxies that
;       have not been fitted (ie, AGN), repair sky-subtraction
;       residuals in the traditional way using ISKYMASK(). 
;
; INPUTS:
;
; OPTIONAL INPUTS:
;       atlasinfo - output from ATLAS_READ_INFO()
;
; KEYWORD PARAMETERS:
;       nuclear    - repair the nuclear spectra (default is to repair
;                    the integrated spectra)
;       repair     - write out new FITS files
;       postscript - write postscript output
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       This routine, of course, assumes that the extracted spectra
;       have been fitted using ISPECLINEFIT() first.  Also, this
;       routine should only be run once.  
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Jul 30, U of A
;       jm08jan21nyu - added RESTORE keyword
;
; Copyright (C) 2005, 2008, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

pro atlas1d_skyrepair, atlas, nuclear=nuclear, repair=repair, $
  maskoi=maskoi, postscript=postscript, restore=restore

    datapath = atlas_path(/atlas1d)
    analysis_path = atlas_path(/analysis)
    backuppath = atlas_path(/atlas1d)+'beforeskyrepair/'

    if keyword_set(nuclear) then atlas = read_nuclear() else atlas = read_integrated()
;   atlas = atlas[speclinefit_locate(atlas,'09425')]
    ngalaxy = n_elements(atlas)

    specfile = strtrim(atlas.specfile,2)
    if keyword_set(nuclear) then begin
       suffix = 'nuclear'
    endif else begin
       suffix = 'integrated'
    endelse

; restore the tarball, if requested

    if keyword_set(restore) then begin
       splog, 'Restore '+backuppath+'atlas1d_beforeskyrepair_'+suffix+'.tar.gz [Y/N]?'
       cc = get_kbrd(1)
       if (strupcase(cc) eq 'Y') then begin
          pushd, datapath
          spawn, 'tar xzvf '+backuppath+'atlas1d_beforeskyrepair_'+suffix+'.tar.gz', /sh
          popd
       endif
       return
    endif

; make a tarball backup copy of every spectrum before overwriting 

    if keyword_set(repair) then begin
       pushd, datapath
       spawn, ['tar czvf '+backuppath+'atlas1d_beforeskyrepair_'+suffix+$
         '.tar.gz '+strjoin(specfile,' ')], /sh
       popd
    endif

    if keyword_set(maskoi) then begin
       splog, 'Masking the [O I] night sky line, which may affect the nebular line at low redshift!'
       oiwave = 6300.32
    endif

; setup some plotting variables

    if keyword_set(repair) then postscript = 1L
    
    if keyword_set(postscript) then begin
       psname = analysis_path+'atlas1d_skyrepair_'+suffix+'.ps'
       dfpsplot, psname, /color, /square
       postthick = 5.0
    endif else begin
       im_window, 0, xratio=0.5, /square
       postthick = 2.0
    endelse
    
    plotsym, 8, 1.0, /fill
    scale = 1E15

; loop on each galaxy in the atlas    
    
;   skywaves = 5577.345
    skywaves = [5577.345,5889.950,6300.32,6363.81]

;   for i = 20, 20 do begin
    for i = 0L, ngalaxy-1L do begin

       skywaves1 = skywaves

; -------------------------
; integrated spectra
; -------------------------

       galaxy = strtrim(atlas[i].galaxy,2)

       scube = rd1dspec(specfile[i],datapath=datapath)
       wave = scube.wave
       flux = scube.spec
       ferr = scube.sigspec
       mask = bytarr(scube.npix)
       header = scube.header

       zobj = atlas[i].z_abs

; mask emission lines and the region around Na D

       emask = emission_mask(wave,z=zobj,width=width,/noskymask)
       if keyword_set(maskoi) then $
         emask[where((wave gt oiwave-width) and (wave lt oiwave+width))] = 1B

       pix = where(emask eq 1B,comp=epix)
       inputmask = mask[pix]
       outmask = mask

; check to see if this object was fitted (AGN were not); if so, then
; restore the best-fitting continuum and compute the residuals;
; otherwise use the observed spectrum
       
       if keyword_set(nuclear) then $
         agnflag = atlas[i].nuclear_agnflag else $
         agnflag = atlas[i].drift_agnflag

       if (agnflag eq 0L) then begin

          specdata = read_atlas_specfit(galaxy,nuclear=nuclear,/silent)
          bestfit = (specdata[*,2]+specdata[*,3])/(1+zobj)
          inputflux = flux - bestfit

; add pixels that are in emission around Na D to the interpolation
; mask since *actual* Na D is always in absorption

          nad = where(((wave gt 5890.0-width) and (wave lt 5890.0+width)) or $
            ((wave gt 5896.0-width) and (wave lt 5896.0+width)) and $
            (flux gt (bestfit+2.0*ferr)),nnad)
          if (nnad ne 0L) then emask[nad] = 1B

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
          modfits, datapath+specfile[i], float(newflux), header, exten_no=0L
       endif
       
; make a diagnostic plot       
       
       djs_plot, wave, scale*flux, xsty=3, ysty=3, ps=10, xr=minmax(skywaves)+[-100,+100], $ ; [5500,5700], $ ;
         xthick=postthick, ythick=postthick, charsize=1.5, charthick=postthick, $
         thick=postthick, xtitle='Wavelength [\AA]', ytitle='Relative Flux'
       djs_oplot, wave, scale*newflux, ps=10, color='orange', thick=postthick
       if (nthesepix ne 0L) then djs_oplot, wave[thesepix], scale*flux[thesepix], ps=8, color='red'
       if (epix[0] ne -1L) then djs_oplot, wave[epix], scale*flux[epix], ps=8, color='dark blue'
       legend, [galaxy,'z = '+strtrim(string(zobj,format='(F12.4)'),2)], /right, $
         /top, box=0, charsize=1.3, charthick=postthick, clear=keyword_set(postscript)
       if (not keyword_set(postscript)) then cc = get_kbrd(1)

       icleanup, scube
       
    endfor

    if keyword_set(postscript) then begin
       dfpsclose
       spawn, 'gzip -f '+psname, /sh
    endif
    
return
end
