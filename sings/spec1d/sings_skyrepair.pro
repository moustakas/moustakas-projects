;+
; NAME:
;       SINGS_SKYREPAIR
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
;       singsinfo - output from SINGS_READ_INFO()
;
; KEYWORD PARAMETERS:
;       nuclear    - repair the nuclear spectra
;       drift20    - repair the drift20 spectra
;       drift56    - repair the drift56 spectra
;       maskoi     - the [O I] 6300 night sky line is particularly
;                    tricky because it overlaps the nebular line in
;                    low-redshift objects; if this keyword is set then
;                    prompt the user to overwrite the default of not
;                    masking sky pixels near [O I]
;       repair     - write out new FITS files
;       postscript - write postscript output
;       restore    - restore the original spectra stored in the
;                    tarballs that are created when REPAIR=1; useful
;                    if you messed up!
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       This routine, of course, assumes that the extracted spectra
;       have been fitted using ISPECLINEFIT() first.  Also, this
;       routine should only be run once.  See also ATLAS1D_SKYREPAIR. 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Aug 25, U of A
;       jm06jan29uofa - added DRIFT56 and MASKOI keywords
;       jm08jan21nyu - added RESTORE keyword
;
; Copyright (C) 2005-2006, 2008, John Moustakas
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

pro sings_skyrepair, sings, nuclear=nuclear, drift20=drift20, $
  drift56=drift56, maskoi=maskoi, repair=repair, postscript=postscript, $
  restore=restore

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

    datapath = sings_path(/spec1d)
    analysis_path = sings_path(/analysis)
    backuppath = sings_path(/spec1d)+'beforeskyrepair/'

    sings = read_sings(nuclear=nuclear,drift20=drift20,drift56=drift56)
    ngalaxy = n_elements(sings)

    specfile = strtrim(sings.specfile,2)
    if keyword_set(nuclear) then begin
       suffix = 'nuclear'
    endif
    if keyword_set(drift20) then begin
       suffix = 'drift20'
    endif
    if keyword_set(drift56) then begin
       suffix = 'drift56'
    endif

; restore the tarball, if requested

    if keyword_set(restore) then begin
       splog, 'Restore '+backuppath+'sings_beforeskyrepair_'+suffix+'.tar.gz [Y/N]?'
       cc = get_kbrd(1)
       if (strupcase(cc) eq 'Y') then begin
          pushd, datapath
          spawn, 'tar xzvf '+backuppath+'sings_beforeskyrepair_'+suffix+'.tar.gz', /sh
          popd
       endif
       return
    endif

; make a tarball backup copy of every spectrum before overwriting 

    if keyword_set(repair) then begin
       pushd, datapath
       spawn, ['tar czvf '+backuppath+'sings_beforeskyrepair_'+suffix+$
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
       psname = analysis_path+'sings_skyrepair_'+suffix+'.ps'
       dfpsplot, psname, /color, /square
       postthick = 5.0
    endif else begin
       im_window, 0, xratio=0.5, /square
       postthick = 2.0
    endelse
    
    plotsym, 8, 1.0, /fill
    scale = 1E15

; loop on each galaxy in the sings    
    
;   skywaves = 5577.345
;   skywaves = 5889.950
    skywaves = [5577.345,5889.950,6300.32,6363.81]

;   for i = 16, 20 do begin
    for i = 0L, ngalaxy-1L do begin

       skywaves1 = skywaves

       galaxy = strtrim(sings[i].galaxy,2)

       scube = rd1dspec(specfile[i],datapath=datapath)
       wave = scube.wave
       flux = scube.spec
       ferr = scube.sigspec
       mask = bytarr(scube.npix)
       header = scube.header

       zobj = sings[i].z

; mask emission lines and the region around Na D

       emask = emission_mask(wave,z=zobj,width=width,/noskymask)
       if keyword_set(maskoi) then $
         emask[where((wave gt oiwave-width) and (wave lt oiwave+width))] = 1B
       
; check to see if this object was fitted (AGN were not); if so, then
; restore the best-fitting continuum and compute the residuals;
; otherwise use the observed spectrum
       
       if keyword_set(nuclear) then agnflag = sings[i].nuclear_agnflag 
       if keyword_set(drift20) then agnflag = sings[i].drift20_agnflag 
       if keyword_set(drift56) then agnflag = sings[i].drift56_agnflag 

       if (agnflag eq 0L) then begin

          specdata = read_sings_specfit(galaxy,nuclear=nuclear,$
            drift20=drift20,drift56=drift56,/silent)
          bestfit = specdata[*,2]/(1+zobj)
;         bestfit = (specdata[*,2]+specdata[*,3])/(1+zobj)
          inputflux = flux - bestfit

; add pixels that are in emission around Na D to the interpolation
; mask since *actual* Na D is always in absorption

          nad = where(((wave gt 5890.0-width) and (wave lt 5890.0+width)) or $
            ((wave gt 5896.0-width) and (wave lt 5896.0+width)) and $
            (flux gt (bestfit+2.0*ferr)),nnad)
          if (nnad ne 0L) then emask[nad] = 1B

          pix = where(emask eq 1B,comp=epix)
          inputmask = mask[pix]
          outmask = mask

          newflux = flux
          tempflux = iskymask(inputflux[pix],ferr[pix],wave[pix],$
            mask=inputmask,skywaves=skywaves1,nsig=1.0,doplot=debug)
          outmask[pix] = inputmask
          thesepix = where(outmask ne 0B,nthesepix,comp=otherpix)

; add some noise to make it look real          

          if (nthesepix ne 0L) then newflux[thesepix] = bestfit[thesepix] + $
            randomn(seed,nthesepix)*djs_median(ferr[otherpix])
          
       endif else begin

          pix = where(emask eq 1B,comp=epix)
          inputmask = mask[pix]
          outmask = mask

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

;      if strmatch(galaxy,'*1097*') then stop

       if keyword_set(repair) then begin
          sxaddhist, "'Sky-subtraction residuals repaired "+im_today()+"'", header
;         sxaddpar, header, 'BITPIX', -32, 'FITS BITS/PIXEL'
          modfits, datapath+specfile[i], float(newflux), header, exten_no=0L
       endif
       
; make a diagnostic plot       
       
       djs_plot, wave, scale*flux, xsty=3, ysty=3, ps=10, xr=minmax(skywaves)+[-100,+100], $; [6200,6400], $ ;[5500,5700], $ ;
         xthick=postthick, ythick=postthick, charsize=1.5, charthick=postthick, $
         thick=postthick, xtitle='Wavelength [\AA]', ytitle='Relative Flux'
       djs_oplot, wave, scale*newflux, ps=10, color='orange', thick=postthick
       if (agnflag eq 0L) then djs_oplot, wave, scale*bestfit, ps=10, color='green'
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
