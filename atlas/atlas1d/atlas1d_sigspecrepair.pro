;+
; NAME:
;       ATLAS_SIGSPECREPAIR
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
;       routine should only be run once.  
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

pro atlas1d_sigspecrepair, atlas, nuclear=nuclear, repair=repair, restore=restore, $
  postscript=postscript

    datapath = atlas_path(/atlas1d)
    analysis_path = atlas_path(/analysis)
    backuppath = atlas_path(/atlas1d)+'beforesigspecrepair/'

    if keyword_set(nuclear) then atlas = read_nuclear() else atlas = read_integrated()
    ngalaxy = n_elements(atlas)

    specfile = strtrim(atlas.specfile,2)
    if keyword_set(nuclear) then begin
       suffix = 'nuclear'
    endif else begin
       suffix = 'integrated'
    endelse

; restore the tarball, if requested

    if keyword_set(restore) then begin
       splog, 'Restore '+backuppath+'atlas1d_beforesigspecrepair_'+suffix+'.tar.gz [Y/N]?'
       cc = get_kbrd(1)
       if (strupcase(cc) eq 'Y') then begin
          pushd, datapath
          spawn, 'tar xzvf '+backuppath+'atlas1d_beforesigspecrepair_'+suffix+'.tar.gz', /sh
          popd
       endif
       return
    endif

; backup the spectra
    
    if keyword_set(repair) then begin
       pushd, datapath
       spawn, ['tar czvf '+backuppath+'atlas1d_beforesigspecrepair_'+suffix+$
         '.tar.gz '+strjoin(specfile,' ')], /sh
       popd
    endif

; setup some plotting variables
    
    if keyword_set(repair) then postscript = 1L

    if keyword_set(postscript) then begin
       psname = analysis_path+'atlas1d_sigspecrepair_'+suffix+'.ps'
       dfpsplot, psname, /color, /square
       postthick1 = 4.0
       postthick2 = 3.0
    endif else begin
       im_window, 0, xratio=0.8
       postthick1 = 2.0
       postthick2 = 2.0
    endelse

    pagemaker, nx=1, ny=2, xspace=0.0, yspace=0.0, xmargin=[1.1,0.2], $
      ymargin=[0.2,1.1], position=pos, /normal

    scale = 1E15
    charsize1 = 1.6
    charsize2 = 1.6

    allfactor = fltarr(ngalaxy)
    
;   for igal = 16, ngalaxy-1 do begin
    for igal = 0L, ngalaxy-1L do begin

       print, format='("Object ",I0,"/",I0,".",A1,$)', igal+1, ngalaxy, string(13b)

       specfit = read_atlas_specfit(atlas[igal].galaxy,nuclear=nuclear,/silent)
       wave = reform(specfit[*,0])
       flux = reform(specfit[*,1])
       cfit = reform(specfit[*,2])
       efit = reform(specfit[*,3])
       mfit = reform(specfit[*,4])
       model = (efit+cfit+mfit)

       z_abs = atlas[igal].z_abs
       thisspecfile = datapath+strtrim(atlas[igal].specfile,2)
       
       header = headfits(thisspecfile,ext=0)
       errheader = headfits(thisspecfile,ext=1)
       ferr = mrdfits(thisspecfile,1,/silent)*(1+z_abs)

; I want the median noise to be given by the rms of the
; model-subtracted (residual) spectrum, but to retan the *shape* of
; the S/N spectrum given by ispec
       
       resid = flux-model
       stats = im_stats(resid,sigrej=4.0)

       histogauss, scale*resid, gterms, xresid, yresid, /noplot       
;      autohist, scale*resid, junk1, junk2, xresid, yresid, /noplot ; use smart binning 
;      yfit = mpfitpeak(xresid,yresid,gterms,/positive,/gaussian)

       noise = gterms[2]/scale
;      noise = stats.sigrej
       
       factor = noise/median(ferr)
;      print, atlas[igal].galaxy, gterms[2], djsig(resid,sigrej=3.0)*scale, factor
       allfactor[igal] = factor
       
       newferr = ferr*factor
       factor_str = strtrim(string(factor,format='(F12.2)'),2)

; write out

       if keyword_set(repair) then begin
          sxaddpar, header, 'EFACTOR', float(factor), ' error spectrum scale factor', before='HISTORY'
          sxaddhist, "'Error spectrum scaled by empirical factor "+im_today()+"'", header
;         sxaddpar, header, 'BITPIX', -32, 'FITS BITS/PIXEL'
          modfits, thisspecfile, 0, header, exten_no=0L
          modfits, thisspecfile, float(newferr), exten_no=1L
       endif
       
; make a diagnostic plot       
       
       djs_plot, wave, scale*flux, ps=10, xsty=3, ysty=3, thick=postthick2, position=pos[*,0], $
         xthick=postthick1, ythick=postthick1, charsize=charsize1, charthick=postthick2, $
         xtitle='', xtickname=replicate(' ',10), ytitle='Flux (10^{-15} '+flam_units()+')', $
         yrange=scale*[stats.minrej,max(model)]
       djs_oplot, wave, scale*model, ps=10, color='red';, thick=postthick2
       djs_oplot, wave, scale*resid, ps=10, color='grey'
       legend, atlas[igal].galaxy, /left, /top, box=0, charsize=charsize1, charthick=postthick2, $
         clear=keyword_set(postscript)

; optional: plot the histogram distribution       
;      histogauss, scale*resid, gterms, xresid, yresid, /noerase, $
;        position=[0.75,0.6,0.95,0.75]

;;       autohist, scale*resid, junk1, junk2, xresid, yresid, /noerase, $
;;;        position=[pos[2,0]-0.5,pos[1,0],pos[2,0]-0.05,pos[3,0]-0.5]
;;         position=[0.75,0.6,0.95,0.75]
;;;      djs_plot, xresid, yresid, ps=10, xsty=3, ysty=3
;;       djs_oplot, xresid, yfit, color='red', thick=2, ps=10

       djs_plot, wave, flux/ferr, ps=10, xsty=3, ysty=3, /noerase, thick=2.0, position=pos[*,1], $
         xthick=postthick1, ythick=postthick1, charsize=charsize1, charthick=postthick2, $
         xtitle='Wavelength (\AA)', ytitle='S/N', color='grey', $
         yrange=[-5.0,(max(flux/ferr)>max(flux/noise))>max(flux/ferr/factor)]
       djs_oplot, wave, flux/noise, ps=10, color='orange'
       djs_oplot, wave, flux/ferr/factor, ps=10, color='dark green'

       legend, textoidl(['S/N_{ispec}','S/N_{empirical}','S/N_{ispec}/'+factor_str]), $
         charsize=charsize2, charthick=postthick2, /left, /top, box=0, $
         clear=keyword_set(postscript), spacing=2.0, margin=0, $
         textcolor=djs_icolor(['grey','orange','dark green'])       
       
       if (not keyword_set(postscript)) then cc = get_kbrd(1)
       
    endfor
    
    if keyword_set(postscript) then begin
       dfpsclose
       spawn, 'gzip -f '+psname, /sh
    endif

    srt = sort(allfactor)
    niceprint, atlas[srt].galaxy, allfactor[srt]
    
return
end


