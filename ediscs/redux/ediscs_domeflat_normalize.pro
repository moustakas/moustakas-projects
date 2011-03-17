;+
; NAME:
;       EDISCS_DOMEFLAT_NORMALIZE
;
; PURPOSE:
;
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; PROCEDURES USED:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Sep 12, U of A
;       jm04nov18uofa - normalize at 6780 Angstroms
;-

pro ediscs_domeflat_normalize, speclist, domeflat, datapath=datapath, $
  domepath=domepath, outpath=outpath, doplot=doplot, wfits=wfits, $
  gzip=gzip, _extra=extra

; SPECLIST - one-dimensional spectra
; DOMEFLAT - one-dimensional dome flat spectrum    
    
    nspec = n_elements(speclist)
    ndome = n_elements(domeflat)

    if (nspec eq 0L) or (ndome eq 0L) then begin
       print, 'Syntax - ediscs_domeflat_normalize, speclist, domeflat'
       return
    endif

    if (n_elements(datapath) eq 0L) then datapath = cwd()
    if (n_elements(domepath) eq 0L) then domepath = datapath

    outname = 'd'+repstr(speclist,'.gz','') ; root output file names

; read the dome flat
    
    if (file_test(domepath+domeflat,/regular) eq 0L) then begin
       splog, 'Unable to find dome flat file '+domepath+domeflat+'.'
       return
    endif

    splog, 'Reading '+domepath+domeflat+'.'

; read and normalize

    domeflux = readfits(domepath+domeflat,domehead,/silent)
    domewave = make_wave(domehead)
    
; loop on each spectrum
    
    for i = 0L, nspec-1L do begin

       spec = rd1dspec(speclist[i],datapath=datapath,/silent)

       flux = spec.spec
       ferr = spec.sigspec
       sky = spec.sky
       mask = spec.mask
       wave = spec.wave
       header = spec.header

       newdomeflux = interpol(domeflux,domewave,wave)
;      newdomeflux = interpol(domeflux,domewave,wave) / interpol(domeflux,domewave,6780.0)
;      newdomeflux = interpol(domeflux,domewave,wave) / (interpol(domeflux,domewave,mean(wave)))[0]
       
       newflux = flux / newdomeflux
       newferr = ferr / newdomeflux

       if keyword_set(doplot) then begin
          yrange = [min(flux)<min(newflux),max(flux)>max(newflux)]
          djs_plot, wave, flux, ps=10, xsty=3, ysty=3, yrange=yrange, $
            charthick=2.0, charsize=2.0, xthick=2.0, ythick=2.0, thick=2.0
          djs_oplot, wave, newflux, ps=10, color='blue', thick=2.0
          if (nspec gt 1L) then cc = get_kbrd(1)
       endif
       
; update the header

       sxaddpar, header, 'DOMEFLAT', domeflat, ' dome flat name', before='HISTORY'
       sxaddhist, "'Normalized by the dome flat "+im_today()+"'", header
       
       if keyword_set(wfits) then begin

          splog, 'Writing '+datapath+outname[i]+'.'
          wrt1dspec, outname[i], newflux, newferr, sky, mask, $
            header, datapath=outpath, gzip=gzip
          
       endif

    endfor

return
end
