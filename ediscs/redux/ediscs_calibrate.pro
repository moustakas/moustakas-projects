;+
; NAME:
;       EDISCS_CALIBRATE
;
; PURPOSE:
;       Flux-calibrate the ediscs spectra.
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
;       J. Moustakas, 2004 Feb 16, U of A - written
;       jm04oct28uofa - updated
;-

pro ediscs_calibrate, speclist, sensname, prefix=prefix, datapath=datapath, $
  senspath=senspath, outpath=outpath, gzip=gzip, wfits=wfits

    nspec = n_elements(speclist)
    nsens = n_elements(sensname)

    if (nspec eq 0L) or (nsens eq 0L) then begin
       print, 'Syntax - ediscs_calibrate'
       return
    endif

    if (n_elements(datapath) eq 0L) then datapath = cwd()
    if (n_elements(senspath) eq 0L) then senspath = datapath
    if (n_elements(outpath) eq 0L) then outpath = datapath
    if (n_elements(extfile) eq 0L) then extfile = 'lasillaextinct.dat'

    if (n_elements(prefix) eq 0L) then prefix = '' ; f

    outname = prefix+repstr(speclist,'.gz','') ; root output file names

; restore the sensitivity function    
    
    if (file_test(senspath+sensname,/regular) eq 0L) then begin
       splog, 'Unable to find sensitivity function '+sensname+'.'
       return
    endif

    splog, 'Reading '+sensname+'.'

    tsens = readfits(senspath+sensname,senshead,/silent)
    sens = 10D^(0.4*tsens)
    senswave = make_wave(senshead)
    
; read in the extinction curve

    extpath = getenv('ISPEC_DIR')+'/etc/'
    if (file_test(extpath+extfile,/regular) eq 0L) then begin
       splog, 'Unable to find extinction file '+extfile+'.'
       return
    endif
          
    splog, 'Reading the extinction file '+extpath+extfile+'.'
    readcol, extpath+extfile, extwave, extvals, format='D,D', /silent

; loop on each spectrum
    
    for i = 0L, nspec-1L do begin

       flux = mrdfits(datapath+speclist[i],0,header,/silent)
       ferr = mrdfits(datapath+speclist[i],1,/silent)
       wave = make_wave(header,cd1_1=dwave)
       exptime = sxpar(header,'EXPTIME')
       airmass= sxpar(header,'AIRMASS')

       splog, 'Extinction correcting and flux calibrating '+speclist[i]+'.'
       extinct = interpol(extvals,extwave,wave)
       eflux = flux * 10D^(0.4*airmass*extinct)
       eferr = ferr * 10D^(0.4*airmass*extinct)

       sensfunc = interpol(sens,senswave,wave)
       sflux = eflux / exptime / dwave / sensfunc
       sferr = eferr / exptime / dwave / sensfunc

; update the header

       sxaddpar, header, 'FLUXCOR', 'T', ' flux-calibrated [T/F]', before='HISTORY'
       sxaddpar, header, 'ZUNITS', 'erg/s/cm2/A', ' flux units', before='HISTORY'
       sxaddpar, header, 'SENSNAME', sensname, ' sensitivity function name', before='HISTORY'
       sxaddhist, "'Flux calibrated "+im_today()+"'", header
       
       if keyword_set(wfits) then begin

          splog, 'Writing '+outpath+outname[i]+'.'
          mwrfits, sflux, outpath+outname[i], header, /create
          mwrfits, sferr, outpath+outname[i]
          if keyword_set(gzip) then spawn, ['gzip -f '+outpath+outname[i]], /sh
          
       endif

    endfor

return
end
