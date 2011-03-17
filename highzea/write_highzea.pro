;+
; NAME:
;       WRITE_HIGHZEA
;
; PURPOSE:
;       Write the final two-dimensional 2D spectra to a common
;       directory (called by WRITE_HIGHZEA_WRAPPER).
;
; CALLING SEQUENCE:
;       write_highzea, specfile, datapath=, outpath=, $
;          outlist=, /wfits, /overwrite, /upcase
;
; INPUTS:
;       specfile - text file containing four columns: (1) iSPEC2d
;                   FITS file name; (2) output FITS file name; (3)
;                   galaxy name; (4) Y or N indicating whether the
;                   data were taken during photometric conditions 
;
; OPTIONAL INPUTS:
;       datapath - path to the input FITS files (column 1, above) 
;       outpath  - output path name (where the files in column 2,
;                  above, are written) 
;
; KEYWORD PARAMETERS:
;       wfits     - write the output FITS file, otherwise just print
;                   the output name
;       overwrite - force overwriting if the output file name already
;                   exists (not recommended); the default is to create
;                   a unique output file name using "_1", "_2", etc.
;       upcase    - use upper case when writing out the FITS file
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       outlist   - output file list
;
; COMMENTS:
;       Also see WRITE_HIGHZEA_WRAPPER and WRITE_ATLAS2D.
;
;       This routine stores the absolute uncertainty in the
;       spectrophotometry in magnitudes in each header.  This
;       uncertainty includes the absolute error in the
;       spectrophotometric standard-star calibrations, the uncertainty
;       in the zero-point shift error, and the median uncertainty in
;       the *observed* sensitivity function, all added in quadrature.  
;
;       If the directory specified by ATLASPATH has old versions of
;       the same images then they need to be deleted manually because
;       this routine will not overwrite images with the same file
;       name, unless the OVERWRITE is set, which in general is not
;       recommended. 
; 
;       STITCHED spectra were wavelength calibrated and *then* sky
;       subtracted (see SINGS2D_STITCH).
;
; PROCEDURES USED:
;       CWD(), ATLAS_PATH(), READCOL, RD2DSPEC(), SXADDHIST, SXADDPAR,
;       WRT2DSPEC, ICLEANUP, SPLOG, DJS_MEAN() 
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2006 Jun 29, U of A - based on WRITE_ATLAS2D
;
; Copyright (C) 2006, John Moustakas
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

pro write_highzea, specfile, datapath=datapath, outpath=outpath1, $
  outlist=outlist, wfits=wfits, overwrite=overwrite, upcase=upcase
    
    if (n_elements(specfile) eq 0L) then begin
       print, 'Syntax - write_highzea, specfile, datapath=, outpath=, $'
       print, '   outlist=, /wfits, /overwrite, /upcase'
       return
    endif

    if (n_elements(datapath) eq 0L) then datapath = cwd()
    if (n_elements(outpath1) eq 0L) then outpath = highzea_path(/spec2d) else outpath = outpath1
    
    if (file_test(datapath+specfile,/regular) eq 0B) then begin
       splog, 'File '+datapath+specfile+' not found...returning.'
       return
    endif

    readcol, datapath+specfile, alist, outnamelist, galaxylist, photflag, $
      /silent, format='A,A,A,A', comment='#', delimiter=' '
    nlist = n_elements(alist)
    outlist = strarr(nlist)

; absolute spectrophotometric calibration error

    caliberr1 = 0.0215 ; (=2%) absolute error in the CALSPEC calibrations [mag]
;   caliberr2 = 0.02   ; uncertainty in the zero point shift [mag]

;   caliberr = sqrt(caliberr1^2 + caliberr2^2) ; total error
    caliberr = caliberr1

    repeatlist = ''
    
    for k = 0L, nlist-1L do begin

       cube = rd2dspec(alist[k],datapath=datapath,wset=wset)
       h = *cube.header

       if tag_exist(cube,'TELLURIC_SPEC') then begin
          telluric_spec = cube.telluric_spec
          telluric_head = cube.telluric_head
       endif

       if keyword_set(upcase) then $
         rootname = strupcase(outnamelist[k]) else $
         rootname = strlowcase(outnamelist[k])
       outname = rootname

; add the calibration error in quadrature with the median absolute
; spectrophotometric error to get the total absolute error

       sxaddpar, h, 'CALIBERR', float(caliberr), format='(F12.4)', $
         ' spectrophotometric calibration error [mag]', before='HISTORY'

       amederr = sxpar(h,'AMEDERR',count=count) ; [mag]
       if (count eq 0L) then amederr = 0.0

       abserror = sqrt(amederr^2.0 + caliberr^2.0)
       sxaddpar, h, 'ABSERROR', float(abserror), format='(F12.4)', $
         ' absolute spectrophotometric error [mag]', before='HISTORY'
       
       sxaddpar, h, 'PHOTFLAG', strtrim(photflag[k],2), ' photometric flag [Y/N]', before='HISTORY'

; add the drift-scan length to the file name

       slit = 'SLIT2.0'
       if (strmatch(outname,'*'+slit+'*',/fold) eq 0L) then count = 0L else count = 1L

; do not overwrite the same file

       if keyword_set(overwrite) then begin

          if (count ne 0L) then outname = outname+'_'+slit

       endif else begin
       
          if (count eq 0L) then checkname = outname+'*' else checkname = outname+'*_'+slit
          files = file_search(outpath+checkname+'.fits',count=fcount)

          if (fcount eq 0L) then $
            if (count eq 0L) then outname = rootname else outname = outname;+'_'+slit

; rename the first existing file to be file number one

          if (strmatch(outname,'*1104*') eq 1B) then stop
       
          if (fcount eq 1L) then begin
             if keyword_set(wfits) then begin
                if (count eq 0L) then begin
                   spawn, ['\mv '+outpath+checkname+'.fits '+outpath+outname+'_1.fits'], /sh
                   repeatlist = [repeatlist,outname+'_1']
                endif else begin
                   spawn, ['\mv '+outpath+checkname+'.fits '+outpath+outname+'_1_'+slit+'.fits'], /sh
                   repeatlist = [repeatlist,outname+'_1_'+slit]
                endelse
             endif
             if (count eq 0L) then outname = outname+'_2' else outname = outname+'_2_'+slit
             repeatlist = [repeatlist,outname]
          endif
          
; do not overwrite any existing files

          if (fcount gt 1L) then begin
             if (count eq 0L) then $
               outname = outname+'_'+string(fcount+1,format='(I0)') else $
               outname = outname+'_'+string(fcount+1,format='(I0)')+'_'+slit
             repeatlist = [repeatlist,outname]
          endif

       endelse

       outname = outname+'.fits'
       outlist[k] = outname

; update the header

       sxaddhist, "'New copy of "+strtrim(alist[k],2)+"'", h
;      sxaddhist, "'New copy of "+datapath+alist[k]+"'", h
       galaxy = strupcase(strcompress(strjoin(strsplit(galaxylist[k],'_',/extract)),/remove))
       sxaddpar, h, 'GALAXY', galaxy, after = 'OBJECT'
       h = h[where(strcompress(h,/remove) ne '')]
       
       if keyword_set(wfits) then begin

          splog, 'Writing '+outpath+outname
          wrt2dspec, outname, cube.image, cube.sigmamap, cube.mask, h, $
            skyimage=cube.sky, wset=wset, telluric_spec=telluric_spec, $
            telluric_head=telluric_head, datapath=outpath, gzip=gzip
          
       endif else splog, 'Image '+outpath+outname
       
       icleanup, cube

    endfor

return
end
