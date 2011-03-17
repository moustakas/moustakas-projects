function read_stelib, silent=silent
; jm03feb18uofa

    datapath = getenv('STELIB_PATH')
    pushd, datapath
    
    speclist = findfile('*.fits',count=nspec)

; read the first spectrum to initialize arrays

    if not keyword_set(silent) then splog, 'Reading '+speclist[0]+'.'
    spec = mrdfits(speclist[0],0,header,/silent)
    if (size(spec,/n_dimension) ne 1L) then begin
       splog, 'Spectrum is not one-dimensional!'
       retall
    endif
    
    specsize = size(spec,/dimension)
    npix = specsize[0]
    
    speccube = {specname: '',           $ ; spectrum name
                datapath: '',           $ ; data path
                star    : '',           $ ; object name
                ra      : '',           $ ; 
                dec     : '',           $ ; 
 ;              header  : strarr(39),   $ ; header
                npix    : 0L,           $ ; number of pixels
                type    : '',           $ ; spectral type
                A_V     : 0.0,          $ ; interstellar A_V
                EBV     : 0.0,          $
                U       : 0.0,          $
                B       : 0.0,          $
                V       : 0.0,          $
                R       : 0.0,          $
;               I       : 0.0,          $
                radvel  : 0.0,          $ ; radial velocity
                teff    : -999.0,       $ ; effective temperature
                logg    : -999.0,       $ ; log gravity
                fe_h    : -999.0,       $ ; metallicity
                spec    : fltarr(npix), $ ; spectrum
                wave    : fltarr(npix)}   ; wavelength vector
    if nspec gt 1L then speccube = replicate(speccube,nspec)

    speccube[0].specname = speclist[0]
    speccube[0].datapath = datapath
    speccube[0].star = strn(sxpar(header,'STARIDEN'))
    speccube[0].ra = im_dec2hms(sxpar(header,'RA')/15.0)
    speccube[0].dec = im_dec2hms(sxpar(header,'DEC'))
;   speccube[0].header = header
    speccube[0].npix = long(sxpar(header,'NAXIS1'))
    speccube[0].type = strn(sxpar(header,'SPECTYPE'))
    speccube[0].A_V = sxpar(header,'ISEXTINC')
    speccube[0].U = sxpar(header,'U')
    speccube[0].B = sxpar(header,'B')
    speccube[0].V = sxpar(header,'V')
    speccube[0].R = sxpar(header,'R')
;   speccube[0].I = sxpar(header,'I')
    speccube[0].radvel = sxpar(header,'RADVEL')

    teff = sxpar(header,'TEFF')
    if strmatch(string(teff),'*-*') ne 1B then speccube[0].teff = teff
    logg = sxpar(header,'LOG_G')
    if strmatch(string(logg),'*-*') ne 1B then speccube[0].logg = logg
    fe_h = sxpar(header,'FE_H')
    if strmatch(string(fe_h),'*-*') ne 1B then speccube[0].fe_h = fe_h

    speccube[0].spec = spec
    speccube[0].wave = make_wave(header)

; read each spectrum

    for i = 1L, nspec-1L do begin
       
       if not keyword_set(silent) then splog, 'Reading ', speclist[i]+'.'
       spec = mrdfits(speclist[i],0,header,/silent)

       szt = size(spec,/dimension)
       if (size(spec,/n_dimension) ne 1L) then begin
          splog, 'Spectrum is not one-dimensional!'
          heap_gc
          retall
       endif else if (specsize[0] ne szt[0]) then begin
          splog, 'Spectra are not the same dimension!'
          heap_gc
          retall
       endif
       
       speccube[i].specname = speclist[i]
       speccube[i].datapath = datapath
       speccube[i].star = strn(sxpar(header,'OBJECT'))
       speccube[i].ra = im_dec2hms(sxpar(header,'RA')/15.0)
       speccube[i].dec = im_dec2hms(sxpar(header,'DEC'))
;      speccube[i].header = header
       speccube[i].npix = long(sxpar(header,'NAXIS1'))
       speccube[i].type = strn(sxpar(header,'SPECTYPE'))
       speccube[i].A_V = sxpar(header,'ISEXTINC')
       speccube[i].U = sxpar(header,'U')
       speccube[i].B = sxpar(header,'B')
       speccube[i].V = sxpar(header,'V')
       speccube[i].R = sxpar(header,'R')
;      speccube[i].I = sxpar(header,'I')
       speccube[i].radvel = sxpar(header,'RADVEL')
       speccube[i].spec = spec
       speccube[i].wave = make_wave(header)

       teff = sxpar(header,'TEFF')
       if strmatch(string(teff),'*-*') ne 1B then speccube[i].teff = teff
       logg = sxpar(header,'LOG_G')
       if strmatch(string(logg),'*-*') ne 1B then speccube[i].logg = logg
       fe_h = sxpar(header,'FE_H')
       if strmatch(string(fe_h),'*-*') ne 1B then speccube[i].fe_h = fe_h

    endfor

    popd
    
return, speccube
end
