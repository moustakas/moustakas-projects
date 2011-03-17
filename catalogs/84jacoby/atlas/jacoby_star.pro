function jacoby_star, starname, jacobypath=jacobypath, print=print
; jm02jan25uofa
; parse the jacoby atlas and retrieve the spectrum of a named star in
; a data structure

    nstar = n_elements(starname)
    if (nstar eq 0L) and (not keyword_set(print)) then begin
       print, 'Syntax - star = jacoby_star(starname,[jacobypath=])'
       return, 0L
    endif

    if n_elements(jacobypath) eq 0L then jacobypath = '/home/ioannis/idl/jacoby_atlas/'

;   openr, lun, jacobypath+'/jacoby_atlas_info.dat', /get_lun
;   readu, lun, name, sp, lum, format='(A10,A4,A3)'
;   readcol, jacobypath+'/jacoby_atlas_info.dat', name, sp, lum, format='(A10,A4,A3)'
;   free_lun, lun
;   reads, info, name, sp, lum, format='(A10,A4,A3)'

    atlasinfo = jacobypath+'/jacoby_atlas_info.dat'
    if file_test(atlasinfo) ne 1L then begin
       print, 'Unable to find '+atlasinfo+'.'
       return, 0L
    endif

    info = djs_readilines(atlasinfo)
    nlines = n_elements(info)

    starinfo = strarr(3,nlines) ; name, spectral type, luminosity class
    for i = 0L, nlines-1L do starinfo[*,i] = [strmid(info[i],0,10),  $
                                              strmid(info[i],10,4), $
                                              strmid(info[i],15,3)]

; print star information if requested
    
    if keyword_set(print) then begin
       print, '     Star     S. Type   L. Class'
       print, '---------------------------------'
       for i = 0L, nlines-1L do print, starinfo[*,i], format='(A12,2x,A6,2x,A5)'
    endif

    if nstar eq 0L then return, 0L else begin

       match = indgen(nstar)
       for j = 0L, nstar-1L do begin
          matchname = match_standard(starname[j],starinfo[0,*],index=index)
          match[j] = index
       endfor

       if nstar eq 1L then match = match[0]

    endelse

; wavelength vector
    
    wave0 = 3510.0D ; starting wavelength [Angstrom]
    dwave = 1.40D   ; dispersion [Angstrom/pixel]
    nwave = 400L*7L-1L
    wave = wave0 + findgen(nwave)*dwave

    atlasdata = jacobypath+'/jacoby_atlas.dat'
    if file_test(atlasdata) ne 1L then begin
       print, 'Unable to find '+atlasdata+'.'
       return, 0L
    endif

    sdata = {star:  '', $
             type:  '', $
             class: '', $
             wave:  wave, $
             flux:  fltarr(nwave)}
    sdata = replicate(sdata,nstar)

    openr, lun1, atlasdata, /get_lun
    stat = fstat(lun1)

    match = match[sort(match)] ; sort

    gotoline = lonarr(nstar)
    gotoline[0] = 400L*match[0]
    gotoline[1L:nstar-1L] = (400L*match - shift(400L*match,1))[1L:nstar-1L]
    
    for k = 0L, nstar-1L do begin

       linenum = gotoline[k]
       dum = strarr(linenum)
       readf, lun1, dum
       dum = ''

       dat = strarr(400) ; read the whole star data as a string
       readf, lun1, dat

       flux = fltarr(7,399)
       reads, dat, flux, format='(10x,7E10.3)' ; read the first 399 lines

       moreflux = fltarr(6)
       reads, dat[399], moreflux, format='(10x,6E10.3)'
       
       spec = [reform(flux,7*399),moreflux] ; spectrum

       sdata[k].star = strcompress(starname[k],/remove)
       sdata[k].type = strcompress(starinfo[1,match[k]],/remove)
       sdata[k].class = strcompress(starinfo[2,match[k]],/remove)
       sdata[k].flux = spec

    endfor
    free_lun, lun1
    
return, sdata
end
