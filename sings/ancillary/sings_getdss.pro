pro sings_getdss, object, dsswidth=dsswidth
; jm05mar29uofa

; retrieve DSS images based on the IRS LL/SL centers, otherwise use
; the NED positions

    dsspath = sings_path()+'DSS/'
    if file_test(dsspath,/directory) ne 0L then begin
       splog, 'Would you like to remove all DSS files from '+dsspath+' [Y/N]?'
       cc = get_kbrd(1)
       if strupcase(cc) eq 'Y' then spawn, ['/bin/rm -f '+dsspath+'*.fits.gz'], /sh
    endif
    
    sings = sings_read_info()
    irs = rsex(sings_path(/observing)+'irs_ll_params.txt')

    if n_elements(object) ne 0L then begin
       doit = match_string(object,sings.galaxy,index=match) ; match by galaxy name
       if match[0] eq -1L then message, 'Object '+object+' not found!'
       sings = sings[match]
    endif
    nsings = n_elements(sings)

    sizes = [5.0,10.0,15.0,20.0,30.0] ; possible image sizes [arcmin]
    
    pushd, dsspath
    for j = 0L, nsings-1L do begin ; loop on each object

       galaxy = strcompress(strlowcase(sings[j].nice_galaxy),/remove)

;      singsgalaxy = strcompress(strlowcase(sings[j].nice_galaxy),/remove)
;      match = where(strmatch(strtrim(irs.nice_galaxy,2),singsgalaxy,/fold) eq 1B,nmatch)
;      if (nmatch gt 1L) then message, 'Problem here!'
;      if (nmatch eq 1L) then begin
;         splog, 'Querying IRS-centered image for '+strupcase(galaxy)+'.'
;         radec = [15.0*im_hms2dec(irs[match].ra),im_hms2dec(irs[match].dec)]
;      endif
;      if (nmatch eq 0L) then begin
;         splog, 'Querying NED-centered image for '+strupcase(galaxy)+'.'
;         radec = [15.0*im_hms2dec(sings[j].ra),im_hms2dec(sings[j].dec)]
;      endif

       radec = [15.0*im_hms2dec(sings[j].ra),im_hms2dec(sings[j].dec)]

; compute the required FITS image size based on the size of the galaxy 

       imsize = ceil(1.2*sings[j].d25_maj) ; [arcmin]
       
;      dmaj = ceil(sings[j].d25_maj) ; [arcmin]
;      if (dmaj/5 eq 0) or (dmaj eq -999) then imsize = sizes[0]
;      if (dmaj/5 eq 1) then imsize = sizes[1]
;      if (dmaj/5 gt 1) and (dmaj/5 lt 3) then imsize = sizes[2]
;      if (dmaj/5 ge 3) and (dmaj/5 lt 5) then imsize = sizes[3]
;      if (dmaj/5 ge 5) then imsize = sizes[4]

       if (n_elements(dsswidth) ne 0L) then imsize = dsswidth
       
       querydss, radec, im, h, imsize=imsize, /ned, /eso

       writefits, galaxy+'.fits', im, h
       spawn, ['gzip -f '+galaxy+'.fits'], /sh
       
    endfor
    popd

return
end
