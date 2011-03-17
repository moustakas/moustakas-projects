pro pilyugin_getdss, object, dsswidth=dsswidth
; jm05mar11uofa

; use the QUERYDSS Goddard routine to grab DSS images from STScI

    root = '04pilyugin'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'

    dsspath = path+'DSS'

    if file_test(dsspath,/directory) ne 0L then begin
       splog, 'Would you like to remove all DSS files from '+dsspath+' [Y/N]?'
       cc = get_kbrd(1)
       if strupcase(cc) eq 'Y' then spawn, ['/bin/rm -f '+dsspath+'/*.fits.gz'], /sh
    endif

; read the catalog

    pilyugin = read_04pilyugin()
    
    if (n_elements(object) ne 0L) then begin
       doit = match_string(object,pilyugin.galaxy,index=match) ; match by galaxy name
       if match[0] eq -1L then message, 'Object '+object+' not found!'
       pilyugin = pilyugin[match]
    endif

    npilyugin = n_elements(pilyugin)

    sizes = [5.0,10.0,15.0,20.0,30.0] ; possible image sizes [arcmin]
    
    for j = 0L, npilyugin-1L do begin ; loop on each object

       galaxy = strcompress(strlowcase(pilyugin[j].galaxy),/remove)
       radec = [15.0*im_hms2dec(pilyugin[j].ra),im_hms2dec(pilyugin[j].dec)]

; compute the required FITS image size based on the size of the galaxy 

       dmaj = ceil(pilyugin[j].dmaj) ; [arcmin]
       if (dmaj/5 eq 0) or (dmaj eq -999) then imsize = sizes[0]
       if (dmaj/5 eq 1) then imsize = sizes[1]
       if (dmaj/5 gt 1) and (dmaj/5 lt 3) then imsize = sizes[2]
       if (dmaj/5 ge 3) and (dmaj/5 lt 5) then imsize = sizes[3]
       if (dmaj/5 ge 5) then imsize = sizes[4]

       if (n_elements(dsswidth) ne 0L) then imsize = dsswidth
       
       splog, 'Querying '+strupcase(galaxy)+'.'
       querydss, radec, im, h, imsize=imsize, survey='2r', /ned
       writefits, dsspath+'/'+galaxy+'.fits', im, h
       spawn, ['gzip -f '+dsspath+'/'+galaxy+'.fits'], /sh
       
    endfor
    
return
end
