pro kenn92_getdss, object
; jm04jan16uofa

; use the QUERYDSS Goddard routine to grab DSS images from STScI;
; write and GZIP the FITS files in
; /home/ioannis/kennicutt/projects/kenn92/analysis/DSS.

    propath = atlas_path(/kenn92)+'webkenn92/'
    fitspath = propath+'DSS'
    if file_test(fitspath,/directory) ne 0L then begin
       splog, 'Would you like to remove all DSS files from '+fitspath+' [Y/N]?'
       cc = get_kbrd(1)
       if strupcase(cc) eq 'Y' then spawn, ['/bin/rm -f '+fitspath+'/*.fits.gz'], /sh
    endif
    
    data = read_kenn92()
    
    if n_elements(object) ne 0L then begin

       doit = match_string(object,data.galaxy,index=match) ; match by galaxy name
       if match[0] eq -1L then message, 'Object '+object+' not found!'
       data = data[match]

    endif
    ndata = n_elements(data)

    sizes = [5.0,10.0,15.0] ; possible image sizes [arcmin]
    
    pushd, fitspath
    for j = 0L, ndata-1L do begin ; loop on each object

       galaxy = strcompress(strlowcase(data[j].galaxy),/remove)
       radec = [15.0*im_hms2dec(data[j].ra),im_hms2dec(data[j].dec)]

; compute the required FITS image size based on the size of the galaxy 

       dmaj = ceil(data[j].d25_maj) ; [arcmin]
       if (dmaj/5 eq 0) or (dmaj eq -999) then imsize = sizes[0]
       if dmaj/5 eq 1 then imsize = sizes[1]
       if dmaj/5 gt 1 then imsize = sizes[2]

       splog, 'Querying '+strupcase(galaxy)+'.'
       querydss, radec, im, h, imsize=imsize, survey='2r', /ned
       writefits, galaxy+'.fits', im, h
       spawn, ['gzip -f '+galaxy+'.fits'], /sh
       
    endfor
    popd

return
end
