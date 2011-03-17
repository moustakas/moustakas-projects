pro read_catalog, catalog, oflux, filters, nsources, catinfo
;+
; NAME:
;	READ_CATALOG
;
; PURPOSE:
;
; CALLING SEQUENCE:
;
; INPUTS:
;	catalog - name of the photometric catalog
;
; OPTIONAL INPUTS:
;	
; OUTPUTS:
;
; COMMENTS:
;
; PROCEDURES USED:
;	READFAST
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 September 12, U of A
;-

    cpath = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='catalogs')

    if not keyword_set((findfile(cpath+catalog))[0]) then begin
       print, 'Catalog '+catalog+' not found!'
       return
    endif

; read the catalog.  skip the first comment line and read in the
; filter names
    
    header = ''
    print, 'Reading '+cpath+catalog
    openr, lun, cpath+catalog, /get_lun
    for i = 0L, 1L do readf, lun, header

    if strmid(header,0,1) ne '#' then $
      message, 'The first line needs to have a comment character.', /info else $
      header = strmid(header,1) ; remove the comment

    filters = strsplit(header,',',/extract) ; parse the filter names
    nbands = n_elements(filters)
    for i = 0L, nbands-1L do filters[i] = strn(filters[i])
    
; need error checking to make sure the filter names match up with the
; known filters

; read the catalog
    
    readfast, cpath+catalog, photocat, nlines=nsources, ncols=ncols, skip=2
    print, format='("There are ",I0," sources in ",A0," observed in ",I0," filters:")', $
      nsources, strlowcase(catalog), nbands
    print, strn(header)
    
; error checking
    
;    if ncols lt nbands then message, 'The number of filters does not match the number of columns!' else $
;      if ncols gt nbands+1L then
;    
;    if ncols gt nbands+1L then begin
;       message, 'There are too many columns in this catalog!'
;       return, -1
;    endif

    oflux = fltarr(2,nbands,nsources)
    oflux[0,*,*] = photocat[0:nbands-1L,*]
    oflux[1,*,*] = photocat[nbands:2L*nbands-1L,*]

; read the catalog information file

    infofile = strmid(catalog,0,strpos(catalog,'.'))+'.info'

    catinfo = {id   : lonarr(nsources), $
               xpix : fltarr(nsources), $
               ypix : fltarr(nsources), $
               ra   : fltarr(nsources), $
               dec  : fltarr(nsources), $
               zspec: fltarr(nsources), $
               type : lonarr(nsources)}
       
    if not keyword_set((findfile(cpath+infofile))[0]) then begin

       print, 'Catalog information file '+infofile+' not found.' 
       catinfo.id = lindgen(nsources)+1L

    endif else begin

       readcol, cpath+infofile, id, xpix, ypix, ra, dec, zspec, type, format='L,F,F,F,F,F,L', /silent
       
       catinfo.id = id
       catinfo.xpix = xpix
       catinfo.ypix = ypix
       catinfo.ra = ra
       catinfo.dec = dec
       catinfo.zspec = zspec
       catinfo.type = type

    endelse

return
end

