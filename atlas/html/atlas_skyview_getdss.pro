pro webatlas_skyview_getdss, object
; jm02feb18uofa
; retrieve FITS DSS images of all the ATLAS2d galaxies using SKYVIEW
; in batch mode

; requires SKYVIEW_BATCH
    
    propath = '/home/ioannis/kennicutt/webatlas/'
    fitspath = propath+'FITS'
    if file_test(fitspath,/directory) ne 0L then begin
       splog, 'Would you like to remove all FITS files from '+fitspath+' [Y/N]?'
       cc = get_kbrd(1)
       if strupcase(cc) eq 'Y' then spawn, ['/bin/rm -f '+fitspath+'/*.fits.gz'], /sh
    endif

    atlaspath = atlas_path(/analysis1d)
    atlas = mrdfits(atlaspath+'atlas1d_info.fits.gz',1,/silent)
    keep = where(atlas.drift eq 'Y',natlas) ; only keep integrated spectra
    if natlas ne 0L then atlas = atlas[keep]
    
    if n_elements(object) ne 0L then begin

       doit = match_string(object,atlas.galaxy,index=match) ; match by galaxy name
       if match[0] eq -1L then message, 'Object '+object+' not found!'
       atlas = atlas[match]

    endif
    natlas = n_elements(atlas)
    
    sizes = [5.0,10.0,15.0] ; possible image sizes [arcmin]

    fitsfiles = strcompress(strlowcase(atlas.galaxy),/remove)+'.fits'

    galaxy = strcompress(strlowcase(atlas.galaxy),/remove)

    remletters = ['se','nw','s','e','w','n']
    
    radec = strarr(natlas)
    sfactr = fltarr(natlas)
    for i = 0L, natlas-1L do begin

       okay = 0
       j = 0L
       
       radec[i] = strjoin(strsplit(strjoin([atlas[i].ra,atlas[i].dec],', '),':',/extract),' ')
       
       dmaj = ceil(atlas[i].d25_maj) ; [arcmin]
       if dmaj lt -900.0 then sfactr[i] = sizes[0]/60.0 else begin
       
          if dmaj/5 eq 0 then sfactr[i] = sizes[0]/60.0
          if dmaj/5 eq 1 then sfactr[i] = sizes[1]/60.0
          if dmaj/5 gt 1 then sfactr[i] = sizes[2]/60.0

       endelse

       last = strlen(galaxy[i])
       while (okay eq 0L) and (j lt n_elements(remletters)) do begin
       
          if strmatch(galaxy[i],'*'+remletters[j],/fold) eq 1B then begin
             galaxy[i] = strmid(galaxy[i],0,strpos(galaxy[i],remletters[j]))
             okay = 1L
          endif

          j = j+1L
          
       endwhile

;      if strmatch(galaxy[i],'*ic2810*') eq 1B then stop

    endfor

; not necessary since we want unique FITS file names (e.g., ic2810se.fits)    
    
;;    srt = sort(galaxy)
;;    galaxy = galaxy[srt]
;;    
;;    unique = uniq(galaxy)
;;    natlas = n_elements(unique)
;;    
;;    galaxy = galaxy[unique]
;;    fitsfiles = fitsfiles[srt[unique]]
;;    sfactr = sfactr[srt[unique]]

;   survey = 'Digitized Sky Survey'
;   survey = '2MASS-J'
    survey = 'DSS2 Red'
;   survey = 'DSS2 Blue'

    t1 = systime(/seconds)
    skyview_batch, fitsfiles, galaxy, survey=survey, sfactr=sfactr, datapath=fitspath+'/'
;   skyview_batch, fitsfiles, radec, survey=survey, sfactr=sfactr, datapath=fitspath+'/'
;   skyview_batch, fitsfiles, radec, survey=survey, datapath=atlas_path(/web)+'dss/'
    print, 'Time elapsed (minutes) = ', (systime(/seconds)-t1)/60.0

; gzip all the data

    pushd, fitspath
    for k = 0L, natlas-1L do spawn, ['gzip '+fitsfiles[k]], /sh
    popd

return
end
