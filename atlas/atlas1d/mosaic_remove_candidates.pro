pro mosaic_remove_candidates, atlas, postscript=postscript
; jm05jul27uofa - make a mosaic plot of the candidates for removal
;                 from the atlas

    analysis_path = atlas_path(/analysis)

    if (n_elements(atlas) eq 0L) then atlas = mrdfits(analysis_path+'atlas1d_info.fits.gz',1,/silent)

    readcol, analysis_path+'remove_candidates.txt', remfile, format='A', /silent
    ngalaxy = n_elements(remfile)

    match, remfile, strtrim(atlas.drift_file,2), indx1, remindx
    remfile = remfile[indx1]
    rematlas = atlas[remindx]

    if keyword_set(postscript) then dfpsplot, analysis_path+'remove_candidates.ps', /landscape, /color
    
;   for i = 0L, 3L do begin
    for i = 0L, ngalaxy-1L do begin
       print, format='("Object ",I0,"/",I0,".",A1,$)', i+1, ngalaxy, string(13b)
       atlas_display_spectrum, rematlas[i], lcharsize=1.2, postscript=postscript
       legend, rematlas[i].galaxy, /left, /top, charsize=1.5, charthick=2.0, box=0
       if (not keyword_set(postscript)) then cc = get_kbrd(1)
    endfor

    if keyword_set(postscript) then dfpsclose

return
end
