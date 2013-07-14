pro hogg_rc3_ned, info, gal=gal
; jm07dec18nyu

;   path = '/Users/ioannis/ay/research/catalogs/hogg_rc3/'
    path = './'
    
    allfile = file_search(path+'*_irg.jpg',count=cc)
    file = file_basename(allfile)

    gal = strarr(cc)
    
    for ii = 0L, cc-1L do begin

       gal[ii] = strcompress(strjoin((strsplit(file[ii],'_',/ext))[0L:1L],''),/remove)
       if strmatch(gal[ii],'A*') then gal[ii] = '[RC2]'+gal[ii]

       if strmatch(file[ii],'Pegasus*') then $
         gal[ii] = strcompress(strjoin((strsplit(file[ii],'_',/ext))[0L:0L],''),/remove)+'Dwarf'
       if strmatch(file[ii],'MCG*') then $
         gal[ii] = repstr(strcompress(strjoin((strsplit(file[ii],'_',/ext))[0L:3L],'-'),/remove),'--','-')
       if strmatch(file[ii],'Ursa*') then $
;        gal[ii] = strcompress(strjoin((strsplit(file[ii],'_',/ext))[0L:1L],''),/remove)
         gal[ii] = 'DDO199'

    endfor

;   niceprint, gal

;   ned_webget_basic, gal, info
    info = struct_addtags(replicate({file: ''},cc),info)
    info.file = file
;   mwrfits, info, path+'hogg_rc3_ned.fits', /create
;   spawn, 'gzip -f '+path+'hogg_rc3_ned.fits', /sh
    
return
end
    
    
