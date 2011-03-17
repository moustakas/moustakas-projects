pro atlas1d_snsort, write=write
; jm03may19uofa - sort by S/N
; jm03sep17uofa - new purpose: compare ispec S/N with
;                 the noise in the continuum 

    atlaspath = atlas_path(/atlas1d)
    pushd, atlaspath

    flist = findfile('*drift*.ms.fits.gz',count=fcount)
    snispec = fltarr(fcount)
    snstat = fltarr(fcount)
    galaxy = strarr(fcount)
    date = strarr(fcount)
    
    wave1 = 5300.0
    wave2 = 6400.0
    
    for i = 0L, fcount-1L do begin

       scube = rd1dspec(flist[i])

       galaxy[i] = sxpar(scube.header,'galaxy')
       date[i] = sxpar(scube.header,'date-obs')
       
       get_element, scube.wave, [wave1,wave2], xx

       snispec[i] = median(scube.spec[xx[0]:xx[1]]/scube.sigspec[xx[0]:xx[1]])
       stats = im_stats(scube.spec[xx[0]:xx[1]],sigrej=3.0)
       snstat[i] = stats.median/stats.sigma_rej
       
    endfor
    
    popd

    diff = snispec-snstat
    ratio = snispec/snstat

;   srt = sort(diff)
;   niceprint, galaxy[srt], date[srt], snispec[srt], snstat[srt], diff[srt]

    srt = sort(ratio)
    for i = 0L, fcount-1L do print, strcompress(galaxy[srt[i]],/remove), $
      strcompress(date[srt[i]],/remove), snispec[srt[i]], snstat[srt[i]], ratio[srt[i]], $
      format='(A20,A13,3F12.5)'

    if keyword_set(write) then begin
       
       openw, lun, 'atlas1d_snsort.txt', /get_lun
       printf, lun, '# Galaxy; Date; S/N [iSPEC]; S/N [Statistical]; Ratio '
       for i = 0L, fcount-1L do printf, lun, strcompress(galaxy[srt[i]],/remove), $
         strcompress(date[srt[i]],/remove), snispec[srt[i]], snstat[srt[i]], ratio[srt[i]], $
         format='(A20,A13,3F12.5)'
       free_lun, lun

    endif
       
stop
    
return
end
    
