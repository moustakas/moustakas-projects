pro compare_jacoby
; jm03feb3uofa

    speclist = findfile('*.ms.fits')
    forage = iforage(speclist)
    
    star = strcompress(forage.object,/remove)
    
; read the JACOBY stellar atlas
    
    jacoby = read_jacoby_atlas()

    doit = match_string(star,jacoby.star,findex=match,/exact)
    good = where(match ne -1L,nstar)
    jacoby = jacoby[match[good]]

    star = star[good]
    speclist = speclist[good]
    forage = forage[good]

    for j = 0L, nstar-1L do begin

;      jebv = sxpar(jacoby[j].header,'Ebv_Obs')

       scube = rd1dspec(speclist[j],/silent)
       wave = scube.wave
       flux = scube.spec
       
       djs_plot, wave, im_normalize(flux,wave), $
         xsty=3, ysty=3, ps=10, yrange=yrange, xrange=xrange, charthick=2.0, $
         charsize=2.0, xtitle='Wavelength ('+angstrom()+')', ytitle='Normalized f_{\lambda}'
       djs_oplot, jacoby[j].wave, im_normalize(jacoby[j].flux,jacoby[j].wave), ps=10, color='green'
       legend, strcompress(star[j],/remove), /right, /top, box=0, charsize=2.0, charthick=2.0
       cc = get_kbrd(1)

       icleanup, scube
       
    endfor
    
return
end    
