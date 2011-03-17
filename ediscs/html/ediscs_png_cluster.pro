pro ediscs_png_cluster, cluster, nsmooth=nsmooth, write=write
; jm04feb17uofa
; write png files of all the spectra and spectral fits
    
    if n_elements(cluster) eq 0L then begin
       print, 'Syntax - '
       return
    endif

    doit = execute('datapath = ediscs_path(/'+repstr(cluster,'-','_')+')')
    webpath = ediscs_path(/html)
    clpath = webpath+cluster+'/'
    pngpath = clpath+'png/'
    
    if keyword_set(write) then splog, 'Writing PNG files for '+cluster+'.'

    if (n_elements(nsmooth) eq 0L) then nsmooth = 5L
    
; grab FITS file headers

    forage = ediscs_headgrab('*.ms.fits.gz',datapath=datapath)
    speclist = forage.galaxy+'.ms.fits.gz'
    galaxy = forage.galaxy
    ngalaxy = n_elements(galaxy)

; default output file names
    
;   pngname = galaxy+'.png'
;   psname = galaxy+'.ps'
    specnamepng = galaxy+'_spec.png'
    specnameps = galaxy+'_spec.ps'

    spec1namepng = galaxy+'_spec1.png'
    spec1nameps = galaxy+'_spec1.ps'
    spec2namepng = galaxy+'_spec2.png'
    spec2nameps = galaxy+'_spec2.ps'
    spec3namepng = galaxy+'_spec3.png'
    spec3nameps = galaxy+'_spec3.ps'
    spec4namepng = galaxy+'_spec4.png'
    spec4nameps = galaxy+'_spec4.ps'

; loop on each galaxy

    if keyword_set(write) then stime0 = systime(1) else $
      window, 0, xs=700, ys=500

    for k = 0L, ngalaxy-1L do begin

       k = k > 0L

       if keyword_set(write) then dfpsplot, pngpath+specnameps[k], /landscape, /color

       scale = 1E17
       ytitle = 'Flux Density [10^{-17} '+flam_units()+']'
       xtitle = 'Observed Wavelength [\AA]'

       scube = rd1dspec(speclist[k],datapath=datapath,/silent)

       flux = smooth(scube.spec,nsmooth)
       ferr = scube.sigspec
       wave = scube.wave

       snrstats = im_stats(flux/ferr)
       snrstr = 'S/N = '+string(snrstats.median,format='(F5.1)')

       stats = im_stats(scale*flux)
       yrange = stats.median+stats.sigma_rej*[-1,4]
       yrange[0] = yrange[0]<0
;      yrange = scale*minmax(flux)*[0.95,1.05]
       djs_plot, wave, scale*flux, xsty=3, ysty=3, ps=10, xthick=2.0, $
         ythick=2.0, xtitle=xtitle, ytitle=ytitle, charsize=1.5, charthick=2.0, $
         yrange=yrange, thick=2.0
       legend, [galaxy[k],'z = '+string(forage[k].z,format='(F6.4)'),snrstr], /left, /top, $
         box=0, charsize=2.0, charthick=2.0
;      legend, ['Observed Spectrum'], /right, /top, box=0, charsize=2.0, $
;        charthick=2.0, color=djs_icolor('default'), line=0, thick=4.0

       if keyword_set(write) then begin
          
          dfpsclose
          splog, 'Writing and compressing '+specnameps[k]+' and '+specnamepng[k]+'.'
          spawn, ['convert -geometry 512 -rotate 270 '+pngpath+specnameps[k]+' '+pngpath+specnamepng[k]], /sh
          spawn, ['gzip -f '+pngpath+specnameps[k]], /sh

       endif

; ---------------------------------------------------------------------------       
; menu options
; ---------------------------------------------------------------------------       

       if (not keyword_set(write)) and (ngalaxy gt 1L) then begin

          print, 'Object '+strn(k,length=3)+' '+galaxy[k]+' [Options: b,g,q]'
          cc = strupcase(get_kbrd(1))

          case strlowcase(strcompress(cc,/remove)) of
             'b': k = (k-2L) ; back
             'g': begin      ; goto 
                number = ''
                read, number, prompt='Goto spectrum number (0-'+strn(ngalaxy-1L)+'): '
                number = 0 > long(number-1L) < (ngalaxy-2L)
                k = number
             end
             'q': return        ; quit
             else: 
          endcase

       endif 

    endfor 

    if keyword_set(write) then $
      splog, format='("Total time for EDISCS_PNG_CLUSTER = ",G0," minutes.")', (systime(1)-stime0)/60.0
    
return
end
