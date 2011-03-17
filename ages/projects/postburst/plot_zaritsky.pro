pro plot_zaritsky, tar=tar, doplot=doplot, nsmooth=nsmooth, postscript=postscript
; jm04sep26uofa
; plot the post-starburst galaxies that Dennis selected    

    if keyword_set(postscript) then doplot = 0L
    if n_elements(nsmooth) eq 0L then nsmooth = 3L
    
    path = ages_path(/analysis)
    tablepath = ages_path(/zaritsky)
    datapath = ages_path(/spec1d)
    
    readcol, tablepath+'postbursts.dat', passfield, cat, format='L,L', /silent

    data = mrdfits(path+'ages_data.fits.gz',1,/silent)
    match, cat, data.catalog_number, indx1, indx2

    passfield = passfield[indx1]
    cat = cat[indx1]
    data = data[indx2]
    plate = string(passfield,format='(I3)')

    ngalaxy = n_elements(data)

    if keyword_set(postscript) then begin
       dfpsplot, tablepath+'postbursts.ps', /landscape, /color
       postthick = 5.0
    endif else postthick = 2.0

    if keyword_set(doplot) then window, 0, xs=650, ys=450

    fitslist = strarr(ngalaxy)
    
    for i = 0L, ngalaxy-1L do begin

       objpath = datapath+plate[i]+'/'
       specfile = strtrim(strlowcase(data[i].galaxy),2)+'.fits.gz'
       fitslist[i] = objpath+specfile
       
       cube = rd1dspec(specfile,datapath=objpath)

       wave = cube.wave
       flux = smooth(1D17*cube.spec,nsmooth)
       
       if keyword_set(doplot) then begin
          
          info = [$
            'Galaxy '+string(cat[i],format='(I0)'),$
            'z = '+strtrim(string(data[i].zosu_z,format='(F12.2)'),2)]
       
;         djs_plot, wave, flux, ps=10, xsty=3, ysty=3, $
;           xthick=postthick, ythick=postthick, charsize=2.0, $
;           charthick=postthick, xtitle='Observed Wavelength [\AA]', $
;           ytitle='Flux [10^{-17} '+flam_units()+']', thick=2.0

          lineid_plot, wave, flux, [3727,4101]*(1+data[i].zosu_z), textoidl(['[O II]','H\delta_{A}']), $
            ps=10, xsty=3, ysty=3, xthick=postthick, ythick=postthick, lcharsize=1.2, lcharthick=postthick, $
            charsize=2.0, charthick=postthick, xtitle='Observed Wavelength [\AA]', $
            ytitle=textoidl('Flux [10^{-17} '+flam_units()+']'), thick=2.0, /extend, position=[0.12,0.15,0.99,0.87]

          legend, info, /right, /top, box=0, charsize=1.5, charthick=postthick
          cc = get_kbrd(1)

       endif
       
    endfor

    if keyword_set(postscript) then dfpsclose

; make a tar-ball

    if keyword_set(tar) then begin

       tarlist = strjoin(fitslist,' ')
       spawn, ['tar czvf '+tablepath+'postbursts.tar.gz '+tarlist], /sh

; write the mapping between catalog number and fiber number

       openw, lun, tablepath+'fibercat.dat', /get_lun
       printf, lun, '# (1) FITS file; (2) Plate number; (3) Catalog number'
       niceprintf, lun, strlowcase(data.galaxy), string(plate,format='(I7)'), string(cat,format='(I7)')
       free_lun, lun

    endif
    

return
end
