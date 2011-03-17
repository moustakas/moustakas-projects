pro view_all_images
; jm02feb19uofa
    
    path = '/home/ioannis/kennicutt/atlas2d/images/'
    pushd, path

    flist = findfile('*.fits.gz',count=fcount)
    ra = strarr(fcount)
    dec = strarr(fcount)
    
    print, 'Retrieving header information for the ATLAS2d galaxies.'
    for i = 0L, fcount-1L do begin
       h = headfits(flist[i])
       ra[i] = im_dec2hms(sxpar(h,'CRVAL1')/15.0)
       dec[i] = im_dec2hms(sxpar(h,'CRVAL2'))
    endfor
       
    srtra = sort(ra)
    flist = flist[srtra]
    ra = ra[srtra]
    dec = dec[srtra]
    
    window, 0, xs=650, ys=650
    print, 'Displaying...'
    for i = 0L, fcount-1L do begin

       if i gt 0L then begin
       
          print, 'Press any key to continue, or (g,b,q).'
          cc = strupcase(get_kbrd(1))
          case strlowcase(strcompress(cc,/remove)) of
             'b': i = i-2L
             'g': begin
                number = ''
                read, number, prompt='Goto image number (0-'+strn(fcount-1L)+'): '
                number = 0 > long(number) < (fcount-2L)
                i = number                
             end
             'q': return
             else: 
          endcase

       endif

       print, strn(i,length=3)+'  '+flist[i]
       im = readfits(flist[i],/silent)

       str = sstretch(im,/silent,npix=1000)
       immax = (str.median + (20*str.sigma)) < max(im)
       immin = (str.median - (2*str.sigma)) < min(im)
       
       display, bytscl(im,min=immin,max=immax,top=!d.table_size-1), /aspect
       galaxy = strupcase(strmid(flist[i],0,strpos(flist[i],'.fits.gz')))
       legend, [galaxy,ra[i],dec[i]], /left, /top, box=0, charsize=2.0, charthick=2.0

;      atv, flist[i], /block  

    endfor

    popd
    
return
end
