pro driftscan_animation, atlas, make_mpeg=make_mpeg
; jm05nov27uofa - 

    datapath = '/home/ioannis/jobs/talk/'
    mpegpath = datapath+'mpegframes/'
    
    if (n_elements(atlas) eq 0L) then atlas = atlas_read_info()

; identify M51 and display it for fun    
    
    indx = speclinefit_locate(atlas,'ngc5194')
    
    atlas_display_image, atlas[indx], /preserve_aspect, labeltype=0, $
      /norc3box, astr=astr, imageinfo=imageinfo, /arcminlabel, $
      barlabelcolor='black', /nolabelbar, imposition=imposition

    postthick = 2.0
    
; read the image    
    
    imagepath = atlas_path(/dss)
    fitsname = imageinfo.fitsname
    image = readfits(imagepath+fitsname,h,/silent)

    image = alog10(image>1)
    zz = zscale_range(image,0.1)
    top = 230L
    img = bytscl(image,min=zz[0],max=zz[1],top=top)
    img = bytscl(top-img,min=-40,top=top)
    
; slit position vector    

    npass = 2L
    
    slitstart = -180.0 ; [arcsec]
    slitend = 180.0    ; [arcsec]
    slitrate = 4.0     ; [arcsec/s]
    slitpos = (findgen((slitend-slitstart)/slitrate+1)*slitrate+slitstart)/60.0
    slitpos = [slitpos,reverse(slitpos)]
    nslitpos = n_elements(slitpos)
    slitpos = reform(rebin(slitpos,nslitpos,npass),nslitpos*npass)
    nslitpos = n_elements(slitpos)

    leftpos = -190.0/60.0
    rightpos = 190.0/60.0

    mpeg_frame_name = 'frame_'+string(lindgen(nslitpos),format='(I4.4)')

    if keyword_set(make_mpeg) then begin
;      mpegid = mpeg_open([479,479])
       set_plot, 'Z'
    endif

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=8.5, height=8.0, $
      xmargin=0.0, ymargin=0.0, xpage=8.5, ypage=8.0, position=position, $
      /normal 

;   if keyword_set(make_mpeg) then set_plot, 'Z'
    for i = 0L, nslitpos-1L do begin

       print, format='("Frame ",I3,"/",I3,".",A5,$)', i, nslitpos, string(13B)
       
       plotimage, img, /normal, margin=0, position=position, $
         imgxrange=minmax(imageinfo.xaxis), imgyrange=minmax(imageinfo.yaxis), $
         charthick=postthick, xthick=postthick, ythick=postthick, $
         /noaxes, noerase=noerase, /preserve_aspect
       oplot, [leftpos,rightpos], slitpos[i]*[1,1], line=0, color=fsc_color('red'), thick=3

; capture the image and write out             

       if keyword_set(make_mpeg) then begin

          x0 = fix((convert_coord(position[0:1],/normal,/to_device))[0])
          nx = fix((convert_coord(position[2:3],/normal,/to_device))[0])-x0
          y0 = fix((convert_coord(position[0:1],/normal,/to_device))[1])
          ny = fix((convert_coord(position[2:3],/normal,/to_device))[1])-y0
          
          thisimage = tvrd(x0+1L,y0+1L-1L,nx-1L,ny-1L,/order)
          if (i eq 0L) then bigimage = thisimage else bigimage = [ [ [bigimage] ], [ [thisimage] ] ]
;         help, bigimage

       endif
          
;      if keyword_set(make_mpeg) then mpeg_put, mpegid, image=thisimage, frame=i
       
;      tvlct, r, g, b, /get
;      write_png, mpegpath+mpeg_frame_name[i], img, r, g, b
          
    endfor

    if keyword_set(make_mpeg) then begin
;      mpeg_save, mpegid, filename=datapath+'driftscan_animation.mpg'
;      mpeg_close, mpegid
       ppmtompeg, bigimage, datapath+'driftscan_animation.mpg', tmpdir=tmpdir
       set_plot, 'X'
    endif

stop
    
return
end
    
    
