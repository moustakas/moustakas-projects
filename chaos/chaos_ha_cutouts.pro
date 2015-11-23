function stretchit, im, scale=scale
    mn = -5.0
    mx = 10.0
    farcsinh = 1.0
    nim = asinh(im*farcsinh/scale)/sqrt(farcsinh)
;   ss = im_stats(nim,/ver)
    img = bytscl(nim,min=mn,max=mx,top=240)
return, img
end

pro chaos_ha_cutouts
; jm15nov21siena - get Ha cutouts of all the HII regions

    gal = 'm101'
    ngal = n_elements(gal)

    path = getenv('IM_PROJECTS_DIR')+'/chaos/'

    ds9 = djs_readlines(path+'m101_obsreg2.reg',nhead=4)
    nhii = n_elements(ds9)
    splog, nhii

    info = replicate({region: '', ra: 0D, dec: 0D, oh: 0.0, $
      ohmodel: 0.0, doh: 0.0},nhii)
    for ii = 0, nhii-1 do begin
       x1 = strpos(ds9[ii],'box(')+4
       x2 = strpos(ds9[ii],'3"')-1-x1
       radec = strsplit(strmid(ds9[ii],x1,x2),',',/extract)
       info[ii].ra = radec[0]*1D
       info[ii].dec = radec[1]*1D

       x1 = strpos(ds9[ii],'xt={')+4
       x2 = strpos(ds9[ii],'}')-x1
       hinfo = strsplit(strmid(ds9[ii],x1,x2),',',/extract)
       info[ii].region = strtrim(hinfo[0],2)
       info[ii].oh = hinfo[1]*1.0
       info[ii].ohmodel = hinfo[2]*1.0
       info[ii].doh = hinfo[3]*1.0
    endfor
    struct_print, info
    srt = reverse(sort(abs(info.doh)))
    info = info[srt]

; get cutouts and pack into a montage
    im = mrdfits('m101_ha.fits',0,hdr)
    extast, hdr, astr
    sz = size(im,/dim)
    pixscale = sqrt(abs(determ(astr.cd)))*3600D ; [arcsec/pixel]
    splog, pixscale
    
    farcsinh = 1.0
;   scale = median(im);+5*djsig(im)
;   ss = im_stats(im,/ver,sigrej=10)
    diam = 40

;   for ii = 0, 5 do begin
    for ii = 0, nhii-1 do begin
       ad2xy, info[ii].ra, info[ii].dec, astr, xx, yy
       x0 = long(xx-diam/2.0)>0L
       x1 = long(xx+diam/2.0)<(sz[0]-1)
       y0 = long(yy-diam/2.0)>0
       y1 = long(yy+diam/2.0)<(sz[1]-1)
       splog, diam, x0, x1, y0, y1
;      if ii eq 10 then stop
       
       hextract, im, hdr, newim, newhdr, x0, x1, y0, y1, /silent

       extast, newhdr, newastr
       ad2xy, info[ii].ra, info[ii].dec, newastr, xcen, ycen
;      splog, xcen, ycen

       obj = strtrim(info[ii].region,2)
       outfile = path+'jpeg/'+string(ii,format='(I2.2)')+'_'+obj+'.png'

;      mm = weighted_quantile(newim,quant=[0.2,0.98])
;      outim = im_imgscl(newim,/neg,/sqrroot)
;      outim = stretchit(newim,scale=scale)
;      ss = im_stats(newim,/ver)
;      newim = newim/scale
;      newim = asinh(newim*farcsinh/scale)/sqrt(farcsinh)
;      outim = cgScaleVector(newim, 0, 255, MINVALUE=-2.0, MAXVALUE=20.0)
;      outim = LogScl(newim, MIN=-1, MAX=20)

       ss = im_stats(newim,sigrej=6)

       delvarx, pos
       set_plot, 'Z'

;      pos = [0.1, 0.1, 0.9, 0.8]
       LoadCT, 0, /silent, bottom=5
       cgImage, newim, /KEEP_ASPECT, POSITION=pos, Erase=1, axes=0, $
         /norm, stretch=5, minvalue=-ss.sigma, maxvalue=ss.maxrej, $
         /save, bottom=5
       im_oplot_box, 10.0/pixscale, 1.0/pixscale, 0.0, $
         xoffset=xcen, yoffset=ycen, thick=4, color=cgcolor('red',2)
;      djs_oplot, [diam/2.0-5/pixscale,diam/2.0+5/pixscale], $
;        diam/2-0.5/pixscale*[1,1], thick=5
;      djs_oplot, [diam/2.0-5/pixscale,diam/2.0+5/pixscale], $
;        diam/2+0.5/pixscale*[1,1], thick=5

;      plotimage, outim, /noaxes, /preserve_aspect, position=pos, /norm
       im_legend, [obj,strtrim(string(info[ii].oh,format='(F12.3)'),2),$
         strtrim(string(info[ii].doh,format='(F12.3)'),2)], $
         /left, /top, box=0, charsize=2.0, $
         charthick=3.0, textcolor=cgcolor('red',2), $
         position=[pos[0]+0.0,pos[3]-0.05], /norm
;      tvcircle, 10.0/pixscale, xcen, ycen, thick=4, color=im_color('red',2)

       for jj = 0, nhii-1 do begin
          ad2xy, info[jj].ra, info[jj].dec, newastr, xh, yh
          if jj ne ii then xyouts, xh, yh, strtrim(info[jj].region,2), $
            color=cgcolor('green',1), /data
       endfor

       LoadCT, 0, /silent, bottom=5

       x0 = fix((convert_coord(pos[0:1],/normal,/to_device))[0])
       nx = fix((convert_coord(pos[2:3],/normal,/to_device))[0])-X0
       y0 = fix((convert_coord(pos[0:1],/normal,/to_device))[1])
       ny = fix((convert_coord(pos[2:3],/normal,/to_device))[1])-y0
       img = tvrd(x0,y0,nx,ny)
       splog, file_basename(outfile), info[ii].ra, info[ii].dec, $
         size(img,/dim), size(newim,/dim)

       tvlct, rr, gg, bb, /get
       write_png, outfile, img, rr, gg, bb
       set_plot, 'X'

;      outim = asinhscl(newim,/negative,min=mm[0],max=mm[1],beta=0)
;      outim = asinhscl(newim,/negative,min=mm[0],max=mm[1],beta=1.0)
;      write_jpeg, outfile, outim
    endfor

    ncol = 10
    nrow = ceil(nhii/(1.0*ncol))

    outfile = path+gal+'_hii.png'
    cmd = 'montage -bordercolor black -borderwidth 1 '+ $
      '-tile '+strtrim(ncol,2)+'x'+strtrim(nrow,2)+' -geometry +0+0 '+$
      '-quality 100 '+path+'jpeg/*.png '+outfile
    spawn, cmd, /sh

stop

return
end
