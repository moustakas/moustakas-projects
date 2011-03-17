pro vimos_thumbnails, cat, band=band, stamp=stamp, $
  nperpage=nperpage, psfile=psfile
; jm08jul25nyu - generate thumbnail postage stamps of objects
; identified in the VIMOS mosaics (based on older code)

; cat [NOBJ] - SE catalog of the sources to display
; stamp - width of the thumbnail (default 20 arcsec) 
; nperpage - number of objects per page (default 16)
; psfile - output postscript file name (if set then generate
;          postscript output)

    common vimos_mosaics, bimg, vimg, rimg

    nobj = n_elements(cat)
    if (nobj eq 0L) then message, 'CAT not defined'

    if (n_elements(band) eq 0L) then band = 'R'
    band = strupcase(band)

    pixscale = 0.205 ; [arcsec/pixel]
    if (n_elements(stamp) eq 0L) then stamp = 20.0/pixscale ; [pixel]
    if (n_elements(nperpage) eq 0L) then nperpage = 16.0

; read the requisite mosaic    
    
    mosaicpath = vimos_path(/mosaics)
    mosaicname = mosaicpath+'sg1120_'+band+'.fits'
    case band of
       'B': begin
          if (n_elements(bimg) eq 0L) then $
            bimg = mrdfits(mosaicname,0,hdr)
          img = bimg
       end
       'V': begin
          if (n_elements(vimg) eq 0L) then $
            vimg = mrdfits(mosaicname,0,hdr)
          img = vimg
       end
       'R': begin
          if (n_elements(rimg) eq 0L) then $
            rimg = mrdfits(mosaicname,0,hdr)
          img = rimg
       end
       else: message, 'VIMOS bandpass not recognized'
    endcase

; some plotting defaults

    npage = ceil(nobj/nperpage)
    ncols = sqrt(nperpage) & nrows = ncols
    xspace = 0.0 & yspace = 0.0
    xmargin = [0.2,0.2]-xspace*(ncols-1.0)*[1,1]/2.0
    ymargin = [0.2,0.2]-yspace*(nrows-1.0)*[1,1]/2.0

    psize = (8.5-total(xmargin)-total(xspace))/ncols
    width = replicate(psize,ncols)
    height = replicate(psize,nrows)

    xpage = total(width)+total(xmargin)+total(xspace)
    ypage = total(height)+total(ymargin)+total(yspace)

    if (n_elements(psfile) ne 0L) then begin
       dfpsplot, psfile, xsize=xpage, ysize=ypage, $
         /color, bits=24
       postthick1 = 4.0
       postthick2 = 3.0
    endif else begin
       postthick1 = 2.0
       postthick2 = 2.0
    endelse

;   for jpage = 0L, 0L do begin
    for jpage = 0L, npage-1L do begin
  
       print, format='("Page ",I0,"/",I0,".",A10,$)', $
         jpage+1, npage, string(13b)
          
       xspace1 = xspace & yspace1 = yspace
       xmargin1 = xmargin & ymargin1 = ymargin
       width1 = width & height1 = height
       xpage1 = xpage & ypage1 = ypage

       arm_plotconfig, /landscape, nx=ncols, ny=nrows, xmargin=xmargin1, $
         ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
         height=height1, coord=pos, xpage=ypage1, ypage=xpage1, bw=0
       cleanplot, /silent
 
       for kperpage = 0L, nperpage-1L do begin

          iobj = fix(jpage*nperpage+kperpage)
          if (iobj lt nobj) then begin

             xcen = cat[iobj].xwin_image
             ycen = cat[iobj].ywin_image

             subimg = img[xcen-stamp/2.0:xcen+stamp/2.0,$
               ycen-stamp/2.0:ycen+stamp/2.0]
             imsize = size(subimg,/dimension)
             xsize = imsize[0] & xcen = xsize/2.0
             ysize = imsize[1] & ycen = ysize/2.0
             xaxis = (findgen(xsize)-xcen);*pixscale ; [arcsec]
             yaxis = (findgen(ysize)-ycen);*pixscale ; [arcsec]

             plotimage, asinhscl(subimg,negative=1,alpha=5.0,$
               beta=16.0,omin=0,omax=250), /normal, $
;            plotimage, logscl(subimg,negative=1,exp=0.5,$
;              mean=1.0,omin=0,omax=255), /normal, $
               position=pos[*,kperpage], margin=0, imgxrange=minmax(xaxis), $
               imgyrange=minmax(yaxis), noerase=(kperpage gt 0L), $
               xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
               xsty=5, ysty=5;, /preserve
;            tvellipse, rcat[indx[iobj]].a_image, rcat[indx[iobj]].b_image, $
;              0.0, 0.0, rcat[indx[iobj]].theta_image, line=0, thick=1.0, /data

             legend, strtrim(string(cat[iobj].number,format='(I0)'),2), $
               /right, /bottom, box=0, charsize=1.2, charthick=2.0, margin=0
             
          endif
       endfor

       if (n_elements(psfile) eq 0L) then cc = get_kbrd(1)

    endfor

    if (n_elements(psfile) ne 0L) then begin
       dfpsclose
       spawn, 'gzip -f '+psfile
    endif

return
end
