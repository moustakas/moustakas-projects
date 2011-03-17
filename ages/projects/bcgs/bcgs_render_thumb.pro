pro bcgs_render_thumb, band=band, all=all
; jm10jul23ucsd - render a thumbnail image of a BCG

    if (n_elements(band) eq 0) then band = 'I'

    bcgspath = ages_path(/projects)+'bcgs/'
    sample = rsex(bcgspath+'bcgs_sample_v3.sex')
    allfits = file_search(bcgspath+'thumbs/*'+band+'*.fits',count=ngal)

    if keyword_set(all) then begin
       psfile = bcgspath+'fortalk/bcgs_allthumbs.ps'
       im_plotconfig, 0, pos, psfile=psfile, bits_per_pixel=24
    endif
;   for ii = 0, 4 do begin
    for ii = 0, ngal-1 do begin
       if (keyword_set(all) eq 0) then begin
          psfile = bcgspath+'fortalk/'+repstr(file_basename(allfits[ii]),'.fits','.ps')
          im_plotconfig, 0, pos, psfile=psfile, bits_per_pixel=24, /keynote
       endif
       image = mrdfits(allfits[ii],0,hdr,/silent)
       extast, hdr, astr
;      imsize = size(image,/dimension)
;      xsize = imsize[0]
;      ysize = imsize[1]
;      xaxis = (findgen(imsize[0])+imsize[0]/2.0)*astr.cd[1,1]/3600.0+astr.crval[0] ; centered on the image [arcsec]
;      yaxis = (findgen(imsize[1])+imsize[1]/2.0)*astr.cd[1,1]/3600.0+astr.crval[1] ; centered on the image [arcsec]
; scale and display       
;      img = logscl(image,/neg,omin=0,omax=254)
       qq = weighted_quantile(image,quant=[0.05,0.99])
       img = asinhscl(image,/neg,alpha=2,beta=10.0,omin=10,omax=250,max=qq[1])
;      img = logscl(image,/neg,omin=10,omax=250,min=qq[0],max=qq[1])
       plotimage, img, /normal, position=pos, /noaxes; /preserve, 
;        imgxrange=minmax(xaxis), imgyrange=minmax(yaxis), $
;        xtitle=textoidl('\alpha_{J2000} (deg)'), $
;        ytitle=textoidl('\delta_{J2000} (deg)')
;      ad2xy, bootes.alpha_j2000, bootes.delta_j2000, astr, xx, yy
;      these = where((xx-astr.crpix[0] gt -astr.naxis[0]) and $
;        (xx-astr.crpix[0] lt astr.naxis[0]) and $
;        (yy-astr.crpix[1] gt -astr.naxis[1]) and $
;        (yy-astr.crpix[1] lt astr.naxis[1]),nthese)
;      djs_oplot, xx[these]-astr.crpix[0], yy[these]-astr.crpix[1], $
;        psym=7, color='red', symsize=1.5
       im_legend, ['BCG '+string(ii+1,format='(I2.2)'),'z='+$
         string(sample[ii].z,format='(F6.4)')], /left, /top, box=0, $
         charsize=2.4, textcolor='white', charthick=4.0, $
         position=[pos[0]+0.02,pos[2]-0.03], /normal
       if (keyword_set(all) eq 0) then $
         im_plotconfig, psfile=psfile, /psclose, /pdf, /keynote
    endfor
    if keyword_set(all) then im_plotconfig, psfile=psfile, /gzip, /psclose

return
end
    
