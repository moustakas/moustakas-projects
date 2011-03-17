pro ldss3_psf, fitpsf=fitpsf, makeplots=makeplots, ps=ps
; jm07may21nyu - measure the PSF in every LDSS3 image

    datapath = ldss3_path(/feb06)+'sg1120/'
    analysis_path = sg1120_path(/analysis)

    imagelist = file_search(datapath+'ra.????_sg1120_?_[g,r].fits')
    weightlist = repstr(imagelist,'.fits','.weight.fits')
    skylist = repstr(imagelist,'.fits','.sky.fits')
    templist1 = repstr(imagelist,'.fits','.c1.fits')
    templist2 = repstr(imagelist,'.fits','.c2.fits')

    if keyword_set(fitpsf) then begin
       
       t0 = systime(1)
       for ii = 0L, n_elements(imagelist)-1L do begin

          print, 'Measuring the PSF for '+imagelist[ii]
          
; read each image and bad pixel map and construct an inverse variance
; map

          im1 = mrdfits(imagelist[ii],1,/silent)
          im2 = mrdfits(imagelist[ii],2,/silent)
          sky1 = mrdfits(skylist[ii],1,/silent)
          sky2 = mrdfits(skylist[ii],2,/silent)
;         w1 = mrdfits(weightlist[ii],1,/silent)
;         w2 = mrdfits(weightlist[ii],2,/silent)
;         sig1 = dsigma(im1)
;         sig2 = dsigma(im2)

          invmap1 = (im1*0.0+1.0/sky1^2)*(im1 gt 0.0)
          invmap2 = (im2*0.0+1.0/sky2^2)*(im2 gt 0.0)
;         invmap1 = (im1*0.0+1.0/sig1^2)*w1
;         invmap2 = (im2*0.0+1.0/sig2^2)*w2
;         invmap1 = 1.0/(im1+(im1 eq 0.0))*w1
;         invmap2 = 1.0/(im2+(im2 eq 0.0))*w2

; now write out each chip separately       
          
          splog, 'Writing temporary file '+templist1[ii]
          mkhdr, hdr1, im1, /extend
          mwrfits, im1, templist1[ii], hdr1, /create
          mwrfits, invmap1, templist1[ii]

          splog, 'Writing temporary file '+templist2[ii]
          mkhdr, hdr2, im2, /extend
          mwrfits, im2, templist2[ii], hdr2, /create
          mwrfits, invmap2, templist2[ii]

          splog, 'Fitting for the PSF in image '+templist1[ii]
          dfitpsf, templist1[ii]
          splog, 'Fitting for the PSF in image '+templist2[ii]
          dfitpsf, templist2[ii]

; remove the temporary FITS files

          rmfile, templist1[ii]
          rmfile, templist2[ii]
          
       endfor
       splog, 'Total time to run = ', (systime(1)-t0)/60.0, ' minutes.'

    endif

    if keyword_set(makeplots) then begin

       fwhm2sig = 2D0*sqrt(2.0*alog(2.0))
       pixscale = 0.188

       psflist_c1 = file_search(datapath+'*.c1-bpsf.fits')
       psflist_c2 = file_search(datapath+'*.c2-bpsf.fits')
       imagelist = repstr(psflist_c1,'.c1-bpsf','')
       nimage = n_elements(imagelist)

       gband = where(strmatch(imagelist,'*_g.fits'))
       rband = where(strmatch(imagelist,'*_r.fits'))
       
       jd_c1 = dblarr(nimage)   & jd_c2 = dblarr(nimage)
       time_c1 = strarr(nimage) & time_c2 = strarr(nimage)
       fwhm_c1 = fltarr(nimage) & fwhm_c2 = fltarr(nimage)

       for it = 0L, nimage-1L do begin
; chip 1
          date = strsplit(sxpar(headfits(imagelist[it],ext=1),'DATE-OBS'),'-',/extract)
          time = strsplit(sxpar(headfits(imagelist[it],ext=1),'UT-TIME'),':',/extract)
          jd_c1[it] = julday(date[1],date[2],date[0],time[0],time[1],time[2])
          time_c1[it] = sxpar(headfits(imagelist[it],ext=1),'UT-TIME')
          fwhm_c1[it] = fwhm2sig*pixscale*sxpar(headfits(psflist_c1[it],ext=1),'PSFSIGMA')
; chip 2
          date = strsplit(sxpar(headfits(imagelist[it],ext=2),'DATE-OBS'),'-',/extract)
          time = strsplit(sxpar(headfits(imagelist[it],ext=2),'UT-TIME'),':',/extract)
          jd_c2[it] = julday(date[1],date[2],date[0],time[0],time[1],time[2])
          time_c2[it] = sxpar(headfits(imagelist[it],ext=2),'UT-TIME')
          fwhm_c2[it] = fwhm2sig*pixscale*sxpar(headfits(psflist_c2[it],ext=1),'PSFSIGMA')
       endfor

       jdoffset = long(min(jd_c1))
       jd_c1 = jd_c1 - jdoffset
       jd_c2 = jd_c2 - jdoffset

       xrange = minmax(jd_c1)
       yrange = [min(fwhm_c1)<min(fwhm_c2),max(fwhm_c1)>max(fwhm_c2)]

       xtitle = 'Julian Date (Relative to '+string(jdoffset,format='(I0)')+')'
       ytitle = 'FWHM Seeing (arcsec)'

       pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
         xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
         position=pos, /normal

       if keyword_set(ps) then begin
          dfpsplot, analysis_path+'ldss3_psf.ps', /square, /color
          postthick = 4.0
       endif else begin
          im_window, 0, xratio=0.4, /square
          postthick = 2.0
       endelse
       
;      djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, charsize=2.0, $
;        charthick=postthick, xthick=postthick, ythick=postthick, xtitle=xtitle, $
;        ytitle=ytitle, xsty=3, ysty=3, position=pos
;      im_symbols, 106, psize=1.5, fill=1, thick=postthick, color=djs_icolor('grey')
;      djs_oplot, jd_c1, fwhm_c1, ps=-8, color=djs_icolor('grey'), line=0, thick=postthick
;      im_symbols, 108, psize=1.5, fill=1, thick=postthick, color=djs_icolor('default')
;      djs_oplot, jd_c2, fwhm_c2, ps=-8, color=djs_icolor('default'), line=2, thick=postthick
;
;      im_legend, ['Chip 1','Chip 2'], /left, /top, psym=[106,108], /fill, $
;        color=djs_icolor(['grey','default']), charsize=2.0, charthick=postthick, $
;        box=0, symsize=1.5, line=[0,2]
       
;      hmlabel = label_date(date_format='%H:%I:%S')
;      plot, time_c1, 2.35*0.188*fwhm_c1, ps=-4, xsty=1, ysty=1, charsize=1.8, $
;        charthick=2.0, xthick=2.0, ythick=2.0, xtickunits='Time', xtickformat='LABEL_DATE'

;      if (not keyword_set(ps)) then cc = get_kbrd(1)

       djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, charsize=2.0, $
         charthick=postthick, xthick=postthick, ythick=postthick, xtitle=xtitle, $
         ytitle=ytitle, xsty=3, ysty=3, position=pos
       im_symbols, 106, psize=1.5, fill=1, thick=postthick, color=djs_icolor('blue')
       djs_oplot, jd_c1[gband], fwhm_c1[gband], ps=-8, color=djs_icolor('blue'), line=0, thick=postthick
       im_symbols, 108, psize=1.5, fill=1, thick=postthick, color=djs_icolor('red')
       djs_oplot, jd_c1[rband], fwhm_c1[rband], ps=-8, color=djs_icolor('red'), line=2, thick=postthick

       im_legend, ['Chip 1, g-band','Chip 1, r-band'], /left, /top, psym=[106,108], /fill, $
         color=djs_icolor(['blue','red']), charsize=2.0, charthick=postthick, $
         box=0, symsize=1.5, line=[0,2]
       
       if (not keyword_set(ps)) then cc = get_kbrd(1)

       djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, charsize=2.0, $
         charthick=postthick, xthick=postthick, ythick=postthick, xtitle=xtitle, $
         ytitle=ytitle, xsty=3, ysty=3, position=pos
       im_symbols, 106, psize=1.5, fill=1, thick=postthick, color=djs_icolor('blue')
       djs_oplot, jd_c2[gband], fwhm_c2[gband], ps=-8, color=djs_icolor('blue'), line=0, thick=postthick
       im_symbols, 108, psize=1.5, fill=1, thick=postthick, color=djs_icolor('red')
       djs_oplot, jd_c2[rband], fwhm_c2[rband], ps=-8, color=djs_icolor('red'), line=2, thick=postthick

       im_legend, ['Chip 2, g-band','Chip 2, r-band'], /left, /top, psym=[106,108], /fill, $
         color=djs_icolor(['blue','red']), charsize=2.0, charthick=postthick, $
         box=0, symsize=1.5, line=[0,2]
       
       if keyword_set(ps) then dfpsclose

       splog, 'Chip 1 stats'
       c1stats = im_stats(100.0*(fwhm_c1[gband]/interpol(fwhm_c1[rband],jd_c1[rband],jd_c1[gband])-1.0),/verbose)
       splog, 'Chip 2 stats'
       c2stats = im_stats(100.0*(fwhm_c2[gband]/interpol(fwhm_c2[rband],jd_c2[rband],jd_c2[gband])-1.0),/verbose)

    endif
    
stop
   
return
end
    
;;    t0 = systime(1)
;;    for ii = 0L, n_elements(imagelist)-1L do begin
;;
;;       print, 'Measuring the PSF for '+imagelist[ii]
;;       
;;; read the images and bad pixel maps and construct an inverse variance
;;; map 
;;
;;       im1 = mrdfits(imagelist[ii],1,/silent)
;;       im2 = mrdfits(imagelist[ii],2,/silent)
;;       w1 = mrdfits(weightlist[ii],1,/silent)
;;       w2 = mrdfits(weightlist[ii],2,/silent)
;;
;;       sig1 = dsigma(im1)
;;       sig2 = dsigma(im2)
;;
;;;      invmap1 = (im1*0.0+1.0/sig1^2)*w1
;;;      invmap2 = (im2*0.0+1.0/sig2^2)*w2
;;       invmap1 = 1.0/(im1+(im1 eq 0.0))*w1
;;       invmap2 = 1.0/(im2+(im2 eq 0.0))*w2
;;
;;; do some simple sky-subtraction
;;
;;       mmm, im1, sky1 & im1 = im1 - sky1
;;       mmm, im2, sky2 & im2 = im2 - sky2
;;
;;; build one large image, crop the edges, and write out 
;;
;;       imsize = size(im1,/dim)
;;       
;;       bigimage = fltarr(imsize[0]*2,imsize[1])
;;       bigimage[0L:imsize[0]-1L,*] = im2
;;       bigimage[imsize[0]:imsize[0]*2L-1L,*] = im1
;;       
;;       biginvmap = fltarr(imsize[0]*2,imsize[1])
;;       biginvmap[0L:imsize[0]-1L,*] = invmap2
;;       biginvmap[imsize[0]:imsize[0]*2L-1L,*] = invmap1
;;
;;       bigimsize = size(bigimage,/dim)
;;
;;       indx = where(bigimage gt 0.0)
;;       x1 = min(indx mod bigimsize[0])
;;       y1 = min(indx/bigimsize[0])
;;       x2 = max(indx mod bigimsize[0])
;;       y2 = max(indx/bigimsize[0])
;;       bigimage = bigimage[x1:x2,y1:y2]
;;       biginvmap = biginvmap[x1:x2,y1:y2]
;;       
;;       mkhdr, hdr, bigimage, /extend
;;
;;       splog, 'Writing temporary file '+templist[ii]
;;       mwrfits, bigimage, templist[ii], hdr, /create
;;       mwrfits, biginvmap, templist[ii]
;;
;;stop
;;       
;;       dfitpsf, templist[ii]
;;
;;; rename the PSF files
;;       
;;       rmfile, templist[ii]
;;       
;;stop       
;;
;;    endfor
;;    splog, 'Total time to run = ', (systime(1)-t0)/60.0, ' minutes.'
