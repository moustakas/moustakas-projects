function stretchit, im, scale=scale
    mn = -1.0
    mx = 10.0
    farcsinh = 1.0
    nim = asinh(im*farcsinh/scale)/sqrt(farcsinh)
    img = bytscl(nim,min=mn,max=mx)
    img = max(img)-img
return, img
end

pro decals_edr3, stars=stars, html_only=html_only
; jm15feb27siena - test the preliminary EDR3 catalogs

    edrpath = '/global/work/decam/cats/edr3/'
    outpath = edrpath+'cutouts/'
    
    cut = 15.0                  ; [arcsec]
    pixscale = 0.262D           ; [arcsec/pixel]
    ncutpix = long(cut/pixscale)
    
    nmosaic = 3600              ; size of the coadds [pixels]
    npad = 60                   ; ignore objects near the edges [pixels]
    beta = 0.01
    scales = [0.0066,0.01,0.025]
    
    alltype = ['S','D','E','C']
    allmytype = ['Star','deVac','Exp','Comp']
    ntype = n_elements(alltype)

    band = ['g','r','z']
    bandindx = [1,2,4]
    nband = n_elements(band)
    
    allbrick = file_basename(file_search(edrpath+'tractor/*',/test_dir,count=nbrick))
    for ii = 0, 0 do begin
;   for ii = 0L, nbrick-1 do begin
       catfile = file_search(edrpath+'tractor/'+allbrick[ii]+'/tractor-*.fits',count=ncat)
       for ic = 0L, 0 do begin
;      for ic = 0L, ncat-1 do begin
          cat = mrdfits(catfile[ic],1)
          brick = strtrim(cat[0].brickname,2)
          coaddpath = edrpath+'coadd/'+allbrick[ii]+'/'+brick+'/'

          html_path = outpath+brick+'/'
          file_mkdir, html_path ; base output path

; read the coadds once
          if keyword_set(html_only) then begin
             hdr = headfits(coaddpath+'decals-'+brick+'-image-'+band[0]+'.fits')
             extast, hdr, astr
             ad2xy, cat.ra, cat.dec, astr, xcat, ycat
          endif else begin
             image = fltarr(nmosaic,nmosaic,nband)
             model = image
             mask = image
             for ib = 0, nband-1 do begin
                imfile = coaddpath+'decals-'+brick+'-image-'+band[ib]+'.fits'
                image[*,*,ib] = mrdfits(imfile,0,hdr)
                model[*,*,ib] = mrdfits(repstr(imfile,'-image','-model'),0)
                mask[*,*,ib] = mrdfits(repstr(imfile,'-image','-invvar'),0)
                extast, hdr, astr
                ad2xy, cat.ra, cat.dec, astr, xcat, ycat
             endfor
             resid = image - model

; parse the photometry             
             decals_to_maggies, cat, maggies, ivarmaggies
             mag = maggies2mag(maggies,ivarmaggies=ivarmaggies,$
               magerr=magerr,magnsigma=magnsigma)
          endelse
          
; one set of webpages per type
          for it = 2, 2 do begin
;         for it = 0, ntype-1 do begin
             html_base = brick+'-'+allmytype[it]
             png_path = html_path+html_base+'/'
             file_mkdir, png_path
             
             these = where(cat.type eq alltype[it] and xcat gt npad and $
               xcat lt nmosaic-npad and ycat gt npad and $
               ycat lt nmosaic-npad and cat.sdss_objid ne 0,nobj)
;            srt = sort(cat[these].decam_flux[1]) ; sort by r-band flux
             srt = reverse(sort(cat[these].decam_flux[2])) ; sort by r-band flux
             these = these[srt]
             these = these[0:14] & nobj = 15

; (re)build just the HTML pages
             if keyword_set(html_only) then begin
                pngfiles = brick+'-'+string(cat[these].objid,format='(I6.6)')+'.png'
                decals_edr2_html, html_base, html_path=html_path, npngcols=5, $
                  title='Brick '+brick+'-'+allmytype[it], /nomenu, $
                  pnglist=pngfiles, ra=cat[these].ra, dec=cat[these].dec
             endif else begin
                xall = xcat[these]
                yall = ycat[these]
;               xall = cat[these].bx
;               yall = cat[these].by

; get image cutouts in all three bands
                for ig = 0L, nobj-1 do begin
                   obj = brick+'-'+string(cat[these[ig]].objid,format='(I6.6)')
                   splog, ig, obj
                   xx = xall[ig]
                   yy = yall[ig]
                   
                   loadct, 0, /silent
                   psfile = png_path+obj+'.ps'
;                  cgps_open, png_path+obj+'.png'
;                  cgdisplay, aspect=1.0 ;/3.0
;                  pos = cglayout([3,nband],aspect=1.0,xgap=0.0,ygap=1,$
;                    oxmargin=[1,1],oymargin=[1,1])
                   im_plotconfig, 0, pos, psfile=psfile
;                  im_plotconfig, 15, pos, psfile=psfile, xmargin=[0.0,0.0], $
;                    ymargin=[0.0,0.0], xspace=0.0, yspace=[0.5,0.5,0.5]
                   pos = im_getposition(nx=4,ny=nband,width=5.0,height=5.0,$
                     xmargin=[0.01,0.01],ymargin=[0.01,0.01],xspace=0.01,yspace=0.02,$
                     xpage=20.05,ypage=15.06)
                   
                   for ib = 0, nband-1 do begin
                      gim = image[long(xx-ncutpix/2.0):long(xx+ncutpix/2.0),$
                        long(yy-ncutpix/2.0):long(yy+ncutpix/2.0),ib]
                      gmask = mask[long(xx-ncutpix/2.0):long(xx+ncutpix/2.0),$
                        long(yy-ncutpix/2.0):long(yy+ncutpix/2.0),ib]
                      gmodel = model[long(xx-ncutpix/2.0):long(xx+ncutpix/2.0),$
                        long(yy-ncutpix/2.0):long(yy+ncutpix/2.0),ib]
                      gresid = resid[long(xx-ncutpix/2.0):long(xx+ncutpix/2.0),$
                        long(yy-ncutpix/2.0):long(yy+ncutpix/2.0),ib]
                      
                      sgim = stretchit(gim,scale=scales[ib])
                      sgmask = stretchit(1.0/sqrt(gmask+(gmask eq 0)),scale=scales[ib]) * (gmask ne 0)
                      sgmodel = stretchit(gmodel,scale=scales[ib])
                      sgresid = stretchit(gresid,scale=scales[ib])
                      
                      plotimage, sgim, /preserve, /noaxes, position=pos[*,4*ib+0], $
                        noerase=ib gt 0, /save
;                     plots, ncutpix/2.0, ncutpix/2.0, psym=cgsymcat(9,thick=5), $
;                       symsize=1.2, color=cgcolor('black')
                      djs_oplot, cat.bx-xx+ncutpix/2.0, cat.by-yy+ncutpix/2.0, psym=cgsymcat(16), $
                        symsize=1.0, color=cgcolor('dodger blue')
                      im_legend, obj, /left, /top, box=0, charsize=1.2, $
                        charthick=4.0, position=[pos[0,4*ib+0]-0.01,pos[3,4*ib+0]-0.06], $
                        /norm, textcolor=cgcolor('white')
                      loadct, 0, /silent
                      if mag[ib,these[ig]] gt 0 then begin
                         thismag = band[ib]+'='+strtrim(string(mag[ib,these[ig]],format='(F12.3)'),2)+$
                           '\pm'+strtrim(string(magerr[ib,these[ig]],format='(F12.3)'),2)
                      endif else begin
                         thismag = band[ib]+'>'+string(magnsigma[ib,these[ig]],format='(F5.2)')
                      endelse
                      im_legend, thismag, /left, /bottom, box=0, charsize=1.2, $
                        charthick=4.0, position=[pos[0,4*ib+0]-0.01,pos[1,4*ib+0]+0.05], /norm, $
                        textcolor=cgcolor('white')
                      loadct, 0, /silent
                      
                      plotimage, sgmask, /preserve, /noaxes, position=pos[*,4*ib+1], /noerase
                      frac = strtrim(string(cat[these[ig]].decam_fracmasked,format='(F12.2)'),2)
                      im_legend, 'fracmasked='+frac, /left, /bottom, box=0, charsize=1.2, $
                        charthick=4.0, position=[pos[0,4*ib+1]-0.01,pos[1,4*ib+1]+0.05], /norm, $
                        textcolor=cgcolor('white')
                      loadct, 0, /silent

                      plotimage, sgmodel, /preserve, /noaxes, position=pos[*,4*ib+2], /noerase, /save
                      chi2 = strtrim(string(cat[these[ig]].decam_rchi2[bandindx[ib]],format='(F12.2)'),2)
                      im_legend, '\chi^{2}_{\nu}='+chi2, /left, /bottom, box=0, charsize=1.2, $
                        charthick=4.0, position=[pos[0,4*ib+2]-0.01,pos[1,4*ib+2]+0.05], /norm, $
                        textcolor=cgcolor('white')
                      loadct, 0, /silent
                      
                      plotimage, sgresid, /preserve, /noaxes, position=pos[*,4*ib+3], /noerase, /save
                   endfor 
;                  cgps_close
                   im_plotconfig, psfile=psfile, /psclose, /png
                endfor 
                
; build the webpage for this brick
                pngfiles = brick+'-'+string(cat[these].objid,format='(I6.6)')+'.png'
                decals_edr2_html, html_base, html_path=html_path, npngcols=4, $
                  title='Brick '+brick+'-'+allmytype[it], /nomenu, ra=cat[these].ra, $
                  dec=cat[these].dec, pnglist=pngfiles
             endelse
stop
          endfor 
       endfor 
    endfor
    
return
end


;; make a scatterplot of r-band magnitudes and build cutouts of some
;; random galaxies
;             djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;               xrange=[14,25], yrange=[-2,2], xtitle='r_{SDSS} (AB mag)', $
;               ytitle='\Delta'+'r (DECaLS - SDSS, AB mag)'
;             djs_oplot, !x.crange, [0,0], line=0, color=cgcolor('grey')
;             djs_oplot, rmag, rmag_decals-rmag, psym=symcat(16), symsize=0.5
;             djs_oplot, rmag[these], rmag_decals[these]-rmag[these], $
;               psym=symcat(15), symsize=1.3, color=cgcolor('orange')
;             im_legend, 'Brick '+brick, /left, /top, box=0
             
;         if keyword_set(stars) then outpath = cutoutpath+brick+'-stars/' else $
;           outpath = cutoutpath+brick+'-galaxies/'
;         file_mkdir, outpath
          

;                     pixscale = sqrt(abs(determ(astr.cd)))*3600D ; [arcsec/pixel]
;                        extast, hdr, astr
;                        ad2xy, cat[these[ig]].ra, cat[these[ig]].dec, astr, xx, yy

;                     hextract, image, hdr, gim, ghdr, x0, x1, y0, y1, /silent
;                     hextract, model, hdr, gmodel, ghdr1, x0, x1, y0, y1, /silent
;                     hextract, resid, hdr, gresid, ghdr1, x0, x1, y0, y1, /silent
                      
