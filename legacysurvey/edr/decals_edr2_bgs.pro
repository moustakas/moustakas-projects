pro decals_edr2_bgs, html_only=html_only
; jm15feb27siena - test the preliminary EDR2 catalogs for the BGS sample

    edrpath = '/global/work/decam/cats/edr2/'
    outpath = edrpath+'bgs/'
    
    cut = 15.0                  ; [arcsec]
    pixscale = 0.262D           ; [arcsec/pixel]
    ncutpix = long(cut/pixscale)
    
    nmosaic = 3600              ; size of the coadds [pixels]
    npad = 60                   ; ignore objects near the edges [pixels]
    beta = 0.01
    
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
       for ic = 0, 0 do begin
;      for ic = 0L, ncat-1 do begin
          cat = mrdfits(catfile[ic],1)
          brick = strtrim(cat[0].brickname,2)
          coaddpath = edrpath+'coadd/'+allbrick[ii]+'/'+brick+'/'

          html_path = outpath;+brick+'/'
          file_mkdir, html_path ; base output path

; read the coadds once
          if keyword_set(html_only) eq 0 then begin
             image = fltarr(nmosaic,nmosaic,nband)
             model = image
             for ib = 0, nband-1 do begin
                imfile = coaddpath+'decals-'+brick+'-image-'+band[ib]+'.fits'
                image[*,*,ib] = mrdfits(imfile,0,hdr)
                model[*,*,ib] = mrdfits(repstr(imfile,'-image','-model'),0)
             endfor
             resid = image - model
          endif

; parse the photometry             
          decals_to_maggies, cat, maggies, ivarmaggies
          mag = maggies2mag(maggies,ivarmaggies=ivarmaggies,$
            magerr=magerr,magnsigma=magnsigma)
          
          html_base = brick+'-BGS'
          png_path = html_path+html_base+'/'
          file_mkdir, png_path

          these = where(cat.type ne 'S' and cat.bx gt npad and cat.bx lt nmosaic-npad and $
            cat.by gt npad and cat.by lt nmosaic-npad and cat.sdss_objid ne 0 and $
            mag[1,*] gt 0 and mag[1,*] lt 19.5,nobj)
          srt = reverse(sort(cat[these].decam_flux[2])) ; sort by r-band flux
          these = these[srt]
;         these = these[0:2] & nobj = 3

; (re)build just the HTML pages
          if keyword_set(html_only) then begin
             pngfiles = brick+'-'+string(cat[these].objid,format='(I6.6)')+'.png'
;            openw, lun, 'junk', /get_lun & niceprintf, lun, pngfiles & free_lun, lun
             decals_edr2_html, html_base, html_path=html_path, npngcols=5, $
               title='Brick '+brick+'-BGS (r<19.5 & type!=Star)', /nomenu, $
               pnglist=pngfiles, ra=cat[these].ra, dec=cat[these].dec
          endif else begin
             xall = cat[these].bx-1
             yall = cat[these].by-1

; get image cutouts in all three bands
             for ig = 0L, nobj-1 do begin
                obj = brick+'-'+string(cat[these[ig]].objid,format='(I6.6)')
                splog, ig, obj
                xx = xall[ig]
                yy = yall[ig]
                
                loadct, 0, /silent
                cgps_open, png_path+obj+'.png'
                cgdisplay, aspect=1.0 ;/3.0
                pos = cglayout([3,nband],aspect=1.0,xgap=0.0,ygap=1,$
                  oxmargin=[1,1],oymargin=[1,1])
                for ib = 0, nband-1 do begin
                   gim = image[long(xx-ncutpix/2.0):long(xx+ncutpix/2.0),$
                     long(yy-ncutpix/2.0):long(yy+ncutpix/2.0),ib]
                   gmodel = model[long(xx-ncutpix/2.0):long(xx+ncutpix/2.0),$
                     long(yy-ncutpix/2.0):long(yy+ncutpix/2.0),ib]
                   gresid = resid[long(xx-ncutpix/2.0):long(xx+ncutpix/2.0),$
                     long(yy-ncutpix/2.0):long(yy+ncutpix/2.0),ib]
                   
                   med = djs_median(gim)
                   sig = djsig(gresid)
                   mn = med - 2.0*sig
                   mx = med + 7.0*sig
                   sgim    = asinhscl(   gim,min=mn,max=mx,/negative,beta=beta)
                   sgmodel = asinhscl(gmodel,min=mn,max=mx,/negative,beta=beta)
                   sgresid = asinhscl(gresid,min=mn,max=mx,/negative,beta=beta)
                   
                   cgimage, sgim, position=pos[*,3*ib+0], noerase=ib gt 0, /save
                   plots, ncutpix/2.0, ncutpix/2.0, psym=cgsymcat(9,thick=5), $
                     symsize=1.2, color=cgcolor('black')
                   djs_oplot, cat.bx-xx+ncutpix/2.0, cat.by-yy+ncutpix/2.0, psym=cgsymcat(16), $
                     symsize=1.0, color=cgcolor('dodger blue')
                   im_legend, obj, /left, /top, box=0, charsize=1.3, $
                     charthick=8.0, margin=0, position=[0.01,0.96], /norm
                   if mag[ib,these[ig]] gt 0 then begin
                      thismag = band[ib]+'='+strtrim(string(mag[ib,these[ig]],format='(F12.3)'),2)+$
                        '\pm'+strtrim(string(magerr[ib,these[ig]],format='(F12.3)'),2)
                   endif else begin
                      thismag = band[ib]+'>'+string(magnsigma[ib,these[ig]],format='(F5.2)')
                   endelse
                   im_legend, thismag, /left, /bottom, box=0, charsize=1.5, $
                     charthick=8.0, position=[pos[0,3*ib+0]-0.01,pos[1,3*ib+0]+0.01], /norm
                   loadct, 0, /silent
                   
                   cgimage, sgmodel, position=pos[*,3*ib+1], /noerase, /save
                   chi2 = strtrim(string(cat[these[ig]].decam_rchi2[bandindx[ib]],format='(F12.2)'),2)
                   im_legend, '\chi^{2}_{\nu}='+chi2, /left, /bottom, box=0, charsize=1.5, $
                     charthick=8.0, position=[pos[0,3*ib+1]-0.01,pos[1,3*ib+1]+0.01], /norm
                   loadct, 0, /silent
                   
                   cgimage, sgresid, position=pos[*,3*ib+2], /noerase, /save
                endfor 
                cgps_close
             endfor 
             
; build the webpage for this brick
             decals_edr2_html, html_base, html_path=html_path, npngcols=4, $
               title='Brick '+brick+'-'+allmytype[it], /nomenu
          endelse 
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
                      
