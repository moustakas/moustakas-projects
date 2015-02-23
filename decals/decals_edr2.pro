function get_type, type, chi2=chi2
    ref = ['','S','D','E','C']
    stype = ['Nothing','Star','de Vac','Exp','Comp']
    this = where(type eq ref)
;   type = stype[this]
    type = stype[this]+' (\chi^{2}_{\nu}='+$
      strtrim(string(chi2[2],format='(F12.2)'),2)+')'
return, type
end

pro decals_edr2, stars=stars, html_only=html_only
; jm15feb13siena - test the preliminary EDR2 catalogs

    edrpath = getenv('HOME')+'/edr2/'
    cutoutpath = getenv('HOME')+'/cutouts/'
    
    cut = 15.0                  ; [arcsec]
    pixscale = 0.262D           ; [arcsec/pixel]
    nmosaic = 3600              ; size of the coadds [pixels]
    npad = 60                   ; ignore objects near the edges [pixels]
    
    alltype = ['S','D','E','C']
    allmytype = ['Star','deVac','Exp','Comp']
    ntype = n_elements(alltype)

    band = ['g','r','z']
    nband = n_elements(band)
    
    allbrick = file_basename(file_search(edrpath+'tractor/*',/test_dir,count=nbrick))
    for ib = 0L, nbrick-1 do begin
       catfile = file_search(edrpath+'tractor/'+allbrick[ib]+'/tractor-*.fits',count=ncat)
       for ic = 0L, ncat-1 do begin
          cat = mrdfits(catfile[ic],1)
          brick = strtrim(cat[0].brickname,2)
          coaddpath = edrpath+'coadd/'+allbrick[ib]+'/'+brick+'/'

          file_mkdir, cutoutpath+brick+'/' ; base output path

; (re)build just the HTML pages
          if keyword_set(html_only) then begin
stop
             if keyword_set(stars) then begin
                decals_edr2_html, brick+'-stars', html_path=cutoutpath, $
                  title='Brick '+brick+' (stars)'
             endif else begin
                decals_edr2_html, brick+'-galaxies', html_path=cutoutpath, $
                  title='Brick '+brick+' (galaxies)'
             endelse
             return
          endif else begin

;; choose N random galaxies for closer analysis
;             nran = 40
;             if keyword_set(stars) then begin
;                inrange = where(rmag gt 15.0 and rmag lt 22.0 and cat.type eq 'S',ngal)
;             endif else begin
;                inrange = where(rmag gt 15.0 and rmag lt 22.0 and cat.type ne 'S',ngal)
;             endelse
;             these = inrange[long(randomu(seed,nran)*ngal)]
;             these = where(cat.objid eq 5330)
;             nran = 1
;;            these = [1933,2140]
;;            nran = n_elements(these)

; one set of webpages per type
             for it = 0, ntype-1 do begin
                these = where(cat.type eq alltype[it] and cat.tx gt npad and $
                  cat.tx lt nmosaic-npad and cat.ty gt npad and $
                  cat.ty lt nmosaic-npad and cat.sdss_objid ne 0,nobj)
                srt = sort(cat[these].decam_flux[1]) ; sort by r-band flux
                these = these[srt]

; parse the photometry             
                decals_to_maggies, cat[these], maggies, ivarmaggies
                mag = maggies2mag(maggies,ivarmaggies=ivarmaggies,magerr=magerr)
                
; make directories
                html_path = cutoutpath+brick+'/'+brick+'-'+allmytype[it]+'/'
                file_mkdir, html_path
                
; get image cutouts in all three bands
                for ig = 0L, nobj-1 do begin
                   obj = string(cat[these[ig]].objid,format='(I6.6)')
                   galaxy = brick+'-'+obj

                   for ib = 0, nband-1 do begin
                      imfile = coaddpath+'decals-'+brick+'-image-'+band[ib]+'.fits'
                      if file_test(imfile) eq 0 then stop
                
                      image = mrdfits(imfile,0,hdr)
                      model = mrdfits(repstr(imfile,'-image','-model'),0)
                      resid = image - model

                      extast, hdr, astr
                      ncutpix = long(cut/pixscale)

;                     type = get_type(cat[these[ig]].type,chi2=cat[these[ig]].decam_rchi2)
                         
                      ad2xy, cat[these[ig]].ra, cat[these[ig]].dec, astr, xx, yy
                      x0 = long(xx-ncutpix/2.0)
                      x1 = long(xx+ncutpix/2.0)
                      y0 = long(yy-ncutpix/2.0)
                      y1 = long(yy+ncutpix/2.0)
stop                      
                      hextract, image, hdr, gim, ghdr, x0, x1, y0, y1, /silent
                      hextract, model, hdr, gmodel, ghdr1, x0, x1, y0, y1, /silent
                      hextract, resid, hdr, gresid, ghdr1, x0, x1, y0, y1, /silent
                         
                      med = median(gim)
                      sig = stdev(gresid)
                      mn = med - 2.0*sig
                      mx = med + 7.0*sig
                      splog, minmax(gim), med, sig, mn, mx
                      beta = 0.01
                      sgim    = gmascl(   gim,min=mn,max=mx,/negative,gamma=gamma)
                      sgmodel = gmascl(gmodel,min=mn,max=mx,/negative,gamma=gamma)
                      sgresid = gmascl(gresid,min=mn,max=mx,/negative,gamma=gamma)
                      sgim    = asinhscl(   gim,min=mn,max=mx,/negative,beta=beta)
                      sgmodel = asinhscl(gmodel,min=mn,max=mx,/negative,beta=beta)
                      sgresid = asinhscl(gresid,min=mn,max=mx,/negative,beta=beta)
                      
                      loadct, 0, /silent
                      psfile = outpath+galaxy+'.png'
                      cgps_open, psfile
                      cgdisplay, aspect=1.0/3.0
                      pos = cglayout([3,1],aspect=1.0,xgap=3.0,ygap=0,$
                        oxmargin=[1,1],oymargin=[1,1])
                      cgimage, sgim, position=pos[*,0]
                      im_legend, [galaxy,'Type='+type,$
                        'r='+string(rmag_decals[these[ig]],format='(F5.2)')], $
                        /left, /top, box=0, charsize=1.3, charthick=4.0, margin=0, $
                        position=[0.01,0.93]
                      loadct, 0, /silent
                      
                      cgimage, sgmodel, position=pos[*,1], /noerase
                      cgimage, sgresid, position=pos[*,2], /noerase
                      cgps_close
                   endfor 
                endfor 

; build the webpage for this brick
                if keyword_set(stars) then begin
                   decals_edr2_html, brick+'-stars', html_path=cutoutpath, $
                     title='Brick '+brick+' (stars)'
                endif else begin
                   decals_edr2_html, brick+'-galaxies', html_path=cutoutpath, $
                     title='Brick '+brick+' (galaxies)'
                endelse
             endfor
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
