pro label_png, image, outfile=outfile, galaxy=galaxy
    sz = size(image,/dim)
    cgps_open, outfile
    cgdisplay, sz[1], sz[2]
    cgimage, image;, /keep_aspect
    if n_elements(galaxy) ne 0 then cgtext, 0.06, 0.92, galaxy, $
      /normal, color='white', charsize=2.5, charthick=1.8, font=1, $
      tt_font='Times'
    cgps_close, resize=100;, /showcmd
return
end

pro make_html, htmlfile, pnglist=pnglist, npngcols=npngcols, $
  title=title, sample=sample1

    if (n_elements(npngcols) eq 0L) then npngcols = 3L

    html_path = file_dirname(htmlfile)
    png_path = html_path+'/png/'

    pngfiles = strtrim(sample1.object,2)+'.png'
    good = where(file_test(png_path+pngfiles) eq 1)
    sample = sample1[good]
    pngfiles = pngfiles[good]
    
    npngfiles = n_elements(pngfiles)
    
;    if (n_elements(pnglist) eq 0L) then begin
;       if keyword_set(or_jpeg) then $
;         pngfiles = [file_search(png_path+'*.png'),$
;         file_search(png_path+'*.jpg')] else $
;           pngfiles = file_search(png_path+'*.png')
;       npngfiles = n_elements(pngfiles)
;       if (npngfiles eq 0L) then begin
;          splog, 'No PNG files found in '+png_path
;          return
;       endif else pngfiles = file_basename(pngfiles)
;    endif else begin
;       pngfiles = pnglist
;       npngfiles = n_elements(pngfiles)
;;      tempfiles = file_search(png_path+pnglist,count=npngfiles)
;;      if (npngfiles ne n_elements(pnglist)) then begin
;;         splog, 'Some files in PNGLIST do not exist in '+png_path
;;         return
;;      endif
;    endelse

    if (not keyword_set(silent)) then splog, 'Found '+$
      string(npngfiles,format='(I0)')+' PNG files'

; generate the web page HTML page
    rootpngfiles = repstr(pngfiles,'.png','')
    nmenucols = npngcols

    if (not keyword_set(only_png)) then begin

       xwidth = string(fix(100/float(npngcols)),format='(I0)')
;      xwidth = string(fix(100/float(npscols)-npscols*3.0-npscols*15.0),format='(I0)')
;      xwidth = string(fix(100/float(npscols)-npscols*5.0),format='(I0)')
;      xwidth = '100'
       
       if (not keyword_set(silent)) then splog, 'Writing '+htmlfile
       
       openw, lun1, htmlfile,/get_lun
       printf, lun1, '<html>'
       printf, lun1, '<head>'
;      if (n_elements(title) eq 0L) then title = $
;        'Postscript Output for '+strupcase(htmlbase)
       printf, lun1, '<title>'+title+'</title>'
       printf, lun1, '<style type="text/css">'
       printf, lun1, '<!--'
       printf, lun1, 'body {'
       printf, lun1, '   color: white;'
       printf, lun1, '   background-color: black;'
       printf, lun1, '   text-align: center;'
       printf, lun1, '   border: 0;'
       printf, lun1, '   font-size: 100%;'
       printf, lun1, '   font-weight: bold;'
       printf, lun1, '}'
       printf, lun1, 'td {'
       printf, lun1, '   font-size: 100%;'
       printf, lun1, '   font-weight: bold;'
       printf, lun1, '   text-align: center;'
       printf, lun1, '}'
       printf, lun1, 'a {'
       printf, lun1, '   color: cyan;'
       printf, lun1, '   visited: cyan;'
       printf, lun1, '   text-decoration: none;'
       printf, lun1, '}'
       printf, lun1, '.center {'
       printf, lun1, '   text-align: center;'
       printf, lun1, '}'
       printf, lun1, '.left {'
       printf, lun1, '   text-align: left;'
       printf, lun1, '}'
       printf, lun1, '.smaller {'
       printf, lun1, '   font-size: 90%;'
       printf, lun1, '}'
       printf, lun1, 'h1 {color: gold;}'
       printf, lun1, '</style>'
       printf, lun1, '-->'
       printf, lun1, ' '
       printf, lun1, '<link href="http://legacysurvey.org/test/assets/css/bootstrap.min.css" rel="stylesheet" type="text/css">'
       printf, lun1, '<link href="http://legacysurvey.org/test/assets/css/rst.css" rel="stylesheet" type="text/css">'
       printf, lun1, '<link href="http://legacysurvey.org/test/assets/css/code.css" rel="stylesheet" type="text/css">'
       printf, lun1, '<link href="http://legacysurvey.org/test/assets/css/colorbox.css" rel="stylesheet" type="text/css">'
       printf, lun1, '<link href="http://legacysurvey.org/test/assets/css/theme.css" rel="stylesheet" type="text/css">'
       printf, lun1, '<link href="http://legacysurvey.org/test/assets/css/custom.css" rel="stylesheet" type="text/css">'
       printf, lun1, ' '
       printf, lun1, '</head>'
       printf, lun1, '<body>'
       printf, lun1, '</br>'
       printf, lun1, '<h1>'+title+'</h1>'
       printf, lun1, '</br>'

;      printf, lun1, '<br />'
;      printf, lun1, '<p class="left">This is text.</p>'
       printf, lun1, '<br />'
       printf, lun1, '<table align="center" cellpadding="5" border="1" cellspacing="5" width="100%">'
       printf, lun1, '<tbody>'

       for i = 0L, ceil(npngfiles/float(npngcols))-1 do begin
; images          
          printf, lun1, '<tr>'
          for j = 0L, npngcols-1 do begin
             indx = j + i*npngcols
             if (indx le npngfiles-1L) then begin
                printf, lun1, '<td width="'+xwidth+'%"><a href="png/'+pngfiles[indx]+'">'+$
                  '<img width="100%" align="center" valign="center" src="png/'+pngfiles[indx]+'"></a></td>'
             endif
          endfor
          printf, lun1, '</tr>'

; image captions          
          printf, lun1, '<tr>'
          for j = 0L, npngcols-1L do begin
             indx = j + i*npngcols
             if (indx le npngfiles-1L) then begin
;               printf, lun1, '<td width="'+xwidth+'%" class="center">'+string(j+1,format='(I0)')+':'+$
;                 string(i+1,format='(I0)')
                printf, lun1, '<td width="'+xwidth+'%" class="left">'+$
                  '<a target="_blank" href="http://legacysurvey.org/viewer/?layer=decals-dr2&ra='+$
                  strtrim(string(sample[indx].ra,format='(F12.4)'),2)+'&dec='+$
                  strtrim(string(sample[indx].dec,format='(F12.5)'),2)+'&zoom=14">'+$
                  strtrim(sample[indx].object,2)+'</a>'
;                 strtrim(repstr(sample[indx].object,'-','/'),2)+'</a>'
;               printf, lun1, rootpngfiles[indx]+' <a href="'+htmlbase+'/'+pngfiles[indx]+'">png</a>'
                printf, lun1, '</td>'
             endif
          endfor
          printf, lun1, '</tr>'
          
       endfor 

       printf, lun1, '</tbody>'
       printf, lun1, '</table>'
       printf, lun1, '</body>'
       printf, lun1, '</html>'
       free_lun, lun1

    endif 

return
end


function get_template, nobj
    template = {id: 0L, object: '', ra: 0D, dec: 0D, $
      radius: 0.0, brickname: '', origin: ''}
return, replicate(template,nobj)
end

function get_ngc
; get and parse the catalog of NGC galaxies
    cat = mrdfits(getenv('CATALOGS_DIR')+'/ngc/ngc2000.fits',1)
    template = get_template(n_elements(cat))
    template.object = 'NGC'+strtrim(cat.ngcnum,2)
    template.ra = cat.ra
    template.dec = cat.dec
    template.radius = (cat.radius*60)>3.0
    template.origin = 'NGC'
; fix some coordinates
    fix = where(strmatch(template.object,'*2783*'))
    template[fix].ra = 15D*hms2dec('09h13m39.4s')
    template[fix].dec = hms2dec('+29d59m35s')
;; add in NGC2783b
;    template1 = get_template(1)
;    template1.object = 'NGC2783B'
;    template1.ra = 
;    template = [template,template1]
return, template
end

function get_abell
; get and parse the catalog of Abell clusters
    cat = mrdfits(getenv('IM_GITREPOS')+'/astrometry.net/catalogs/abell-all.fits',1)
    template = get_template(n_elements(cat))
    template.object = 'NGC'+strtrim(cat.ngcnum,2)
    template.ra = cat.ra
    template.dec = cat.dec
    template.radius = (cat.radius*60)>1.5
    template.origin = 'NGC'
return, template
end

function get_messier
; get and parse the catalog of Messier objects
    cat = rsex(getenv('CATALOGS_DIR')+'/ngc/messier.cat')
    template = get_template(n_elements(cat))
    template.object = 'M'+strtrim(cat.messier,2)
    template.ra = 15D*hms2dec(cat.ra)
    template.dec = hms2dec(cat.dec)
    template.radius = 3.0
    template.origin = 'Messier'
return, template
end

function get_rc3, mindiam=mindiam
; get and parse the catalog of RC3 objects

    if n_elements(mindiam) eq 0 then mindiam = 2.0 ; [arcmin]
    
    cat = mrdfits(getenv('CATALOGS_DIR')+'/rc3/rc3_catalog.fits',1)
;   cat = read_rc3()

    d25_maj = 0.1*10.0^cat.logd_25 ; [arcmin]

    keep = where(d25_maj gt mindiam)
    cat = cat[keep]
    d25_maj = d25_maj[keep]
    
    template = get_template(n_elements(cat))

    fix = where(strmatch(strtrim(cat.name1,2),'A*'))
    cat[fix].name1 = strtrim(cat[fix].name2) ; replace A names with UGC 

    template.object = strcompress(cat.name1,/remove)
    need = where(strcompress(template.object,/remove) eq '')
    if need[0] ne -1 then template[need].object = $
      strcompress(cat[need].name2,/remove)
    need = where(strcompress(template.object,/remove) eq '')
    if need[0] ne -1 then template[need].object = $
      strcompress(cat[need].name3,/remove)
    need = where(strcompress(template.object,/remove) eq '')
    if need[0] ne -1 then template[need].object = $
      strcompress(cat[need].pgc,/remove)    
    
    template.object = repstr(template.object,'?','')    
    template.ra = cat.ra
    template.dec = cat.dec
    template.radius = d25_maj ; [arcmin]
    template.origin = 'RC3'

; fix some issues here
    fix = where(strtrim(template.object,2) eq 'PGC36594')
    if fix[0] ne 0 then begin
       template[fix].ra = hms2dec('11h44m54.3s')*15D
       template[fix].dec = hms2dec('+02d09m47s')
    endif
    
return, template
end
;Pegasus dIrr, ~800 kpc
;M_B = -11.2, 1.2d8 Msun in stars
;352.14208       14.746667 
;
;Leo A, ~800 kpc
;M_B = -11.6, 7d7 Msun in stars
;149.86000       30.746389

function resolve_duplicates, sample, debug=debug
; resolve duplicate objects
    grp = spheregroup(sample.ra,sample.dec,10D/3600D)
    ngrp = histogram(grp,reverse_indices=ri)
    ivec = ri[0:n_elements(ngrp)-1]
    
    oneobj = where(ngrp eq 1)   ; only one object - keep it!
    out_sample = sample[ri[ivec[oneobj]]]
    
; more than 1 object
    multobj = where(ngrp gt 1,nmult)
    for ii = 0L, nmult-1 do begin
       print, format='("Resolving duplicates ",I0,"/",I0, A10,$)', ii+1, nmult, string(13b)
       match = where(grp eq multobj[ii],nmatch)
       this = where(strtrim(sample[match].origin,2) eq 'RC3')
       if this[0] eq -1 then this = where(strtrim(sample[match].origin,2) eq 'Messier')
       if this[0] eq -1 then this = where(strtrim(sample[match].origin,2) eq 'NGC')
       out_sample = [out_sample,sample[match[this]]]
       if keyword_set(debug) then begin
          struct_print, sample[match] & struct_print, sample[match[this[0]]], /no_head
          print & cc = get_kbrd(1)
       endif
    endfor          
return, out_sample
end

function group_objects, sample, debug=debug
; group distinct objects which are near one another in the sky
    grp = spheregroup(sample.ra,sample.dec,1.5D/60D)
    ngrp = histogram(grp,reverse_indices=ri)
    ivec = ri[0:n_elements(ngrp)-1]
    
    oneobj = where(ngrp eq 1)   ; only one object - keep it!
    out_sample = sample[ri[ivec[oneobj]]]
    
; more than 1 object
    sample1 = get_template(1)
    multobj = where(ngrp gt 1,nmult)
    for ii = 0L, nmult-1 do begin
       print, format='("Grouping objects ",I0,"/",I0, A10,$)', ii+1, nmult, string(13b)
       match = where(grp eq multobj[ii],nmatch)
       obj = strtrim(sample[match].object,2)
       obj = obj[uniq(obj,sort(obj))]
       sample1.object = strjoin(obj,'-')
       sample1.ra = djs_mean(sample[match].ra)
       sample1.dec = djs_mean(sample[match].dec)
       sample1.origin = sample[match[0]].origin
       dra = (max(sample[match].ra)-min(sample[match].ra))/cos(sample[match[0]].dec*!dtor)
       ddec = max(sample[match].dec)-min(sample[match].dec)
       sample1.radius = ((dra>ddec)*60*1.5)>3.0
       out_sample = [out_sample,sample1]
       if keyword_set(debug) then begin
          struct_print, sample[match] & struct_print, sample1, /no_head
          print & cc = get_kbrd(1)
       endif
    endfor          
return, out_sample
end

pro decals_dr2_gallery, build_sample=build_sample, runbrick=runbrick, $
  make_png=make_png, html=html, debug=debug, mindiam=mindiam, clobber=clobber
; jm15mar19siena - build some pretty pictures for DR2

    dr = 'dr2'

    drdir = getenv('DECALS_DR2_DIR')+'/'
;   drdir = decals_path(dr=dr)
;   gallerydir = decals_path(dr=dr,/gallery)
;   gallerydir = '/project/projectdirs/cosmo/work/legacysurvey/gallery-dr1/'
    gallerydir = getenv('SCRATCH')+'/gallery-dr2/'
    if file_test(gallerydir,/dir) eq 0 then file_mkdir, gallerydir
    if file_test(gallerydir+'png',/dir) eq 0 then file_mkdir, gallerydir+'png'
    if file_test(gallerydir+'fits',/dir) eq 0 then file_mkdir, gallerydir+'fits'
    if file_test(gallerydir+'logs',/dir) eq 0 then file_mkdir, gallerydir+'logs'
    
    band = ['g','r','z']
    nband = n_elements(band)

    pixscale = 0.26D            ; [arcsec/pixel]
;   pixscale = 0.262D           ; [arcsec/pixel]
    barlen30 = fix(30.0/pixscale) ; [pixels]
    nmosaic = 3600

    if n_elements(mindiam) eq 0 then mindiam = 0.75 ; [arcmin]
    
; --------------------------------------------------
; build a sample of objects for the DR1 gallery
    if keyword_set(build_sample) then begin
       bricks = mrdfits(drdir+'decals-bricks-dr2.fits',1)
       keep = where(bricks.nobs_med_g gt 1 and bricks.nobs_med_r gt 1 and bricks.nobs_med_z gt 1)
       bricks = bricks[keep]

;      ngc = get_ngc()
;      messier = get_messier()
;      rc3 = get_rc3()
;      cat = [ngc,messier,rc3]
       cat = get_rc3(mindiam=mindiam)
       sample = resolve_duplicates(cat,debug=debug)
       splog, n_elements(sample), n_elements(cat)

;; throw out some junk
;       junk = where($
;         strtrim(sample.object,2) eq 'NGC7670'$ ; no object at this position (together with NGC7668,7669
;         ,comp=good)            ; 
;       sample = sample[good]

; group the matches on 1 arcminute scales and get cutouts centered on
; the center of the group (to take care of multiple close systems)
       sample = group_objects(sample,debug=debug)
       nobj = n_elements(sample)

; find matching bricks
       spherematch, sample.ra, sample.dec, bricks.ra, $
         bricks.dec, 0.25, icat, ibrick ;, maxmatch=0
       sample[icat].brickname = bricks[icat].brickname
       icat = icat[uniq(icat,sort(icat))]
       out_sample = sample[icat]
       out_sample = out_sample[sort(strtrim(out_sample.object,2))]
       ngal = n_elements(out_sample)

       out_sample.id = lindgen(ngal)
       struct_print, out_sample
       splog, ngal
       
;      djs_plot, sample.ra, sample.dec, psym=6, xsty=3, ysty=3, symsize=0.5
;      djs_oplot, bricks.ra, bricks.dec, psym=6, color='orange'
;      djs_oplot, out_sample.ra, out_sample.dec, psym=7, color='blue'
       im_mwrfits, out_sample, gallerydir+'gallery_sample.fits', /clobber
    endif 

; --------------------------------------------------
; run legacypipe centered on each galaxy
    if keyword_set(runbrick) then begin
       sample = mrdfits(gallerydir+'gallery_sample.fits.gz',1)
       gal = strcompress(sample.object,/remove)
       nobj = n_elements(sample)

;      for ii = 23, nobj-1 do begin
       for ii = 0L, nobj-1 do begin
          width = floor(sample[ii].radius*1.5*60.0/pixscale)
          logfile = gallerydir+'logs/'+gal[ii]+'.log'
          splog, logfile

          if file_test(gallerydir+'fits/'+gal[ii]+'-r.fits') eq 0 or keyword_set(clobber) then begin
             cmd = 'python -u ${LEGACYPIPE}/py/legacypipe/runbrick.py --no-write --splinesky --pixpsf '+$
               '--no-sdss --stage image_coadds --radec '+strtrim(sample[ii].ra,2)+' '+$
               strtrim(sample[ii].dec,2)+' --width '+strtrim(width,2)+' --height '+$
               strtrim(width,2)+' --pixscale '+strtrim(string(pixscale,format='(G0.0)'),2)+$
               ' > & '+logfile+' &'
             splog, cmd
             spawn, cmd
          endif
       endfor
    endif
       
; --------------------------------------------------
; make the three-color images
    if keyword_set(make_png) then begin
       sample = mrdfits(gallerydir+'gallery_sample.fits.gz',1)
       gal = strcompress(sample.object,/remove)
       nobj = n_elements(sample)

;      for ii = 56, 56 do begin
;      for ii = 26, nobj-1 do begin
       for ii = 0, nobj-1 do begin

          if sample[ii].dec lt 0.0 then pm = 'm' else pm = 'p'
          cusdir = 'custom-'+string(1000*sample[ii].ra,format='(I6.6)')+pm+$
            string(1000*abs(sample[ii].dec),format='(I5.5)')
          fitsfile = strarr(3)
          for jj = 0, 2 do fitsfile[jj] = gallerydir+'coadd/cus/'+cusdir+$
            '/decals-'+cusdir+'-image-'+band[jj]+'.fits'

          if total(file_test(fitsfile)) ne 3 then begin
             splog, gal[ii]
             if n_elements(prob) eq 0L then prob = ii else prob = [prob,ii]
          endif else begin
             if file_test(gallerydir+'png/'+gal[ii]+'.png') eq 0 or keyword_set(clobber) then begin
                for jj = 0, 2 do begin
                   spawn, 'cp -f '+gallerydir+'coadd/cus/'+cusdir+'/decals-'+$
                     cusdir+'-image-'+band[jj]+'.fits '+$
                     gallerydir+'fits/'+gal[ii]+'-'+band[jj]+'.fits', /sh
                endfor
          
                infile = gallerydir+'png/'+gal[ii]+'.in'
                openw, lun, infile, /get_lun
                printf, lun, 'B'
                printf, lun, file_search(gallerydir+'fits/'+gal[ii]+'-g.fits')
                printf, lun, ''
                printf, lun, 'G'
                printf, lun, file_search(gallerydir+'fits/'+gal[ii]+'-r.fits')
                printf, lun, ''
                printf, lun, 'R'
                printf, lun, file_search(gallerydir+'fits/'+gal[ii]+'-z.fits')
                printf, lun, ''
                printf, lun, 'indir '+gallerydir+'fits/'
                printf, lun, 'outname '+gallerydir+'png/'+gal[ii]
                printf, lun, 'noiselum 0.2'
;               printf, lun, 'scaling  '+gallerydir+'levels.txt'
                printf, lun, 'show 0'
                printf, lun, 'legend 0'
                printf, lun, 'testfirst 0'
                free_lun, lun
                
                cmd = 'python '+getenv('IMPY_DIR')+'/image/trilogy.py '+infile
                splog, cmd
                spawn, cmd, /sh
                
; clean up the parameter files          
                file_delete, gallerydir+'png/'+['levels.txt', 'trilogyfilterlog.txt',$
                  gal[ii]+'.in',gal[ii]+'_filters.txt'], /quiet
                file_delete, gallerydir+['levels.txt', 'trilogyfilterlog.txt'], /quiet
                
; read the image back in, label it, and put a scale bar
                im = read_png(gallerydir+'png/'+gal[ii]+'.png')
                sz = size(im,/dim)
                s0 = fix(sz[1]*0.05)
                ww = fix(sz[1]*0.0035)>2
                im[*,s0:s0+barlen30,s0:s0+ww] = cgcolor('white')
                label_png, im, galaxy=gal[ii], outfile=gallerydir+'png/'+gal[ii]+'.png'
             endif
          endelse 
       endfor
       if n_elements(prob) ne 0L then begin
          splog, 'Problem running runbrick '
          struct_print, sample[prob]
       endif
    endif

; --------------------------------------------------
; build a simple web page
    if keyword_set(html) then begin
       sample = mrdfits(gallerydir+'gallery_sample.fits.gz',1)
       gal = strcompress(sample.object,/remove)
       these = where(file_test(gallerydir+'png/'+gal+'.png'),nthese)

; missing data on these       
       skip = ['IC1256','IC2327','IC534','MCG0-32-4','NGC2576',$
         'NGC2743','NGC3323','NGC4382','NGC4536','NGC5534',$
         'NGC5656','NGC5806','NGC5838','NGC5850','NGC5956',$
         'NGC7046','NGC7081','UGC10176','UGC10547','UGC10831',$
         'UGC10837','UGC10905','UGC11859','UGC4285','UGC4340',$
         'UGC6440','UGC6821','UGC7170','UGC9117','UGC9760',$
         'UGC9900','UGC9977','IC1516','MCG-1-59-27','CGCG119-119',$
         'MCG2-55-4','NGC5983','UGC10023','UGC11807','UGC4571',$
         'UGC4884','UGC6379','UGC9782','UGC9949','MCG0-29-28']
       keep = these
       match, gal[these], skip, m1, m2
       remove, m1, keep

       title = 'DECaLS/DR2 Image Gallery'
       
       htmlfile = gallerydir+'index.html'
       make_html, htmlfile, npngcols=5, title=title, $
         pnglist=pngfiles, sample=sample[keep]
    endif
    
stop    
    
return
end

