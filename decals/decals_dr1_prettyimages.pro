pro decals_dr1_prettyimages, make_dr1bricks=make_dr1bricks, get_cutouts=get_cutouts, $
  make_png=make_png
; jm15mar19siena - build some pretty pictures for DR1

    dr1dir = getenv('DECALS_DIR')+'/'
    outdir = getenv('HOME')+'/dr1-color/'
    if file_test(outdir,/dir) eq 0 then file_mkdir, outdir
    if file_test(outdir+'fits',/dir) eq 0 then file_mkdir, outdir+'fits'
    
    band = ['g','r','z']
    nband = n_elements(band)

    pixscale = 0.262D           ; [arcsec/pixel]
    nmosaic = 3600
    
    width = 1.5D                ; cutout diameter [arcmin]
    link = 3D                   ; linking length [arcmin]
    
; --------------------------------------------------
; temporary hack to make the dr1-bricks.fits file
    if keyword_set(make_dr1bricks) then begin
       all = mrdfits(dr1dir+'decals-bricks.fits',1)
       allbrick = file_basename(file_search(dr1dir+'tractor/*',/test_dir,count=nbrick))
       for ii = 0L, nbrick-1 do begin
          print, format='("Brick ",I0,"/",I0, A10,$)', ii+1, nbrick, string(13b)
          catfile = file_search(dr1dir+'tractor/'+allbrick[ii]+'/tractor-*.fits')
          brickname1 = repstr(repstr(file_basename(catfile),'tractor-',''),'.fits','')
          if n_elements(brickname) eq 0L then brickname = brickname1 else $
            brickname = [brickname,brickname1]
; check for bricks with 3-band coverage
          delvarx, cat
          for ic = 0, n_elements(catfile)-1 do begin
             cat1 = mrdfits(catfile[ic],1,/silent,rows=0)
             if n_elements(cat) eq 0L then cat = cat1 else cat = [cat,cat1]
          endfor
          grz = where(total(cat.decam_nobs gt 0,1) ge 3)
          if grz[0] ne -1 then begin
             brickname1_grz = brickname1[grz]
             if n_elements(brickname_grz) eq 0L then brickname_grz = brickname1_grz else $
               brickname_grz = [brickname_grz,brickname1_grz]
          endif
       endfor
; write out the list of grz dr1 bricks       
       match, all.brickname, brickname_grz, these_grz
       mwrfits, all[these_grz], dr1dir+'dr1-grz-bricks.fits', /create
; write out the list of all dr1 bricks       
       match, all.brickname, brickname, these
       mwrfits, all[these], dr1dir+'dr1-bricks.fits', /create
    endif

; --------------------------------------------------
; find matching objects
    if keyword_set(get_cutouts) then begin
       bricks = mrdfits(dr1dir+'dr1-grz-bricks.fits',1)
;      allbigcat = read_ngc()

; put together a list of "big" objects we want pretty pictures of and
; an estimate of their diameters       
       allbigcat = mrdfits(getenv('CATALOGS_DIR')+'/ngc/ngc2000.fits',1)
       allbigcat = struct_addtags(replicate({name: ''},$
         n_elements(allbigcat)),allbigcat)
       allbigcat.name = 'NGC'+string(allbigcat.ngcnum,format='(I4.4)')
;      allbigcat = read_rc3()
;      big = where(allbigcat.d25_maj gt 0.75)
;      allbigcat = allbigcat[big]

; get a list of all the bricks that contain a "big" RC3 object       
       spherematch, allbigcat.ra, allbigcat.dec, bricks.ra, $
         bricks.dec, 0.25/2, ibigcat, ibrick, maxmatch=0
       unq = uniq(ibrick,sort(ibrick))
       ibrick = ibrick[unq]

;      djs_plot, bricks.ra, bricks.dec, psym=6
;      djs_oplot, allbigcat[ibigcat].ra, allbigcat[ibigcat].dec, psym=6, color='orange'
       
; loop on each matching (unique) brick and demand 3-color photometry
       nbrick = n_elements(ibrick)
       for jj = 0L, nbrick-1 do begin
          thisbrick = bricks[ibrick[jj]].brickname
          subdir = strmid(thisbrick,0,3)
          tractordir = dr1dir+'tractor/'+subdir+'/'
          coadddir = dr1dir+'coadd/'+subdir+'/'+thisbrick+'/'

          imfile = file_search(coadddir+'decals-'+thisbrick+'-image-'+band+'.fits',count=nfile)
          if nfile ne nband then stop
          if nfile eq nband then begin
             splog, subdir+'-'+thisbrick

; save time by reading the images in all three bands
             image = fltarr(nmosaic,nmosaic,nband)
             for ib = 0, nband-1 do image[*,*,ib] = mrdfits(imfile[ib],0,hdr)
             
; match back to the catalog             
             spherematch, allbigcat.ra, allbigcat.dec, bricks[ibrick[jj]].ra, $
               bricks[ibrick[jj]].dec, 0.25D/2, m1, m2, maxmatch=0
             if m1[0] eq -1 then message, 'This should not happen!'
             bigcat = allbigcat[m1]

;            allcat = mrdfits(tractordir+'tractor-'+thisbrick+'.fits',1)
;            keep = where(allcat.brick_primary eq 'T' and $
;              total(allcat[0].decam_nobs,1) ge 3)
             
; group the matches on 1 arcminute scales and get cutouts centered on
; the center of the group (to take care of multiple close systems)
             grp = spheregroup(bigcat.ra,bigcat.dec,link/60D)
             ugrp = grp[uniq(grp,sort(grp))]
             for ig = 0, n_elements(ugrp)-1 do begin
                these = where(ugrp[ig] eq grp)
                outname = strcompress(strjoin(bigcat[these].name,'-'),/remove)
                if outname eq '' then outname = strcompress(strjoin(bigcat[these].altname,'-'),/remove)
                if outname eq '' then stop
                
                racen = djs_mean(bigcat[these].ra)
                deccen = djs_mean(bigcat[these].dec)
                adxy, hdr, racen, deccen, xx, yy
                x0 = (xx-width*60D/pixscale)>0
                x1 = (xx+width*60D/pixscale)<(nmosaic-1)
                y0 = (yy-width*60D/pixscale)>0
                y1 = (yy+width*60D/pixscale)<(nmosaic-1)
                
                for ib = 0, nband-1 do begin
                   hextract, reform(image[*,*,ib]), hdr, newim, newhdr, $
                     x0, x1, y0, y1, /silent                   
                   outfile = outdir+'fits/'+outname+'-'+band[ib]+'.fits'
                   splog, 'Writing '+outfile
                   mwrfits, newim, outfile, newhdr, /create
                endfor 
             endfor 
          endif 
       endfor 
    endif 
       
; --------------------------------------------------
; make pretty pictures cutouts
    if keyword_set(make_png) then begin
       allfile = file_search(outdir+'fits/*-r.fits',count=nobj)
       for ii = 0, nobj-1 do begin
          name = repstr(file_basename(allfile[ii]),'-r.fits','')
          pushd, outdir

; clean up previously created parameter files          
          file_delete, outdir+['levels.txt', 'trilogyfilterlog.txt',$
            name+'-image.in',name+'-image_filters.txt'], /quiet
          
          infile = outdir+name+'-image.in'
          openw, lun, infile, /get_lun
          printf, lun, 'B'
          printf, lun, file_search(outdir+'fits/'+name+'-g.fits')
          printf, lun, ''
          printf, lun, 'G'
          printf, lun, file_search(outdir+'fits/'+name+'-r.fits')
          printf, lun, ''
          printf, lun, 'R'
          printf, lun, file_search(outdir+'fits/'+name+'-z.fits')
          printf, lun, ''
          printf, lun, 'indir '+outdir+'fits/'
          printf, lun, 'outname '+outdir+name+'-image'
          printf, lun, 'noiselum 0.15'
;         printf, lun, 'scaling  '+outdir+'levels.txt'
;         printf, lun, 'scaling  '+getenv('CLASH_ARCHIVE')+'/color_images/scaling/levels_ir.txt'
          printf, lun, 'show 0'
          printf, lun, 'legend 0'
          printf, lun, 'testfirst 0'
          free_lun, lun

          cmd = 'python '+getenv('IMPY_DIR')+'/image/trilogy.py '+infile
          splog, cmd
          spawn, cmd, /sh
          popd

; clean up the parameter files          
          file_delete, outdir+['levels.txt', 'trilogyfilterlog.txt',$
            name+'-image.in',name+'-image_filters.txt'], /quiet
       endfor

    endif

;    jj = mrdfits('dr1/dr1-grz-bricks.fits',1)
;    openw, lun, 'decals-grz.txt', /get_lun
;    niceprintf, lun, 'decals-'+jj.brickname+'-image-*.fits'
;    free_lun, lun
;    
;    rsync -auvPn --delete --include='*/' --include-from='decals-grz.txt' --exclude='*' nersc:'/project/projectdirs/cosmo/work/decam/release/dr1/coadd' '/global/work/decam/release/dr1/'
    
return
end
