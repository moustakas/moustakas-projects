pro decals_dr1_gallery, build_sample=build_sample, get_cutouts=get_cutouts, $
  make_png=make_png
; jm15mar19siena - build some pretty pictures for DR1

    dr = 'dr1'

    dr1dir = decals_path(dr=dr)
;   outdir = getenv('HOME')+'/dr1-color/'
;   if file_test(outdir,/dir) eq 0 then file_mkdir, outdir
;   if file_test(outdir+'fits',/dir) eq 0 then file_mkdir, outdir+'fits'
    
    band = ['g','r','z']
    nband = n_elements(band)

    pixscale = 0.262D           ; [arcsec/pixel]
    nmosaic = 3600
    
    width = 1.5D                ; cutout diameter [arcmin]
    link = 3D                   ; linking length [arcmin]

    bricks = mrdfits(dr1dir+'decals-bricks-in-dr1-grz.fits',1)

; --------------------------------------------------
; build a sample of objects for the DR1 gallery
    sample_template = {object: '', ra: 0D, dec: 0D, $
      radius: 0.0, brickname: ''}
    if keyword_set(build_sample) then begin

;; RC3
;       cat = read_rc3()
;       cat = cat[where(cat.d25_maj gt 1.0)]
;
;       spherematch, cat.ra, cat.dec, bricks.ra, $
;         bricks.dec, 0.25/2, icat, ibrick, maxmatch=0
;       icat = icat[uniq(icat,sort(icat))]
;       sample1 = replicate(sample_template,n_elements(icat))
;
;       sample1.object = strcompress(cat[icat].name,/remove)
;       need = where(strcompress(sample1.object,/remove) eq '')
;       if need[0] ne -1 then sample1[need].object = $
;         strcompress(cat[icat[need]].altname,/remove)
;       need = where(strcompress(sample1.object,/remove) eq '')
;       if need[0] ne -1 then sample1[need].object = $
;         strcompress(cat[icat[need]].pgc,/remove)
;
;       sample1.ra = cat[icat].ra
;       sample1.dec = cat[icat].dec
;       sample1.radius = cat[icat].d25_maj*3
;
;       sample = sample1

; NGC objects
       cat = mrdfits(getenv('CATALOGS_DIR')+'/ngc/ngc2000.fits',1)

       spherematch, cat.ra, cat.dec, bricks.ra, $
         bricks.dec, 0.25/2, icat, ibrick, maxmatch=0
       icat = icat[uniq(icat,sort(icat))]
       sample1 = replicate(sample_template,n_elements(icat))

       sample1.object = 'NGC'+strtrim(cat[icat].ngcnum,2)
       sample1.ra = cat[icat].ra
       sample1.dec = cat[icat].dec
       sample1.radius = (cat[icat].radius*60)>1.5

       sample = sample1
;      sample = [sample,sample1]

       djs_plot, bricks.ra, bricks.dec, psym=6, xsty=3, ysty=3
       djs_oplot, sample.ra, sample.dec, psym=6, color='orange'
       djs_oplot, sample1.ra, sample1.dec, psym=6, color='blue'

;; Abell clusters
;       cat = mrdfits(getenv('IM_GITREPOS')+'/astrometry.net/catalogs/abell-all.fits',1)
;
;       spherematch, cat.ra, cat.dec, bricks.ra, $
;         bricks.dec, 0.25/2, icat, ibrick, maxmatch=0
;       icat = icat[uniq(icat,sort(icat))]
;       sample1 = replicate(sample_template,n_elements(icat))
;
;       sample1.object = 'NGC'+strtrim(cat[icat].ngcnum,2)
;       sample1.ra = cat[icat].ra
;       sample1.dec = cat[icat].dec
;       sample1.radius = (cat[icat].radius*60)>1.0
;       sample = [sample,sample1]

; group the matches on 1 arcminute scales and get cutouts centered on
; the center of the group (to take care of multiple close systems)
       grp = spheregroup(sample.ra,sample.dec,1.5D/60D)
       ngrp = histogram(grp,reverse_indices=ri)
       ivec = ri[0:n_elements(ngrp)-1]
       
       oneobj = where(ngrp eq 1) ; only one object - keep it!
       out_sample = sample[ri[ivec[oneobj]]]
    
; more than 1 object
       multobj = where(ngrp gt 1,nmult)
       for ii = 0L, nmult-1 do begin
          match = where(grp eq multobj[ii],nmatch)
          struct_print, sample[match]
          sample1 = sample_template
          sample1.object = strjoin(sample[match].object,'/')
          sample1.ra = djs_mean(sample[match].ra)
          sample1.dec = djs_mean(sample[match].dec)
          sample1.radius = 3.0
          out_sample = [out_sample,sample1]
       endfor          
       nobj = n_elements(out_sample)

; find the matching brick
       for ii = 0, nobj-1 do begin
          this = where(out_sample[ii].ra gt bricks.ra1 and $
            out_sample[ii].ra lt bricks.ra2 and $
            out_sample[ii].dec gt bricks.dec1 and $
            out_sample[ii].dec lt bricks.dec2,nthis)
          if nthis ne 1 then message, 'Multiple matches!'
          out_sample[ii].brickname = bricks[this].brickname
       endfor
       struct_print, out_sample[sort(out_sample.ra)]

       im_mwrfits, out_sample, dr1dir+'gallery_sample.fits', /clobber
    endif 

; --------------------------------------------------
; find matching objects
    if keyword_set(get_cutouts) then begin
       sample = mrdfits(dr1dir+'gallery_sample.fits.gz',1)
       nobj = n_elements(sample)

; loop on each object (and demand 3-color photometry)
       for jj = 0L, nobj-1 do begin
          thisbrick = sample[jj].brickname
          splog, sample[jj].object, thisbrick

          subdir = strmid(thisbrick,0,3)
          tractordir = dr1dir+'tractor/'+subdir+'/'
          coadddir = dr1dir+'coadd/'+subdir+'/'+thisbrick+'/'

          imfile = file_search(coadddir+'decals-'+thisbrick+'-image-'+band+'.fits',count=nfile)
          if nfile ne nband then message, 'Missing some images!'

; save time by reading the images in all three bands
          image = fltarr(nmosaic,nmosaic,nband)
          for ib = 0, nband-1 do image[*,*,ib] = mrdfits(imfile[ib],0,hdr)
             
          racen = djs_mean(sample[jj].ra)
          deccen = djs_mean(sample[jj].dec)
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
       
; --------------------------------------------------
; make the three-color images
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


;; --------------------------------------------------
;; temporary hack to make the dr1-bricks.fits file
;    if keyword_set(make_dr1bricks) then begin
;       all = mrdfits(dr1dir+'decals-bricks.fits',1)
;       allbrick = file_basename(file_search(dr1dir+'tractor/*',/test_dir,count=nbrick))
;       for ii = 0L, nbrick-1 do begin
;          print, format='("Brick ",I0,"/",I0, A10,$)', ii+1, nbrick, string(13b)
;          catfile = file_search(dr1dir+'tractor/'+allbrick[ii]+'/tractor-*.fits')
;          brickname1 = repstr(repstr(file_basename(catfile),'tractor-',''),'.fits','')
;          if n_elements(brickname) eq 0L then brickname = brickname1 else $
;            brickname = [brickname,brickname1]
;; check for bricks with 3-band coverage
;          delvarx, cat
;          for ic = 0, n_elements(catfile)-1 do begin
;             cat1 = mrdfits(catfile[ic],1,/silent,rows=0)
;             if n_elements(cat) eq 0L then cat = cat1 else cat = [cat,cat1]
;          endfor
;          grz = where(total(cat.decam_nobs gt 0,1) ge 3)
;          if grz[0] ne -1 then begin
;             brickname1_grz = brickname1[grz]
;             if n_elements(brickname_grz) eq 0L then brickname_grz = brickname1_grz else $
;               brickname_grz = [brickname_grz,brickname1_grz]
;          endif
;       endfor
;; write out the list of grz dr1 bricks       
;       match, all.brickname, brickname_grz, these_grz
;       mwrfits, all[these_grz], dr1dir+'dr1-grz-bricks.fits', /create
;; write out the list of all dr1 bricks       
;       match, all.brickname, brickname, these
;       mwrfits, all[these], dr1dir+'dr1-bricks.fits', /create
;    endif
;
