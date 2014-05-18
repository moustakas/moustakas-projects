pro get_colorcutout, image, outfile=outfile, cluster=cluster

    ps_start, outfile
    cgdisplay, 1000, 1000
    cgimage, image, /keep_aspect
    if n_elements(cluster) ne 0 then cgtext, 0.08, 0.85, cluster, $
      /normal, color='white', charsize=4.0, charthick=2.0
    ps_end, /png
    
;    delvarx, pos
;    set_plot, 'Z'
;    cgimage, image, /keep_aspect, position=pos;, /norm
;;   plotimage, image, /noaxes, /preserve_aspect, position=pos, /norm
;    if n_elements(cluster) ne 0 then begin
;       im_legend, strupcase(cluster), /left, /top, box=0, charsize=1.5, $
;         charthick=1.0, margin=0, position=[pos[0]+0.0,pos[3]-0.05], /norm, $
;         textcolor=cgcolor('white',0)
;    endif
;    x0 = fix((convert_coord(pos[0:1],/normal,/to_device))[0])
;    nx = fix((convert_coord(pos[2:3],/normal,/to_device))[0])-X0
;    y0 = fix((convert_coord(pos[0:1],/normal,/to_device))[1])
;    ny = fix((convert_coord(pos[2:3],/normal,/to_device))[1])-y0
;
;    img = tvrd(x0,y0,nx,ny)
;    tvlct, rr, gg, bb, /get
;    write_png, outfile, img, rr, gg, bb
;    set_plot, 'X'

return
end

pro plotbcgmstar_montage, cfirst, clast, dobcg=dobcg, donobcg=donobcg, $
  doimage=doimage, onecluster_montage=onecluster_montage, final_montage=final_montage, $
  doitall=doitall
; jm13nov01siena - build a color montage of all the BCGs 

; note that /DOIMAGE should be run first so that we can build a good
; 'levels.txt' file!

;   echo "plotbcgmstar_montage, /doitall" | idl > & plotbcgmstar_montage.log &
    
    if keyword_set(doitall) then begin
       doimage = 1
       dobcg = 1
       donobcg = 1
;      onecluster_montage = 1
    endif
    
    coloroutpath = bcgmstar_path(/colormosaics)
    paperpath = bcgmstar_path(/paper)

    sample = read_bcgmstar_sample()
    ncl = n_elements(sample)

; choose the red, blue, and green filters
    blue = ['f105w','f110w']
    green = ['f125w','f140w']
    red = 'f160w'

    if n_elements(cfirst) eq 0 then cfirst = 0
    if n_elements(clast) eq 0 then clast = ncl-1
    
; ---------------------------------------------------------------------------
; build the original-image trilogy mosaics
    if keyword_set(doimage) then begin
       for ic = cfirst, clast do begin
          cluster = strtrim(sample[ic].shortname,2)
          outpath = coloroutpath+cluster+'/'
          if file_test(outpath,/dir) eq 0 then file_mkdir, outpath
          pushd, outpath

; clean up previously created parameter files          
          file_delete, outpath+['levels.txt', 'trilogyfilterlog.txt',cluster+'-image_filters.txt'], /quiet
          if cluster eq 'macs1149' then begin ; cleaned mosaics are missing!
             subdir = ''
             suffix = ''
          endif else begin
             subdir = 'clnpix_cosmetic/'
             suffix = '_cln'
          endelse
          
          hstpath = getenv('CLASH_ARCHIVE')+strtrim(sample[ic].dirname,2)+'/HST/'
          imagepath = hstpath+'images/mosaicdrizzle_image_pipeline/scale_65mas/'+subdir
          fitsfile = file_search(imagepath+cluster+'_mosaic_065mas_*_'+red[0]+suffix+'_drz_*.fits*',count=nfits)
          if nfits ne 1 then message, 'File '+fitsfile+' not found!'
          mwrfits, 0, outpath+cluster+'-hdr.fits', headfits(fitsfile), /create, /silent
          
          infile = outpath+cluster+'-image.in'
          openw, lun, infile, /get_lun
          printf, lun, 'B'
          for ii = 0, n_elements(blue)-1 do printf, lun, $
            file_search(imagepath+cluster+'_mosaic_065mas_*_'+blue[ii]+suffix+'_drz_*.fits*')
          printf, lun, ''
          printf, lun, 'G'
          for ii = 0, n_elements(green)-1 do printf, lun, $
            file_search(imagepath+cluster+'_mosaic_065mas_*_'+green[ii]+suffix+'_drz_*.fits*')
          printf, lun, ''
          printf, lun, 'R'
          for ii = 0, n_elements(red)-1 do printf, lun, $
            file_search(imagepath+cluster+'_mosaic_065mas_*_'+red[ii]+suffix+'_drz_*.fits*')
          printf, lun, ''
          printf, lun, 'indir '+imagepath
          printf, lun, 'outname '+outpath+cluster+'-image'
          printf, lun, 'noiselum 0.15'
;         printf, lun, 'scaling  '+outpath+'levels.txt'
;         printf, lun, 'scaling  '+getenv('CLASH_ARCHIVE')+'/color_images/scaling/levels_ir.txt'
          printf, lun, 'show 0'
          printf, lun, 'testfirst 0'
          free_lun, lun

          spawn, 'python '+getenv('IM_RESEARCH_DIR')+'/mybin/trilogy.py '+infile, /sh
          popd
       endfor
    endif

; ---------------------------------------------------------------------------
; build the BCG-only trilogy mosaics; there's a problem with
; how trilogy reads the 
    if keyword_set(dobcg) then begin
       for ic = cfirst, clast do begin
          cluster = strtrim(sample[ic].shortname,2)
          outpath = coloroutpath+cluster+'/'
          if file_test(outpath,/dir) eq 0 then file_mkdir, outpath
          pushd, outpath ; pushd so that we can get the levels.txt file in the right spot
          
          file_delete, outpath+cluster+'-bcg_filters.txt', /quiet

          hstpath = getenv('CLASH_ARCHIVE')+strtrim(sample[ic].dirname,2)+'/HST/'
          bcgpath = hstpath+'galaxy_subtracted_images/marc/'

          infile = outpath+cluster+'-bcg.in'
          openw, lun, infile, /get_lun
          printf, lun, 'B'
          for ii = 0, n_elements(blue)-1 do printf, lun, $
            file_search(bcgpath+cluster+'_mosaic_065mas_*_'+blue[ii]+'_drz_*_BCG.fits.gz')
          printf, lun, ''
          printf, lun, 'G'
          for ii = 0, n_elements(green)-1 do printf, lun, $
            file_search(bcgpath+cluster+'_mosaic_065mas_*_'+green[ii]+'_drz_*_BCG.fits.gz')
          printf, lun, ''
          printf, lun, 'R'
          for ii = 0, n_elements(red)-1 do printf, lun, $
            file_search(bcgpath+cluster+'_mosaic_065mas_*_'+red[ii]+'_drz_*_BCG.fits.gz')
          printf, lun, ''
          printf, lun, 'indir '+bcgpath
          printf, lun, 'outname '+outpath+cluster+'-bcg'
          printf, lun, 'noiselum 0.15'
          printf, lun, 'scaling  '+outpath+'levels.txt'
;         printf, lun, 'scaling  '+getenv('CLASH_ARCHIVE')+'/color_images/scaling/levels_ir.txt'
          printf, lun, 'show 0'
          printf, lun, 'testfirst 0'
          free_lun, lun

          spawn, 'python '+getenv('IM_RESEARCH_DIR')+'/mybin/trilogy.py '+infile, /sh
          popd
       endfor
    endif

; ---------------------------------------------------------------------------
; build the BCG-subtracted trilogy mosaics
    if keyword_set(donobcg) then begin
       tmppath = getenv('HOME')+'/tmp/'
       for ic = cfirst, clast do begin
          cluster = strtrim(sample[ic].shortname,2)
          outpath = coloroutpath+cluster+'/'
          if file_test(outpath,/dir) eq 0 then file_mkdir, outpath
          pushd, outpath

          file_delete, outpath+cluster+'-nobcg_filters.txt', /quiet
          if cluster eq 'macs1149' then begin ; cleaned mosaics are missing!
             subdir = ''
             suffix = ''
          endif else begin
             subdir = 'clnpix_cosmetic/'
             suffix = '_cln'
          endelse

          hstpath = getenv('CLASH_ARCHIVE')+strtrim(sample[ic].dirname,2)+'/HST/'
          imagepath = hstpath+'images/mosaicdrizzle_image_pipeline/scale_65mas/'+subdir
          bcgpath = hstpath+'galaxy_subtracted_images/marc/'

; need to do the subtraction myself
          allfilt = [blue,green,red]
          for ii = 0, n_elements(allfilt)-1 do begin
             imagefile = file_search(imagepath+cluster+'_mosaic_065mas_*_'+allfilt[ii]+suffix+'_drz_*.fits*')
             bcgfile = file_search(bcgpath+cluster+'_mosaic_065mas_*_'+allfilt[ii]+'_drz_*_BCG.fits.gz')
             
             data = mrdfits(imagefile,0,hdr,/silent)
             bcg = mrdfits(bcgfile,0,/silent)

             im_mwrfits, data-bcg, tmppath+cluster+'-'+allfilt[ii]+'.fits', hdr, /silent, /clobber
          endfor 

          infile = outpath+cluster+'-nobcg.in'
          openw, lun, infile, /get_lun
          printf, lun, 'B'
          for ii = 0, n_elements(blue)-1 do printf, lun, tmppath+cluster+'-'+blue[ii]+'.fits.gz'
;         for ii = 0, n_elements(blue)-1 do printf, lun, $
;           file_search(bcgpath+cluster+'_mosaic_065mas_*_'+blue[ii]+'_drz_*_BCG.fits.gz')
          printf, lun, ''
          printf, lun, 'G'
          for ii = 0, n_elements(green)-1 do printf, lun, tmppath+cluster+'-'+green[ii]+'.fits.gz'
;         for ii = 0, n_elements(green)-1 do printf, lun, $
;           file_search(bcgpath+cluster+'_mosaic_065mas_*_'+green[ii]+'_drz_*_BCG.fits.gz')
          printf, lun, ''
          printf, lun, 'R'
          for ii = 0, n_elements(red)-1 do printf, lun, tmppath+cluster+'-'+red[ii]+'.fits.gz'
;         for ii = 0, n_elements(red)-1 do printf, lun, $
;           file_search(bcgpath+cluster+'_mosaic_065mas_*_'+red[ii]+'_drz_*_BCG.fits.gz')
          printf, lun, ''
          printf, lun, 'indir '+bcgpath
          printf, lun, 'outname '+outpath+cluster+'-nobcg'
          printf, lun, 'noiselum 0.15'
          printf, lun, 'scaling  '+outpath+'levels.txt'
;         printf, lun, 'scaling  '+getenv('CLASH_ARCHIVE')+'/color_images/scaling/levels_ir.txt'
          printf, lun, 'show 0'
          printf, lun, 'testfirst 0'
          free_lun, lun

          spawn, 'python '+getenv('IM_RESEARCH_DIR')+'/mybin/trilogy.py '+infile, /sh
          popd
          
; clean up          
          for ii = 0, n_elements(allfilt)-1 do rmfile, tmppath+cluster+'-'+allfilt[ii]+'.fits.gz'
          
       endfor
    endif

; ---------------------------------------------------------------------------
; get color cutouts of the original image, BCG model, and
; BCG-subtracted images
    if keyword_set(onecluster_montage) then begin
       dd = 300D                ; extraction diameter [kpc]
       
       for ic = cfirst, clast do begin
          cluster = strtrim(sample[ic].shortname,2)
          title = repstr(repstr(strupcase(sample[ic].dirname),'_',' '),'ABELL','Abell')
          arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
       
; extract the center of the mosaic centered on the BCG
          extast, headfits(coloroutpath+cluster+'/'+cluster+'-hdr.fits'), astr
          pixscale = sqrt(abs(determ(astr.cd)))*3600D ; [arcsec/pixel]

          pixscale = sqrt(abs(determ(astr.cd)))*3600D ; [arcsec/pixel]
          ad2xy, 15D*hms2dec(sample[ic].ra), hms2dec(sample[ic].dec), astr, xcen, ycen
          xcen = long(xcen) & ycen = long(ycen)
          
          diam_deg = dd/arcsec2kpc/3600D      ; [DD in degrees]
          diam = long(dd/arcsec2kpc/pixscale) ; [DD in pixels]
          
          x0 = xcen-diam/2
          x1 = xcen+diam/2
          y0 = ycen-diam/2
          y1 = ycen+diam/2

          im = read_png(coloroutpath+cluster+'/'+cluster+'-image.png')
          bcg = read_png(coloroutpath+cluster+'/'+cluster+'-bcg.png')
          nobcg = read_png(coloroutpath+cluster+'/'+cluster+'-nobcg.png')

          get_colorcutout, im[*,x0:x1,y0:y1], outfile=coloroutpath+cluster+'/'+cluster+'-image-cutout.png', cluster=title
          get_colorcutout, bcg[*,x0:x1,y0:y1], outfile=coloroutpath+cluster+'/'+cluster+'-bcg-cutout.png'
          get_colorcutout, nobcg[*,x0:x1,y0:y1], outfile=coloroutpath+cluster+'/'+cluster+'-nobcg-cutout.png'

          outfile = coloroutpath+cluster+'-montage.png'
          cmd = 'montage -bordercolor black -borderwidth 1 '+ $
            '-tile 3x1 -geometry +0+0 -quality 100 '+$ ; -resize 200x200 '+$
            strjoin(coloroutpath+cluster+'/'+cluster+'-'+['image','bcg','nobcg']+'-cutout.png',' ')+' '+outfile
          splog, cmd
          spawn, cmd, /sh
       endfor
    endif

; ---------------------------------------------------------------------------
; build the final montage!  write a full-resolution 
    if keyword_set(final_montage) then begin
       outfile = bcgmstar_path()+'bcg-finalmontage.png'
       infile = strjoin(file_search(coloroutpath+'*-montage.png'),' ')
       
       cmd = 'montage -bordercolor white -borderwidth 1 '+ $
         '-tile 3x5 -geometry +0+0 -quality 100 '+$ ; -resize 1024x1024 '+$
         infile+' '+outfile
       splog, cmd
       spawn, cmd, /sh

       spawn, 'convert -scale 50% '+outfile+' '+paperpath+$
         'bcg-finalmontage-rescaled.png', /sh
    endif
    
return
end
    

