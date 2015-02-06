pro bcgmstar_get_bcg, debug=debug
; jm13sep04siena - cut out the BCG

; note! images in units of [10^-12 erg/s/cm^2/Hz] (pico-maggies)
    getbcgpath = bcgmstar_path(/bcg)
    skypath = bcgmstar_path(/skysub)
    skyinfopath = bcgmstar_path()+'skyinfo/'
    qapath = bcgmstar_path()+'qaplots-getbcg/'

;   splog, 'Update BCGMODELPATH! to the archive!'
;   bcgmodelpath = '/Users/ioannis/archive/bcg_models/'

; read the sample
    sample = read_bcgmstar_sample()
    ncl = n_elements(sample)

    pixscale = 0.065D ; [arcsec/pixel]
    rmaxkpc = 200D     ; [kpc]

; wrap on each cluster    
;   for ic = 11, 11 do begin
    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       splog, 'Working on cluster '+cluster
       outpath = getbcgpath+cluster+'/'
       if file_test(outpath,/dir) eq 0 then file_mkdir, outpath

       bcgmodelpath = getenv('CLASH_ARCHIVE')+'/'+strtrim(sample[ic].dirname,2)+$
         '/HST/galaxy_subtracted_images/marc/'
       bcgqafile = qapath+'qa_getbcg_'+cluster+'.ps'

; get a fixed RMAXKPC cutout centered on the BCG
       ebv = sample[ic].ebv
       arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
       rmax = ceil(rmaxkpc/arcsec2kpc/pixscale)         ; [pixels]

; read the skyinfo structure and figure out which bands have
; Marc's BCG model photometry
       skyinfo = mrdfits(skyinfopath+'skyinfo-'+cluster+'.fits.gz',1,/silent)
       short = strtrim(skyinfo.band,2)
       these = where(file_test(bcgmodelpath+cluster+'_mosaic_065mas_*_'+$
         short+'_drz_*_BCG.fits.gz'),nfilt,comp=missing)
       if missing[0] ne -1 then begin
          misslist1 = strupcase(cluster)+': '+strjoin(strupcase(short[missing]),', ')
          if n_elements(misslist) eq 0 then misslist = misslist1 else $
            misslist = [misslist,misslist1]
          splog, 'Missing BCG models in '+strupcase(cluster)+': '+$
            strjoin(strupcase(short[missing]),', ')
       endif

; the F435W model for MACS1149 is not reliable (the galaxy is
; basically invisible and the model fit is like a point source); get
; rid of it here 
       if cluster eq 'macs1149' then begin
          these = these[where(short[these] ne 'f435w',nfilt)]
          splog, 'Unreliable model: MACS1149/F435W'
       endif

; temporarily drop the F390W image of RXJ2129       
       if cluster eq 'rxj2129' then begin
          these = these[where(short[these] ne 'f390w',nfilt)]
          splog, 'Hack!!! Dropping RXJ2129/F390W!!'
       endif

       outinfo = struct_addtags(skyinfo[these],replicate({ra: 0D, dec: 0D, $
         mge_ellipticity: 0.0, mge_posangle: 0.0, mge_size: 0.0, mge_size_kpc: 0.0},nfilt))
       weff = outinfo.weff
       short = strtrim(outinfo.band,2)
       reffilt = where(short eq 'f160w') ; reference filter

; read the model in the reference band to find the BCG and get its
; basic properties
       bcgmodelfile = file_search(bcgmodelpath+cluster+'_mosaic_065mas_*_'+$
         short[reffilt]+'_drz_*_BCG.fits.gz',count=nn)
       if nn eq 0 then message, 'This should not happen!'
       splog, 'Reading '+bcgmodelfile

       model = mrdfits(bcgmodelfile,0,hdr,/silent)
       find_galaxy, model, size, ellipticity, posangle, xcen, $
         ycen, xcen_lum, ycen_lum, fraction=0.05, /quiet
       xyad, hdr, xcen_lum, ycen_lum, ra, dec
       outinfo[reffilt].ra = ra
       outinfo[reffilt].dec = dec
       outinfo[reffilt].mge_size = size*pixscale
       outinfo[reffilt].mge_size_kpc = size*pixscale*arcsec2kpc
       outinfo[reffilt].mge_posangle = posangle
       outinfo[reffilt].mge_ellipticity = ellipticity
       help, outinfo[reffilt], /st

; build a 4-panel QAplot in all the bands
       im_plotconfig, 5, pos, psfile=bcgqafile, xspace=0.02, $
         yspace=0.02, xmargin=[0.4,0.4]

; now loop through each band, cut out the BCG, measure its properties,
; and write out
;      for ib = nfilt-1, nfilt-2, -1 do begin
       for ib = nfilt-1, 0, -1 do begin

; read the full-cluster images and isolate the BCG          
          imfile = skypath+cluster+'/'+cluster+'-'+short[ib]+'.fits.gz'
          splog, 'Reading '+imfile
          image = gz_mrdfits(imfile,0,hdr,/silent)
          invvar = gz_mrdfits(imfile,1,ivarhdr,/silent)
          mask = gz_mrdfits(imfile,2,maskhdr,/silent)
          sz = size(image,/dim)

          adxy, hdr, ra, dec, xcen, ycen
          x0 = (xcen-rmax)>0
          y0 = (ycen-rmax)>0
          x1 = (xcen+rmax)<(sz[0]-1)
          y1 = (ycen+rmax)<(sz[1]-1)

          hextract, image, hdr, cutimage, cuthdr, x0, x1, y0, y1, /silent
          hextract, invvar, ivarhdr, cutinvvar, cutivarhdr, x0, x1, y0, y1, /silent
          hextract, mask, maskhdr, cutmask, cutmaskhdr, x0, x1, y0, y1, /silent

; read the BCG model, sky-subtract, scale to picomaggies (which
; includes dust correction: see BCGMSTAR_SKYSUBTRACT), and finally crop
; to the same region
          bcgmodelfile = file_search(bcgmodelpath+cluster+$
            '_mosaic_065mas_*_'+short[ib]+'_drz_*_BCG.fits.gz')
          if file_test(bcgmodelfile) eq 0 then message, 'This should not happen!'
          splog, 'Reading '+bcgmodelfile

          fullmodel = mrdfits(bcgmodelfile,0,fullmodelhdr,/silent)
          fullmodel = fullmodel*outinfo[ib].factor - outinfo[ib].mode

          hline = hdr[where(strmatch(hdr,'*Extracted Image*'))]
          hline = strmid(hline,strpos(hline,'['),strpos(hline,']'))
          xxyy = long((strsplit(hline,'[]:,',/extract))[0:3])
          hextract, fullmodel, fullmodelhdr, model1, modelhdr1, $
            xxyy[0], xxyy[1], xxyy[2], xxyy[3], /silent
          if total(size(model1,/dim) eq sz) ne 2 then message, 'Mismatched dimensions!'

          hextract, model1, modelhdr1, model, modelhdr, x0, x1, y0, y1, /silent

; get some info in this band
          find_galaxy, model, size, ellipticity, posangle, xcen, $
            ycen, xcen_lum, ycen_lum, fraction=0.05, /quiet
          xyad, cuthdr, xcen_lum, ycen_lum, ra, dec
          outinfo[ib].ra = ra
          outinfo[ib].dec = dec
          outinfo[ib].mge_size = size*pixscale
          outinfo[ib].mge_size_kpc = size*pixscale*arcsec2kpc
          outinfo[ib].mge_posangle = posangle
          outinfo[ib].mge_ellipticity = ellipticity
          help, outinfo[ib], /str

;; create an object mask so that we can do photometry on the data
;; itself in BCGMSTAR_ELLIPSE; a smoothing scale larger than 3 is
;; usually too aggressive in the core of the BCG and increasing PLIM
;; leaves too many faint sources undetected
;          if cluster eq 'clj1226' then domask = 0 else domask = 1 ; come back to the mask
;          if domask then begin
;             simage = dsmooth(cutimage-model,3L)
;             plim = 30.0
;             delvarx, obj
;             dobjects, simage, objects=obj, plim=plim, nlevel=1L
;             mask = obj eq -1 ; 0 = masked, 1 = unmasked
;          endif else mask = cutimage*0+1 

;; find residual significant peaks and mask those, too, using a fixed
;; 4-pixel rectangular
;          dpeaks, (cutimage-model)*(mask eq 0), xc=xc, yc=yc, /refine, /check
;          cgimage, cutimage, /negative, /keep, minvalue=0.0, maxvalue=5.0, stretch=5, /save
;          cgoplot, xc, yc, color=cgcolor('red'), psym=7

;         cgimage, cutimage*(mask eq 0), /negative, /keep, minvalue=0.0, maxvalue=5.0, stretch=5
          
;; here's some *good* code to smooth & grow the mask, but it may not be
;; needed for most clusters
;          sdev = 1.0
;          for ii = 0, 0 do mask = convol(mask*1.0,psf_gaussian($
;            npix=long(sdev*5),st_dev=sdev,/norm)) gt 0.0

;         cgimage, cutimage*(mask eq 0), /negative, /keep, minvalue=0.0, maxvalue=5.0, stretch=5
;         cgimage, cutimage, /negative, /keep, minvalue=0.0, maxvalue=5.0, stretch=5
;         cgimage, mask, /negative, /keep, minvalue=0, maxvalue=1, stretch=5
;         cgimage, (cutimage-model)*(mask eq 0), /negative, $
;           /keep, minvalue=0.0, maxvalue=5.0, stretch=5

; write out
          outfile = outpath+cluster+'-'+short[ib]+'.fits'
          splog, 'Writing '+outfile

          mwrfits, cutimage, outfile, cuthdr, /create, /silent
          mwrfits, cutinvvar, outfile, cutivahdr, /silent
          mwrfits, cutmask, outfile, cutmaskhdr, /silent
          mwrfits, model, outfile, cuthdr, /silent
          spawn, 'gzip -f '+outfile
          
          cgloadct, 0, /silent
          mx1 = weighted_quantile(cutimage,quant=0.9)
          cgimage, cutimage, clip=3, /negative, stretch=5, minvalue=0.0, maxvalue=mx1, $
            margin=0, /keep_aspect, position=pos[*,0], /save
          xyouts, pos[0,0]+0.02, pos[3,0]-0.05, strupcase(cluster)+'/'+$
            strupcase(short[ib]), /norm, charsize=1.5
          cgimage, cutimage*cutmask, clip=5, /negative, /noerase, $
            stretch=5, margin=0, /keep_aspect, position=pos[*,1], $
            minvalue=0.0, maxvalue=mx1
;         cgimage, cutinvvar, clip=3, /negative, /noerase, $
;           stretch=5, margin=0, /keep_aspect, position=pos[*,1]
          cgimage, model, clip=3, /negative, /noerase, minvalue=0.0, $
            maxvalue=mx1, stretch=5, margin=0, /keep_aspect, position=pos[*,2]
          cgimage, cutimage-model, clip=3, /negative, /noerase, minvalue=0.0, $
            maxvalue=mx1, stretch=5, margin=0, /keep_aspect, position=pos[*,3]
       endfor                   ; close filter loop
       im_plotconfig, /psclose, psfile=bcgqafile, /pdf

; write out the info file       
       outfile = skyinfopath+cluster+'-mgeskyinfo.fits'
       im_mwrfits, outinfo, outfile, /clobber
    endfor                      ; close cluster loop

    if n_elements(misslist) ne 0 then begin
       print
       splog, 'Missing clusters/filters :'
       niceprint, misslist
    endif
    
return
end
    
