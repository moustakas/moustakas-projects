pro streams_bcg_postman, debug=debug
; jm13sep04siena - do the preliminary BCG subtraction

; note! images in units of [10^-12 erg/s/cm^2/Hz] (pico-maggies)

    qapath = streams_path(/postman_bcg)+'qaplots/'
    bcgmodelpath = '/Users/ioannis/archive/bcg_models/' ; FIX THIS!!!
    
; read the sample
    sample = rsex(streams_path(/propath)+'streams_sample.sex')
    splog, 'IGNORING A2261!!!'
    keep = where(strtrim(sample.shortname,2) ne 'a2261')
    sample = sample[keep]
;   struct_print, sample
    ncl = n_elements(sample)

    pixscale = 0.065D ; [arcsec/pixel]
    rmaxkpc = 200D     ; [kpc]

; wrap on each cluster    
;   for ic = 7, 7 do begin
    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       splog, 'Working on cluster '+cluster
       skypath = streams_path(/skysub)+cluster+'/'
       outpath = streams_path(/postman_bcg)+cluster+'/'
       if file_test(outpath,/dir) eq 0 then file_mkdir, outpath

       bcgqafile = qapath+cluster+'_bcg.ps'

; get a fixed RMAXKPC cutout centered on the BCG
       arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
       rmax = ceil(rmaxkpc/arcsec2kpc/pixscale)         ; [pixels]

; read the skyinfo structure and figure out which bands have
; Marc's BCG model photometry
       skyinfo = mrdfits(streams_path(/skysub)+'skyinfo-'+$
         cluster+'.fits.gz',1,/silent)
       short = strtrim(skyinfo.band,2)
       these = where(file_test(bcgmodelpath+cluster+'_mosaic_065mas_*_'+$
         short+'_drz_*_BCG.fits.gz'),nfilt,comp=missing)
       if missing[0] ne -1 then splog, 'Missing BCG models in '+$
         strupcase(cluster)+': '+strjoin(strupcase(short[missing]),', ')

       outinfo = struct_addtags(skyinfo[these],replicate({ra: 0D, dec: 0D, $
         ellipticity: 0.0, posangle: 0.0, size: 0.0, size_kpc: 0.0},nfilt))
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
       outinfo[reffilt].size = size*pixscale
       outinfo[reffilt].size_kpc = size*pixscale*arcsec2kpc
       outinfo[reffilt].posangle = posangle
       outinfo[reffilt].ellipticity = ellipticity
       help, outinfo[reffilt], /st

; build a 4-panel QAplot in all the bands
       im_plotconfig, 5, pos, psfile=bcgqafile, xspace=0.02, $
         yspace=0.02, xmargin=[0.4,0.4]

; now loop through each band, cut out the BCG, measure its properties,
; and write out
;      for ib = nfilt-1, nfilt-3, -1 do begin
       for ib = nfilt-1, 0, -1 do begin

; read the full-cluster images and isolate the BCG          
          imfile = skypath+cluster+'-'+short[ib]+'.fits.gz'
          splog, 'Reading '+imfile
          image = gz_mrdfits(imfile,0,hdr,/silent)
          invvar = gz_mrdfits(imfile,1,ivarhdr,/silent)
          sz = size(image,/dim)

          adxy, hdr, ra, dec, xcen, ycen
          x0 = (xcen-rmax)>0
          y0 = (ycen-rmax)>0
          x1 = (xcen+rmax)<(sz[0]-1)
          y1 = (ycen+rmax)<(sz[1]-1)

          hextract, image, hdr, cutimage, cuthdr, x0, x1, y0, y1, /silent
          hextract, invvar, ivarhdr, cutinvvar, cutivarhdr, x0, x1, y0, y1, /silent

; read the BCG model, sky-subtract and scale to picomaggies, and
; finally crop to the same region; need to check where the model
; becomes unreliable!
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
          outinfo[ib].size = size*pixscale
          outinfo[ib].size_kpc = size*pixscale*arcsec2kpc
          outinfo[ib].posangle = posangle
          outinfo[ib].ellipticity = ellipticity
          help, outinfo[ib], /str

; create an object mask so that we can do photometry on the data
; itself in STREAMS_BCG_APPHOT; a smoothing scale larger than 3 is
; usually too aggressive in the core of the BCG and increasing PLIM
; leaves too many faint sources undetected
          domask = 0 ; come back to the mask
          if domask then begin
             simage = dsmooth(cutimage-model,3L)
             dobjects, simage, objects=obj, plim=10, nlevel=1L
             mask = obj eq -1 ; 0 = masked, 1 = unmasked
          endif else mask = cutimage*0

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
          mwrfits, model, outfile, cuthdr, /silent
          mwrfits, cutinvvar, outfile, cutivahdr, /silent
          mwrfits, mask, outfile, cuthdr, /silent
          spawn, 'gzip -f '+outfile
       
          cgloadct, 0, /silent
          cgimage, cutimage, clip=3, /negative, stretch=5, minvalue=0.0, maxvalue=1.0, $
            margin=0, /keep_aspect, position=pos[*,0], /save
          xyouts, pos[0,0]+0.02, pos[3,0]-0.05, strupcase(cluster)+'/'+$
            strupcase(short[ib]), /norm, charsize=1.5
          cgimage, cutinvvar, clip=3, /negative, /noerase, $
            stretch=5, margin=0, /keep_aspect, position=pos[*,1]
          cgimage, model, clip=3, /negative, /noerase, minvalue=0.0, $
            maxvalue=1.0, stretch=5, margin=0, /keep_aspect, position=pos[*,2]
          cgimage, cutimage-model, clip=3, /negative, /noerase, minvalue=0.0, $
            maxvalue=1.0, stretch=5, margin=0, /keep_aspect, position=pos[*,3]
       endfor                   ; close filter loop
       im_plotconfig, /psclose, psfile=bcgqafile, /pdf

; write out the info file       
       outfile = outpath+cluster+'-mgeskyinfo.fits'
       im_mwrfits, outinfo, outfile, /clobber
    endfor                      ; close cluster loop

return
end
    
