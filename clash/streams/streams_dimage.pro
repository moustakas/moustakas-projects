pro streams_dimage, build_links=build_links, parents=parents, $
  clean_parents=clean_parents, stargal=stargal, children=children, $
  measure=measure, rebuild_mosaic=rebuild_mosaic, clobber=clobber, $
  debug=debug
; jm13sep09siena - process the STREAMS clusters through dimage 

; note! images in units of [10^-12 erg/s/cm^2/Hz] (pico-maggies)

    sample = rsex(streams_path(/propath)+'streams_sample.sex')
    struct_print, sample 
    ncl = n_elements(sample) 

    pixscale = 0.065D ; [arcsec/pixel]

; wrap on each cluster    
    for ic = 6, 6 do begin
;   for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       splog, 'Working on cluster '+cluster

       outpath = streams_path(/dimage)+cluster+'/'
       if file_test(outpath,/dir) eq 0 then file_mkdir, outpath
          
; read the skyinfo structure to get the filters
       skyinfo = mrdfits(streams_path(/skysub)+'skyinfo-'+$
         cluster+'.fits.gz',1,/silent)
       short = strtrim(skyinfo.band,2)
       reffilt = (where(short eq 'f160w'))[0] ; reference filter
       nfilt = n_elements(skyinfo)

       psffiles = file_search(streams_path(/psfs)+short+'-bpsf.fits*')

; --------------------------------------------------
; build links to the mosaics (either original sky-subtracted, or sky-
; and BCG-subtracted)
       if keyword_set(build_links) then begin
          for ii = 0, nfilt-1 do spawn, 'ln -sf ../../skysub/'+cluster+'/'+$
            cluster+'-'+short[ii]+'.fits.gz '+outpath+cluster+'-'+short[ii]+'.fits.gz'
;         for ii = 0, nfilt-1 do spawn, 'ln -sf ../../nobcg/'+cluster+'/'+$
;           cluster+'-'+short[ii]+'.fits.gz '+outpath+cluster+'-'+short[ii]+'.fits.gz'
       endif

; --------------------------------------------------
; find parents; creates pcat, pimage, and parents files 
       if keyword_set(parents) then begin
          splog, 'Finding parents'
          pushd, outpath

; clean up first?
          if keyword_set(clobber) eq 0 then begin
             splog, 'Delete all files from a previous call to STREAMS_DIMAGE, /PARENTS [Y/N]?'
             cc = get_kbrd()
             if strupcase(cc) eq 'Y' then begin
                if file_test('parents',/dir) then file_delete, 'parents/', /recursive, /quiet
                delfiles = file_search([cluster+'-pcat*.fits*',cluster+'-*-pimage.fits*',$
                  cluster+'-pset.fits*'],count=ndel)
                if ndel ne 0 then file_delete, delfiles, /quiet
             endif
          endif

          imfiles = cluster+'-'+short+'.fits.gz'
          pset = {$
            nlevel:             1L,$
            plim:             25.0,$ ; detection significance level
            pbuffer:           0.2,$ ; 10% buffer
            base:          cluster,$
            imfiles:       imfiles,$
            ref:           reffilt,$
            puse:    intarr(nfilt),$
            dopsf:  intarr(nfilt)+1}
          pset.puse[reffilt] = 1 ; use just F160W for detection

          im_mwrfits, pset, outpath+cluster+'-pset.fits', /clobber

          dparents, cluster, outpath+imfiles, ref=pset.ref, plim=pset.plim, $
            pbuffer=pset.pbuffer, nlevel=pset.nlevel, puse=pset.puse
          heap_gc
          popd
       endif

;; --------------------------------------------------
;; clean up parents by tossing out objects that have more than 1%
;; masked pixels (edges!)
;       if keyword_set(clean_parents) then begin
;          splog, 'Cleaning up parents'
;          cat = mrdfits(outpath+base+'-pcat.fits',1)
;          ngal = n_elements(cat)
;          fracmask = fltarr(ngal)
;          for ii = 0L, ngal-1 do begin
;             im = mrdfits(outpath+'parents/'+base+'-parent-'+$
;               strtrim(ii,2)+'.fits',reffilt*2-1,/silent)
;             fracmask[ii] = total(im eq 0)/cmproduct(size(im,/dim))
;          endfor
;
;          keep = where(fracmask lt 0.01,nkeep,comp=toss)
;          splog, 'Keeping '+strtrim(nkeep,2)+'/'+strtrim(ngal,2)+' parents'
;          djs_plot, cat.xst[0]+cat.xc[0], cat.yst[0]+cat.yc[0], psym=6, title=cluster
;          djs_oplot, cat[toss].xst[0]+cat[toss].xc[0], $
;            cat[toss].yst[0]+cat[toss].yc[0], psym=6, color='orange'
;          
;          file_copy, outpath+base+'-pcat.fits', outpath+base+'-pcat-original.fits', /overwrite
;          mwrfits, cat[keep], outpath+base+'-pcat.fits', /create
;
;; rebuild the PARENTS directory
;          file_move, outpath+'parents', outpath+'oldparents', /overwrite
;          file_mkdir, outpath+'parents'
;          for ii = 0, nkeep-1 do file_copy, outpath+'oldparents/'+$
;            base+'-parent-'+strtrim(keep[ii],2)+'.fits', $
;            outpath+'parents/'+base+'-parent-'+strtrim(ii,2)+'.fits'
;          
;; test code
;;         pimage = mrdfits(outpath+base+'-pimage.fits',0)
;;         nstamp = lonarr(ngal)
;;         for ii = 0L, ngal-1 do nstamp[ii] = total(pimage eq ii)
;       endif
       
; --------------------------------------------------
; classify objects into stars and galaxies; also subtract foreground
; PSF stars
       if keyword_set(stargal) then begin
          splog, 'Performing star-galaxy separation.'
          pushd, outpath

; clean up first?
          if keyword_set(clobber) eq 0 then begin
             splog, 'Delete all files from a previous call to STREAMS_DIMAGE, /STARGAL [Y/N]?'
             cc = get_kbrd()
             if strupcase(cc) eq 'Y' then begin
                if file_test('atlases',/dir) then file_delete, 'atlases/', /recursive, /quiet
             endif
          endif
          
          glim = 15.0   ; minimum significance above the noise
          gsmooth = 3.0 ; 3.0 ; default is 10, which is too aggressive
          gsaddle = 3.0 ; 3.0
          maxnstar = 50L ; 100L

;         glim = 20.0           ; minimum significance above the noise
;         gsmooth = 10.0 ; 3.0 ; default is 10, which is too aggressive
;         gsaddle = 10.0 ; 3.0
;         maxnstar = 50L ; 100L
          
          streams_stargal, cluster, plot=plot, gsmooth=gsmooth, glim=glim, $
            gsaddle=gsaddle, nsigma=nsigma, noclobber=noclobber, $
            maxnstar=maxnstar, psffiles=psffiles
          heap_gc
          popd
          
; make a QAplot?!?          
          
       endif

; --------------------------------------------------
; deblend into children
       if keyword_set(children) then begin
          splog, 'Deblending parents into children'
          pushd, outpath
          streams_children, cluster, noclobber=noclobber, debug=debug, /sersic
          heap_gc
          popd
       endif

; --------------------------------------------------
; model each object
       if keyword_set(measure) then begin
          splog, 'Modeling!'
          pushd, outpath

; need to deal with the fact that streams_measure uses the same sky
; value from the reference image for all the other images

          
          
          
; test code          
;         cutimage = mrdfits('atlases/303/macs1206-303-templates-0.fits',0,hdr)
;         mge1 = streams_mge(cutimage,badpixels=badpixels,$
;           pixscale=pixscale,model=cutmodel,twist=twist)
;         mge1 = streams_mge(cutimage,badpixels=badpixels,pixscale=pixscale,$
;           model=cutmodel,/twist,minlevel=3*ss[12].sigma)
;         help, mge1, /str

;         bulgedisk = 1 ; do two-component modeling
          streams_measure, cluster, psffiles=psffiles, noclobber=noclobber, $
            bulgedisk=bulgedisk, mgefit=mgefit
          heap_gc
          popd
       endif

; --------------------------------------------------
; rebuild the mosaic with the galaxies subtracted
       if keyword_set(rebuild_mosaic) then begin
          splog, 'Rebuilding the mosaic!'

          sub = 'atlases'
          postfix = ''
          
          for ib = nfilt-1, nfilt-1 do begin
;         for ib = 0, nfilt-1 do begin
             splog, 'Working on filter '+short[ib]

             imfile = outpath+cluster+'-'+short[ib]+'.fits.gz'
             outfile = outpath+cluster+'-'+short[ib]+'-model.fits'
             image = gz_mrdfits(imfile,0,hdr,/silent)
             invvar = gz_mrdfits(imfile,1,ivarhdr,/silent)
             model = image*0
             
             pcat = gz_mrdfits(outpath+cluster+'-pcat.fits',1)
;            for iparent = 266, 275 do begin
             for iparent = 0L, n_elements(pcat)-1 do begin
                splog, 'Parent ', iparent
    
                pstr = strtrim(string(iparent),2)
                mfile = outpath+'/'+sub+'/'+pstr+'/'+cluster+'-'+pstr+ $
                  '-measure'+postfix+'.fits.gz'
                sfile = outpath+'/'+sub+'/'+pstr+'/'+cluster+'-'+pstr+ $
                  '-sersic'+postfix+'.fits.gz'
                if file_test(sfile) then begin
;                  if file_test(mfile) eq 0 then stop
;                  mall = mrdfits(mfile,1,/silent)
                   model1 = gz_mrdfits(sfile,ib,hdr1,/silent)
;                  model1 = model1-mall.sky[ib] ; subtract the sky
                   embed_stamp, model, model1, pcat[iparent].xst[ib], pcat[iparent].yst[ib]
                endif
             endfor
             im_mwrfits, model, outfile, hdr, /nogzip, /clobber
             im_mwrfits, image-model, outfile, hdr, /gzip, /append
          endfor 
       endif 
    endfor ; close cluster loop

return
end
    
