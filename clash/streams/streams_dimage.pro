pro streams_dimage, skysubtract=skysubtract, build_psf=build_psf, $
  parents=parents
; jm13may26siena - sky-subtract the CLASH clusters containing tidal
; streams 

    propath = getenv('CLASH_DIR')+'/streams/'
    datapath = getenv('IM_ARCHIVE_DIR')+'/projects/streams/'
;   datapath = getenv('IM_PROJECTS_DIR')+'/clash/streams/'

    filt = bcgimf_filterlist(short=short,instr=instr,weff=weff,zpt=zpt)
    these = [9,12] ; = [F110W, F160W]
    filt = filt[these]
    short = short[these]
    instr = instr[these]
    weff = weff[these]
    zpt = zpt[these]
    nfilt = n_elements(filt)

    reffilt = (where(short eq 'f160w'))[0] ; reference filter
    kl = k_lambda(weff,/odon)

    clash = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    cluster = 'macs1206' ; start with one cluster
    ncl = n_elements(cluster) 

    pixscale = 0.065D ; [arcsec/pixel]
    
; wrap on each cluster    
    for ic = 0, ncl-1 do begin
       splog, 'Working on cluster '+cluster[ic]
       this =  where(cluster[ic] eq strtrim(clash.cluster_short,2))
       ebv = clash[this].ebv

       if file_test(datapath+cluster[ic],/dir) eq 0 then $
         spawn, 'mkdir -p '+datapath+cluster[ic]
       mosaicpath = getenv('IM_ARCHIVE_DIR')+'/'+cluster[ic]+'/images/'
;      mosaicpath = getenv('IM_ARCHIVE_DIR')+'/'+cluster[ic]+'/patched_images/'
       catpath = getenv('IM_ARCHIVE_DIR')+'/'+cluster[ic]+'/catalogs/'

; --------------------------------------------------
; do the sky subtraction
       if keyword_set(skysubtract) then begin
          skyinfo = replicate({file: '', band: '', sblimit: 0.0, $
            mode: 0.0, sigma: 0.0, skew: 0.0, nsky: 0L},nfilt)
          skyinfo.band = short
          
; read the original mosaics and inverse variance maps; crop to the
; maximum range of the reference filter, F160W
          data = im_fits_cube(file_search(mosaicpath+cluster[ic]+$
            '_mosaic_065mas_'+instr+'_'+short+'_drz_????????.fits*'))
          ivardata = im_fits_cube(file_search(mosaicpath+cluster[ic]+$
            '_mosaic_065mas_'+instr+'_'+short+'_wht_????????.fits*'))
          xcrop = ceil(minmax(where(total(data[reffilt].image,1) gt 0))/100.0)*100
          ycrop = ceil(minmax(where(total(data[reffilt].image,2) gt 0))/100.0)*100

; aggresively build an object mask
          dobjects, data.image[xcrop[0]:xcrop[1],ycrop[0]:ycrop[1]], $
            object=oimage, plim=5.0, fobject=fobj, dpsf=dpsf
          sz = size(fobj,/dim)

; sky-subtract each band independently
          skyinfo.file = file_basename(data.fname)
          for ib = 0, nfilt-1 do begin
             splog, 'Sky-subtracting band '+short[ib]
             good = where(ivardata[ib].image[xcrop[0]:xcrop[1],ycrop[0]:ycrop[1]] gt 0 and $
               fobj eq -1,ngood)
             mmm, (data[ib].image[xcrop[0]:xcrop[1],ycrop[0]:ycrop[1]])[good], $
               mode1, sigma1, skew1, nsky=nsky1, /silent;, /debug
             skyinfo[ib].mode = mode1
             skyinfo[ib].sigma = sigma1
             skyinfo[ib].skew = skew1
             skyinfo[ib].nsky = nsky1
          endfor

; compute the 3-sigma surface brightness limit
          skyinfo.sblimit = -2.5*alog10(3*skyinfo.sigma)+5*alog10(pixscale)+zpt-kl*ebv
          
;; fit a plane          
;          xx = findgen(sz[0])#replicate(1,sz[1])
;          yy = replicate(1,sz[1])#findgen(sz[0])
;          for ib = 0, nfilt-1 do begin
;             acoeff = djs_sfit(data[ib].image[xcrop[0]:xcrop[1],ycrop[0]:ycrop[1]],$
;               xx,yy,1,1,yfit=yfit,$
;               sqivar=ivardata[ib].image[xcrop[0]:xcrop[1],ycrop[0]:ycrop[1]]*(fobj eq -1))
;          endfor

          struct_print, skyinfo
          im_mwrfits, skyinfo, datapath+cluster[ic]+'/'+cluster[ic]+'_skyinfo.fits', /clobber
          
; write out the sky-subtracted images
          for ib = 0, nfilt-1 do begin
             outfile = datapath+cluster[ic]+'/'+cluster[ic]+'-'+short[ib]+'.fits'
             splog, 'Writing '+outfile
             hextract, data[ib].image-skyinfo[ib].mode, data[ib].header[0:data[ib].nhead-1], $
               newim, newhdr, xcrop[0], xcrop[1], ycrop[0], ycrop[1], /silent
             hextract, ivardata[ib].image, ivardata[ib].header[0:ivardata[ib].nhead-1], $
               newivar, newivarhdr, xcrop[0], xcrop[1], ycrop[0], ycrop[1], /silent
             mwrfits, newim, outfile, newhdr, /create
             mwrfits, newivar, outfile, newivarhdr
             spawn, 'gzip -f '+outfile
          endfor 
       endif 

; --------------------------------------------------
; build the PSF in each band; for now use Molino's PSFs
       if keyword_set(build_psf) then begin
          splog, 'Building the PSFs'
          
          psfpath = getenv('CLASH_DATA')+'/psfs/indivfilters/'
          instr1 = repstr(instr,'acs','acswfc')
          for ib = 0, nfilt-1 do begin
             bpsf = mrdfits(psfpath+instr1[ib]+'_'+short[ib]+'.psfmos.fits',0,hdr,/silent)
             bpsf = bpsf/total(bpsf) ; normalize
             dfit_mult_gauss, bpsf, 2, amp, psfsig, model=model, /quiet
             sxaddpar, hdr, 'PSFSIGMA', float(psfsig[0]), ' Gaussian sigma [pixel]'

             outfile = datapath+cluster[ic]+'/'+cluster[ic]+'-'+short[ib]+'-bpsf.fits'
             splog, 'Writing '+outfile
             mwrfits, float(bpsf), outfile, hdr, /create, /silent
             mwrfits, float(model), outfile, hdr, /silent

; variable PSF; just use a constant here
             sz = size(bpsf,/dim)
             outfile = datapath+cluster[ic]+'/'+cluster[ic]+'-'+short[ib]+'-vpsf.fits'
             mwrfits, float(bpsf), outfile, hdr, /create, /silent
             mwrfits, bpsf*0+1, outfile, /silent
             mwrfits, fltarr(1,1), outfile, /silent
          endfor
          
; read the standard SE catalog; should replace this more generally
; with my own call to SE or a kludge call to dimage
;         imfiles = datapath+cluster[ic]+'/'+cluster[ic]+'-'+short+'.fits.gz'
;         for ib = 0, nfilt-1 do begin
;            streams_fitpsf, imfiles[ib], base=cluster[ic], maxnstar=50
;         endfor
             
;          cat = rsex(getenv('IM_ARCHIVE_DIR')+'/'+cluster[ic]+'/'+cluster[ic]+'_IR.cat')
;          for ib = 0, nfilt-1 do begin
;             image = gz_mrdfits(imfiles[ib],0,hdr)
;             extast, hdr, astr
;
;             cat = rsex(catpath+cluster[ic]+'_'+short[ib]+'.cat')
;;            wstar = 
;;            magindx = tag_indx(cat,short[ib]+'_mag')
;;            wstar = where(cat.stel gt 0.9 and cat.(magindx) gt 0 and cat.(magindx) lt 90,nw)
;;            ad2xy, cat[wstar].ra, cat[wstar].dec, astr, xx, yy
;          endfor

       endif

; --------------------------------------------------
; find parents; creates pcat, pimage, and parents files 
       if keyword_set(parents) then begin
          splog, 'Finding parents'

; write the PSET structure
          base = cluster[ic]
          imfiles = base+'/'+cluster[ic]+'-'+short+'.fits.gz'
          
          pset = {$
            base: base, $
            imfiles: imfiles, $
            ref: reffilt, $
            puse: intarr(nfilt)+1, $
            dopsf: intarr(nfilt)+1}
          im_mwrfits, pset, datapath+cluster[ic]+'/'+cluster[ic]+'-pset.fits', /clobber
          
          plim = 10.0
          dparents, cluster[ic], datapath+imfiles, ref=reffilt, $
            noclobber=noclobber, plim=plim, pbuffer=pbuffer
       endif

; --------------------------------------------------
; do stuff       

       if keyword_set(doit) then begin
; find galaxies
          skyinfo = mrdfits(datapath+cluster[ic]+'/'+cluster[ic]+'_skyinfo.fits.gz',1,/silent)

          image = mrdfits('parents/macs1206-parent-832.fits',2*reffilt,hdr)
          ivar = mrdfits('parents/macs1206-parent-832.fits',2*reffilt+1,ihdr)
          extast, hdr, astr
          
          simage = dsmooth(image,5)
          ssig = skyinfo[reffilt].sigma
          gsaddle = 1.0
          minpeak = 2.0*ssig
          dpeaks, simage, xc=xc, yc=yc, sigma=ssig, minpeak=minpeak, $
            npeaks=ngals, /check, saddle=gsaddle;, /refine
          
          atv, image
          atvplot, xc, yc, psym=8, color='red'

          xy2ad, xc, yc, astr, ra1, dec1
;         adxy, thdr, acat.racen, acat.deccen, xgals, ygals
          dtemplates, image, xc, yc, templates=curr_templates, $
            sersic=sersic, ikept=ikept

;         psf = mrdfits('macs1206-f160w-bpsf.fits')
          
stop          
          streams_stargal, /plot


          
      

          
          
stop          



; read the PSFs          
          for ib = 0, nfilt-1 do begin
             psffile = datapath+cluster[ic]+'/'+cluster[ic]+'-'+short[ib]+'-bpsf.fits'
             bpsf1 = mrdfits(psffile,0,/silent)
             if ib eq 0 then bpsf = bpsf1 else bpsf = [bpsf,bpsf1]
          endfor
          
;         pbuffer = 0.5

stop          
          
          dchildren, cluster[ic], 370, ref=reffilt, psfs=bpsf
          
          
stop          
       endif       
       
    endfor

return
end
    
