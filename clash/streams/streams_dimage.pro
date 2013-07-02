pro streams_dimage, skysubtract=skysubtract, build_psf=build_psf, bcg_prelim=bcg_prelim, $
  parents=parents, stargal=stargal, children=children, measure=measure, $
  rebuild_mosaic=rebuild_mosaic
; jm13may26siena - sky-subtract the CLASH clusters containing tidal
; streams 

    datapath = streams_path()

    filt = bcgimf_filterlist(short=short,instr=instr,weff=weff,zpt=zpt)
;   niceprint, lindgen(n_elements(filt)), short, instr
;   these = [9,12] ; = [F110W, F160W]
;   these = [8,9,10,11,12]
;   these = [6,9,12]
    these = lindgen(n_elements(filt))
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

; region file that defines by the WFC3 footprint (created by
; J. Rogers)     
    wfc3regfile = datapath+'macs1206_f160w.reg'
    
; wrap on each cluster    
    for ic = 0, ncl-1 do begin
       splog, 'Working on cluster '+cluster[ic]
       path = datapath+cluster[ic]+'/'

       base = cluster[ic]
;      base = cluster[ic]+'-nobcg'
       
       this =  where(cluster[ic] eq strtrim(clash.shortname,2))
       arcsec2kpc = dangular(clash[this].z,/kpc)/206265D ; [kpc/arcsec]
       ebv = clash[this].ebv

       if file_test(path,/dir) eq 0 then spawn, 'mkdir -p '+datapath+cluster[ic]
       mosaicpath = getenv('IM_ARCHIVE_DIR')+'/'+cluster[ic]+'/HST/images/'+$
         'mosaicdrizzle_image_pipeline/scale_65mas/'
;      mosaicpath = getenv('IM_ARCHIVE_DIR')+'/'+cluster[ic]+'/patched_images/'

; --------------------------------------------------
; subtract the sky
       if keyword_set(skysubtract) then begin
          skyinfo = replicate({file: '', band: '', sblimit: 0.0, $
            mode: 0.0, sigma: 0.0, skew: 0.0, nsky: 0L},nfilt)
          skyinfo.band = short
          
; read the original mosaics and inverse variance maps; crop all the
; images to match the WFC3 footprint (since the reference filter is
; F160W); also convert to physical units
          splog, 'Reading data'
          data = im_fits_cube(file_search(mosaicpath+cluster[ic]+$
            '_mosaic_065mas_'+instr+'_'+short+'_drz_????????.fits*'))
          ivardata = im_fits_cube(file_search(mosaicpath+cluster[ic]+$
            '_mosaic_065mas_'+instr+'_'+short+'_wht_????????.fits*'))

          wfc3mask = ds9polygon_indices(wfc3regfile,$
            header=data[reffilt].header,/inverse)

; [counts/s] --> [10^-12 erg/s/cm^2/Hz] (pico-maggies)          
          factor = 10^(-0.4*(zpt-kl*ebv))*1D12 
          for ib = 0, nfilt-1 do begin
             data[ib].image = data[ib].image*factor[ib]
             ivardata[ib].image = ivardata[ib].image/factor[ib]^2
             
             data[ib].image[wfc3mask] = 0.0
             ivardata[ib].image[wfc3mask] = 0.0
          endfor

; round to the nearest hundred pixels          
          xcrop = ceil(minmax(where(total(data[reffilt].image,1) gt 0))/100.0)*100
          ycrop = ceil(minmax(where(total(data[reffilt].image,2) gt 0))/100.0)*100

; aggresively build an object mask
          splog, 'Building object mask'
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

; compute the 1-sigma surface brightness limit
;         skyinfo.sblimit = -2.5*alog10(skyinfo.sigma)+5*alog10(pixscale) ; +zpt-kl*ebv
          skyinfo.sblimit = -2.5*alog10(skyinfo.sigma)+5*alog10(pixscale)-2.5*alog10(1D-12)

;; fit a plane          
;          xx = findgen(sz[0])#replicate(1,sz[1])
;          yy = replicate(1,sz[1])#findgen(sz[0])
;          for ib = 0, nfilt-1 do begin
;             acoeff = djs_sfit(data[ib].image[xcrop[0]:xcrop[1],ycrop[0]:ycrop[1]],$
;               xx,yy,1,1,yfit=yfit,$
;               sqivar=ivardata[ib].image[xcrop[0]:xcrop[1],ycrop[0]:ycrop[1]]*(fobj eq -1))
;          endfor

          struct_print, skyinfo
          im_mwrfits, skyinfo, path+cluster[ic]+'-skyinfo.fits', /clobber
          
; write out the sky-subtracted images
          for ib = 0, nfilt-1 do begin
             outfile = path+cluster[ic]+'-'+short[ib]+'.fits'
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
          
          psfpath = getenv('CLASH_DATA')+'/psfs-molino/indivfilters/'
          instr1 = repstr(instr,'acs','acswfc')
          for ib = 0, nfilt-1 do begin
             bpsf = mrdfits(psfpath+instr1[ib]+'_'+short[ib]+'.psfmos.fits',0,hdr,/silent)
             bpsf = bpsf/total(bpsf) ; normalize
             dfit_mult_gauss, bpsf, 2, amp, psfsig, model=model, /quiet
             sxaddpar, hdr, 'PSFSIGMA', float(psfsig[0]), ' Gaussian sigma [pixel]'

             outfile = path+cluster[ic]+'-'+short[ib]+'-bpsf.fits'
             splog, 'Writing '+outfile
             mwrfits, float(bpsf), outfile, hdr, /create, /silent
             mwrfits, float(model), outfile, hdr, /silent

;; variable PSF; just use a constant here
;             sz = size(bpsf,/dim)
;             outfile = path+cluster[ic]+'-'+short[ib]+'-vpsf.fits'
;             mwrfits, float(bpsf), outfile, hdr, /create, /silent
;             mwrfits, bpsf*0+1, outfile, /silent
;             mwrfits, fltarr(1,1), outfile, /silent
          endfor
       endif

; --------------------------------------------------
; do a preliminary subtraction of the BCG in each band separately  
       if keyword_set(bcg_prelim) then begin
          splog, 'Preliminary BCG subtraction'

          rmaxkpc = 60D                            ; [kpc]
          rmax = ceil(rmaxkpc/arcsec2kpc/pixscale) ; [pixels]

          for ib = nfilt-1, 0, -1 do begin
;         for ib = 0, nfilt-1 do begin
             imfile = cluster[ic]+'-'+short[ib]+'.fits.gz'
             splog, 'Reading '+imfile
             image = mrdfits(imfile,0,hdr,/silent)
             invvar = mrdfits(imfile,1,ivarhdr,/silent)
             
; get a RMAXKPC by RMAXKPC cutout centered on the BCG
             extast, hdr, astr
             ad2xy, clash[this].ra_bcg, clash[this].dec_bcg, astr, xcen, ycen
             xcen = fix(xcen)
             ycen = fix(ycen)
             hextract, image, hdr, cutimage, cuthdr, xcen-rmax, $
               xcen+rmax, ycen-rmax, ycen+rmax, /silent
             hextract, invvar, ivarhdr, cutinvvar, cutivarhdr, xcen-rmax, $
               xcen+rmax, ycen-rmax, ycen+rmax, /silent
             badpixels = where(cutinvvar le 0,nbad)
             if nbad eq 0L then delvarx, badpixels

             mge1 = streams_mge(cutimage,badpixels=badpixels,$
               pixscale=pixscale,model=cutmodel,twist=twist)
             help, mge1, /str

; rebuild the original image with the BCG subtracted and write out 
             model = image*0.0
             embed_stamp, model, cutmodel, xcen-rmax, ycen-rmax
;            atv, image-model, /bl

             outfile = path+cluster[ic]+'-nobcg-'+short[ib]+'.fits'
             splog, 'Writing '+file_basename(outfile)

             mwrfits, image-model, outfile, hdr, /create
             mwrfits, invvar, outfile, ivarhdr, /silent
             spawn, 'gzip -f '+outfile
          endfor
       endif 

; find and subtract the stars in the reference image and in all the
; other images (see streams_stargal)
;      skyinfo = mrdfits(path+cluster[ic]+'-skyinfo.fits.gz',1,/silent)
;      image = mrdfits(cluster[ic]+'-'+short[reffilt]+'.fits.gz',0,hdr)

       reffilt = 12
       skyinfo = mrdfits('macs1206-skyinfo.fits.gz',1,/silent)
       glim = 3.0
       gsmooth = 5.0
       gsaddle = 10.0
       ssig = skyinfo[reffilt].sigma
       saddle = gsaddle*ssig

       psf = mrdfits('macs1206-f160w-bpsf.fits',0)
       image = mrdfits('macs1206-f160w.fits.gz',0,hdr)
       invvar = mrdfits('macs1206-f160w.fits.gz',1,ihdr)
       simage = dsmooth(image,gsmooth)
       dpeaks, simage, xc=xgals, yc=ygals, sigma=ssig, minpeak=glim*ssig, $
         npeaks=ngals, saddle=gsaddle, /refine, /check

;; build the templates and deblend into children
;       ww = where(xgals gt 535 and xgals lt 825 and ygals gt 1000 and ygals lt 1190)
;       dtemplates, image, xgals[ww], ygals[ww], templates=temp, $
;         sersic=sersic, ikept=ikept, sigma=ssig
;
;       dweights, image, invvar, temp, weights=weights, /nonneg
;       dfluxes, image, temp, weights, xgals[ww], ygals[ww], $
;         children=children, sigma=ssig

; the above didn't work well; try just getting a cutout of the
; largest galaxies
       x0 = 420 & x1 = 720 & y0 = 900 & y1 = 1200
       hextract, image, hdr, im, hd, x0, x1, y0, y1
       hextract, invvar, ihdr, ivar, ihd, x0, x1, y0, y1
       dpeaks, im, xc=xc, yc=yc, sigma=ssig, minpeak=glim*ssig, $
         npeaks=ngals, saddle=gsaddle, /refine, /check, maxnpeaks=1
       dmeasure, im, ivar, xcen=xc, ycen=yc, measure=measure, $
         check=check, cpetrorad=cpetrorad, /fixcen, faper=faper

stop       
       
       dsersic2, im, ivar, xcen=xc, ycen=yc, sersic=ser, psf=psf, /fixcen, $
         fixsky=0, model=model, bulge=bulge, disk=disk, /reinit

       mosaic = image*0 & embed_stamp, mosaic, model-ser.sky, x0, y0

       
       
stop       
       
; --------------------------------------------------
; find parents; creates pcat, pimage, and parents files 
       if keyword_set(parents) then begin
          splog, 'Finding parents'

          imfiles = base+'-'+short+'.fits.gz'
          pset = {$
            base:    base, $
            imfiles: imfiles, $
            ref: reffilt, $
            puse: intarr(nfilt)+1, $
            dopsf: intarr(nfilt)+1}
          im_mwrfits, pset, path+base+'-pset.fits', /clobber

          plim = 8.0
          pbuffer = 0.1
          dparents, base, path+imfiles, ref=reffilt, /cenonly, $
            noclobber=noclobber, plim=plim, pbuffer=pbuffer
          heap_gc

; now clean up the parents!
       endif

; --------------------------------------------------
; classify objects into stars and galaxies; also subtract foreground
; PSF stars
       if keyword_set(stargal) then begin
          splog, 'Performing star-galaxy separation.'

          skyinfo = mrdfits(path+cluster[ic]+'-skyinfo.fits.gz',1,/silent)
          niceprint, skyinfo.band, skyinfo.sigma
          psffiles = file_search(path+cluster+'-'+short+'-bpsf.fits*')

          streams_stargal, base, plot=plot, gsmooth=gsmooth, glim=glim, $
            gsaddle=gsaddle, nsigma=nsigma, noclobber=noclobber, $
            psffiles=psffiles
          heap_gc
       endif

; --------------------------------------------------
; deblend into children
       if keyword_set(children) then begin
          splog, 'Deblending parents into children'
          streams_children, base, sersic=sersic, noclobber=noclobber
          heap_gc
       endif

; --------------------------------------------------
; model each object
       if keyword_set(measure) then begin

; test code          
          cutimage = mrdfits('atlases/303/macs1206-303-templates-0.fits',0,hdr)

          mge1 = streams_mge(cutimage,badpixels=badpixels,$
            pixscale=pixscale,model=cutmodel,twist=twist)
          mge1 = streams_mge(cutimage,badpixels=badpixels,pixscale=pixscale,$
            model=cutmodel,/twist,minlevel=3*ss[12].sigma)
          help, mge1, /str
          


stop
          splog, 'Modeling!'
          psffiles = file_search(path+cluster+'-'+short+'-bpsf.fits*')
          streams_measure, base, psffiles=psffiles, noclobber=noclobber
          heap_gc
       endif

; --------------------------------------------------
; rebuild the mosaic with the galaxies subtracted
       if keyword_set(rebuild_mosaic) then begin
          splog, 'Rebuilding the mosaic!'
          base = cluster[ic]+'-nobcg'
          base2 = cluster[ic]

          sub = 'atlases'
          postfix = ''
          
          for ib = 0, nfilt-1 do begin
             splog, 'Working on filter '+short[ib]

             imfile = path+base2+'-'+short[ib]+'.fits.gz'
             outfile = path+base2+'-'+short[ib]+'-model.fits'
             image = gz_mrdfits(imfile,0,hdr,/silent)
             invvar = gz_mrdfits(imfile,1,ivarhdr,/silent)
             model = image*0
             
             pcat = gz_mrdfits(path+base+'-pcat.fits',1)
             for iparent = 0L, n_elements(pcat)-1 do begin
                splog, 'Parent ', iparent
    
                pstr = strtrim(string(iparent),2)
                mfile = path+'/'+sub+'/'+pstr+'/'+base+'-'+pstr+ $
                  '-measure'+postfix+'.fits.gz'
                sfile = path+'/'+sub+'/'+pstr+'/'+base+'-'+pstr+ $
                  '-sersic'+postfix+'.fits.gz'
                if file_test(sfile) then begin
                   model1 = gz_mrdfits(sfile,ib,hdr1,/silent)
                   embed_stamp, model, model1, pcat[iparent].xst[ib], pcat[iparent].yst[ib]
                endif
             endfor
             im_mwrfits, model, outfile, hdr, /nogzip, /clobber
             im_mwrfits, image-model, outfile, hdr, /gzip, /append
          endfor 
stop
       endif 
       
    endfor ; close cluster loop

return
end
    
