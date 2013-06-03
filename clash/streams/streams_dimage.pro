function streams_mge, image, badpixels=badpixels, pixscale=pixscale, $
  model=model, twist=twist
; fit the BCG using MGE; ToDo: fit the PSF using an MGE

    forward_function multi_gauss, multi_gauss_twist
    resolve_routine, 'mge_print_contours', /compile_full_file, /no_recompile
    resolve_routine, 'mge_print_contours_twist', /compile_full_file, /no_recompile

; find the galaxy              
    find_galaxy, image, majoraxis, ellipticity, posangle, xcen, ycen, $
      xcen_lum, ycen_lum, fraction=fraction, index=index, level=level, $
      nblob=nblob, plot=plot, /quiet
    
; get the photometry
    if keyword_set(twist) then begin
       sectors_photometry_twist, image, ellipticity, posangle, xcen, ycen, $
         radius, phi, counts, n_sectors=n_sectors, sector_width=sector_width, $
         badpixels=badpixels, minlevel=minlevel
    endif else begin
       sectors_photometry, image, ellipticity, posangle, xcen, ycen, $
         radius, phi, counts, n_sectors=n_sectors, sector_width=sector_width, $
         badpixels=badpixels, minlevel=minlevel
    endelse

; do the MGE fitting
    if keyword_set(twist) then begin
       mge_fit_sectors_twist, radius, phi, counts, ellipticity, ngauss=ngauss, $
         scale=pixscale, bulge_disk=bulge_disk, fastnorm=fastnorm, $
         linear=linear, negative=negative, normpsf=normpsf, print=print, $
         qbounds=qbounds, rbounds=rbounds, /quiet, sigmapsf=sigmapsf, $
         outer_slope=outer_slope, sol=sol, absdev=absdev
;      mge_print_contours_twist, image, ang, xpeak, ypeak, sol, model=model
       model = multi_gauss_twist(sol,image,sigmapsf,normpsf,xpeak,ypeak,posangle)
    endif else begin
       mge_fit_sectors, radius, phi, counts, ellipticity, ngauss=ngauss, $
         scale=pixscale, bulge_disk=bulge_disk, fastnorm=fastnorm, $
         linear=linear, negative=negative, normpsf=normpsf, print=print, $
         qbounds=qbounds, rbounds=rbounds, /quiet, sigmapsf=sigmapsf, $
         outer_slope=outer_slope, sol=sol, absdev=absdev
;      mge_print_contours, cutimage, ang, xpeak, ypeak, sol, model=model
;      model = multi_gauss(sol,image,sigmapsf,normpsf,xcen,ycen,posangle)
       model = float(multi_gauss(sol,image,0.0,1.0,xcen,ycen,posangle))
       model1 = multi_gauss(sol,image,0.0,1.0,xcen,ycen,90) ; aligned
    endelse
    
; get the surface brightness profile of the model along the semimajor
; axis 
    sz = size(image,/dim)
    sma = range(1.0,200,200,/log)
;   sma = range(0.0,200,201)
    sb_sma = interpolate(model1[xcen:sz[0]-1,ycen],$ ;missing=0.0,$
      findex(range(xcen,sz[0]-1,sz[0]-xcen),xcen+sma))
    
    mge = {$
      majoraxis:   majoraxis,$
      ellipticity: ellipticity,$
      posangle:    posangle,$
      xcen:        xcen,$
      ycen:        ycen,$
      xcen_lum:    xcen_lum,$
      ycenlum:     ycen_lum,$
      absdev:      absdev,$
      sol:         sol,$
      sma:         sma,$
      sb_sma:      sb_sma}
;     radius:      radius,$
;     phi:         phi,$
;     counts:      counts,$
;     sbprofile:   radius*0.0}
;   mge.sbprofile = sol[0,*]/(2.0*!dpi*sol[1,*]^2*sol[2,*])

return, mge
end
    
pro streams_dimage, skysubtract=skysubtract, build_psf=build_psf, bcg_prelim=bcg_prelim, $
  parents=parents, stargal=stargal, children=children, measure=measure, $
  rebuild_mosaic=rebuild_mosaic
; jm13may26siena - sky-subtract the CLASH clusters containing tidal
; streams 

    propath = getenv('CLASH_DIR')+'/streams/'
    datapath = getenv('IM_ARCHIVE_DIR')+'/projects/clash/streams/'
;   datapath = getenv('IM_PROJECTS_DIR')+'/clash/streams/'

    filt = bcgimf_filterlist(short=short,instr=instr,weff=weff,zpt=zpt)
;   these = [9,12] ; = [F110W, F160W]
;   these = [8,9,10,11,12]
    these = [8,9,12]
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
       path = datapath+cluster[ic]+'/'
       
       this =  where(cluster[ic] eq strtrim(clash.shortname,2))
       arcsec2kpc = dangular(clash[this].z,/kpc)/206265D ; [kpc/arcsec]
       ebv = clash[this].ebv

       if file_test(path,/dir) eq 0 then spawn, 'mkdir -p '+datapath+cluster[ic]
       mosaicpath = getenv('IM_ARCHIVE_DIR')+'/'+cluster[ic]+'/images/'
;      mosaicpath = getenv('IM_ARCHIVE_DIR')+'/'+cluster[ic]+'/patched_images/'
       catpath = getenv('IM_ARCHIVE_DIR')+'/'+cluster[ic]+'/catalogs/'

; --------------------------------------------------
; subtract the sky
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
; do a preliminary subtraction of the BCG in each band separately  
       if keyword_set(bcg_prelim) then begin
          splog, 'Preliminary BCG subtraction'

          rmaxkpc = 60D                            ; [kpc]
          rmax = ceil(rmaxkpc/arcsec2kpc/pixscale) ; [pixels]

          for ii = 0, nfilt-1 do begin
             imfile = cluster[ic]+'-'+short[ii]+'.fits.gz'
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

             outfile = path+cluster[ic]+'-nobcg-'+short[ii]+'.fits'
             splog, 'Writing '+file_basename(outfile)
             mwrfits, image-model, outfile, hdr, /create
             mwrfits, invvar, outfile, ivarhdr, /silent
             spawn, 'gzip -f '+outfile
          endfor
       endif 
       
; --------------------------------------------------
; find parents; creates pcat, pimage, and parents files 
       if keyword_set(parents) then begin
          splog, 'Finding parents'

; write the PSET structure
          base = cluster[ic]+'-nobcg'
          imfiles = base+'-'+short+'.fits.gz'
          
          pset = {$
            base: base, $
            imfiles: imfiles, $
            ref: reffilt, $
            puse: intarr(nfilt)+1, $
            dopsf: intarr(nfilt)+1}
          im_mwrfits, pset, path+cluster[ic]+'-nobcg-pset.fits', $
            /clobber, /nogzip

          plim = 20.0
          dparents, cluster[ic]+'-nobcg', path+imfiles, ref=reffilt, $
            noclobber=noclobber, plim=plim, pbuffer=pbuffer
       endif

; --------------------------------------------------
; classify objects into stars and galaxies; also subtract foreground
; PSF stars
       if keyword_set(stargal) then begin
          splog, 'Performing star-galaxy separation.'
          base = cluster[ic]+'-nobcg'

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
          base = cluster[ic]+'-nobcg'
          streams_children, base, noclobber=noclobber
          heap_gc
       endif

; --------------------------------------------------
; model each object
       if keyword_set(measure) then begin
          splog, 'Modeling!'
          base = cluster[ic]+'-nobcg'
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
          
          for ii = 0, nfilt-1 do begin
             splog, 'Working on filter '+short[ii]

             imfile = path+base2+'-'+short[ii]+'.fits.gz'
             outfile = path+base2+'-'+short[ii]+'-model.fits'
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
                   model1 = gz_mrdfits(sfile,ii,hdr1,/silent)
                   embed_stamp, model, model1, pcat[iparent].xst[ii], pcat[iparent].yst[ii]
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
    
