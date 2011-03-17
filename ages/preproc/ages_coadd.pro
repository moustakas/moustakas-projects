pro ages_coadd, writeout=writeout
; jm07feb23nyu - excised from HS_EXTRACT

    rootpath = ages_path(/rawdata)
    rerun = ['0301','0051']

    nfib = 300L

    badplates = [$
      '106',$ ; not fluxed
      '110',$ ; not fluxed
      '209',$ ; not fluxed
      '310',$ ; not fluxed
      '311' $ ; not fluxed
      ]

;   for irun = 1L, n_elements(rerun)-1L do begin
    for irun = 0L, n_elements(rerun)-1L do begin

       datapath = rootpath+rerun[irun]+'/'
       passfield = strmid(file_basename(file_search(datapath+'spHect-???-'+rerun[irun]+'.fits.gz',count=npass)),7,3)

       splog, filename=datapath+'splog_ages_coadd.txt'
       t0 = systime(1)

;      for ipass = 40L, 40L do begin
       for ipass = 0L, npass-1L do begin

          splog, 'Rerun '+rerun[irun]+', pass/field '+passfield[ipass]
          filelist = file_search(datapath+'spObs-'+passfield[ipass]+'-'+rerun[irun]+'-????.fits.gz',count=nobs)
          outname = datapath+'spObs-'+passfield[ipass]+'-'+rerun[irun]+'.fits'

          lam = mrdfits(filelist[0],0,/silent)
          
          lamout = dblarr(nobs,n_elements(lam[*,0]),nfib)
          fluxout = dblarr(nobs,n_elements(lam[*,0]),nfib)
          ivarout =  dblarr(nobs,n_elements(lam[*,0]),nfib)
          pixmaskout =  dblarr(nobs,n_elements(lam[*,0]),nfib)

          for ii = 0L, nobs-1L do begin
             splog, 'Reading '+filelist[ii]
             lamout[ii,*,*] = mrdfits(filelist[ii],0,/silent)
             fluxout[ii,*,*] = mrdfits(filelist[ii],1,/silent)
             ivarout[ii,*,*] = mrdfits(filelist[ii],2,/silent)
             pixmaskout[ii,*,*] = mrdfits(filelist[ii],3,/silent)
             if (ii eq 0L) then plugmapout = mrdfits(filelist[ii],5,/silent) ; correct?
;            if (ii eq 0L) then plugmapout = mrdfits(filelist[ii],5,/silent) else $
;              plugmapout = [plugmapout,mrdfits(filelist[ii],5,/silent)]
          endfor

          if (total(strmatch(badplates,passfield[ipass])) eq 0.0) then ftweak = 1L else ftweak = 0L
          
          if (ftweak eq 1L) then begin
             
             sn = fltarr(nobs)
             for im = 0, nobs-1 do begin
                hs_snr,  transpose(transpose(lamout[im,*,*])), $
                  transpose(transpose(fluxout[im,*,*])), plugmapout, $
                  sn2=sn21, /noplot, /nomag
                sn(im) =sn21  
             endfor
             
             expid = filelist
             sn_exp = sn
             maxval = max(sn_exp/median(sn_exp), iframe_best)
             if not keyword_set(best_exposure) then $
               best_exposure = expid[iframe_best]
             splog, 'Best Exposure is: ' + best_exposure 

; Compute the exposure-to-exposure corrections
             
             lam = lamout
             flux = fluxout
             ivar = ivarout
             mask = pixmaskout
             
             ndim = size(lam,/n_dimensions)
             lamorig = lam
             airtovac, lamorig
             
             if (ndim eq 3L) then begin
                
                size1 = size(lam)
                nim = size1(1)
                npix = size1(2)
                nfib = size1(3)
                newlam = dblarr(npix, nfib*nim)
                
                newflux = newlam
                newivar = newlam
                newmask = newlam
                
                count = 0L
                for i = 0L, nim -1L do begin
                   
                   newlam(*,count:count+(nfib-1L))  =  $
                     transpose(transpose(lam(i,  *, *)))
                   newflux(*,count:count+(nfib-1L)) = $
                     transpose(transpose(flux(i, *, *)))
                   newivar(*,count:count+(nfib-1L)) =  $
                     transpose(transpose(ivar(i, *, *)))
                   newmask(*,count:count+(nfib-1L)) =  $
                     transpose(transpose(mask(i, *, *)))
                   plugmapout.frames = i 
                   plugmapout.expid = filelist[i]
                   print, filelist[i]
                   
                   if (i eq 0L) then newplug = plugmapout
                   if (i ne 0L) then newplug = [newplug, plugmapout]
                   count = nfib*(i+1L)
                   
                endfor

                lam = newlam
                flux = newflux
                ivar = newivar
                mask = newmask
                plugmap = newplug
                
                corrfiles = datapath+'spFluxcorr-'+file_basename(filelist)
;               corrfiles = 'reduction/' + rerun +'/spFluxcorr-' + junk1

                if (ftweak eq 1L) then begin
                   
;                  hs_frame_flux_tweak, alog10(lam), flux, ivar,$ ; correct?
;                    best_exposure, plugmap, corrfiles, /diag
                   
; Read back the corrections and apply
        
                   for iexp = 0L, nobs-1L do begin

                      corrfile = corrfiles[iexp]
                      splog, 'Reading '+corrfile
                      corrset = mrdfits(corrfile,1,/silent)

                      indx = where(plugmap.expid eq expid[iexp])
                      traceset2xy, corrset, alog10(lam[*,indx]), corrimg

; Don't let the flux correction be more than a factor of 10

                      invertcorr = 1.0 / corrimg
                      tempflux = flux[*,indx]
                      tempivar = ivar[*,indx]
                      
                      divideflat, tempflux, invvar=tempivar, $
                        invertcorr, minval=0.1

                      splog, mean(corrimg), median(corrimg), stddev(corrimg), minmax(corrimg)
                      
                      tempmask = mask[*, indx] OR $
                        (corrimg GE 10) * pixelmask_bits('BADFLUXFACTOR')
                      tempmask = mask[*, indx] OR $
                        (corrimg LE 0.1) * pixelmask_bits('BADFLUXFACTOR') ;
                      
                      fluxout[iexp,*,*] = tempflux
                      lamout[iexp,*,*] = lam[*,indx]
                      ivarout[iexp,*,*] = tempivar
                      pixmaskout[iexp,*,*] = tempmask
                      
                   endfor
                   
                endif else begin

                   for iexp = 0L, nobs-1L do begin
                      
                      indx = where(plugmap.expid eq expid[iexp])
                      tempflux = flux[*,indx]
                      tempivar = ivar[*,indx]
                      tempmask = mask[*, indx]
                      fluxout[iexp,*,*] = tempflux
                      lamout[iexp,*,*] = lam[*,indx]
                      ivarout[iexp,*,*] = tempivar
                      pixmaskout[iexp,*,*] = tempmask
                      
                   endfor

                endelse 

             endif 

          endif 

; now combine each of the fibers into a composite spectrum

          newlam = transpose(transpose(lamout[0,*,*]))
          newflux = transpose(transpose((fluxout[0,*,*]))) * 0.0
          newfluxivar = transpose(transpose(ivarout[0,*,*])) * 0.0
          newmask = transpose(transpose(ivarout[0,*,*])) * 0.0
          ormask = newmask
          
          for i = 0L, nfib-1L do begin
             
             lamin = transpose(lamout(*,*,i))
             
             print, format='("Co-adding ",i4," of ",i4,a1,$)', $
               i, nfib, string(13b)
             
             hs_combine1fiber,(lamin), $
               transpose(fluxout(*,*,i)), transpose(ivarout(*,*,i)), $
               newlam=(newlam(*,i)), $
               newflux=newflux1, newivar=newivar1, $
               finalmask = transpose(pixmaskout(*,*,i)), $
               andmask = newmask1, ormask=ormask1
             newflux[*,i] = newflux1
             newfluxivar[*,i] = newivar1
             newmask[*,i] = newmask1
             ormask[*,i] = ormask1
          endfor
          
          fluxout = newflux
          ivarout = newfluxivar
          lamout = newlam
          pixmaskout = newmask
          ormaskout = ormask

          if keyword_set(writeout) then begin
             
             splog, 'Writing coadded spObs file ' + outname
             
             sxaddpar, objhdr, 'NAXIS2', n_elements(lamout[*,0]), after='NAXIS'
             sxaddpar, objhdr, 'EXTEND', 'T', after='NAXIS2'

             mwrfits, float(lamout), outname, /create, objhdr
             mwrfits, float(fluxout), outname
             mwrfits, float(ivarout), outname
             mwrfits, float(pixmaskout), outname
             mwrfits, float(ormaskout), outname
             mwrfits, plugmapout, outname
             spawn, 'gzip -f '+outname

          endif 

       endfor
       splog, 'Total time to unpack rerun '+rerun[irun]+' = '+string((systime(1)-t0)/60.0,format='(G0.0)')+' minutes.'
       splog, /close

    endfor 
       
return
end

;pro ages_coadd, writeout=writeout
;; jm07feb23nyu - excised from HS_EXTRACT
;    
;    rootpath = ages_path(/rawdata)
;
;    rerun = ['0nfib','0051']
;    for irun = 0L, n_elements(rerun)-1L do begin
;
;       datapath = rootpath+rerun[irun]+'/'
;       passfield = strmid(file_basename(file_search(datapath+'spHect-???-'+rerun[irun]+'.fits.gz',count=npass)),7,3)
;
;;      for ipass = 39L, 40L do begin
;       for ipass = 0L, npass-1L do begin
;
;          splog, 'Rerun '+rerun[irun]+', pass/field '+passfield[ipass]
;          filelist = file_search(datapath+'spObs-'+passfield[ipass]+'-'+rerun[irun]+'-????.fits.gz',count=nobs)
;          outname = datapath+'spObs-'+passfield[ipass]+'-'+rerun[irun]+'.fits'
;
;          lam = mrdfits(filelist[0],0,/silent)
;          
;          lamout = dblarr(nobs,n_elements(lam[*,0]),nfib)
;          fluxout = dblarr(nobs,n_elements(lam[*,0]),nfib)
;          ivarout =  dblarr(nobs,n_elements(lam[*,0]),nfib)
;          pixmaskout =  dblarr(nobs,n_elements(lam[*,0]),nfib)
;
;          for ii = 0L, nobs-1L do begin
;             splog, 'Reading '+filelist[ii]
;             lamout[ii,*,*] = mrdfits(filelist[ii],0,/silent)
;             fluxout[ii,*,*] = mrdfits(filelist[ii],1,/silent)
;             ivarout[ii,*,*] = mrdfits(filelist[ii],2,/silent)
;             pixmaskout[ii,*,*] = mrdfits(filelist[ii],3,/silent)
;             if (ii eq 0L) then plugmapout = mrdfits(filelist[ii],5,/silent)
;;            if (ii eq 0L) then plugmapout = mrdfits(filelist[ii],5,/silent) else $
;;              plugmapout = [plugmapout,mrdfits(filelist[ii],5,/silent)]
;          endfor
;
;; now combine each of the fibers into a composite spectrum; also
;; correct for Milky Way dust
;
;          newlam = transpose(transpose(lamout[0,*,*]))
;          newflux = transpose(transpose((fluxout[0,*,*]))) * 0.0
;          newfluxivar = transpose(transpose(ivarout[0,*,*])) * 0.0
;          newmask = transpose(transpose(ivarout[0,*,*])) * 0.0
;          ormask = newmask
;          
;          newplug_struct = create_struct(plugmapout[0],$
;            'tsobjid', lonarr(5),'tsobj_mag',fltarr(5), 'ebv_sfd',0.0)
;          newplug = make_array(val=newplug_struct,dim=nfib)
;          struct_assign, plugmapout, newplug
;          newplug.tsobj_mag = plugmapout.mag
;          plugmapout = newplug
;
;          glactc, plugmapout.ra, plugmapout.dec, 2000, gall, galb, 1, /degre
;          ebv_sfd = dust_getval(gall,galb,/interp)
;          plugmapout.ebv_sfd = ebv_sfd
;
;          for ifib = 0L, nfib-1L do begin
;             
;             lamin = transpose(lamout[*,*,ifib])
;             
;             print, format='("Co-adding ",i4," of ",i4,a1,$)', $
;               ifib, nfib, string(13b)
;             
;             hs_combine1fiber,(lamin), $
;               transpose(fluxout[*,*,ifib]),transpose(ivarout[*,*,ifib]),newlam=(newlam[*,ifib]),$
;               newflux=newflux1,newivar=newivar1,finalmask=transpose(pixmaskout[*,*,ifib]),$
;               andmask=newmask1,ormask=ormask1
;
;             dustcor = 10.0^(0.4*k_lambda(reform(newlam[*,ifib]),/odonnell,r_v=3.1)*ebv_sfd[ifib])
;             
;             newflux[*,ifib] = newflux1*dustcor
;             newfluxivar[*,ifib] = newivar1/dustcor^2.0
;             newmask[*,ifib] = newmask1
;             ormask[*,ifib] = ormask1
;
;          endfor
;          
;          fluxout = newflux
;          ivarout = newfluxivar
;          lamout = newlam
;          pixmaskout = newmask
;          ormaskout = ormask
;
;          if keyword_set(writeout) then begin
;             
;             splog, 'Writing coadded spObs file ' + outname
;             
;             sxaddpar, objhdr, 'NAXIS2', n_elements(lamout[*,0]), after='NAXIS'
;             sxaddpar, objhdr, 'EXTEND', 'T', after='NAXIS2'
;
;             mwrfits, float(lamout), outname, /create, objhdr
;             mwrfits, float(fluxout), outname
;             mwrfits, float(ivarout), outname
;             mwrfits, float(pixmaskout), outname
;             mwrfits, float(ormaskout), outname
;             mwrfits, plugmapout, outname
;             spawn, 'gzip -f '+outname
;
;          endif
;
;       endfor
;
;    endfor
;       
;return
;end
