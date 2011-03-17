pro sc1120_write_sigspec, doplot=doplot, wfits=wfits
; jm05jan17uofa
; this routine should be run before all others    

    datapath = sc1120_path(/unfluxed)
    outpath = sc1120_path(/sigspec)
    zpath = sc1120_path(/analysis)

; read the redshift data

    zfile = 'sc1120-z.cat.sex'
    splog, 'Reading '+zpath+zfile+'.'
    zinfo = rsex(zpath+zfile)

    pagemaker, nx=1, ny=3, ymargin=[0.2,1.0], yspace=0.0, $
      xmargin=[1.3,0.2], position=pos, /normal

; construct the object list

    splog, 'Pushing into '+datapath & print

    pushd, datapath
    biglist = file_search('*.fits',count=nbiglist)
    bigindx = lindgen(nbiglist)
    
    skyindx = where(strmatch(biglist,'*sky*',/fold) eq 1B,nsky,comp=objindx,ncomp=nobj)

    objlist = biglist[objindx]
    skylist = biglist[skyindx]

    popd

; loop on each object

;   for iobj = 107, 107 do begin
    for iobj = 0L, nobj-1L do begin

       if (not keyword_set(wfits)) then $
         print, format='("Processing object ",I0,"/",I0,".",A1,$)', iobj+1, nobj, string(13b)

; extract the root object name

       galaxy = strupcase(repstr(objlist[iobj],'.fits',''))
       rootgalaxy = strupcase(strmid(objlist[iobj],0,strpos(objlist[iobj],'.')))
;      galaxy = strupcase(strmid(objlist[iobj],0,strpos(objlist[iobj],'.')))
       rootname = strmid(objlist[iobj],0,strpos(objlist[iobj],'.',/reverse_search))

       if strmatch(rootname,'*s*') then $
         skymatch = (where(strmatch(skylist,repstr(rootname,'s','skys')+'.fits',/fold) eq 1B,nskymatch))[0] else $
         skymatch = (where(strmatch(skylist,rootname+'sky.fits',/fold) eq 1B,nskymatch))[0]

       if (nskymatch ne 1L) then begin

          splog, 'No sky spectrum for '+objlist[iobj]+'.' ;message, 'Houston, we have a problem.'
          continue
          
       endif

;      print, objlist[iobj]+' '+skylist[skymatch]
       
; match this object to the redshift data file

       match = where(strmatch(zinfo.galid,rootgalaxy,/fold) eq 1B,nmatch)

       if (nmatch eq 0L) then begin
          print
          splog, '**** No REDSHIFT match for '+objlist[iobj]+' ****'
          print
       endif
       if (nmatch gt 1L) then begin
          splog, 'Multiple matches for '+objlist[iobj]+'.'
          stop
       endif

       if (nmatch eq 1L) then begin

;         print & splog, 'Match found:'

          zmatch = zinfo[match]

; read the object and sky spectra

          objflux = readfits(datapath+objlist[iobj],hobj,/silent) ; object spectrum [ADU]
          skyflux = readfits(datapath+skylist[skymatch],hsky,/silent) ; sky spectrum [ADU]
          objwave = make_wave(hobj)
          info = vlt_header_forage(hobj,/silent)

          gain = float(info.det_out1_conad)  ; gain [e/ADU]
          rdnoise = float(info.det_out1_ron) ; read-noise [e]
          ncombine = float(info.ncombine)    ; number of images combined
          
; generate the error spectrum; the 1D science spectra are the sum of
; the 5 central rows and the sky spectra = 5x(average 1D sky spectrum
; over slit); consequently, if n_pix = 5 and S_obj and S_sky are in
; ADU, then the noise spectrum is just:

;   N = sqrt( S_obj/G + S_sky/G + n_pix * (r/G)^2 ) ; [ADU]

          npix = 5L ; fixed!

          term1 = sqrt(objflux^2) / gain + sqrt(skyflux^2) / gain ; [ADU^2]
          term2 = npix * (rdnoise/gain)^2                         ; [ADU^2]
          
          objferr = sqrt(term1 + term2) ; [ADU]
          
          snrstats = im_stats(objflux/objferr,/no_head);,/verbose)
          snr_median = string(snrstats.median,format='(F4.1)') ; median S/N

; crop the wavelength range

          case objlist[iobj] of
             'g10359.3.fits': minwave = 5025.0
             'g10400.3.fits': minwave = 4950.0
             'g10407.3.fits': minwave = 4950.0
             'g10447.3.fits': minwave = 4900.0
             'g10991.2.fits': minwave = 5000.0
             'g20220.1.fits': minwave = 4950.0
             'g20298.3.fits': minwave = 5030.0
             'g20322.1.fits': minwave = 5500.0
             'g20332.1.fits': minwave = 5450.0
             'g20342.3.fits': minwave = 4950.0
             'g20374.1.fits': minwave = 5180.0
             'g20402.3.fits': minwave = 4875.0
             'g20409.1.fits': minwave = 5250.0
             'g20416.1.fits': minwave = 4980.0
             'g20424.1.fits': minwave = 5300.0
             'g20448.1.fits': minwave = 5010.0
             'g20492.1.fits': minwave = 5100.0
             'g20497.1.fits': minwave = 4960.0
             'g20526.1.fits': minwave = 5000.0
             'g20542.1.fits': minwave = 5100.0
             'g20661.1.fits': minwave = 4890.0
             'g22809.2.fits': minwave = 5700.0
             'g30082.1.fits': minwave = 5800.0
             'g30224.1.fits': minwave = 5450.0
             'g30236.1.fits': minwave = 5300.0
             'g30295.1.fits': minwave = 5250.0
             'g32149.2.fits': minwave = 5010.0
             'g32245.2.fits': minwave = 5160.0
             'g32687.1.fits': minwave = 5750.0
             'g40024.2.fits': minwave = 5250.0
             'g40029.2.fits': minwave = 5200.0
             'g40097.2.fits': minwave = 5600.0
             'g40138.2.fits': minwave = 5110.0
             'g40331.1.fits': minwave = 5450.0
             'g40335.3.fits': minwave = 4950.0
             'g40384.1.fits': minwave = 5120.0
             'g40440.3.fits': minwave = 5020.0
             'g40504.1.fits': minwave = 5350.0
             'g40534.1.fits': minwave = 5300.0
             'g43191.2.fits': minwave = 4920.0
             else: minwave = 4850.0
          endcase
             
          maxwave = 9200.0
          dwave = info.cd1_1

          newwave = findgen((maxwave-minwave)/dwave+1)*dwave+minwave
          objflux = interpol(objflux,objwave,newwave)
          objferr = sqrt(interpol(objferr^2.0,objwave,newwave))
          skyflux = interpol(skyflux,objwave,newwave)
          objwave = newwave

; interpolate over negative pixels

          neg = where(objflux le 0.0,nneg)
          if (nneg ne 0L) then begin

             while (nneg ne 0L) do begin

                mflux = dkboxstats(objflux,boxstat='mean',xwidth=30)
                objflux[neg] = mflux[neg]
                neg = where(objflux le 0.0,nneg)
                
             endwhile
             
          endif
          
;         good = where((objwave ge minwave) and (objwave le maxwave) and (objflux ne 0.0))
;         objflux = objflux[good]
;         objwave = objwave[good]
;         skyflux = skyflux[good]
          
; now, fringing is a pain in the ass redward of 7700 Angstroms; so, to
; enhance the statistical error in that regime, enhance the error by
; the rms variation in the spectrum

          fringes = where(objwave gt 7700.0,nfringes)
          if (nfringes ne 0L) then begin

;            result = im_medxbin(objwave[fringes],objflux[fringes],150.0,minpts=5.0)
;            sflux = interpol(result.medy,result.binctr,objwave[fringes])

             sflux = smooth(objflux[fringes],100,/edge_truncate)
             fsigma = djsig(objflux[fringes]-sflux,sigrej=3.0)

;            result = im_medxbin(objwave[fringes],objflux[fringes]-sflux,50.0,minpts=5.0)
;            fsigma = interpol(result.sigy,result.binctr,objwave[fringes])
             
;            ploterror, objwave[fringes], objflux[fringes]-sflux, replicate(fsigma,nfringes), xsty=3, ysty=3, ps=10
;            plot, objwave[fringes], objflux[fringes]-sflux, xsty=3, ysty=3, ps=10
;            djs_oplot, objwave[fringes], sflux, ps=10, color='yellow'
;            cc = get_kbrd(1)

;            objferr[fringes] = objferr[fringes] > fsigma
             
          endif
          
; correct for foreground Galactic extinction       
          
          glactc, info.ra, info.dec, 2000.0, gl, gb, 1, /degree
          ebv_mw = dust_getval(gl,gb,/interp)

          kl = k_lambda(objwave,/odonnell,R_V=3.1)

          objflux = objflux * 10^(0.4*kl*ebv_mw)
          objferr = objferr * 10^(0.4*kl*ebv_mw)

          if keyword_set(doplot) then begin

             xrange = minmax(objwave)
;            xrange = [7000,max(objwave)]
;            xrange = [min(objwave),5200]
             
             ploterror, objwave, objflux, objferr, ps=10, xsty=3, ysty=3, position=pos[*,0], /normal, $
;            djs_plot, objwave, objflux, ps=10, xsty=3, ysty=3, position=pos[*,0], /normal, $
               thick=2.0, xthick=2.0, ythick=2.0, xtickname=replicate(' ',10), $
               ytitle='Object [Counts]', charsize=1.5, charthick=2.0, xrange=xrange
             legend, objlist[iobj], /left, /top, box=0, charsize=1.5, charthick=2.0
             if (nneg ne 0L) then djs_oplot, objwave[neg], objflux[neg], ps=4, color='red'
             
             djs_plot, objwave, skyflux, ps=10, xsty=3, ysty=3, position=pos[*,1], /normal, $
               /noerase, xthick=2.0, ythick=2.0, xtickname=replicate(' ',10), $
               ytitle='Sky [Counts]', charsize=1.5, charthick=2.0, xrange=xrange
             djs_plot, objwave, objflux/objferr, ps=10, xsty=3, ysty=3, position=pos[*,2], /normal, $
               /noerase, xthick=2.0, ythick=2, ytitle='S/N', charsize=1.5, charthick=2.0, $
               xtitle='Observed Wavelength [\AA]', xrange=xrange

             cc = get_kbrd(1)
             
          endif

; update the object header; in some object headers EXTEND appears
; twice: remove the second occurance

          mask = strmid(zmatch.mask,3,1,/reverse)
          quadrant = strmid(zmatch.mask,1,2,/reverse)
          
          airmass = (float(info.tel_airm_start) + float(info.tel_airm_end)) / 2.0
          
;         print, objlist[iobj]+' '+galaxy+' '+mask+' '+quadrant, gain, rdnoise
          
;         mkhdr, header, objflux, /extend
;         sxdelpar, header, 'COMMENT'

          header = hobj[where(strcompress(hobj,/remove) ne '')]
          keep = where(strmatch(header,'*EXTEND*',/fold) eq 0B)
          header = header[keep]

          sxaddpar, header, 'EXTEND', 'T', ' File may contain extensions', after='NAXIS1'
          sxaddpar, header, 'RA', repstr(strtrim(im_dec2hms(info.ra/15.0D),2),' ',':')
          sxaddpar, header, 'DEC', repstr(strtrim(im_dec2hms(info.dec),2),' ',':')
          sxaddpar, header, 'GALAXY', galaxy
          sxaddpar, header, 'MASK', mask, ' mask number'
          sxaddpar, header, 'QUADRANT', quadrant, ' quadrant number'
          sxaddpar, header, 'SLIT', zmatch.slit, ' slit number'
          sxaddpar, header, 'Z', float(zmatch.z), ' redshift'
          sxaddpar, header, 'Z_ERR', float(0.0), ' redshift error'
;         sxaddpar, header, 'Z_ERR', float(zmatch.zerr), ' redshift error' ; not available for all objects
          sxaddpar, header, 'Q', zmatch.q, ' redshift quality flag'
          sxaddpar, header, 'YOBJ', zmatch.yobj, ' object Y position '
          sxaddpar, header, 'YSKY', zmatch.ysky, ' sky Y position '
          sxaddpar, header, 'RA_VLT', zmatch.ra_vlt, ' VLT RA '
          sxaddpar, header, 'DEC_VLT', zmatch.dec_vlt, ' VLT DEC '
          sxaddpar, header, 'GAIN', gain, ' Gain [electron/ADU]'
          sxaddpar, header, 'RDNOISE', rdnoise, ' Read Noise [electron]'
          sxaddpar, header, 'AIRMASS', float(airmass), ' Mean airmass'
          sxaddpar, header, 'CRVAL1', minwave, ' wavelength at CRPIX1'
          sxaddpar, header, 'CRPIX1', 1.0, ' reference pixel number'
          sxaddpar, header, 'CDELT1', dwave, ' dispersion [Angstrom/pixel]'
          sxaddpar, header, 'CTYPE1', 'LINEAR', ' projection type'
          sxaddpar, header, 'CD1_1', dwave, ' dispersion [Angstrom/pixel]'

          if keyword_set(wfits) and (zmatch.z gt 0.0) and (zmatch.q ne 1L) then begin

             outfile = strtrim(objlist[iobj],2)
;            outfile = strtrim(repstr(objlist[iobj],'.fits','_sigspec.fits'),2)

             splog, 'Writing '+outpath+outfile+'.'
             mwrfits, objflux, outpath+outfile, header, /create
             mwrfits, objferr, outpath+outfile
;            spawn, ['gzip -f '+outpath+outfile], /sh
                
          endif

       endif
          
    endfor 

stop
    
return
end
    
