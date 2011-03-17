pro ediscs_coadd, inpath, outpath, gzip=gzip, wfits=wfits
; jm04feb16uofa

    pushd, inpath
;   flist = file_search('*.ms.fits.gz',count=fcount)
    flist = file_search('f*.ms.fits.gz',count=fcount)
    templist = flist
    forage = ediscs_headgrab(templist)
    popd

    unique = uniq(forage.galaxy)
    nunique = n_elements(unique)

    ugalaxy = forage[unique].galaxy
    outname = ugalaxy+'.ms.fits'

    for i = 0L, nunique-1L do begin

       match = where(ugalaxy[i] eq forage.galaxy,nmatch)

; read the base (first) spectrum       

       cube = rd1dspec(flist[match[0]],datapath=inpath,/silent)
       flux = cube.spec
       ferr = cube.sigspec
       ivar = 1.0/ferr^2.0
       header = cube.header
       wave = make_wave(header)
       
       if (nmatch gt 1L) then begin

          bigoutflux = fltarr(cube.npix,nmatch)
          bigoutivar = fltarr(cube.npix,nmatch)

          bigoutflux[*,0] = flux
          bigoutivar[*,0] = ivar
          
          for imatch = 1L, nmatch-1L do begin

; read the spectrum             
             
             nextcube = rd1dspec(flist[match[imatch]],datapath=inpath,/silent)
             nextflux = cube.spec
             nextferr = cube.sigspec
             nextivar = cube.sigspec
             nextheader = cube.header
             nextwave = make_wave(nextheader)

; interpolate onto the wavelength spacing of the zeroth spectrum  
             
             combine1fiber, alog10(nextwave), nextflux, nextivar, $
               newloglam=alog10(wave), newflux=newflux, newivar=newivar

             bigoutflux[*,imatch] = newflux
             bigoutivar[*,imatch] = newivar

          endfor

; coadd the spectra optimally

          combine1fiber, alog10(wave)#(fltarr(nmatch)+1), bigoutflux, bigoutivar, $
            newloglam=alog10(wave), newflux=outflux, newivar=outivar

; interpolate over zeros in the output variance spectrum

          mask = bytarr(cube.npix)
          zero = where(outivar eq 0.0,nzero)
          if (nzero ne 0L) then begin
             mask[zero] = 1B
             outivar = djs_maskinterp(outivar,mask)
          endif

          outferr = 1.0/sqrt(outivar)
          outhead = header
          sxaddpar, outhead, 'NCOMB', nmatch, ' number of co-added spectra'

       endif else begin

          outflux = flux
          outferr = ferr
          outhead = header
          
       endelse
       
       if keyword_set(wfits) then begin

          splog, 'Writing '+outpath+outname[i]+'.'
          mwrfits, outflux, outpath+outname[i], outhead, /create
          mwrfits, outferr, outpath+outname[i]
          if keyword_set(gzip) then spawn, ['gzip -f '+outpath+outname[i]], /sh

       endif
       
    endfor
    
return
end    
