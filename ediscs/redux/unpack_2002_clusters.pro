pro unpack_2002_clusters, cluster, doplot=doplot, wfits=wfits
; jm04sep23uofa
; this routine assumes the following directory structure:

;   cluster/
;      domeflats/  
;         cluster_mask1/
;         cluster_mask2/
;      fluxed/ ; <-- this directory contains the fluxed 1D spectra
;      object_spectra_1d/  
;         cluster_mask1/
;         cluster_mask2/
;      sky_spectra_1d/
;         cluster_mask1/
;         cluster_mask2/
;      spec1d/ ; <-- this directory contains the telluriced/fluxed
;                    (final) 1D spectra

;   cluster = 'cl1040-1155'
;   cluster = 'cl1054-1146'
;   cluster = 'cl1054-1245'
    cluster = 'cl1216-1201'
;   cluster = 'cl1232-1250'
    
    ncluster = n_elements(cluster)
    if (ncluster eq 0L) then begin
       print, 'Syntax - '
       return
    endif

    splog, 'Cluster '+cluster & print
    
    rootpath = ediscs_path(/d2002)
    datapath = rootpath+cluster+'/'
    objpath = datapath+'object_spectra_1d/'
    skypath = datapath+'sky_spectra_1d/'
    domepath = datapath+'domeflats/'
    flatpath = datapath+'flatdivided/'
    fluxedpath = datapath+'fluxed/'
    spec1dpath = datapath+'spec1d/'
    zpath = datapath+'redshifts/'

; read the redshift data

    splog, 'Reading '+zpath+cluster+'.'
    zinfo = rsex(zpath+cluster+'.dat')

    splog, 'Reading '+zpath+cluster+'.'

    pagemaker, nx=1, ny=3, ymargin=[0.2,1.0], yspace=0.0, $
      xmargin=[1.3,0.2], position=pos, /normal

; move into the object directory to figure out how many masks there
; are; we assume that the appropriate directories also exist in
; SKYPATH and DOMEPATH

    splog, 'Pushing into '+objpath & print

    pushd, objpath
    masklist = file_search(cluster+'_mask?',count=nmask)
    
    splog, 'Found '+string(nmask,format='(I0)')+' masks:'
    niceprint, masklist & print

    popd
    
    if keyword_set(wfits) then begin
       splog, 'Delete all FITS files from: [Y/N]?'
       splog, '   '+flatpath
       splog, '   '+fluxedpath
       splog, '   '+spec1dpath
       cc = get_kbrd(1)
       if strupcase(cc) eq 'Y' then begin
          spawn, ['/bin/rm -f '+flatpath+'*.fits'], /sh
          spawn, ['/bin/rm -f '+fluxedpath+'*.fits'], /sh
          spawn, ['/bin/rm -f '+spec1dpath+'*.fits'], /sh
       endif
       print
    endif

; ###########################################################################
; PART 1: Parse the individual masks; generate the error spectra;
; normalize by the dome flat; write out to FLATDIVIDED
; ###########################################################################

    for imask = 0L, nmask-1L do begin

       pushd, objpath+masklist[imask]
       objlist = file_search('*.ms.*fits',count=objcount)
       popd

       pushd, domepath+masklist[imask]
       domelist = file_search('*.fits',count=domecount)
       popd

;      skylist = file_search(skypath+masklist[imask]+'/*.ms.*fits',count=skycount)

       if (objcount ne domecount) then begin
          splog, 'Dome flats missing.'
          return
       endif

       for iobj = 0L, objcount-1L do begin
       
; read the object, sky, and dome spectra

          cube = ediscs_rd1dspec(objlist[iobj],datapath=objpath+masklist[imask]+'/')
          hobj = headfits(objpath+masklist[imask]+'/'+objlist[iobj])
          
; match this spectrum to an object in the redshift data file generated
; by Bianca & company

          match = where((zinfo.cluster eq cube.cluster) and (string(zinfo.mask,format='(I1)') eq cube.mask) and $
            (string(zinfo.slit,format='(I2.2)') eq cube.slit) and (zinfo.subslit eq cube.subslit),nmatch)

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

             print
             splog, 'Match found:'

             zmatch = zinfo[match]
             
             matchinfo = replicate({cluster: '', mask: '', slit: '', subslit: ''},2)
             matchinfo.cluster = [cube.cluster,zmatch.cluster]
             matchinfo.mask = [cube.mask,zmatch.mask]
             matchinfo.slit = [cube.slit,zmatch.slit]
             matchinfo.subslit = [cube.subslit,zmatch.subslit]

             struct_print, matchinfo
             print

; read the data             

             objflux = cube.spec ; object spectrum [ADU]
             objwave = cube.wave
             gain = cube.out1gain ; gain [e/ADU]
             rdnoise = cube.out1ron ; read-noise [e]

             skyflux = total(readfits(skypath+masklist[imask]+'/'+cube.skyname,hsky,/silent),2) ; total sky spectrum [ADU]
             skywave = make_wave(hsky)

             domeflux = readfits(domepath+masklist[imask]+'/'+domelist[iobj],hdome,/silent)
             domewave = make_wave(hdome)

; compute the width of the object, sky, and dome flat apertures using
; the header information

             objapnum1 = sxpar(hobj,'APNUM1',count=apcount)
             if (apcount eq 0L) then begin
                splog, 'Improper object header!'
                return
             endif

             objaps = float(strsplit(objapnum1,' ',/extract))
             nobj = objaps[3]-objaps[2] ; object aperture [pixels]

             domeapnum1 = sxpar(hdome,'APNUM1',count=apcount)
             if (apcount eq 0L) then begin
                splog, 'Improper dome header!'
                return
             endif

             domeaps = float(strsplit(domeapnum1,' ',/extract))
             ndome = domeaps[3]-domeaps[2] ; dome flat aperture [pixels]

             skyapnum1 = sxpar(hsky,'APNUM1',count=apcount1)
             skyapnum2 = sxpar(hsky,'APNUM2',count=apcount2)
             if (apcount1 eq 0L) or (apcount2 eq 0L) then begin
                splog, 'Improper sky header!'
                return
             endif

             skyaps1 = float(strsplit(skyapnum1,' ',/extract))
             skyaps2 = float(strsplit(skyapnum2,' ',/extract))

             nsky1 = skyaps1[3]-skyaps1[2]
             nsky2 = skyaps2[3]-skyaps2[2]
             nsky = nsky1 + nsky2 ; total sky aperture [pixels]

; generate the error spectrum
             
             Nskyflux = skyflux / nsky ; sky flux per pixel [ADU/pixel]

             objfvar = sqrt(objflux^2) / gain + nobj * (1 + nobj/nsky) * $
               (sqrt(Nskyflux^2) / gain + (rdnoise/gain)^2) ; [ADU^2]

             objferr = sqrt(objfvar) ; [ADU]

             snrstats = im_stats(objferr/objflux)
             snr_median = string(snrstats.median,format='(F4.1)') ; median S/N

; generate the dome flat error spectrum

             domeferr = sqrt(domeflux/gain + ndome*(rdnoise/gain)^2) ; read noise is negligible [ADU]
             
; divide the object by the dome flat spectrum and propagate the error 

             interpdomeflux = interpol(domeflux,domewave,objwave)
             interpdomeferr = sqrt(interpol(domeferr^2,domewave,objwave))

             objdomeflux = objflux / interpdomeflux
             objdomeferr = sqrt( (objferr/interpdomeflux)^2 + interpdomeferr^2*(objflux/interpdomeflux^2)^2 )

             if keyword_set(doplot) then begin

                djs_plot, objwave, objflux, ps=10, xsty=3, ysty=3, position=pos[*,0], /normal, $
                  thick=2.0, xthick=2.0, ythick=2.0, xtickname=replicate(' ',10), $
                  ytitle='Object [Counts]', charsize=1.5, charthick=2.0
                legend, objlist[iobj], /left, /top, box=0, charsize=1.5, charthick=2.0
                
                djs_plot, objwave, skyflux, ps=10, xsty=3, ysty=3, position=pos[*,1], /normal, $
                  /noerase, xthick=2.0, ythick=2.0, xtickname=replicate(' ',10), $
                  ytitle='Sky [Counts]', charsize=1.5, charthick=2.0
                djs_plot, objwave, objflux/objferr, ps=10, xsty=3, ysty=3, position=pos[*,2], /normal, $
                  /noerase, xthick=2.0, ythick=2, ytitle='S/N', charsize=1.5, charthick=2.0, $
                  xtitle='Observed Wavelength [\AA]'

                cc = get_kbrd(1)
                
             endif

; update the object header

             header = hobj[where(strcompress(hobj,/remove) ne '')]

             sxaddpar, header, 'EXTEND', 'T'
             sxaddpar, header, 'CLUSTER', zmatch.cluster, ' galaxy cluster'
             sxaddpar, header, 'MASK', zmatch.mask, ' mask number'
             sxaddpar, header, 'SLIT', zmatch.slit, ' slit number'
             sxaddpar, header, 'SUBSLIT', zmatch.subslit, ' extraction number'
             sxaddpar, header, 'MEMBER', zmatch.member, ' cluster membership flag'
             sxaddpar, header, 'ID', zmatch.id, ' V2.1 photometric catalog ID'
             sxaddpar, header, 'Z', float(zmatch.z_final), ' redshift'
             sxaddpar, header, 'ZFLAG', zmatch.z_flag, ' redshift flag'
             sxaddpar, header, 'TYPE', zmatch.type, ' spectral type'
             sxaddpar, header, 'TYPEFLAG', zmatch.type_flag, ' type flag'

             sxaddpar, header, 'DOMEFLAT', domelist[iobj], 'Dome Flat File', before='HISTORY'
             sxaddhist, "'Dome flat normalization applied "+im_today()+"'", header

; write out the FLATDIVIDED spectra; use the original file names but
; write them to a common directory; the commented code below can be
; used to write out files with their ID numbers, without overwriting
; repeated objects
             
             if keyword_set(wfits) then begin

                outfile = strtrim(objlist[iobj],2)

;               rootfile = strtrim(zmatch.id,2)
;               if (rootfile eq '?') then begin 
;
;                  outfile = strtrim(objlist[iobj],2)
;
;               endif else begin
;                  
;                  files = file_search(flatpath+rootfile+'*',count=fcount)
;                  
;                  if (fcount eq 0L) then outfile = rootfile+'.fits'
;                  if (fcount eq 1L) then begin
;
;                     oldname = repstr(files[0],flatpath,'')
;                     newname = repstr(oldname,'.fits','_1.fits')
;                     spawn, ['/bin/mv '+flatpath+oldname+' '+flatpath+newname], /sh
;                     outfile = rootfile+'_2.fits'
;
;                  endif
;                  
;                  if (fcount gt 1L) then outfile = rootfile+'_'+string(fcount+1,format='(I0)')+'.fits'
;
;               endelse
                
                splog, 'Writing '+flatpath+outfile+'.'
                mwrfits, objdomeflux, flatpath+outfile, header, /create
                mwrfits, objdomeferr, flatpath+outfile
;               spawn, ['gzip -f '+flatpath+outfile], /sh
                print
                
             endif

          endif

       endfor

    endfor 

stop
    
; ###########################################################################
; PART 2: flux-calibrate
; ###########################################################################

; sensitivity function and telluric spectrum

    senspath = datapath
    sensname = 'sensitivity_2002.fits'
    extfile = 'extlasil.dat'

    pushd, flatpath
    speclist = file_search('*.fits',count=count)
    popd
    
    splog, 'Flux-calibrating '+string(count,format='(I0)')+' objects.'
    ediscs_calibrate, speclist, sensname, datapath=flatpath, $
      senspath=rootpath, outpath=fluxedpath, wfits=wfits

stop
    
; ###########################################################################
; PART 3: telluric-correct
; ###########################################################################

    tellfile = 'telluric_2002.fits'
    
    
             

stop    
    
return
end
    
