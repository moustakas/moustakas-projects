pro alpha_develop_linelist, build=build, doit=doit, clobber=clobber, $
  solve=solve, stats=stats
; jm10jan28ucsd - develop the custom linelist

    setup = 1
    side = 2
    niter = 4 ; number of iterations
    
    alphapath = getenv('DEEP2_ALPHA_DIR')+'/'
    devpath = getenv('DEEP2_ALPHA_DIR')+'/develop_linelist/'
    lindir = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/'

    allnight = [$
      'ut080414',$
      'ut080415',$
      'ut080416',$
      'ut080417']
;     'ut080918',$
;     'ut080919',$
;     'ut080920',$
;     'ut080921']
    nnight = n_elements(allnight)

; build the files
    if keyword_set(build) then begin
       for iter = 0, niter-1 do begin
          iterpath = devpath+'iter'+string(iter+1,format='(I2.2)')+'/'
          if (file_test(iterpath,/dir) eq 0) then $
            spawn, 'mkdir -p '+iterpath, /sh
          if (file_test(iterpath+'qaplots/',/dir) eq 0) then $
            spawn, 'mkdir -p '+iterpath+'qaplots/', /sh
          for ii = 0, nnight-1 do begin
             inpath = alphapath+allnight[ii]+'/'
             outpath = iterpath+allnight[ii]+'/'
; base paths
             if (file_test(outpath,/dir) eq 0) then $
               spawn, 'mkdir -p '+outpath, /sh
             if (file_test(outpath+'Raw/',/dir) eq 0) then $
               spawn, 'mkdir -p '+outpath+'Raw/', /sh
             if (file_test(outpath+'Flats/',/dir) eq 0) then $
               spawn, 'mkdir -p '+outpath+'Flats/', /sh
             if (file_test(outpath+'Arcs/',/dir) eq 0) then $
               spawn, 'mkdir -p '+outpath+'Arcs/', /sh
             if (file_test(outpath+'Arcs/Fits/',/dir) eq 0) then $
               spawn, 'mkdir -p '+outpath+'Arcs/Fits/', /sh
             if (file_test(outpath+'Arcs/TRC/',/dir) eq 0) then $
               spawn, 'mkdir -p '+outpath+'Arcs/TRC/', /sh
             if (file_test(outpath+'Arcs/TRC/',/dir) eq 0) then $
               spawn, 'mkdir -p '+outpath+'Arcs/TRC/', /sh
; arcs
             spawn, 'cp -f '+inpath+'/mike.fits '+outpath+'/mike.fits', /sh
             mike = mrdfits(inpath+'/mike.fits',1)
             arc = where(strtrim(mike.type,2) eq 'ARC',narc)
             for jj = 0, narc-1 do spawn, 'ln -sf ../../../../'+allnight[ii]+'/Raw/'+$
               mike[arc[jj]].img_root+'.gz '+outpath+$
               '/Raw/'+mike[arc[jj]].img_root+'.gz', /sh
; flats
             allflats = file_search(inpath+'Flats/*',count=nflats)
             for jj = 0, nflats-1 do spawn, 'ln -sf ../../../../'+allnight[ii]+$
               '/Flats/'+file_basename(allflats[jj])+' '+outpath+'Flats/'+$
               file_basename(allflats[jj]), /sh
             
          endfor 
       endfor
    endif

; --------------------------------------------------    
; do it!
;   for iter = 1, 1 do begin
    for iter = 0, niter-1 do begin
       iterpath = devpath+'iter'+string(iter+1,format='(I2.2)')+'/'
       qaplotpath = iterpath+'qaplots/'
       allstatsfile = qaplotpath+'arcfit_08apr.fits'
       statsfile = repstr(allstatsfile,'.fits','_stats.fits.gz')
       qafile = qaplotpath+'qaplot_arcfit_08apr.ps'
; input and output linelists       
       if (iter eq 0) then spawn, '/bin/cp -f '+lindir+'alpha_thar_murphy.lst '+$
         lindir+'alpha_thar_custom_'+string(iter+1,format='(I2.2)')+'.lst'
       linlist = lindir+'alpha_thar_custom_'+string(iter+1,format='(I2.2)')+'.lst'
       outlinlist = lindir+'alpha_thar_custom_'+string(iter+2,format='(I2.2)')+'.lst'
; process the data if requested
;      if keyword_set(solve) then begin
;      for ii = 0, 0 do begin
       for ii = 0, nnight-1 do begin
          datapath = iterpath+allnight[ii]+'/'
          pushd, datapath
          mike = mike_ar('mike.fits')
          mike_allarc, mike, setup, side, fits=datapath+'mike.fits', $
            clobber=clobber, chk=chk, linlist=linlist, nycoeff=nycoeff, $
            nocoeff=nocoeff
          popd
       endfor
;      endif
; collate the results and write out statistics
;      if keyword_set(stats) then begin
       get_alpha_arcfit, allnight, setup=setup, side=side, $
         linlist=linlist, outfile=allstatsfile, $
         alphapath=iterpath
       qaplot_alpha_arcfit, allstatsfile+'.gz', qafile=qafile
;      endif
; read the statistics file and build a new linelist
       if (iter lt niter-1) then begin
          splog, 'Reading '+file_basename(statsfile)+$
            ' and '+file_basename(linlist)
          x_arclist, linlist, lines
       
          stats = mrdfits(statsfile,1,/silent)
          rej = where((stats.nbad/float(stats.nused) gt 0.5),$
            nrej,comp=keep,ncomp=nkeep)
;         niceprint, stats[rej].wave, lines[stats[rej].lineid].wave
          splog, 'Rejecting '+string(nrej,format='(I0)')+'/'+$
            string(n_elements(lines),format='(I0)')+' lines'
       
          good = lindgen(n_elements(lines))
          remove, stats[rej].lineid, good
          outlines = lines[good]

          txt = djs_readlines(linlist)
          hdr = txt[where(strmatch(txt,'*#*'))]
          openw, lun, outlinlist, /get_lun
          for jj = 0, n_elements(hdr)-1 do printf, lun, hdr[jj]
          for jj = 0, n_elements(outlines)-1 do printf, lun, $
            outlines[jj].wave, outlines[jj].flg_qual, $
            outlines[jj].name, format='(F12.6,1x,I1,1x,A0)'
          free_lun, lun
       end
;      if iter eq 1 then stop
    endfor ; close ITER

return
end
    
