pro alpha_develop_linelist, build=build, solve=solve, clobber=clobber
; jm10jan28ucsd - develop the custom linelist

    setup = 1
    side = 2
    niter = 6    ; number of iterations
    rejcut = 0.7 ; fraction of times bad to be rejected
    sigrej_2darc = 5.0 ; for 2D arc fitting
;   fordr = 41 ; final order to analyze (~8250 A)
    
    alphapath = getenv('DEEP2_ALPHA_DIR')+'/'
    devpath = getenv('DEEP2_ALPHA_DIR')+'/develop_linelist/'
    xidllindir = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/'
    lindir = devpath+'linelist/'

; initial and final linelists    
    initlinlist = xidllindir+'mike_thar_murphy.lst'    
    finallinlist = xidllindir+'mike_thar_murphy_custom.lst'
    
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
       if (file_test(devpath+'qaplots/',/dir) eq 0) then $
         spawn, 'mkdir -p '+devpath+'qaplots/', /sh
       for iter = 5, niter do begin
;      for iter = 1, niter do begin
          iterpath = devpath+'iter'+string(iter,format='(I2.2)')+'/'
          if (file_test(iterpath,/dir) eq 0) then $
            spawn, 'mkdir -p '+iterpath, /sh
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
       return 
    endif

; --------------------------------------------------    
; do it!
; start with the Murphy linelist
    x_arclist, initlinlist, lines
    openw, lun, lindir+'alpha_thar_custom_01.lst', /get_lun
    printf, lun, '# Fine-structure constant project by J. Moustakas'
    printf, lun, '# Initial line-list: mike_thar_murphy.lst'
    printf, lun, '# Iteration 01'
    for jj = 0, n_elements(lines)-1 do printf, lun, lines[jj].wave, $
      lines[jj].flg_qual, lines[jj].name, format='(F12.6,1x,I1,1x,A0)'
    free_lun, lun

    for iter = 1, niter do begin
       splog, '##################################################'
       iterstr = 'iter'+string(iter,format='(I2.2)')
       splog, strupcase(iterstr)
       iterpath = devpath+iterstr+'/'
       qaplotpath = devpath+'qaplots/'
       qafile = qaplotpath+'qaplot_arcfit_08apr_'+iterstr+'.ps'
       allstatsfile = qaplotpath+'arcfit_08apr_'+iterstr+'.fits'
       statsfile = repstr(allstatsfile,'.fits','_stats.fits')
; input and output linelists
       linlist = lindir+'alpha_thar_custom_'+string(iter,format='(I2.2)')+'.lst'
       outlinlist = lindir+'alpha_thar_custom_'+string(iter+1,format='(I2.2)')+'.lst'
       splog, 'Input linelist '+linlist
; process the data if requested
       if keyword_set(solve) then begin
          for ii = 0, nnight-1 do begin
             splog, '#########################'
             splog, 'Working on NIGHT '+allnight[ii]
             datapath = iterpath+allnight[ii]+'/'
             pushd, datapath
             mike = mike_ar('mike.fits')
             mike_allarc, mike, setup, side, fits=datapath+'mike.fits', $
               clobber=clobber, chk=chk, linlist=linlist, nycoeff=nycoeff, $
               nocoeff=nocoeff, sig2drej=sigrej_2darc
             popd
          endfor
       endif

; collate the results and write out statistics
       get_alpha_arcfit, allnight, setup=setup, side=side, linlist=linlist, $
         outfile=allstatsfile, alphapath=iterpath
       qaplot_alpha_arcfit, allstatsfile+'.gz', rejcut=rejcut, $
         iter=iter, qafile=qafile
       
; read the statistics file and build a new linelist
       if (iter lt niter) then begin
          splog, 'Reading '+file_basename(statsfile)+$
            ' and '+file_basename(linlist)
          x_arclist, linlist, lines
          
          stats = mrdfits(statsfile+'.gz',1,/silent)
          rej = where((stats.nbad/float(stats.nused) gt rejcut),$ ; reject
            nrej,comp=keep,ncomp=nkeep)
          
          if (nrej gt 0) then begin
;            niceprint, stats[rej].wave, lines[stats[rej].lineid].wave
             splog, 'Rejecting '+string(nrej,format='(I0)')+'/'+$
               string(n_elements(lines),format='(I0)')+' lines'
             good = lindgen(n_elements(lines))
             remove, stats[rej].lineid, good
             outlines = lines[good]
             
             splog, 'Writing '+outlinlist
             openw, lun, outlinlist, /get_lun
             txt = djs_readlines(linlist)
             printf, lun, '# Fine-structure constant project by J. Moustakas'
             printf, lun, '# Initial line-list: mike_thar_murphy.lst'
             printf, lun, '# Iteration '+string(iter+1,format='(I2.2)')
             for jj = 0, n_elements(outlines)-1 do printf, lun, $
               outlines[jj].wave, outlines[jj].flg_qual, $
               outlines[jj].name, format='(F12.6,1x,I1,1x,A0)'
             free_lun, lun
          endif else begin
             splog, 'No lines rejected! - the line-list has converged!!'
             spawn, '/bin/cp -f '+linlist+' '+finallinlist, /sh
             return
          endelse
       endif 
    endfor ; close ITER

return
end
    
