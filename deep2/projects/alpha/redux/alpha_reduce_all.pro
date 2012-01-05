pro alpha_reduce_all, night, preproc=preproc, blue=blue, fixheaders=fixheaders, $
  makelinks=makelinks, init=init, clobber=clobber, flat=flat, arc=arc, $
  slitflat=slitflat, proc=proc, emlines=emlines, dotrace=dotrace, skysub=skysub, $
  extract=extract, calibrate=calibrate, coadd=coadd, dostandards=dostandards, $
  makesens=makesens, stage1=stage1, stage2=stage2, stage3=stage3, $
  final_spec1d=final_spec1d, tarballs=tarballs, qaplot_arcfit=qaplot_arcfit, $
  combine=combine
; jm09jan06nyu - uber-wrapper to reduce *all* the alpha data; the
; default is to just reduce the RED data, unless /BLUE is set 

; ToDo [as of 09jan08]:
;  1) Fix all the sensitivity function stuff!  How can I combine the
;     standards observed over multiple nights?
;  2) Incorporate Dennis' preprocessing steps: median-filtering
;     plus reject cosmic rays
;  3) Include more intelligent coadding of the orders; for example, if
;     a line is near the edge of the order, then don't combine
;     it with the adjacent order
;  4) Include more intelligent coadding of repeat observations; in
;     particular, don't include crappy spectra when coadding!
;  5) Derive more accurate wavelength maps; build a more optimal
;     line-list 
;  6) Optimize the sky-subtraction; verify that the scattered-light
;     model is reasonable
;  7) Verify that the sky-line wavelength shifts are reasonable
    
;   mike_rslvall
    
    alphapath = getenv('DEEP2_ALPHA_DIR')+'/'
    spec1dpath = alphapath+'spec1d/'
    qaplotpath = alphapath+'qaplots/'
    if (n_elements(night) ne 0) then allnight = night else $
      allnight = [$
      'ut080414',$
      'ut080415',$
      'ut080416',$
      'ut080417',$
      'ut080918',$
      'ut080919',$
      'ut080920',$
      'ut080921']
    nnight = n_elements(allnight)

; custom line-list    
    linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/mike_thar_murphy_custom.lst'
;   linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/hires_thar.lst'
;   linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/mike_thar_alpha_custom.lst'
;   nycoeff = 4
;   nocoeff = 6

    sigrej_2darc = 5.0 ; for 2D arc fitting

;; --------------------------------------------------    
;; test code
;    ordrs = reverse(fix(im_array(37,70,1)))
;    qafile = qaplotpath+'qaplot_08apr14_arc0001.ps'
;    templfil = 'templ_arc_2x2R.idl'
;    mike_tweakarc, 'Arcs/Fits/mr0001_fit.idl', ordrs, $
;      templfil, linlist=linlist, qafil=qafile
;;stop    
;; --------------------------------------------------    
    
    setup = 1                   ; default setup
    if keyword_set(blue) then side = 1 else side = 2 ; red

; initial reductions
    if keyword_set(stage1) then begin
       clobber = 1 & flat = 1 & arc = 1 & slitflat = 1
    endif
    if keyword_set(stage2) then begin
       clobber = 1 & proc = 1 & emlines = 1
       dotrace = 1 & skysub = 1 & extract = 1
    endif
    if keyword_set(dostandards) then begin
       clobber = 1 & proc = 1 & dotrace = 1
       skysub = 1 & extract = 1
    endif
    if keyword_set(stage3) then begin
       calibrate = 1 & coadd = 1
    endif

; ##################################################
; data reduction stuff    
    for ii = 0, nnight-1 do begin
       datapath = alphapath+allnight[ii]+'/'
; ##################################################
; fix headers
       if keyword_set(fixheaders) then begin
          case allnight[ii] of
             'ut080417': begin ; #108 (sdss) was observed, not #128 (which is a deep2 galaxy)
                these = datapath+['rawdata/r006'+['0','2'],$
                  'rawdata/r0061']+'.fits.gz' 
                for tt = 0L, n_elements(these)-1L do begin
                   hh = headfits(these[tt])
                   obj = sxpar(hh,'OBJECT')
                   obj = repstr(obj,'128','108')
                   sxaddpar, hh, 'OBJECT', obj
                   spawn, 'gunzip -f '+these[tt], /sh
                   modfits, repstr(these[tt],'.gz',''), 0, hh
                   spawn, 'gzip -f '+repstr(these[tt],'.gz',''), /sh
                endfor
             end 
             'ut080919': begin ; [r,b]0065.fits mislabeled as a 'comp'
                these = datapath+'rawdata/r0065.fits.gz' 
                for tt = 0L, n_elements(these)-1L do begin
                   hh = headfits(these[tt])
                   sxaddpar, hh, 'EXPTYPE', 'Object'
                   spawn, 'gunzip -f '+these[tt], /sh
                   modfits, repstr(these[tt],'.gz',''), 0, hh
                   spawn, 'gzip -f '+repstr(these[tt],'.gz',''), /sh
                endfor
             end 
             else:
          endcase
          continue
       endif
; ##################################################
; do the pre-processing
       if keyword_set(preproc) then begin
          alpha_preproc, datapath
          continue
       endif
; ##################################################
; build the requisite soft links
       if keyword_set(makelinks) then begin
          alpha_make_softlinks, datapath
          continue
       endif
; ##################################################
; initialize all the directories and run the preprocessing required by
; the MIKE pipeline; this only should be run once!!
       if keyword_set(init) then begin
; now run the auto-id routine
          if file_test(datapath+'rawmike.fits',/regular) then $
            rawmike = mike_ar(datapath+'rawmike.fits') else begin
             pushd, datapath & mike_strct, rawmike, /noedit, $
               outfil=datapath+'rawmike.fits' & popd
          endelse
; each night needs a customized routine to edit the initial MIKE
; structure; for example, to fix headers, reject crappy spectra, etc.
          case allnight[ii] of
             'ut080414': mike = edit_mike_ut080414(rawmike)
             'ut080415': mike = edit_mike_ut080415(rawmike)
             'ut080416': mike = edit_mike_ut080416(rawmike)
             'ut080417': mike = edit_mike_ut080417(rawmike)
             'ut080918': mike = edit_mike_ut080918(rawmike)
             'ut080919': mike = edit_mike_ut080919(rawmike)
             'ut080920': mike = edit_mike_ut080920(rawmike)
             'ut080921': mike = edit_mike_ut080921(rawmike)
             else: stop
          endcase
; run some final setup steps and then write out 
          mike_setup, mike, outfil=datapath+'mike_summ.txt'
          mike_wrstrct, mike, fits=datapath+'mike.fits', $
            outfil=datapath+'mike.list'
          continue
       endif 
;; ##################################################
;; build a QAplot of the arc-line fitting
;       if keyword_set(qaplot_arcfit) then begin
;          pushd, datapath
;          mike = mike_ar('mike.fits')
;          qaplot_alpha_arcfit_night, mike, setup, side, linlist=linlist, $
;            qafile=qaplotpath+'qaplot_arcfit_'+allnight[ii]+'.ps'
;          popd
;       endif
       
; ##################################################
; now reduce the night
       if (keyword_set(final_spec1d) eq 0) and (keyword_set(tarballs) eq 0) and $
         (keyword_set(qaplot_arcfit) eq 0) then begin
          splog, 'Working directory '+datapath
          alpha_reduce_night, datapath, setup=setup, side=side, clobber=clobber, $
            flat=flat, arc=arc, slitflat=slitflat, proc=proc, emlines=emlines, $
            dotrace=dotrace, skysub=skysub, extract=extract, calibrate=calibrate, $
            coadd=coadd, dostandards=dostandards, makesens=makesens, linlist=linlist, $
            nycoeff=nycoeff, nocoeff=nocoeff, sigrej_2darc=sigrej_2darc, combine=combine
       endif
    endfor ; close NIGHT

;; ##################################################
;; build a QAplot of the arc-line fitting for each run separately (see
;; ALPHA_DEVELOP_LINELIST for proper use of these routines)
;    if keyword_set(qaplot_arcfit) then begin
;       arcfile = qaplotpath+'arcfit_08apr.fits'
;       get_alpha_arcfit, ['ut080414','ut080415','ut080416','ut080417'], $
;         outfile=arcfile, linlist=linlist, setup=setup, side=side
;
;       qafile = qaplotpath+'qaplot_arcfit_08apr.ps'
;       qaplot_alpha_arcfit, arcfile+'.gz', qafile=qafile
;    endif

stop    
    
; ##################################################
; build the final 1D spectra
    if keyword_set(final_spec1d) then begin
       alpha_final_spec1d, alphapath+allnight+'/spec1d/', spec1dpath, fluxed=0
;      alpha_final_spec1d, alphapath+allnight+'/spec1d/', spec1dpath, fluxed=1
    endif
       
; ##################################################
; make tarballs
    if keyword_set(tarballs) then begin
; 1d spectra from the individual nights
       for ii = 0L, nnight-1L do begin
          tarname = alphapath+'tar/'+allnight[ii]+'_spec1d.tar.gz'
          pushd, alphapath+allnight[ii]+'/spec1d/'
          flist = [file_search('*.fits'),file_search('*.ps.gz')]
          spawn, 'tar czvf '+tarname+' '+strjoin(flist,' '), /sh
          popd
       endfor
; final coadded spectra
       tarname = alphapath+'tar/final_spec1d.tar.gz'
       pushd, alphapath+'spec1d/'
       flist = [file_search('*.fits'),file_search('*.ps.gz')]
       spawn, 'tar czvf '+tarname+' '+strjoin(flist,' '), /sh
       popd
    endif
       
stop    

return
end
    
