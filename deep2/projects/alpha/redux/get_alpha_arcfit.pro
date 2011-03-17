;+
; NAME:
;   GET_ALPHA_ARCFIT
;
; PURPOSE:
;   Gather information on the MIKE arc-line fitting.
;
; INPUTS: 
;   night - night of the observations (e.g., 'ut080414')
;
; OPTIONAL INPUTS: 
;   setup - see MIKE reduction code
;   side - see MIKE reduction code
;   outfile - output file name
;   linlist - linelist used during arc-line fitting
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Jan 04, UCSD
;-

pro get_alpha_arcfit, night, setup=setup, side=side, $
  outfile=outfile, linlist=linlist, alphapath=alphapath

    if (n_elements(alphapath) eq 0) then $
      alphapath = getenv('DEEP2_ALPHA_DIR')+'/'
    if (n_elements(setup) eq 0) then setup = 1
    if (n_elements(side) eq 0) then side = 2

; read the line-list used in X_ARCFIT    
    if (n_elements(linlist) eq 0) then linlist = $
      getenv('XIDL_DIR')+'/Spec/Arcs/Lists/hires_thar.lst'
    x_arclist, linlist, lines
    lines = lines[sort(lines.wave)]
    lineid = lindgen(n_elements(lines)) ; unique ID number (zero-indexed)

; loop on all the nights in this run    
    nnight = n_elements(night)
    for inight = 0, nnight-1 do begin
       datapath = alphapath+night[inight]+'/'
       splog, 'Working path ', datapath

       pushd, datapath
       mike = mike_ar('mike.fits')

       ordr_str = mike_getfil('ordr_str',setup,side=side)
       arcs = where(mike.type EQ 'ARC' AND mike.flg_anly NE 0 $
         AND mike.setup EQ setup AND mike.side EQ side,narc)
       thesearcs = strtrim(mike[arcs].img_root,2)

;      for iarc = 0, 0 do begin
       for iarc = 0, narc-1 do begin
; generate the arc name and the IDL save set name from the raw outfile
          arc_root = 'Raw/'+thesearcs[iarc]
          arc_fil = 'Arcs/Arc_m'+thesearcs[iarc]
          arc_saveset = 'Arcs/Fits/m'+repstr(thesearcs[iarc],'.fits','_fit.idl')
; get MJD
          hdr = xheadfits(arc_fil)
          date = strtrim(sxpar(hdr,'UT-DATE'),2)
          ut = sxpar(hdr,'UT-START')
          mjd = x_setjdate(date,ut)
; restore the IDL save set
          restore, arc_saveset
; loop on each order
          for ii = 0L, n_elements(ordr_str)-1L do begin
             if rejstr[ii].ngdf EQ 0L then continue
             these = lindgen(rejstr[ii].ngdf) ; lines used
             gdfit = rejstr[ii].gdfpt[these]
             if (rejstr[ii].nrej NE 0L) then $
               rejpt = rejstr[ii].rejpt[0:rejstr[ii].nrej-1] else $
                 rejpt = -1L
             fit = 10^x_calcfit(double(rejstr[ii].gdfpx[these]),FITSTR=all_arcfit[ii])
             goodbad = these*0B+1B ; 1=good, 0=bad
             if (rejpt[0] ne -1L) then goodbad[rejpt] = 0B
; pack into a structure
             info1 = {night: night[inight], mjd: mjd, arc: thesearcs[iarc], $
               order: -1, lineid: -1, pixel: 0.0, $
               lambda_true: 0.0D, lambda_fit: 0.0D, good: 0B}
             info1 = replicate(info1,n_elements(gdfit))
             info1.order = ordr_str[ii].order
             info1.pixel = rejstr[ii].gdfpx[these]
             info1.lambda_true = lines[gdfit].wave
             info1.lambda_fit = fit
             info1.lineid = lineid[gdfit]
             info1.good = goodbad
             if (inight eq 0) and (iarc eq 0) and (ii eq 0) then $
               info = info1 else info = [info,info1]
          endfor                ; close order
       endfor                   ; close arc lamp
       popd
    endfor                      ; close NIGHT

; write out
    im_mwrfits, info, outfile, /clobber

; now get some statistics on each unique line and write out
    allid = info.lineid
    uindx = uniq(allid,sort(allid))
    id = allid[uindx]
    nid = n_elements(id)
    stats = replicate({lineid: 0L, wave: 0.0D, nused: 0L, $
      nbad: 0L, ngood: 0L, mean: -1.0D, sigma: -1.0D, $
      order: 0},nid)
    for jj = 0L, nid-1L do begin
       these = where(id[jj] eq allid,nthese)
       stats[jj].nused = nthese
       stats[jj].lineid = id[jj]
       stats[jj].wave = info[these[0]].lambda_true
       stats[jj].order = info[these[0]].order

       good = where(info[these].good eq 1,ngood,ncomp=nbad)
       stats[jj].ngood = ngood
       stats[jj].nbad = nbad
       if (ngood ne 0) then begin
          diff = info[these[good]].lambda_fit-info[these[good]].lambda_true
          stats[jj].mean = djs_mean(diff)
          stats[jj].sigma = djsig(diff)
       endif
    endfor

; write out    
    statsfile = repstr(outfile,'.fits','_stats.fits')
    im_mwrfits, stats, statsfile, /clobber

return
end
    
