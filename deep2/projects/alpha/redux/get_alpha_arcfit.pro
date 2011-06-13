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
;   jm11jun08ucsd - updated to use the 2D fitting results
;-

pro get_alpha_arcfit, night, setup=setup, side=side, nycoeff=nycoeff, $
  nocoeff=nocoeff, outfile=outfile, linlist=linlist, alphapath=alphapath

    if (n_elements(alphapath) eq 0) then $
      alphapath = getenv('DEEP2_ALPHA_DIR')+'/'
    if (n_elements(setup) eq 0) then setup = 1
    if (n_elements(side) eq 0) then side = 2
    if (n_elements(sigrej) eq 0) then sigrej = 6.0

; read the line-list used in X_ARCFIT    
    if (n_elements(linlist) eq 0) then linlist = $
      getenv('XIDL_DIR')+'/Spec/Arcs/Lists/mike_thar_murphy.lst'
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
          arc_fil = 'Arcs/Arc_m'+thesearcs[iarc]
          arc_fit = 'Arcs/Fits/Arc_m'+repstr(thesearcs[iarc],'.fits','_fit2D.fits')
          arc_saveset = 'Arcs/Fits/m'+repstr(thesearcs[iarc],'.fits','_fit.idl')
; get MJD
          hdr = xheadfits(arc_fil)
          date = strtrim(sxpar(hdr,'UT-DATE'),2)
          ut = sxpar(hdr,'UT-START')
          mjd = x_setjdate(date,ut)
; restore the fitting results and the IDL save set
          restore, arc_saveset
          fit = mrdfits(arc_fit,1,/silent)
          nrm = fit.nrm
          nrmt = fit.nrmt
          nycoeff = fit.ny
          nocoeff = fit.no
          res = fit.res

; identify the problematic arc lines; the following code is all in
; x_fit2darc
          gd_ordr = where(sv_lines.nlin NE 0 AND ordr_str.flg_anly NE (-1), ngd_ordr)
          npix = round(total(sv_lines[gd_ordr].nlin))
          t = dblarr(npix)

          cnt = 0L
          for j=0L,ngd_ordr-1 do begin
             ;; Dummy indx
             i = gd_ordr[j]
             ;; PIX
             if j EQ 0 then all_pix = [sv_lines[i].pix[0:sv_lines[i].nlin-1]] $
             else all_pix = [all_pix,sv_lines[i].pix[0:sv_lines[i].nlin-1]]
             ;; WV
             if j EQ 0 then all_wv = $
               [sv_lines[i].wv[0:sv_lines[i].nlin-1]]*ordr_str[i].order $
             else all_wv = [all_wv,sv_lines[i].wv[0:sv_lines[i].nlin-1] $
               *ordr_str[i].order]
             ;; t
             ;; Order #
             t[cnt:cnt+sv_lines[i].nlin-1] = ordr_str[i].order
             cnt = cnt + sv_lines[i].nlin
          endfor

          ;; NORMALIZE PIX then ORDER
          pix_nrm = 2. * (all_pix - nrm[0])/nrm[1]
          t_nrm = 2. * (t - nrmt[0])/nrmt[1]
          invvar = replicate(1., npix)

;  Setup the Functions
          work2d = dblarr(npix,nycoeff*nocoeff)
          worky = flegendre(pix_nrm[*], nycoeff)
          workt = flegendre(t_nrm[*], nocoeff)
          
          for i=0,nocoeff-1 do for j=0,nycoeff-1 do $
            work2d[*,j*nocoeff+i] = worky[*, j] * workt[*,i]
          wv_mod = work2d # res

; find outliers          
          djs_iterstat, (wv_mod-all_wv), sigrej=sigrej, mask=msk
          gd = where(msk EQ 1B, complement=bad, ncomplement=nbad)

; pack into a structure
          info1 = replicate({night: night[inight], mjd: mjd, arc: thesearcs[iarc], $
            order: -1, lineid: -1, pixel: 0.0, lambda_true: 0.0D, lambda_fit: 0.0D, $
            good: 0B, allgood: 0B},npix)
          info1.order = t
          info1.pixel = all_pix
          info1.lambda_true = all_wv/t
          info1.lambda_fit = wv_mod/t
          info1.lineid = round(findex(lines.wave,info1.lambda_true))
          info1.good = msk ; 1=good, 0=bad (rejected)
          if (inight eq 0) and (iarc eq 0) then $
            info = info1 else info = [info,info1]

;; --------------------------------------------------          
;; the code below is all correct, but it's wrong in that it uses
;; the 1D fitting results from x_fitarc, not the final 2D fitting
;; results from x_fittrcarc
;
;; generate the arc name and the IDL save set name from the raw outfile
;          arc_root = 'Raw/'+thesearcs[iarc]
;          arc_fil = 'Arcs/Arc_m'+thesearcs[iarc]
;          arc_saveset = 'Arcs/Fits/m'+repstr(thesearcs[iarc],'.fits','_fit.idl')
;; get MJD
;          hdr = xheadfits(arc_fil)
;          date = strtrim(sxpar(hdr,'UT-DATE'),2)
;          ut = sxpar(hdr,'UT-START')
;          mjd = x_setjdate(date,ut)
;; restore the IDL save set
;          restore, arc_saveset
;; loop on each order
;          for ii = 0L, n_elements(ordr_str)-1L do begin
;             if rejstr[ii].ngdf EQ 0L then continue
;             these = lindgen(rejstr[ii].ngdf) ; lines used
;             gdfit = rejstr[ii].gdfpt[these]
;             if (rejstr[ii].nrej NE 0L) then $
;               rejpt = rejstr[ii].rejpt[0:rejstr[ii].nrej-1] else $
;                 rejpt = -1L
;             fit = 10^x_calcfit(double(rejstr[ii].gdfpx[these]),FITSTR=all_arcfit[ii])
;             goodbad = these*0B+1B ; 1=good, 0=bad
;             if (rejpt[0] ne -1L) then goodbad[rejpt] = 0B
;; pack into a structure
;             info1 = {night: night[inight], mjd: mjd, arc: thesearcs[iarc], $
;               order: -1, lineid: -1, pixel: 0.0, $
;               lambda_true: 0.0D, lambda_fit: 0.0D, good: 0B}
;             info1 = replicate(info1,n_elements(gdfit))
;             info1.order = ordr_str[ii].order
;             info1.pixel = rejstr[ii].gdfpx[these]
;             info1.lambda_true = lines[gdfit].wave
;             info1.lambda_fit = fit
;             info1.lineid = lineid[gdfit]
;             info1.good = goodbad
;             if (inight eq 0) and (iarc eq 0) and (ii eq 0) then $
;               info = info1 else info = [info,info1]
;          endfor                ; close order
;; --------------------------------------------------          
       endfor                   ; close arc lamp
       popd
    endfor                      ; close NIGHT

; compute GOOD using *all* the available nights and arcs
    djs_iterstat, info.lambda_true-info.lambda_fit, mask=msk, sigrej=sigrej
    info.allgood = msk
    
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

       good = where(info[these].allgood eq 1,ngood,ncomp=nbad)
;      good = where(info[these].good eq 1,ngood,ncomp=nbad)
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
    
