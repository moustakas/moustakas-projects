;+
; NAME:
;   deep_vmax
; PURPOSE:
;   Calculate VMAX values for DEEP2
; CALLING SEQUENCE:
;   deep_vmax
; COMMENTS:
;   Returns R^3 in equation V = omega R^3/ 3
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------

; echo "deep2_oiilf_vmax" | idl > & deep2_oiilf_vmax.log &

pro deep2_oiilf_vmax, outstr, sdss=sdss, nocolor=nocolor, $
  nooii=nooii, test=test

    red, h100=0.7, omega0=0.3, omega_lambda=0.7
    h100 = redh100()

    numran = 12000L

    analysis_path = deep2_path(/analysis)
    outpath = deep2_path(/projects)+'oiilf/'
    
    fattempt = mrdfits(analysis_path+'deep2_completeness.fits',0,/silent)
    fgot = mrdfits(analysis_path+'deep2_completeness.fits',1,/silent)
    num = 15L
    rmilim = [0.0,2.0]
    rlim = [19.0,24.1]

    zbinsize = 0.01 ; for the completeness functions

    if (not keyword_set(sdss)) then begin
       alldeep1 = read_deep2(/kcorr)
       alldeep2 = read_deep2(/specdata)
       alldeep = struct_addtags(alldeep1,alldeep2)
       if keyword_set(test) then alldeep = alldeep[lindgen(500)+15000L]
;      if keyword_set(test) then $
;        alldeep = (read_deep2_zcat(/good))[10000:10050] else $
;        alldeep = read_deep2_zcat(/good)
       splog, 'Computing K-corrections for the full sample.'
       kcorig = deep_kcorrect(alldeep.zhelio,zcat=alldeep,$
         coeffs=allcoeffs,absmag=allabsmag,/sdss,/silent)
       nallgal=n_elements(alldeep)
;; compute the continuum flux at 3727 A
;       k_load_vmatrix, vmatrix, lambda, vname=vname
;       restwave = k_lambda_to_centers(lambda)
;       cflux_3727 = fltarr(nallgal)
;       for ii = 0L, nallgal-1L do cflux_3727[ii] = $
;         interpol(reform(vmatrix#allcoeffs[*,ii]),restwave,3727.0) ; [erg/s/cm2/A]
    endif else begin
       message, 'Not supported yet!'
       post = mrdfits(getenv('EVOLUTION_DIR')+'/data/deep2noev.fits',1)
       cat = mrdfits(getenv('EVOLUTION_DIR')+'/data/deep2noev.fits',2)
       bri = mrdfits(getenv('EVOLUTION_DIR')+'/data/deep2noev.fits',3)
       err = replicate(0.1,3,n_elements(bri))
       kcorig = deep_kcorrect(bri.z,mag=bri.bri,$
         err=err,coeffs=coeffs,absmag=absmag,/sdss)
       ngal=n_elements(bri)
       deep=bri
    endelse

    zlo = [0.7,0.9,1.1,1.3,0.7,0.8,1.0,0.8]
    zhi = [0.9,1.1,1.3,1.5,1.5,1.0,1.2,1.2]
;   zlo = [0.7] & zhi = [0.9]

    for iz = 0L, n_elements(zlo)-1L do begin

       zlo1 = zlo[iz]
       zhi1 = zhi[iz]

       zrange = where((alldeep.z gt zlo1) and $
         (alldeep.z lt zhi1),ngal)
       deep = alldeep[zrange]
       coeffs = allcoeffs[*,zrange]
       absmag = allabsmag[*,zrange]

       vmax = fltarr(ngal)

; initialize the completeness function array

       nzbins = fix((zhi1-zlo1)/zbinsize)+1L
       completeness = fltarr(nzbins,ngal)

; store NQUANT (quantile) statistics on the model colors and
; magnitudes 
       
       quantile = [0.05,0.16,0.25,0.5,0.75,0.84,0.95]
       nquant = n_elements(quantile)
       bmr_stats  = fltarr(nquant,ngal)-999.0 
       rmi_stats  = fltarr(nquant,ngal)-999.0
       rmag_stats = fltarr(nquant,ngal)-999.0
       
; now loop through each object       
       
;      for igal = 9L, ngal-1L do begin
       for igal = 0L, ngal-1L do begin

          if ((igal mod 10) eq 0) then splog, igal

; uniformly distribute NUMRAN realizations of this objects in volume
; between ZLO1 and ZHI1
          
          vinner = (lf_comvol(zlo1))[0]
          vouter = (lf_comvol(zhi1))[0]
          vol = vinner+(vouter-vinner)*randomu(seed,numran)
          mcz = lf_vtoz(vol) ; Volume-->z

; compute K-corrections from observed BRI to NUV,u,g
          
          curr_coeffs=fltarr(5,n_elements(vol))
          for j = 0L, 4L do curr_coeffs[j,*]=coeffs[j,igal]

          kc = deep_kcorrect(mcz, coeffs=curr_coeffs, /sdss, /silent)
          dm = lf_distmod(mcz)
          appm = fltarr(3,n_elements(vol))
          for j = 0L, 2L do appm[j,*] = absmag[j,igal]+kc[j,*]+dm ; apparent BRI
          bmr = transpose(appm[0,*]-appm[1,*]) ; B-R color
          rmi = transpose(appm[1,*]-appm[2,*]) ; R-I color
          rmag = transpose(appm[1,*])          ; R mag
          
; compute the color/magnitude weights for each object and check to see
; if the galaxy passes the selection criteria
          
          irmi = long((rmi-rmilim[0])/(rmilim[1]-rmilim[0])*float(num))
          ir = long((rmag-rlim[0])/(rlim[1]-rlim[0])*float(num))
          chances = fattempt[ir,irmi]*fgot[ir,irmi]
          iout = where((irmi ge num) or (irmi lt 0) or (ir ge num) or (ir lt 0),nout)
          if (nout gt 0) then chances[iout] = 1.0 ; throw out

          if (not keyword_set(nocolor)) then begin
             ikeep = where(deep2_color_cuts(appm[0,*],appm[1,*],appm[2,*]) gt 0 and $
               randomu(seed,n_elements(vol)) lt chances,nkeep,comp=toss)
          endif else begin
             message, 'Not supported yet!'
             ikeep=where(appm[1,*] lt 24.1 and $
               randomu(seed,n_elements(vol)) lt chances,nkeep,comp=toss)
          endelse

; now check which of the mock galaxies also pass the [O II] selection
                   
          if (deep[igal].oii_3727_1_ew[1] eq -2.0) or $
            (deep[igal].oii_3727_2_ew[1] eq -2.0) then begin

             ikeep = -1L
             nkeep = 0L

          endif else begin

             if (nkeep gt 0L) and (not keyword_set(nooii)) then begin
                
                oii_3727_1_limit = deep[igal].oii_3727_1_ew_limit
                oii_3727_2_limit = deep[igal].oii_3727_2_ew_limit

                oii_3727_1_mock = deep[igal].oii_3727_1_ew[0]/(1.0+mcz[ikeep])
                oii_3727_2_mock = deep[igal].oii_3727_2_ew[0]/(1.0+mcz[ikeep])
                
                oii_keep = where((oii_3727_1_mock gt oii_3727_1_limit) and $
                  (oii_3727_2_mock gt oii_3727_2_limit),noii_keep)
                if (noii_keep ne 0L) then begin
                   ikeep_cmc = ikeep ; after color/magnitude cuts
                   nkeep_cmc = nkeep
                   ikeep = ikeep[oii_keep]
                endif
                nkeep = noii_keep

;               yrange_1 = [min(oii_3727_1_mock)<oii_3727_1_limit,max(oii_3727_1_mock)>oii_3727_1_limit]
;               yrange_2 = [min(oii_3727_2_mock)<oii_3727_2_limit,max(oii_3727_2_mock)>oii_3727_2_limit]
;               yrange = [yrange_1[0]<yrange_2[0],yrange_1[1]>yrange_2[1]]
;
;               pagemaker, nx=1, ny=2, position=pos, /normal
;               djs_plot, [0], [0], /nodata, xrange=[zlo1,zhi1], yrange=yrange, $
;                 xsty=3, ysty=3, charsize=2, charthick=2, position=pos[*,0], $
;                 xtickname=replicate(' ',10)
;               djs_oplot, mcz[ikeep], oii_3727_1_mock[ikeep], ps=4, color='cyan'
;               djs_oplot, mcz[toss], oii_3727_1_mock[toss], ps=3, color='green'
;               djs_oplot, !x.crange, oii_3727_1_limit*[1,1], line=2, thick=2
;
;               djs_plot, [0], [0], /nodata, xrange=[zlo1,zhi1], yrange=yrange, $
;                 xsty=3, ysty=3, charsize=2, charthick=2, position=pos[*,1], /noerase, $
;                 xtitle='Redshift'
;               djs_oplot, mcz[ikeep], oii_3727_2_mock[ikeep], ps=4, color='cyan'
;               djs_oplot, mcz[toss], oii_3727_2_mock[toss], ps=3, color='green'
;               djs_oplot, !x.crange, oii_3727_2_limit*[1,1], line=2, thick=2
;               cc = get_kbrd(1)

             endif else noii_keep = 0L
                
          endelse 

; finally compute vmax          
          
          vmax[igal] = (vouter-vinner)*float(nkeep)/float(n_elements(vol)) ; compute VMAX
;         print, vmax[igal], nkeep, noii_keep
          
; calculate the redshift completeness function for this object

          if (nkeep gt 0L) then begin
             
             yall = im_hist1d(mcz,binsize=zbinsize,obin=xall,histmin=zlo1,histmax=zhi1,binedge=-1)
             ykeep = im_hist1d(mcz[ikeep],binsize=zbinsize,obin=xkeep,histmin=zlo1,histmax=zhi1,binedge=-1)
             if (n_elements(xall) ne nzbins) then message, 'Problem here!'
             
             cgood = where(yall gt 0.0,ncgood)
             if (ncgood ne 0L) then completeness[cgood,igal] = float(ykeep[cgood])/float(yall[cgood])

;             djs_plot, xall, yall/float(numran), xsty=1, ysty=1, ps=10, thick=2.0, xthick=2.0, ythick=2.0, $
;               charsize=2.0, charthick=2.0, xtitle='Redshift', ytitle='Fraction'
;             djs_oplot, xkeep, ykeep/float(numran), ps=10, thick=2.0, line=2
;             
;;            niceprint, xall[good], xkeep[good], ykeep[good]/yall[good]
;             djs_plot, xall[good], ykeep[good]/yall[good], line=0, thick=2, ps=10, xthick=2.0, ythick=2.0, xsty=1, ysty=1, $
;               charsize=2.0, charthick=2.0, xtitle='Redshift', ytitle='Completeness', yrange=[0,1.05], xrange=[zlo1,zhi1]
;             cc = get_kbrd(1)

; store statistics on the model colors

             bmr_stats[*,igal]  = weighted_quantile(bmr[ikeep],quant=quantile)
             rmi_stats[*,igal]  = weighted_quantile(rmi[ikeep],quant=quantile)
             rmag_stats[*,igal] = weighted_quantile(rmag[ikeep],quant=quantile)
             
          endif
             
       endfor 

       outstr1 = create_struct(alldeep[0],'vmax', -999.0, 'absmag', fltarr(3)-999.0, $
         'completeness', fltarr(nzbins), 'bmr_stats', fltarr(nquant)-999.0, $
         'rmi_stats', fltarr(nquant)-999.0, 'rmag_stats', fltarr(nquant)-999.0)
       outstr = replicate(outstr1,nallgal)
       struct_assign, alldeep, outstr

       outstr[zrange].vmax = vmax/h100^3.0 ; h=1-->h=0.7
       outstr[zrange].absmag = absmag + 5.0*alog10(h100) ; [NUV,u,g] h=1-->h=0.7
       outstr[zrange].completeness = completeness
       outstr[zrange].bmr_stats = bmr_stats
       outstr[zrange].rmi_stats = rmi_stats
       outstr[zrange].rmag_stats = rmag_stats

       if (NOT keyword_set(sdss)) then begin
          outname = 'deep2_vmax_'+string(zlo1,format='(F3.1)')+'_'+$
            string(zhi1,format='(F3.1)')+'.fits'
       endif else begin
          message, 'Not supported yet!'
          outname = 'deep2noev_vmax.fits'
       endelse

       if (keyword_set(nocolor)) then outname='noc_'+outname

;      splog, 'Writing '+outpath+outname
;      mwrfits, outstr, outpath+outname, /create

    endfor
      
return
end
