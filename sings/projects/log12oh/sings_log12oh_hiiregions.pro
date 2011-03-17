function log12oh_gradient, hiiinfo, hiigalaxy, method=method, debug=debug
; compute the abundance gradient for a given set of abundances

;   oh_minerr = 0.0 ; abundance error added in quadrature to the observed error
    oh_minerr = 0.05 ; abundance error added in quadrature to the observed error
    rr25_char = 0.4  ; metallicity at R/R25=0.4 ("characteristic" metallicity)
    rr25_cent = 0.1  ; metallicity at R/R25=0.1
    
    branchtag = tag_indx(hiigalaxy,'r23_branch_'+method)
    ohtag = tag_indx(hiigalaxy,'log12oh_'+method)
    rr25avgtag = tag_indx(hiiinfo,'hii_'+method+'_rr25_avg')
    
    good = where(((hiigalaxy.(ohtag))[0,*] gt -900.0) and $
      (hiigalaxy.rc3_rr25 gt -900.0),ngood)
    if (ngood eq 0) then message, 'This should not happen'

    branch = hiigalaxy[good].(branchtag)
    oh = (hiigalaxy[good].(ohtag))[0,*]
    oh_err = (hiigalaxy[good].(ohtag))[1,*]
    rr25 = hiigalaxy[good].rc3_rr25
    rr25_err = rr25*0.0
    p = hiigalaxy[good].p[0]

; do not fit HII regions with ambiguous branches
    these = where((branch ne 'A') and (branch ne 'X'),nthese)
    rej = where((branch eq 'A') or (branch eq 'X'))
;   if (method eq 'pt05') then begin
;      these = where((branch ne 'A') and (branch ne 'X') and (p gt 0.1),nthese,comp=rej)
;      rej = where((branch eq 'A') or (branch eq 'X') or (p lt 0.1))
;   endif else begin
;      these = where((branch ne 'A') and (branch ne 'X'),nthese)
;      rej = where((branch eq 'A') or (branch eq 'X'))
;   endelse
    oh_err[these] = sqrt(oh_err[these]^2 + oh_minerr^2)
;   if strmatch(hiiinfo.galaxy,'*0925*',/fold) then stop
;   ww = where(p lt 0.1)
;   if (ww[0] ne -1) then stop

    coeff = linfit(rr25[these],oh[these],sigma=coeff_err,$
      chisq=chisq,measure_errors=oh_err[these],yfit=ohmodel,$
      covar=covar)
    chisq = chisq/sqrt(nthese-2) ; reduced chi^2
    corr = covar[0,1]/sqrt(covar[0,0]*covar[1,1]) ; correlation coefficient

    hiiinfo.(tag_indx(hiiinfo,'hii_'+method+'_slope'))           = [coeff[1],coeff_err[1]]
    hiiinfo.(tag_indx(hiiinfo,'hii_'+method+'_log12oh_nuclear')) = [coeff[0],coeff_err[0]]
    hiiinfo.(tag_indx(hiiinfo,'hii_'+method+'_gradient_rms'))    = djsig(oh[these]-ohmodel)
    hiiinfo.(tag_indx(hiiinfo,'hii_'+method+'_gradient_chisq'))  = chisq
    hiiinfo.(tag_indx(hiiinfo,'hii_'+method+'_gradient_covar'))  = covar

; the central and characteristic abundances are defined by the
; best-fitting model, but to get the errors we have to use a Monte
; Carlo method; get NRAND realizations of the best-fitting model using
; the full covariance matrix and then compute the standard deviation
; at the desired radius
    nrand = 1000
    rand = mrandomn(seed,covar,nrand)
    rand[*,0] = rand[*,0] + coeff[0]
    rand[*,1] = rand[*,1] + coeff[1]

    ell = covar2ellipse(covar,nsigma=1.0) ; get models within 1-sigma
    indx = get_ellipse_indices(rand[*,0],rand[*,1],$
      major=ell.major,minor=ell.minor,angle=ell.angle, $
      xcenter=coeff[0],ycenter=coeff[1])
    nindx = n_elements(indx)

    char = fltarr(nindx)
    cent = fltarr(nindx)
    for ii = 0, nindx-1 do begin
       char[ii] = poly(rr25_char,rand[indx[ii],*])
       cent[ii] = poly(rr25_cent,rand[indx[ii],*])
    endfor

    hiiinfo.(tag_indx(hiiinfo,'hii_'+method+'_log12oh_central')) = $
      [poly(rr25_cent,coeff),djsig(cent)]
    hiiinfo.(tag_indx(hiiinfo,'hii_'+method+'_log12oh_char'))    = $
      [poly(rr25_char,coeff),djsig(char)]

;; some debugging plots    
;    window, 0
;    djs_plot, rand[*,0], rand[*,1], ps=6, xsty=3, ysty=3
;    tvellipse, ell.major, ell.minor, coeff[0], coeff[1], $
;      ell.angle, color=djs_icolor('blue'), /data
;    djs_oplot, rand[indx,0], rand[indx,1], ps=6, color='cyan'
;
;    window, 2
;    rr25axis = findgen(501)/100.0
;    ploterror, rr25[these], oh[these], oh_err[these], psym=6, xsty=3, ysty=3, $
;      title=hiiinfo.galaxy, xrange=[0,1], yrange=[8.6,9.5]
;    for ii = 0, nindx-1 do djs_oplot, rr25axis, poly(rr25axis,rand[indx[ii],*])
;    oploterror, rr25[these], oh[these], oh_err[these], psym=6, $
;      color=djs_icolor('orange'), errcolor=djs_icolor('orange')
;    djs_oplot, rr25axis, poly(rr25axis,coeff), line=0, color='blue'

; OLD, OBSOLETE LINEAR-FITTING CODE!    
;   res = sings_fit_gradient(rr25[these],rr25_err[these],$
;     oh[these],oh_err[these],nmonte=nmonte,xchar=rr25_char,$
;     xcent=rr25_cent,yminerr=oh_minerr)
;   hiiinfo.(tag_indx(hiiinfo,'hii_'+method+'_slope'))           = res.slope
;   hiiinfo.(tag_indx(hiiinfo,'hii_'+method+'_gradient_rms'))    = res.rms
;   hiiinfo.(tag_indx(hiiinfo,'hii_'+method+'_gradient_chisq'))  = res.chisq
;   hiiinfo.(tag_indx(hiiinfo,'hii_'+method+'_log12oh_nuclear')) = res.int
;   hiiinfo.(tag_indx(hiiinfo,'hii_'+method+'_log12oh_central')) = res.ycent
;   hiiinfo.(tag_indx(hiiinfo,'hii_'+method+'_log12oh_char'))    = res.ychar

; store the average R/R25 radius    
    hiiinfo.(rr25avgtag) = [djs_mean(rr25),djsig(rr25)]
    
    if keyword_set(debug) then begin
       rr25axis = findgen(501)/100.0
       djs_plot, [0], [0], /nodata, xsty=3, ysty=3, $
         title=hiiinfo.galaxy, xrange=[0,1], yrange=[7.6,9.6]
       for ii = 0, nindx-1 do djs_oplot, rr25axis, poly(rr25axis,rand[indx[ii],*])
       oploterror, rr25[these], oh[these], oh_err[these], psym=6, $
         symsize=1.2, color=djs_icolor('green'), errcolor=djs_icolor('green')
       if (rej[0] ne -1) then djs_oplot, rr25[rej], oh[rej], color='red', $
         psym=6, symsize=1.2
       djs_oplot, rr25axis, poly(rr25axis,coeff), line=0, color='orange'
;      ww = where(strmatch(strtrim(hiigalaxy[good].region,2),'-033-118',/fold))
;      ww = where(strmatch(strtrim(hiigalaxy[good].region,2),'+062-170',/fold))
;      ww = where(strmatch(strtrim(hiigalaxy[good].region,2),'S3A1',/fold))
;      ww = where(strmatch(strtrim(hiigalaxy[good].region,2),'S3A2',/fold))
;      ww = where(strmatch(strtrim(hiigalaxy[good].region,2),'-085+173',/fold))
;      ww = where(strmatch(strtrim(hiigalaxy[good].region,2),'-077+167',/fold))
;      plots, rr25[ww], oh[ww], psym=7, sym=2
       if (method eq 'pt05') then niceprint, rr25, oh, oh_err, hiigalaxy[good].region, $
         hiigalaxy[good].p, hiigalaxy[good].r23_branch_pt05
       cc = get_kbrd(1)
    endif
       
return, hiiinfo
end

function log12oh_avg, hiiinfo, hiigalaxy, thesehii, $
  method=method, weighted_avg=weighted_avg, lun=lun
; compute the average abundance, for the specified method

    branchtag = tag_indx(hiigalaxy,'r23_branch_'+method)
    uptag = tag_indx(hiigalaxy,'log12oh_'+method+'_upper')
    lotag = tag_indx(hiigalaxy,'log12oh_'+method+'_lower')
    ohtag = tag_indx(hiigalaxy,'log12oh_'+method)
    reftag = tag_indx(hiiinfo,'hii_'+method+'_texrefs')
    avgtag = tag_indx(hiiinfo,'hii_'+method+'_log12oh_avg')
    usetag = tag_indx(hiiinfo,'hii_'+method+'_nhii_used')

; assign metallicities according to the desired R23 branch    
    kk04 = 0
    pt05 = 0
    case method of
       'kk04': kk04 = 1
       'pt05': pt05 = 1
       else: message, 'Update me'
    endcase
    branch = sings_assign_r23branch(thesehii,$
      r23branch=strtrim(hiigalaxy.r23_branch,2),$
      /justflux,kk04=kk04,pt05=pt05,silent=0,$
      /debug,title=strtrim(hiiinfo.galaxy,2)+'/'+$
      strupcase(method))
    
    btag1 = tag_indx(branch,'r23branch_'+method)
    ohtag1 = tag_indx(branch,'zstrong_12oh_'+method)
    errtag1 = tag_indx(branch,'zstrong_12oh_'+method+'_err')

    hiigalaxy.(branchtag) = strtrim(branch.(btag1),2)
    hiigalaxy.(ohtag) = transpose([[branch.(ohtag1)],[branch.(errtag1)]])
    
; deal with references
    texrefs = hiigalaxy.texref
    texrefs = strtrim(texrefs[uniq(texrefs,sort(texrefs))],2)
    hiiinfo.(reftag) = strjoin(texrefs,',')

; compute the average metallicity, excluding HII regions with
; ambiguous branches and other predetermined regions; note that these
; HII regions are also not used when fitting the abundance gradients;
; only exclude regions if they haven't already been excluded by
; SINGS_ASSIGN_R23BRANCH (e.g., S3A2 in NGC3621)
    good = where($
      (strmatch(hiigalaxy.(branchtag),'*rej*',/fold) eq 0) and $
      (strtrim(hiigalaxy.(branchtag),2) ne 'A'),ngood)

    rej = -1
    case strupcase(strcompress(hiiinfo.galaxy,/remove)) of
; outer regions +062-170 and -033-118 have a really strange ionization
; level, but keep them for now; region +013+097 has a *very* low
; excitation (P~0.05), so the PT05 abundance may be wonky; -010-025
; yields a very low abundance for its R/R25 position and makes the
; gradient positive; the excitation P~0.13, so it's not crazy
; to exclude it 
       'NGC3521': begin
          rej = where((strmatch(strtrim(hiigalaxy[good].region,2),'*-033-118*',/fold) eq 1))
          if (method eq 'pt05') then rej = where((strmatch(strtrim(hiigalaxy[good].region,2),'*-010-025*',/fold) eq 1))
;         rej = where((strmatch(strtrim(hiigalaxy[good].region,2),'*+062-170*',/fold) eq 1))
;         rej = where((strmatch(strtrim(hiigalaxy[good].region,2),'*+013+097*',/fold) eq 1))
;        (strmatch(strtrim(hiigalaxy[good].region,2),'+062-170',/fold) eq 1) or $
;        (strmatch(strtrim(hiigalaxy[good].region,2),'-033-118',/fold) eq 1))
       end
; the gradient is the same whether or not we reject these two outer
; HII regions from Ryder+95, so just keep them
;     'NGC3621': rej = where($
;       (strmatch(hiigalaxy[good].region,'*S3A2*',/fold)) or $
;       (strmatch(hiigalaxy[good].region,'*S3A1*',/fold)))
       else: rej = -1
    endcase 
    if (rej[0] ne -1) then hiigalaxy[good[rej]].(branchtag) = 'X'

; compute the average metallicity    
    good = where($
      (strmatch(hiigalaxy.(branchtag),'*rej*',/fold) eq 0) and $
      (strtrim(hiigalaxy.(branchtag),2) ne 'A') and $
      (strtrim(hiigalaxy.(branchtag),2) ne 'X'),ngood)
    if (ngood eq 0) then begin
       splog, 'No unambiguous branches for method '+$
         strupcase(method)+' - tread carefully!!'
       good = where($
         (strmatch(hiigalaxy.(branchtag),'*rej*',/fold) eq 0) and $
         (strtrim(hiigalaxy.(branchtag),2) ne 'X'),ngood)
    endif
    if (ngood eq 0) then message, 'This should not happen'
    hiiinfo.(usetag) = ngood

    alloh = hiigalaxy[good].(ohtag)
    oh = reform(alloh[0,*])
    oh_err = reform(alloh[1,*])

; compute the weighted and unweighted averages, for comparison
    avg = djs_mean(oh)
    sig = djsig(oh)
    meansig = djs_mean(oh_err)
    wavg = sings_weighted_mean(oh,oh_err,wsigma=wsig)
    if keyword_set(weighted_avg) then $
      hiiinfo.(avgtag) = [wavg,wsig>meansig>0.01] else $
      hiiinfo.(avgtag) = [avg,sig>meansig>0.01]

; write to the output file
    tags = tag_names(hiigalaxy)
    printf, lun, '--------------------------------------------------'
    printf, lun, '### '+hiiinfo.galaxy
    printf, lun, '### '+strupcase(method)+' calibration'
    struct_print, struct_trimtags(hiigalaxy,select=['region',$
      tags[ohtag],tags[branchtag],'texref','p']), /no_head, lun=lun
    printf, lun, '  N(HII) = '+strtrim(ngood,2)
    printf, lun, '  Unweighted Average = '+string([avg,sig],format='(F7.5,1x,F7.5)')
    printf, lun, '  Weighted   Average = '+string([wavg,wsig],format='(F7.5,1x,F7.5)')
    printf, lun, '  Mean error = '+string(meansig,format='(F7.5)')
    printf, lun, ' '
    
return, hiiinfo
end
    
pro sings_log12oh_hiiregions, weighted_avg=weighted_avg, $
  debug=debug, clobber=clobber
; jm10mar10ucsd - derive the radial gradients and mean HII-region
; abundances for the SINGS galaxies 

; read the data    
    version = sings_log12oh_version()
    outpath = sings_path(/projects)+'log12oh/'

    outfile = outpath+'sings_log12oh_hiiregions_'+version+'.fits'
    if file_test(outfile+'.gz') and (keyword_set(clobber) eq 0) then begin
       splog, 'Output file '+outfile+' exists; use /CLOBBER'
       return
    endif

; read the data       
    sings = sings_read_info()
    hii = read_sings_log12oh_samples(/parent_hii)
    branch = mrdfits(outpath+'sings_r23_branch_'+$
      version+'.fits.gz',1,/silent)
    ngal = n_elements(sings)

; initialize the HII-region and output data structures    
    hiigalaxy1 = {$
      number:                                0, $
      galaxy:                               '', $
      hii_galaxy:                           '', $
      ned_galaxy:                           '', $
      nice_galaxy:                          '', $
      region:                               '', $
      raoffset:                         -999.0, $
      deoffset:                         -999.0, $

      r23:                     [-999.0,-999.0], $
      o32:                     [-999.0,-999.0], $
      p:                       [-999.0,-999.0], $
      
      log12oh_kk04:            [-999.0,-999.0], $
      log12oh_pt05:            [-999.0,-999.0], $
;     log12oh_o3n2:            [-999.0,-999.0], $
      log12oh_te:              [-999.0,-999.0], $ ; electron temperature
      log12oh_te_lit:          [-999.0,-999.0], $ ; literature electron temperature
      
      log12oh_pt05_upper:      [-999.0,-999.0], $
      log12oh_kk04_upper:      [-999.0,-999.0], $
      log12oh_pt05_lower:      [-999.0,-999.0], $
      log12oh_kk04_lower:      [-999.0,-999.0], $

;     hii_kk04_r23_flag:                    0L, $ ; lower branch solution greater than upper branch
;     hii_pt05_r23_flag:                    0L, $ ; lower branch solution greater than upper branch

      n2:                      [-999.0,-999.0], $
      n2o2:                    [-999.0,-999.0], $
      r23_branch:                          '?', $
      r23_branch_kk04:                     '?', $
      r23_branch_pt05:                     '?', $

      rc3_radius:                       -999.0, $
      rc3_rr25:                         -999.0, $
      twomass_radius:                   -999.0, $
      twomass_rr25:                     -999.0, $
      texref:                               '', $
      reference:                            ''}

    hiiinfo = {$
      sings_id:                                 0, $
      galaxy:                                  '', $
      ned_galaxy:                              '', $
      nice_galaxy:                             '', $
      type:                                    '', $
      t:                                       '', $

      gradient_flag:                            0, $ ; abundance gradient? (0 or 1)

      hii_nhii:                                 0, $ ; number of HII regions
      hii_r25frac:                            0.0, $ ; fraction of the R25 radius spanned by the HII regions

      hii_kk04_nhii_used:                       0, $ ; number of HII regions used
      hii_kk04_slope:             [-999.0,-999.0], $ ; KK04 abundance gradient slope
;     hii_kk04_slope_flag:                      0, $ ; formal positive slope
      hii_kk04_slope_max:                  -999.0, $ ; see PLOTSINGS_LOG12OH_HIIREGIONS
      hii_kk04_slope_min:                  -999.0, $ ; see PLOTSINGS_LOG12OH_HIIREGIONS
      hii_kk04_log12oh_nuclear:   [-999.0,-999.0], $ ; KK04 abundance gradient intercept
      hii_kk04_log12oh_central:   [-999.0,-999.0], $ ; KK04 metallicity at R=0.1*R25
      hii_kk04_log12oh_char:      [-999.0,-999.0], $ ; KK04 metallicity at R=0.4*R25
      hii_kk04_log12oh_avg:       [-999.0,-999.0], $ ; average KK04 metallicity
      hii_kk04_gradient_rms:               -999.0, $ ; rms scatter about the best-fit line 
      hii_kk04_gradient_chisq:             -999.0, $ ; chi-squared statistic
      hii_kk04_gradient_covar:        fltarr(2,2), $ ; covariance matrix
      hii_kk04_rr25_avg:          [-999.0,-999.0], $ ; average RR25 position of the KK04 metallicity
      hii_kk04_texrefs:                        '', $ ; references

      hii_pt05_nhii_used:                       0, $ ; number of HII regions used
      hii_pt05_slope:             [-999.0,-999.0], $ ; PT05 abundance gradient slope
;     hii_pt05_slope_flag:                      0, $ ; formal positive slope
      hii_pt05_log12oh_nuclear:   [-999.0,-999.0], $ ; PT05 abundance gradient intercept
      hii_pt05_log12oh_central:   [-999.0,-999.0], $ ; PT05 metallicity at R=0.1*R25
      hii_pt05_log12oh_char:      [-999.0,-999.0], $ ; PT05 metallicity at R=0.4*R25
      hii_pt05_log12oh_avg:       [-999.0,-999.0], $ ; average PT05 metallicity
      hii_pt05_gradient_rms:               -999.0, $ ; rms scatter about the best-fit line 
      hii_pt05_gradient_chisq:             -999.0, $ ; chi-squared statistic
      hii_pt05_gradient_covar:        fltarr(2,2), $ ; covariance matrix
;     hii_pt05_p_avg:             [-999.0,-999.0], $ ; average excitation P-parameter
      hii_pt05_rr25_avg:          [-999.0,-999.0], $ ; average RR25 position of the PT05 metallicity
      hii_pt05_texrefs:                        '', $ ; references
                                  
      hii_te_nhii_used:                         0, $ ; number of HII regions used
      hii_te_log12oh_avg:         [-999.0,-999.0], $ ; average Te metallicity
      hii_te_texrefs:                          '', $ ; references
      
      r23_branch:                             '?'}
    hiiinfo = replicate(hiiinfo,ngal)

    hiiinfo.sings_id    = sings.sings_id
    hiiinfo.galaxy      = strtrim(sings.galaxy,2)
    hiiinfo.ned_galaxy  = strtrim(sings.ned_galaxy,2)
    hiiinfo.nice_galaxy = strtrim(sings.nice_galaxy,2)
    hiiinfo.type        = strtrim(sings.lit_type,2)
    hiiinfo.t           = sings.lit_t

    hiiinfo.r23_branch = strtrim(branch.r23_branch,2)

; define the pre-determined set of galaxies for which we will compute
; abundance gradients
    gradient_galaxies = ['NGC0628','NGC0925','NGC1097','NGC2403','NGC2841','NGC3031',$
      'NGC3184','NGC3198','NGC3351','NGC3521','NGC3621','NGC4254','NGC4321','NGC4559',$
;     'NGC4725',$
      'NGC4736','NGC5033','NGC5055','NGC5194','NGC6946','NGC7331',$
      'NGC7793'] ; ,'NGC6822'

; list of HII regions at R>R25 that were rejected in
; BUILD_SINGS_LOG12OH_SAMPLES:
; 
;   NGC0628: FGW628C, FGW628D, FGW628E, FGW628F
;   NGC0925: -149+177, -022+227
;   NGC3031: MUENCH1
;   NGC6946: FGW6946C

; loop on each object
;   for igal = 37, 37 do begin ; n3521
;   for igal = 21, 21 do begin ; n2841
;   for igal = 49, 49 do begin ; n4559
;   for igal = 38, 38 do begin ; n3621
    for igal = 0, ngal-1 do begin
       if (igal eq 0) then begin
          openw, lun1, repstr(outfile,'.fits','_avgstats.txt'), /get_lun
          psfile1 = repstr(outfile,'.fits','.ps')
          im_plotconfig, 0, pos, psfile=psfile1, charsize=1.5 ; sings_assign_r23branch QAplot
       endif
; store HII-region abundances, if available
       indx = where(strtrim(strupcase(hiiinfo[igal].ned_galaxy),2) eq $
         strtrim(strupcase(hii.ned_galaxy),2),nindx)
       if (nindx eq  0) then begin
          if (keyword_set(silent) eq 0) then splog, 'No HII-region data for '+$
            strtrim(hiiinfo[igal].galaxy,2)
          printf, lun1, '--------------------------------------------------'
          printf, lun1, '### '+strtrim(hiiinfo[igal].galaxy,2)
          printf, lun1, 'No HII regions'
          printf, lun1, ' '
       endif else begin
          splog, 'Working on '+strtrim(hiiinfo[igal].galaxy,2)          
          
          thesehii = hii[indx]
          hiiinfo[igal].hii_nhii = nindx

          hiigalaxy                = replicate(hiigalaxy1,nindx)
          hiigalaxy.number         = lindgen(nindx)+1L
          hiigalaxy.galaxy         = strtrim(hiiinfo[igal].galaxy,2)
          hiigalaxy.nice_galaxy    = strtrim(hiiinfo[igal].nice_galaxy,2)
          hiigalaxy.hii_galaxy     = strtrim(thesehii.hii_galaxy,2)
          hiigalaxy.ned_galaxy     = strtrim(thesehii.ned_galaxy,2)
          hiigalaxy.region         = strtrim(thesehii.hii_region,2)
          hiigalaxy.raoffset       = thesehii.hii_raoffset
          hiigalaxy.deoffset       = thesehii.hii_deoffset

          hiigalaxy.reference      = thesehii.reference
          hiigalaxy.texref         = thesehii.texref
          hiigalaxy.rc3_radius     = thesehii.hii_rc3_radius
          hiigalaxy.rc3_rr25       = thesehii.hii_rc3_rr25
          hiigalaxy.twomass_radius = thesehii.hii_twomass_radius
          hiigalaxy.twomass_rr25   = thesehii.hii_twomass_rr25

          hiigalaxy.n2   = transpose([ [thesehii.zstrong_niiha], [thesehii.zstrong_niiha_err] ])
          hiigalaxy.n2o2 = transpose([ [thesehii.zstrong_niioii], [thesehii.zstrong_niioii_err] ])

          hiigalaxy.r23 = transpose([ [thesehii.zstrong_r23], [thesehii.zstrong_r23_err] ])
          hiigalaxy.o32 = transpose([ [thesehii.zstrong_o32], [thesehii.zstrong_o32_err] ])
          hiigalaxy.p   = transpose([ [thesehii.zstrong_p], [thesehii.zstrong_p_err] ])

          hiigalaxy.log12oh_kk04_upper = transpose([ [thesehii.zstrong_12oh_kk04_upper], [thesehii.zstrong_12oh_kk04_upper_err] ])
          hiigalaxy.log12oh_kk04_lower = transpose([ [thesehii.zstrong_12oh_kk04_lower], [thesehii.zstrong_12oh_kk04_lower_err] ])
          hiigalaxy.log12oh_pt05_upper = transpose([ [thesehii.zstrong_12oh_pt05_upper], [thesehii.zstrong_12oh_pt05_upper_err] ])
          hiigalaxy.log12oh_pt05_lower = transpose([ [thesehii.zstrong_12oh_pt05_lower], [thesehii.zstrong_12oh_pt05_lower_err] ])

;         hiigalaxy.log12oh_o3n2   = transpose([ [thesehii.zstrong_12oh_oiiinii_pettini], [thesehii.zstrong_12oh_oiiinii_pettini_err] ])
          hiigalaxy.log12oh_te     = transpose([ [thesehii.zt_log12oh_te], [thesehii.zt_log12oh_te_err] ])
          hiigalaxy.log12oh_te_lit = transpose([ [thesehii.lit_log12oh_te], [thesehii.lit_log12oh_te_err] ])

; by default, assign all the HII regions the R23 branch of the galaxy
          hiigalaxy.r23_branch = strtrim(branch[igal].r23_branch,2)

; compute the average abundances
;         if strmatch(hiiinfo[igal].galaxy,'*0925*',/fold) then stop
          hiiinfo[igal] = log12oh_avg(hiiinfo[igal],hiigalaxy,thesehii,$
            method='pt05',weighted_avg=weighted_avg,lun=lun1)
          hiiinfo[igal] = log12oh_avg(hiiinfo[igal],hiigalaxy,thesehii,$
            method='kk04',weighted_avg=weighted_avg,lun=lun1)

; compute the range of the disk diameter spanned by all the regions 
          rr25_good = where((hiigalaxy.rc3_radius gt -900.0),nrr25_good)
          if (nrr25_good ne 0) then hiiinfo[igal].hii_r25frac = max(hiigalaxy[rr25_good].rc3_rr25) - $
            min(hiigalaxy[rr25_good].rc3_rr25)

; compute the abundance gradients
          hiiinfo[igal].gradient_flag = fix(total(strmatch(gradient_galaxies,'*'+$
            strtrim(hiiinfo[igal].galaxy,2)+'*')))

          if hiiinfo[igal].gradient_flag then begin
             hiiinfo[igal] = log12oh_gradient(hiiinfo[igal],hiigalaxy,method='kk04',debug=debug)
             hiiinfo[igal] = log12oh_gradient(hiiinfo[igal],hiigalaxy,method='pt05',debug=debug)
          endif

; store the Te abundances, if any
;
;; average electron temperature abundances, except the abundance
;; gradient galaxies (i.e., NGC2903, NGC5236=M83), which have a handful
;; of (non-representative) HII regions with measured electron
;; temperatures; there's also a temporary hack here to adopt the
;; published (literature value) of the electron temperature and oxygen
;; abundance, mostly from Izotov et 2006 and Kniazev et 2004 of SDSS
;; galaxies without [OII] 3727 detections (here, the Te was derived
;; from the [OII] 7375 doublet)
;
;          good_te = where(((hiigalaxy.log12oh_te[0] gt -900.0) or (hiigalaxy.log12oh_te_lit[0] gt -900.0)),ngood_te)
;          if (ngood_te ne  0) then begin
;
;;            if (max(abs(hiigalaxy[good_te].log12oh_te_lit[0]-djs_median(hiigalaxy[good_te].log12oh_te_lit[0]))) gt 0.01) then begin
;;               struct_print, im_struct_trimtags(hiigalaxy[good_te],select=['galaxy','ned_galaxy','region',$
;;                 'r23_branch','log12oh_pt05','log12oh_o3n2','log12oh_te','log12oh_te_lit','reference'],newtags=['galaxy','ned_galaxy',$
;;                 'region','branch','oh_pt05','oh_o3n2','oh_te','oh_te_lit','ref']) ;, no_head=(k gt  0)
;;               cc = get_kbrd(1)
;;            endif
;
;             flag = where(((hiigalaxy[good_te].log12oh_te_lit[0] gt -900.0) and (hiigalaxy[good_te].log12oh_te[0] lt -900.0)),nflag)
;;            flag = where(((hiigalaxy[good_te].log12oh_te_lit[0] gt -900.0) and (hiigalaxy[good_te].log12oh_te[0] lt -900.0)) or $
;;              ((hiigalaxy[good_te].log12oh_te_lit[0] lt -900.0) and (hiigalaxy[good_te].log12oh_te[0] gt -900.0)),nflag)
;             if (nflag gt  0) then begin
;
;;               struct_print, im_struct_trimtags(hiigalaxy[good_te],select=['galaxy','ned_galaxy','region',$
;;                 'r23_branch','log12oh_pt05','log12oh_o3n2','log12oh_te','log12oh_te_lit','reference'],newtags=['galaxy','ned_galaxy',$
;;                 'region','branch','oh_pt05','oh_o3n2','oh_te','oh_te_lit','ref']) ;, no_head=(k gt  0)
;;               cc = get_kbrd(1)             
;
;                oh = [-999.0]
;                oh_err = [-999.0]
;                te_texrefs = ['']
;                
;                oh_lit = where((hiigalaxy[good_te].log12oh_te[0] lt -900.0) and (hiigalaxy[good_te].log12oh_te_lit[0] gt -900.0),$
;                  noh_lit,comp=oh_me,ncomp=noh_me)
;                if (noh_lit ne  0) then begin
;                   oh         = [oh,hiigalaxy[good_te[oh_lit]].log12oh_te_lit[0]]
;                   oh_err     = [oh_err,hiigalaxy[good_te[oh_lit]].log12oh_te_lit[1]]
;                   te_texrefs = [te_texrefs,hiigalaxy[good_te[oh_lit]].texref]
;                endif
;                if (noh_me ne  0) then begin
;                   oh         = [oh,hiigalaxy[good_te[oh_me]].log12oh_te[0]]
;                   oh_err     = [oh_err,hiigalaxy[good_te[oh_me]].log12oh_te[1]]
;                   te_texrefs = [te_texrefs,hiigalaxy[good_te[oh_me]].texref]
;                endif
;
;                ngood_te = n_elements(oh)-1L ; offset
;                hiiinfo[igal].hii_te_nhii_used = ngood_te
;                
;                oh         = oh[1L:ngood_te]
;                oh_err     = oh_err[1L:ngood_te]
;
;                te_texrefs = te_texrefs[1L:ngood_te]
;                te_texrefs = strtrim(te_texrefs[uniq(te_texrefs,sort(te_texrefs))],2)
;                hiiinfo[igal].hii_te_texrefs = strjoin(te_texrefs,',')
;                
;             endif else begin
;             
;                hiiinfo[igal].hii_te_nhii_used = ngood_te
;
;                te_texrefs = hiigalaxy[good_te].texref
;                te_texrefs = strtrim(te_texrefs[uniq(te_texrefs,sort(te_texrefs))],2)
;                hiiinfo[igal].hii_te_texrefs = strjoin(te_texrefs,',')
;
;                oh     = hiigalaxy[good_te].log12oh_te[0]
;                oh_err = hiigalaxy[good_te].log12oh_te[1]
;
;             endelse             
;
;             if keyword_set(weighted_avg) then begin
;                hiiinfo[igal].hii_te_log12oh_avg[0] = total(oh/oh_err^2.0)/total(1.0/oh_err^2.0)
;                hiiinfo[igal].hii_te_log12oh_avg[1] = 1.0/sqrt(total(1.0/oh_err^2.0))
;             endif else begin
;                hiiinfo[igal].hii_te_log12oh_avg[0] = im_mean(oh,sigrej=sigrej1)
;                if (ngood_te gt 1L) then $
;                  hiiinfo[igal].hii_te_log12oh_avg[1] = djsig(oh,sigrej=sigrej1)>djs_median(oh_err) else $ ; NOTE!!
;;                 hiiinfo[igal].hii_te_log12oh_avg[1] = djsig(oh,sigrej=sigrej1)/sqrt(ngood_te) else $
;                    hiiinfo[igal].hii_te_log12oh_avg[1] = oh_err
;             endelse
;
;          endif
;          
;;         good_te = where((hiigalaxy.log12oh_te_lit[0] gt -900.0) or (hiigalaxy.log12oh_te[0] gt -900.0),ngood_te)
;;         if (ngood_te gt  0) then begin
;;
;;            flag = where(((hiigalaxy[good_te].log12oh_te_lit[0] gt -900.0) and (hiigalaxy[good_te].log12oh_te[0] lt -900.0)) or $
;;              ((hiigalaxy[good_te].log12oh_te_lit[0] lt -900.0) and (hiigalaxy[good_te].log12oh_te[0] gt -900.0)),nflag)
;;            if (nflag gt  0) then begin
;;               struct_print, im_struct_trimtags(hiigalaxy[good_te],select=['galaxy','ned_galaxy','region',$
;;                 'r23_branch','log12oh_pt05','log12oh_o3n2','log12oh_te','log12oh_te_lit','reference'],newtags=['galaxy','ned_galaxy',$
;;                 'region','branch','oh_pt05','oh_o3n2','oh_te','oh_te_lit','ref']) ;, no_head=(k gt  0)
;;               cc = get_kbrd(1)             
;;            endif
;;            
;;         endif

; concatenate the HIIGALAXY structure          
          if (n_elements(allhiigalaxy) eq  0) then allhiigalaxy = hiigalaxy else $
            allhiigalaxy = [allhiigalaxy,hiigalaxy]
       endelse
       if (igal eq ngal-1) then begin
          free_lun, lun1
          im_plotconfig, psfile=psfile1, /psclose, /gzip
       endif
    endfor

    allhiigalaxy.number = lindgen(n_elements(allhiigalaxy))+1 ; unique ID number

;; give some statistics on the number of abundances
;
;    nucoh = where((hiiinfo.nuclear_log12oh_kk04[0] gt -900.0) or $
;      (hiiinfo.drift20_log12oh_kk04[0] gt -900.0) or (hiiinfo.hii_kk04_log12oh_central[0] gt -900.0),nnucoh)
;    charoh = where((hiiinfo.drift56_log12oh_kk04[0] gt -900.0) or (hiiinfo.hii_kk04_log12oh_char[0] gt -900.0),ncharoh)
;    print
;    splog, 'Nuclear 12+log(O/H): '+string(nnucoh,format='(I0)')+'/75'+'='+string(100.0*nnucoh/75.0,format='(I0)')+'%'
;    splog, 'Characteristic 12+log(O/H): '+string(ncharoh,format='(I0)')+'/75'+'='+string(100.0*ncharoh/75.0,format='(I0)')+'%'
;    print

; write out
    splog, 'Writing '+outfile
    mwrfits, hiiinfo, outfile, /create
    mwrfits, allhiigalaxy, outfile
    spawn, 'gzip -f '+outfile

return
end
    
