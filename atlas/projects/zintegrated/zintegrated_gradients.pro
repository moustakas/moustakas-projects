;+
; NAME:
;       ZINTEGRATED_GRADIENTS
;
; PURPOSE:
;       Compute abundance gradients for the ZINTEGRATED sample.
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Nov 08, U of A - written
;       jm06jan19uofa - code excised from SINGS_LOG12OH
;       jm06mar26uofa - updated/re-written
;-

pro zintegrated_gradients, intdust, intnodust, hii, result, plotavg=plotavg, $
  errorbars=errorbars, weightedfit=weightedfit, silent=silent, debug=debug, $
  blackwhite=blackwhite, paper=paper, write=write

; NOTE: We use the 2MASS de-projected radius (hii_twomass_radius) but
; normalized to the RC3 R25
    
    plot_m91 = 1L
;   plot_kk04 = 1L
;   plot_zkh = 1L
    plot_pt05 = 1L
;   plot_o3n2 = 1L
    
    outpath = atlas_path(/projects)+'zintegrated/'
    pspath = atlas_path(/papers)+'zintegrated/FIG_ZINTEGRATED/'

    if (n_elements(intdust) eq 0L) then intdust = read_zintegrated_sample(intnodust=intnodust,hii=hii)
    ngalaxy = n_elements(intdust)
    nhii = n_elements(hii)

    if keyword_set(paper) then write = 1L
    if keyword_set(write) then nmonte = 500L else nmonte = 3L

    rr25axis = findgen(51)/10.0
    rr25_char = 0.4
    nhii_min = 4L ; 6L
    r25frac_min = 0.1 ; 0.3 ; minimum HII-region disk covering fraction
    r23_cut = 1.1
    
    ohrange = [7.6,9.3]
    ohrange2 = [6.8,9.9]
    ohrange3 = [7.5,9.6]

    hiigalaxy1 = {$
      galaxy:                               '', $
      region:                               '', $

      zstrong_12oh_o3n2:       [-999.0,-999.0], $
      zstrong_12oh_pt05:       [-999.0,-999.0], $
      zstrong_12oh_kk04:       [-999.0,-999.0], $
      zstrong_12oh_zkh:        [-999.0,-999.0], $
      zstrong_12oh_m91:        [-999.0,-999.0], $

      zstrong_12oh_pt05_upper: [-999.0,-999.0], $
      zstrong_12oh_kk04_upper: [-999.0,-999.0], $
      zstrong_12oh_m91_upper:  [-999.0,-999.0], $
      zstrong_12oh_pt05_lower: [-999.0,-999.0], $
      zstrong_12oh_kk04_lower: [-999.0,-999.0], $
      zstrong_12oh_m91_lower:  [-999.0,-999.0], $

      hii_pt05_r23_flag:                    0L, $ ; lower branch solution greater than upper branch
      hii_kk04_r23_flag:                    0L, $ ; lower branch solution greater than upper branch
      hii_m91_r23_flag:                     0L, $ ; lower branch solution greater than upper branch

      n2:                      [-999.0,-999.0], $
      n2o2:                    [-999.0,-999.0], $
      r23_branch:                          'U', $ ; assume all on the upper branch
      rc3_radius:                       -999.0, $
      rc3_rr25:                         -999.0, $
      twomass_radius:                   -999.0, $
      twomass_rr25:                     -999.0, $
      texref:                               '', $
      reference:                            ''}
    
; initialize the output data structure    
    
    result = {$

      atlas_id:                                0L, $
      galaxy:                                  '', $
      ned_galaxy:                              '', $
      nice_galaxy:                             '', $
      m_b:                                 -999.0, $
      type:                                    '', $
      r25:                                 -999.0, $
      incl:                                -999.0, $
      pa:                                  -999.0, $
      gradient_flag:                           1L, $ ; abundance gradient? (0 or 1)
                                       
      hii_nhii:                                0L, $ ; number of HII regions
      hii_r25frac:                            0.0, $ ; fraction of the R25 radius spanned by the HII regions

      hii_o3n2_nhii_used:                      0L, $ ; number of HII regions used
      hii_o3n2_slope:             [-999.0,-999.0], $ ; O3N2 abundance gradient slope
      hii_o3n2_slope_flag:                     0L, $ ; formal positive slope
      hii_o3n2_log12oh_central:   [-999.0,-999.0], $ ; O3N2 abundance gradient intercept
      hii_o3n2_log12oh_char:      [-999.0,-999.0], $ ; characteristic O3N2 metallicity (R=0.4R25)
      hii_o3n2_log12oh_avg:       [-999.0,-999.0], $ ; average O3N2 metallicity
      hii_o3n2_rr25_avg:          [-999.0,-999.0], $ ; average RR25 position of the O3N2 metallicity
      hii_o3n2_texrefs:                        '', $ ; references
                                      
      hii_zkh_nhii_used:                       0L, $ ; number of HII regions used
      hii_zkh_slope:              [-999.0,-999.0], $ ; M91 abundance gradient slope
      hii_zkh_slope_flag:                      0L, $ ; formal positive slope
      hii_zkh_log12oh_central:    [-999.0,-999.0], $ ; M91 abundance gradient intercept
      hii_zkh_log12oh_char:       [-999.0,-999.0], $ ; characteristic M91 metallicity (R=0.4R25)
      hii_zkh_log12oh_avg:        [-999.0,-999.0], $ ; average M91 metallicity
      hii_zkh_rr25_avg:           [-999.0,-999.0], $ ; average RR25 position of the M91 metallicity
      hii_zkh_texrefs:                         '', $ ; references
                                  
      hii_m91_nhii_used:                       0L, $ ; number of HII regions used
      hii_m91_slope:              [-999.0,-999.0], $ ; M91 abundance gradient slope
      hii_m91_slope_flag:                      0L, $ ; formal positive slope
      hii_m91_log12oh_central:    [-999.0,-999.0], $ ; M91 abundance gradient intercept
      hii_m91_log12oh_char:       [-999.0,-999.0], $ ; characteristic M91 metallicity (R=0.4R25)
      hii_m91_log12oh_avg:        [-999.0,-999.0], $ ; average M91 metallicity
      hii_m91_rr25_avg:           [-999.0,-999.0], $ ; average RR25 position of the M91 metallicity
      hii_m91_texrefs:                         '', $ ; references
                                  
      hii_pt05_nhii_used:                      0L, $ ; number of HII regions used
      hii_pt05_slope:             [-999.0,-999.0], $ ; PT05 abundance gradient slope
      hii_pt05_slope_flag:                     0L, $ ; formal positive slope
      hii_pt05_log12oh_central:   [-999.0,-999.0], $ ; PT05 abundance gradient intercept
      hii_pt05_log12oh_char:      [-999.0,-999.0], $ ; characteristic PT05 metallicity (R=0.4R25)
      hii_pt05_log12oh_avg:       [-999.0,-999.0], $ ; average PT05 metallicity
      hii_pt05_rr25_avg:          [-999.0,-999.0], $ ; average RR25 position of the PT05 metallicity
      hii_pt05_texrefs:                        '', $ ; references
                                  
      hii_kk04_nhii_used:                      0L, $ ; number of HII regions used
      hii_kk04_slope:             [-999.0,-999.0], $ ; KK04 abundance gradient slope
      hii_kk04_slope_flag:                     0L, $ ; formal positive slope
      hii_kk04_log12oh_central:   [-999.0,-999.0], $ ; KK04 abundance gradient intercept
      hii_kk04_log12oh_char:      [-999.0,-999.0], $ ; characteristic KK04 metallicity (R=0.4R25)
      hii_kk04_log12oh_avg:       [-999.0,-999.0], $ ; average KK04 metallicity
      hii_kk04_rr25_avg:          [-999.0,-999.0], $ ; average RR25 position of the KK04 metallicity
      hii_kk04_texrefs:                        '', $ ; references

      r23_branch:                             'U', $ ; assume all on the upper branch
      
      int_class:                              '?', $
      int_ebv:                    [-999.0,-999.0], $
      int_ewhb:                            -999.0, $

      int_log12oh_lower_limit:                  0, $

      int_obs_log12oh_o3n2:       [-999.0,-999.0], $
      int_obs_log12oh_pt05:       [-999.0,-999.0], $
      int_obs_log12oh_kk04:       [-999.0,-999.0], $
      int_obs_log12oh_zkh:        [-999.0,-999.0], $
      int_obs_log12oh_m91:        [-999.0,-999.0], $
      int_obs_r23:                [-999.0,-999.0], $

      int_obs_rr25_o3n2:                   -999.0, $ ; radius where the abundance intercepts the gradient
      int_obs_rr25_pt05:                   -999.0, $
      int_obs_rr25_kk04:                   -999.0, $
      int_obs_rr25_zkh:                    -999.0, $
      int_obs_rr25_m91:                    -999.0, $

      int_cor_log12oh_pt05:       [-999.0,-999.0], $
      int_cor_log12oh_kk04:       [-999.0,-999.0], $
      int_cor_log12oh_zkh:        [-999.0,-999.0], $
      int_cor_log12oh_m91:        [-999.0,-999.0], $
      int_cor_r23:                [-999.0,-999.0], $

      int_cor_rr25_pt05:                   -999.0, $
      int_cor_rr25_kk04:                   -999.0, $
      int_cor_rr25_zkh:                    -999.0, $
      int_cor_rr25_m91:                    -999.0, $

      int_ew_log12oh_pt05:        [-999.0,-999.0], $
      int_ew_log12oh_kk04:        [-999.0,-999.0], $
      int_ew_log12oh_zkh:         [-999.0,-999.0], $
      int_ew_log12oh_m91:         [-999.0,-999.0], $
      int_ew_r23:                 [-999.0,-999.0], $

      int_ew_rr25_pt05:                    -999.0, $
      int_ew_rr25_kk04:                    -999.0, $
      int_ew_rr25_zkh:                     -999.0, $
      int_ew_rr25_m91:                     -999.0}

    result = replicate(result,ngalaxy)

    result_stats = {$

      mean_int_obs_char_pt05:    [-999.0,-999.0], $
      median_int_obs_char_pt05:           -999.0, $
      int_obs_char_pt05_ebv_coeff:        -999.0, $ ; spearman rank correlation coefficient
      int_obs_char_pt05_ebv_prob:         -999.0, $ ; probability of above
      int_obs_char_pt05_slope_coeff:      -999.0, $
      int_obs_char_pt05_slope_prob:       -999.0, $
      int_obs_char_pt05_incl_coeff:       -999.0, $
      int_obs_char_pt05_incl_prob:        -999.0, $
      int_obs_char_pt05_hasb_coeff:       -999.0, $
      int_obs_char_pt05_hasb_prob:        -999.0, $

      mean_int_cor_char_pt05:    [-999.0,-999.0], $
      median_int_cor_char_pt05:           -999.0, $
      int_cor_char_pt05_ebv_coeff:        -999.0, $ ; spearman rank correlation coefficient
      int_cor_char_pt05_ebv_prob:         -999.0, $ ; probability of above
      int_cor_char_pt05_slope_coeff:      -999.0, $
      int_cor_char_pt05_slope_prob:       -999.0, $
      int_cor_char_pt05_incl_coeff:       -999.0, $
      int_cor_char_pt05_incl_prob:        -999.0, $
      int_cor_char_pt05_hasb_coeff:       -999.0, $
      int_cor_char_pt05_hasb_prob:        -999.0, $

      mean_int_ew_char_pt05:    [-999.0,-999.0], $
      median_int_ew_char_pt05:           -999.0, $
      int_ew_char_pt05_ebv_coeff:        -999.0, $ ; spearman rank correlation coefficient
      int_ew_char_pt05_ebv_prob:         -999.0, $ ; probability of above
      int_ew_char_pt05_slope_coeff:      -999.0, $
      int_ew_char_pt05_slope_prob:       -999.0, $
      int_ew_char_pt05_incl_coeff:       -999.0, $
      int_ew_char_pt05_incl_prob:        -999.0, $
      int_ew_char_pt05_hasb_coeff:       -999.0, $
      int_ew_char_pt05_hasb_prob:        -999.0, $

      mean_int_obs_char_m91:     [-999.0,-999.0], $
      median_int_obs_char_m91:            -999.0, $
      int_obs_char_m91_ebv_coeff:         -999.0, $ ; spearman rank correlation coefficient
      int_obs_char_m91_ebv_prob:          -999.0, $ ; probability of above
      int_obs_char_m91_slope_coeff:       -999.0, $
      int_obs_char_m91_slope_prob:        -999.0, $
      int_obs_char_m91_incl_coeff:        -999.0, $
      int_obs_char_m91_incl_prob:         -999.0, $
      int_obs_char_m91_hasb_coeff:        -999.0, $
      int_obs_char_m91_hasb_prob:         -999.0, $

      mean_int_cor_char_m91:     [-999.0,-999.0], $
      median_int_cor_char_m91:            -999.0, $
      int_cor_char_m91_ebv_coeff:         -999.0, $ ; spearman rank correlation coefficient
      int_cor_char_m91_ebv_prob:          -999.0, $ ; probability of above
      int_cor_char_m91_slope_coeff:       -999.0, $
      int_cor_char_m91_slope_prob:        -999.0, $
      int_cor_char_m91_incl_coeff:        -999.0, $
      int_cor_char_m91_incl_prob:         -999.0, $
      int_cor_char_m91_hasb_coeff:        -999.0, $
      int_cor_char_m91_hasb_prob:         -999.0, $

      mean_int_ew_char_m91:     [-999.0,-999.0], $
      median_int_ew_char_m91:            -999.0, $
      int_ew_char_m91_ebv_coeff:         -999.0, $ ; spearman rank correlation coefficient
      int_ew_char_m91_ebv_prob:          -999.0, $ ; probability of above
      int_ew_char_m91_slope_coeff:       -999.0, $
      int_ew_char_m91_slope_prob:        -999.0, $
      int_ew_char_m91_incl_coeff:        -999.0, $
      int_ew_char_m91_incl_prob:         -999.0, $
      int_ew_char_m91_hasb_coeff:        -999.0, $
      int_ew_char_m91_hasb_prob:         -999.0}

    result.atlas_id    = intdust.atlas_id
    result.galaxy      = intdust.galaxy
    result.ned_galaxy  = intdust.ned_galaxy
    result.nice_galaxy = intdust.nice_galaxy
    result.type        = strtrim(intdust.lit_type,2)
    result.r25         = intdust.d25_maj/2.0
    result.m_b         = intdust.rc3_m_b

; take inclination and position angles from WRITE_HII_REGIONS (cf,
; especially NGC5194, where the inclination was set to be 20)

    match, strtrim(intdust.ned_galaxy,2), strtrim(hii.ned_galaxy,2), indx1, indx2
;   niceprint, intdust[indx1].ned_galaxy, hii[indx2].ned_galaxy, intdust[indx1].twomass_inclination, $
;     hii[indx2].galaxy_twomass_incl, intdust[indx1].twomass_posangle, hii[indx2].galaxy_twomass_pa, $
;     intdust[indx1].twomass_origin, intdust[indx1].d25_origin
    
    result[indx1].incl = hii[indx2].galaxy_twomass_incl
    result[indx1].pa   = hii[indx2].galaxy_twomass_pa

; compute the H-alpha surface brightness

    lee = rsex(atlas_path(/projects)+'zintegrated/zintegrated_halpha.sex')
    hanii = intdust.h_alpha[0]/(intdust.nii_6548[0]+intdust.nii_6584[0]); & niceprint, intdust.galaxy, 1.0/hanii
    lee_ha = lee.ha_flux - alog10(1+hanii)
    lee_ha_cor = 10^lee_ha*10^(0.4*k_lambda(6563.0,/odonnell)*(intnodust.ebv_hahb+intdust.ebv_mw))
;   lee_ha_cor = 10^(lee_ha + 0.44)
;   niceprint, alog10(intnodust.h_alpha[0]), alog10(intdust.h_alpha[0]), lee_ha, lee_ha_cor
    
;   plot, lee_ha, alog10(intdust.h_alpha[0]), ps=4, xsty=3, ysty=3, xr=[-12.5,-10.5], yr=[-12.5,-10.5]
;   oplot, !x.crange, !y.crange, thick=2
;   stats = im_stats(lee_ha-alog10(intdust.h_alpha[0]),/verbose)

;   hasb = alog10((4.0*intnodust.h_alpha[0]/(result.r25/60.0*!dtor)*3.086D18)^2.0)
;   lee_hasb = alog10((4.0*lee_ha_cor/(result.r25/60.0*!dtor)*3.086D18)^2.0)
;   niceprint, hasb, lee_hasb, lee.log_sigma_sfr+34.65

    result.hasb = alog10(4.0*intnodust.h_alpha[0]*3.086D18^2/(result.r25/60.0*!dtor)^2) ; [erg/s/pc2]
    
;   lee_hasb = alog10(4.0*lee_ha_cor*3.086D18^2/(result.r25/60.0*!dtor)^2) ; [erg/s/pc2]
;   niceprint, result.hasb-34.65, lee_hasb-34.65, lee.log_sigma_sfr

;   w = where(lee.log_sigma_sfr gt -900.0)
;   stats = im_stats(result[w].hasb-34.65-lee[w].log_sigma_sfr,/verbose)
;   stats = im_stats(result.hasb-lee_hasb,/verbose)
    
; store the observed, reddening-corrected, and EW abundances

    result.int_class = intdust.bpt_nii_mixture_class
    result.int_ebv   = transpose([[intnodust.ebv_hahb], [intnodust.ebv_hahb_err]])
    result.int_ewhb  = intdust.h_beta_ew[0]

    result.int_log12oh_lower_limit  = intdust.oh_lower_limit
    
    result.int_obs_log12oh_o3n2 = transpose([[intdust.zstrong_12oh_oiiinii_niiha], [intdust.zstrong_12oh_oiiinii_niiha_err]])
    result.int_obs_log12oh_pt05 = transpose([[intdust.zstrong_12oh_pt05_upper],    [intdust.zstrong_12oh_pt05_upper_err]])
    result.int_obs_log12oh_kk04 = transpose([[intdust.zstrong_12oh_kk04_r23_upper],[intdust.zstrong_12oh_kk04_r23_upper_err]])
    result.int_obs_log12oh_zkh  = transpose([[intdust.zstrong_12oh_zkh94],         [intdust.zstrong_12oh_zkh94_err]])
    result.int_obs_log12oh_m91  = transpose([[intdust.zstrong_12oh_m91_upper],     [intdust.zstrong_12oh_m91_upper_err]])
    result.int_obs_r23          = transpose([[intdust.zstrong_r23],                [intdust.zstrong_r23_err]])

    result.int_cor_log12oh_pt05 = transpose([[intnodust.zstrong_12oh_pt05_upper],    [intnodust.zstrong_12oh_pt05_upper_err]])
    result.int_cor_log12oh_kk04 = transpose([[intnodust.zstrong_12oh_kk04_r23_upper],[intnodust.zstrong_12oh_kk04_r23_upper_err]])
    result.int_cor_log12oh_zkh  = transpose([[intnodust.zstrong_12oh_zkh94],         [intnodust.zstrong_12oh_zkh94_err]])
    result.int_cor_log12oh_m91  = transpose([[intnodust.zstrong_12oh_m91_upper],     [intnodust.zstrong_12oh_m91_upper_err]])
    result.int_cor_r23          = transpose([[intnodust.zstrong_r23],                [intnodust.zstrong_r23_err]])

    result.int_ew_log12oh_pt05 = transpose([[intdust.zstrong_ew_12oh_pt05_upper],    [intdust.zstrong_ew_12oh_pt05_upper_err]])
    result.int_ew_log12oh_kk04 = transpose([[intdust.zstrong_ew_12oh_kk04_r23_upper],[intdust.zstrong_ew_12oh_kk04_r23_upper_err]])
    result.int_ew_log12oh_zkh  = transpose([[intdust.zstrong_ew_12oh_zkh94],         [intdust.zstrong_ew_12oh_zkh94_err]])
    result.int_ew_log12oh_m91  = transpose([[intdust.zstrong_ew_12oh_m91_upper],     [intdust.zstrong_ew_12oh_m91_upper_err]])
    result.int_ew_r23          = transpose([[intdust.zstrong_ew_r23],                [intdust.zstrong_ew_r23_err]])

; by how much would the EW abundances change if we assume alpha=0.9?

;   zz = im_abundance(intdust,ewalpha=0.9,/silent,snrcut=1.0)
;   stats = im_stats(zz.zstrong_ew_12oh_m91_upper-intdust.zstrong_ew_12oh_m91_upper,/verbose)
;   stats = im_stats(zz.zstrong_ew_12oh_pt05_upper-intdust.zstrong_ew_12oh_pt05_upper,/verbose)
    
;   struct_print, struct_trimtags(result,select=['GALAXY','INT_*LOG12OH_KK04*'])
    
    galaxy = strtrim(result.galaxy,2)

    if keyword_set(write) then debug = 0L
    if keyword_set(debug) then im_window, 0, xratio=0.5, /square
    
; ---------------------------------------------------------------------------        
; initialize abundance gradient plotting variables
; ---------------------------------------------------------------------------    

;   xrange = [0.1,3.0]
    xrange = [-0.2,2.1]
    yrange = ohrange
    xtitle = '\rho / \rho_{25}'       
    ytitle = '12 + log (O/H)'
;   plotsym, 0, 1.0, /fill

    o3n2_psym = 8L
    pt05_psym = 8L
    kk04_psym = 0L
    zkh_psym  = 0L
    m91_psym  = 0L

    if keyword_set(blackwhite) then o3n2_color = 'black' else o3n2_color = 'blue'
    if keyword_set(blackwhite) then pt05_color = 'black' else pt05_color = 'blue'
    if keyword_set(blackwhite) then kk04_color = 'black' else kk04_color = 'red'
    if keyword_set(blackwhite) then zkh_color  = 'black' else zkh_color  = 'red'
    if keyword_set(blackwhite) then m91_color  = 'black' else m91_color  = 'red'

    avg_psym = 4 & avg_color = 'grey'
    
    if keyword_set(write) or keyword_set(write) then begin

       postthick = 8.0 
       postthick2 = 5.0 
       pcharsize = 1.3
       lcharsize = 1.0
       psize = 0.8
       psize2 = 0.8

    endif else begin

       postthick = 2.0 
       postthick2 = 2.0 
       pcharsize = 1.5
       lcharsize = 1.5
       psize = 1.3
       psize2 = 3.0

    endelse

; open the ZINTEGRATED_GRADIENTS postscript file
       
    if keyword_set(write) then begin

       if keyword_set(blackwhite) then suffix = '_blackwhite' else suffix = ''

       if keyword_set(paper) then $
         psname = 'zintegrated_gradients'+suffix+'.eps' else $
         psname = 'zintegrated_gradients'+suffix+'.ps'

       ncols = 3L & nrows = 5L
       xspace = 0.0 & yspace = 0.0
       xmargin = [1.0,0.3] & ymargin = [0.3,0.8]
       width = 2.4*(lonarr(ncols)+1) & height = 2.0*(lonarr(nrows)+1)
;      width = 1.85*(lonarr(ncols)+1) & height = 1.65*(lonarr(nrows)+1)
       xpage = total(width)+total(xmargin)+total(xspace)
       ypage = total(height)+total(ymargin)+total(yspace)

       arm_plotconfig, nx=ncols, ny=nrows, xmargin=xmargin, $
         ymargin=ymargin, xspace=xspace, yspace=yspace, width=width, $
         height=height, coord=pos, xpage=xpage, ypage=ypage, $
         psfile=pspath+psname, /writeover, bw=blackwhite
       cleanplot, /silent

;      dfpsplot, pspath+psname, xsize=6.6, ysize=6.6, /color, /encap
;      pagemaker, nx=nx, ny=ny, xspace=xspace, yspace=yspace, width=width1, height=height1, $
;        xmargin=xmargin, ymargin=ymargin, xpage=xpage, ypage=ypage, $
;        position=pos, /normal
 
    endif

; ---------------------------------------------------------------------------    
; loop on each object
; ---------------------------------------------------------------------------    

    gradient_counter = 0L

;   for k = 8, ngalaxy-1L do begin
    for k = 0L, ngalaxy-1L do begin

;      print, k, galaxy[k]

; ###########################################################################       
; find the matching HII regions for this galaxy, and fill the
; corresponding structure
; ###########################################################################       

       indx = where(strtrim(strupcase(result[k].ned_galaxy),2) eq strtrim(strupcase(hii.ned_galaxy),2),nindx)
       if (nindx eq 0L) then message, 'Problem here.'

       result[k].hii_nhii = nindx
       
       hiigalaxy                = replicate(hiigalaxy1,nindx)
       hiigalaxy.galaxy         = strtrim(hii[indx].hii_galaxy,2)
       hiigalaxy.region         = strtrim(hii[indx].hii_region,2)
       hiigalaxy.reference      = hii[indx].reference
       hiigalaxy.texref         = hii[indx].texref
       hiigalaxy.rc3_radius     = hii[indx].hii_rc3_radius
       hiigalaxy.rc3_rr25       = hii[indx].hii_rc3_rr25
       hiigalaxy.twomass_radius = hii[indx].hii_twomass_radius
       hiigalaxy.twomass_rr25   = hii[indx].hii_twomass_radius/(hii[indx].galaxy_rc3_r25*60.0)
;      niceprint, hii[indx].hii_twomass_radius/(hii[indx].galaxy_rc3_r25*60.0), hii[indx].hii_twomass_rr25

;      hiigalaxy.twomass_rr25   = hii[indx].hii_twomass_rr25

       hiigalaxy.n2   = transpose([ [hii[indx].zstrong_niiha], [hii[indx].zstrong_niiha_err] ])
       hiigalaxy.n2o2 = transpose([ [hii[indx].zstrong_niioii], [hii[indx].zstrong_niioii_err] ])

       hiigalaxy.zstrong_12oh_o3n2 = transpose([ [hii[indx].zstrong_12oh_oiiinii_niiha], [hii[indx].zstrong_12oh_oiiinii_niiha_err] ])
       hiigalaxy.zstrong_12oh_zkh  = transpose([ [hii[indx].zstrong_12oh_zkh94], [hii[indx].zstrong_12oh_zkh94_err] ])

       hiigalaxy.zstrong_12oh_pt05_upper = transpose([ [hii[indx].zstrong_12oh_pt05_upper], [hii[indx].zstrong_12oh_pt05_upper_err] ])
       hiigalaxy.zstrong_12oh_pt05_lower = transpose([ [hii[indx].zstrong_12oh_pt05_lower], [hii[indx].zstrong_12oh_pt05_lower_err] ])
       hiigalaxy.zstrong_12oh_kk04_upper = transpose([ [hii[indx].zstrong_12oh_kk04_r23_upper], [hii[indx].zstrong_12oh_kk04_r23_upper_err] ])
       hiigalaxy.zstrong_12oh_kk04_lower = transpose([ [hii[indx].zstrong_12oh_kk04_r23_lower], [hii[indx].zstrong_12oh_kk04_r23_lower_err] ])
       hiigalaxy.zstrong_12oh_m91_upper  = transpose([ [hii[indx].zstrong_12oh_m91_upper], [hii[indx].zstrong_12oh_m91_upper_err] ])
       hiigalaxy.zstrong_12oh_m91_lower  = transpose([ [hii[indx].zstrong_12oh_m91_lower], [hii[indx].zstrong_12oh_m91_lower_err] ])

       radsort = sort(hiigalaxy.rc3_rr25)

; ###########################################################################       
; figure out the appropriate branch for each HII region
; ###########################################################################       

; figure out the appropriate R23 branch
       
;      lo = where((hiigalaxy.n2[0] gt -900.0) and (hiigalaxy.n2[0] lt -1.0),nlo)
;      if (nlo ne 0L) then hiigalaxy[lo].r23_branch = 'L'
;      up = where((hiigalaxy.n2[0] gt -900.0) and (hiigalaxy.n2[0] gt -1.0),nup)
;      if (nup ne 0L) then hiigalaxy[up].r23_branch = 'U'

       lo = where((hiigalaxy.n2[0] gt -900.0) and (hiigalaxy.n2[0] lt -1.0) and $
         (hiigalaxy.n2o2[0] gt -900.0) and (hiigalaxy.n2o2[0] lt -1.05),nlo)
       if (nlo ne 0L) then hiigalaxy[lo].r23_branch = 'L'
       up = where((hiigalaxy.n2[0] gt -900.0) and (hiigalaxy.n2[0] gt -1.0) and $
         (hiigalaxy.n2o2[0] gt -900.0) and (hiigalaxy.n2o2[0] gt -0.8),nup)
       if (nup ne 0L) then hiigalaxy[up].r23_branch = 'U'
       ambig = where((hiigalaxy.n2o2[0] gt -900.0) and (hiigalaxy.n2o2[0] gt -1.05) and (hiigalaxy.n2o2[0] lt -0.8),nambig)
       if (nambig ne 0L) then hiigalaxy[ambig].r23_branch = 'A'

       w = where((hiigalaxy[radsort].r23_branch eq '?') or (hiigalaxy[radsort].r23_branch eq 'A'),nw)
       if (not keyword_set(silent)) then $
         struct_print, struct_trimtags(hiigalaxy[radsort],select=['GALAXY','REGION','N2','N2O2',$
           'ZSTRONG_12OH_ZKH','R23_BRANCH','RC3_RR25','REFERENCE'])
;          'ZSTRONG_12OH_KK04_UPPER','ZSTRONG_12OH_KK04_LOWER','R23_BRANCH','RC3_RR25','REFERENCE'])

; assign ambigious branches by hand for galaxies with abundance gradients

;      unknown = where((hiigalaxy.r23_branch eq '?') or (hiigalaxy.r23_branch eq 'A'),nunknown)
;      if (nunknown ne 0L) then hiigalaxy[unknown].r23_branch = 'U' ; put everything on the upper branch

       case galaxy[k] of
          'NGC1058': begin
             hiigalaxy.r23_branch = 'U'
             hiigalaxy[where($
               (strmatch(strtrim(hiigalaxy.region,2),'FGW1058H',/fold) eq 1B) or $
               (strmatch(strtrim(hiigalaxy.region,2),'FGW1058G',/fold) eq 1B))].r23_branch = 'A'
          end
          'NGC1569': begin
          end
          'NGC2541': begin
          end
          'NGC2903': begin
;            hiigalaxy[where($
;              (strmatch(strtrim(hiigalaxy.region,2),'-062-085',/fold) eq 1B) or $
;              (strmatch(strtrim(hiigalaxy.region,2),'-065-073',/fold) eq 1B) or $
;              (strmatch(strtrim(hiigalaxy.region,2),'-067-061',/fold) eq 1B) or $
;              (strmatch(strtrim(hiigalaxy.region,2),'+171+243',/fold) eq 1B))].r23_branch = 'A'
          end
          'NGC3198': begin
          end
          'NGC3344': begin
             hiigalaxy.r23_branch = 'U' ; places +007+157 from Vilchez et al. (1988) on the upper branch (also measured by McCall et al. 1985)
          end
          'NGC3351': begin
          end
          'NGC3521': begin
             hiigalaxy[where($
               (strmatch(strtrim(hiigalaxy.region,2),'-033-118',/fold) eq 1B))].r23_branch = 'A'
          end
          'NGC4254': begin
          end
          'NGC4303': begin
             hiigalaxy.r23_branch = 'U' ; places H278 and H234 from Henry et al. (1992) on the upper branch
          end
          'NGC4321': begin
          end
          'NGC4651': begin
             hiigalaxy.r23_branch = 'U' ; places -077-043 from Skillman et al. (1996) on the upper branch
          end
          'NGC4713': begin
             hiigalaxy.r23_branch = 'U'
          end
          'NGC4736': begin
          end
          'NGC4861': begin
          end
          'NGC5194': begin
          end
          else: 
       endcase 

       if (not keyword_set(silent)) then struct_print, struct_trimtags(hiigalaxy[radsort],$
;        select=['GALAXY','REGION','ZSTRONG_12OH_PT05_UPPER','ZSTRONG_12OH_PT05_LOWER','R23_BRANCH','RC3_RR25','REFERENCE'])
;        select=['GALAXY','REGION','ZSTRONG_12OH_ZKH','R23_BRANCH','RC3_RR25','REFERENCE'])
         select=['GALAXY','REGION','ZSTRONG_12OH_M91_UPPER','ZSTRONG_12OH_M91_LOWER','R23_BRANCH','RC3_RR25','REFERENCE'])
;           select=['GALAXY','REGION','ZSTRONG_12OH_KK04_UPPER','ZSTRONG_12OH_KK04_LOWER','R23_BRANCH','RC3_RR25','REFERENCE'])
;           select=['GALAXY','REGION','N2','N2O2','ZSTRONG_12OH_KK04_UPPER','ZSTRONG_12OH_KK04_LOWER',$
;           'R23_BRANCH','RC3_RR25','REFERENCE'])

; do some R23 checking and flag these HII regions as "ambigious"
       
       flag = where((hiigalaxy.zstrong_12oh_pt05_lower[0] gt -900.0) and (hiigalaxy.zstrong_12oh_pt05_upper[0] gt -900.0) and $
         (hiigalaxy.zstrong_12oh_pt05_lower[0] gt hiigalaxy.zstrong_12oh_pt05_upper[0]) and $
         ((hiigalaxy.zstrong_12oh_pt05_upper[0]+hiigalaxy.zstrong_12oh_pt05_upper[1]) lt $
         (hiigalaxy.zstrong_12oh_pt05_lower[0]-hiigalaxy.zstrong_12oh_pt05_lower[1])),nflag)
       if (nflag ne 0L) then hiigalaxy[flag].hii_pt05_r23_flag = 1L

       flag = where((hiigalaxy.zstrong_12oh_m91_lower[0] gt -900.0) and (hiigalaxy.zstrong_12oh_m91_upper[0] gt -900.0) and $
         (hiigalaxy.zstrong_12oh_m91_lower[0] gt hiigalaxy.zstrong_12oh_m91_upper[0]) and $
         ((hiigalaxy.zstrong_12oh_m91_upper[0]+hiigalaxy.zstrong_12oh_m91_upper[1]) lt $
         (hiigalaxy.zstrong_12oh_m91_lower[0]-hiigalaxy.zstrong_12oh_m91_lower[1])),nflag)
       if (nflag ne 0L) then hiigalaxy[flag].hii_m91_r23_flag = 1L

       flag = where((hiigalaxy.zstrong_12oh_kk04_lower[0] gt -900.0) and (hiigalaxy.zstrong_12oh_kk04_upper[0] gt -900.0) and $
         (hiigalaxy.zstrong_12oh_kk04_lower[0] gt hiigalaxy.zstrong_12oh_kk04_upper[0]) and $
         ((hiigalaxy.zstrong_12oh_kk04_upper[0]+hiigalaxy.zstrong_12oh_kk04_upper[1]) lt $
         (hiigalaxy.zstrong_12oh_kk04_lower[0]-hiigalaxy.zstrong_12oh_kk04_lower[1])),nflag)
       if (nflag ne 0L) then hiigalaxy[flag].hii_kk04_r23_flag = 1L

;      if (nflag ne 0L) then hiigalaxy[flag].hii_pt05_r23_flag = 1L
;      flag = where((hiigalaxy.zstrong_12oh_kk04_lower[0] gt hiigalaxy.zstrong_12oh_kk04_upper[0]),nflag)
;      if (nflag ne 0L) then hiigalaxy[flag].hii_kk04_r23_flag = 1L
;      flag = where((hiigalaxy.zstrong_12oh_m91_lower[0] gt hiigalaxy.zstrong_12oh_m91_upper[0]),nflag)
;      if (nflag ne 0L) then hiigalaxy[flag].hii_m91_r23_flag = 1L

;         flag = where((hiigalaxy.zstrong_12oh_pt05_lower[0] gt hiigalaxy.zstrong_12oh_pt05_upper[0]) or $
;           (hiigalaxy.zstrong_12oh_kk04_lower[0] gt hiigalaxy.zstrong_12oh_kk04_upper[0]) or $
;           (hiigalaxy.zstrong_12oh_m91_lower[0] gt hiigalaxy.zstrong_12oh_m91_upper[0]),nflag)
;         if (nflag ne 0L) then begin
;            niceprint, replicate('Check: ',nflag), hiigalaxy[flag].galaxy, hiigalaxy[flag].region, $
;              hiigalaxy[flag].zstrong_12oh_kk04_upper[0], hiigalaxy[flag].zstrong_12oh_kk04_lower[0], $
;              hiigalaxy[flag].zstrong_12oh_pt05_upper[0], hiigalaxy[flag].zstrong_12oh_pt05_lower[0]
;         endif
       
; ###########################################################################       
; compute the average HII-region abundances
; ###########################################################################       
       
       good_o3n2 = where((hiigalaxy.zstrong_12oh_o3n2[0] gt -900.0) and (hiigalaxy.r23_branch ne 'A') and $
         (hiigalaxy.r23_branch ne '?'),ngood_o3n2)
       if (ngood_o3n2 ne 0L) then begin

          result[k].hii_o3n2_nhii_used = ngood_o3n2

          o3n2_texrefs = hiigalaxy[good_o3n2].texref
          o3n2_texrefs = strtrim(o3n2_texrefs[uniq(o3n2_texrefs,sort(o3n2_texrefs))],2)
          result[k].hii_o3n2_texrefs = strjoin(o3n2_texrefs,',')

          oh     = hiigalaxy[good_o3n2].zstrong_12oh_o3n2[0]
          oh_err = hiigalaxy[good_o3n2].zstrong_12oh_o3n2[1]
          result[k].hii_o3n2_log12oh_avg[0] = djs_mean(oh)
          if (ngood_o3n2 gt 1L) then $
            result[k].hii_o3n2_log12oh_avg[1] = djsig(oh) else $
            result[k].hii_o3n2_log12oh_avg[1] = oh_err
;            result[k].hii_o3n2_log12oh_avg[0] = total(oh/oh_err^2)/total(1.0/oh_err^2)
;            result[k].hii_o3n2_log12oh_avg[1] = 1.0/sqrt(total(1.0/oh_err^2))

          rad = where(hiigalaxy[good_o3n2].rc3_rr25 gt -900.0,nrad)
          if (nrad ne 0L) then $
            result[k].hii_o3n2_rr25_avg = [djs_mean(hiigalaxy[good_o3n2[rad]].rc3_rr25),$
              djsig(hiigalaxy[good_o3n2[rad]].rc3_rr25)]

       endif

       good_kk04 = where((hiigalaxy.zstrong_12oh_kk04_upper[0] gt -900.0),ngood_kk04) ; upper or lower would work here
       if (ngood_kk04 ne 0L) then begin

          up = where((hiigalaxy[good_kk04].r23_branch eq 'U') and (hiigalaxy[good_kk04].hii_kk04_r23_flag eq 0L),nup)
          if (nup ne 0L) then hiigalaxy[good_kk04[up]].zstrong_12oh_kk04 = hiigalaxy[good_kk04[up]].zstrong_12oh_kk04_upper
          lo = where((hiigalaxy[good_kk04].r23_branch eq 'L') and (hiigalaxy[good_kk04].hii_kk04_r23_flag eq 0L),nlo)
          if (nlo ne 0L) then hiigalaxy[good_kk04[lo]].zstrong_12oh_kk04 = hiigalaxy[good_kk04[lo]].zstrong_12oh_kk04_lower

          kk04_avg = where((hiigalaxy[good_kk04].zstrong_12oh_kk04[0] gt -900.0),nkk04_avg)
          if (nkk04_avg ne 0L) then begin

             result[k].hii_kk04_nhii_used = nkk04_avg ; exclude ambigious objects

             kk04_texrefs = hiigalaxy[good_kk04[kk04_avg]].texref

             kk04_texrefs = strtrim(kk04_texrefs[uniq(kk04_texrefs,sort(kk04_texrefs))],2)
             result[k].hii_kk04_texrefs = strjoin(kk04_texrefs,',')

             oh     = hiigalaxy[good_kk04[kk04_avg]].zstrong_12oh_kk04[0]
             oh_err = hiigalaxy[good_kk04[kk04_avg]].zstrong_12oh_kk04[1]
             result[k].hii_kk04_log12oh_avg[0] = djs_mean(oh)
             if (ngood_kk04 gt 1L) then $
               result[k].hii_kk04_log12oh_avg[1] = djsig(oh) else $
               result[k].hii_kk04_log12oh_avg[1] = oh_err
;               result[k].hii_kk04_log12oh_avg[0] = total(oh/oh_err^2)/total(1.0/oh_err^2)
;               result[k].hii_kk04_log12oh_avg[1] = 1.0/sqrt(total(1.0/oh_err^2))

             rad = where(hiigalaxy[good_kk04[kk04_avg]].rc3_rr25 gt -900.0,nrad)
             if (nrad ne 0L) then $
               result[k].hii_kk04_rr25_avg = [djs_mean(hiigalaxy[good_kk04[kk04_avg[rad]]].rc3_rr25),$
                 djsig(hiigalaxy[good_kk04[kk04_avg[rad]]].rc3_rr25)]
          endif
          
       endif

       good_pt05 = where((hiigalaxy.zstrong_12oh_pt05_upper[0] gt -900.0),ngood_pt05) ; upper or lower would work here
       if (ngood_pt05 ne 0L) then begin

          up = where((hiigalaxy[good_pt05].r23_branch eq 'U') and (hiigalaxy[good_pt05].hii_pt05_r23_flag eq 0L),nup)
          if (nup ne 0L) then hiigalaxy[good_pt05[up]].zstrong_12oh_pt05 = hiigalaxy[good_pt05[up]].zstrong_12oh_pt05_upper
          lo = where((hiigalaxy[good_pt05].r23_branch eq 'L') and (hiigalaxy[good_pt05].hii_pt05_r23_flag eq 0L),nlo)
          if (nlo ne 0L) then hiigalaxy[good_pt05[lo]].zstrong_12oh_pt05 = hiigalaxy[good_pt05[lo]].zstrong_12oh_pt05_lower

          pt05_avg = where((hiigalaxy[good_pt05].zstrong_12oh_pt05[0] gt -900.0),npt05_avg)
          if (npt05_avg ne 0L) then begin

             result[k].hii_pt05_nhii_used = npt05_avg ; exclude ambigious objects

             pt05_texrefs = hiigalaxy[good_pt05[pt05_avg]].texref
             pt05_texrefs = strtrim(pt05_texrefs[uniq(pt05_texrefs,sort(pt05_texrefs))],2)
             result[k].hii_pt05_texrefs = strjoin(pt05_texrefs,',')

             oh     = hiigalaxy[good_pt05[pt05_avg]].zstrong_12oh_pt05[0]
             oh_err = hiigalaxy[good_pt05[pt05_avg]].zstrong_12oh_pt05[1]
             result[k].hii_pt05_log12oh_avg[0] = djs_mean(oh)
             if (ngood_pt05 gt 1L) then $
               result[k].hii_pt05_log12oh_avg[1] = djsig(oh) else $
               result[k].hii_pt05_log12oh_avg[1] = oh_err
;               result[k].hii_pt05_log12oh_avg[0] = total(oh/oh_err^2)/total(1.0/oh_err^2)
;               result[k].hii_pt05_log12oh_avg[1] = 1.0/sqrt(total(1.0/oh_err^2))

             rad = where(hiigalaxy[good_pt05[pt05_avg]].rc3_rr25 gt -900.0,nrad)
             if (nrad ne 0L) then $
               result[k].hii_pt05_rr25_avg = [djs_mean(hiigalaxy[good_pt05[pt05_avg[rad]]].rc3_rr25),$
                 djsig(hiigalaxy[good_pt05[pt05_avg[rad]]].rc3_rr25)]
          endif

       endif

       good_zkh = where((hiigalaxy.zstrong_12oh_zkh[0] gt -900.0) and (hiigalaxy.r23_branch ne 'A') and $
         (hiigalaxy.r23_branch ne '?'),ngood_zkh)
;      good_zkh = where(hiigalaxy.zstrong_12oh_zkh[0] gt -900.0,ngood_zkh)
       if (ngood_zkh ne 0L) then begin

          result[k].hii_zkh_nhii_used = ngood_zkh

          zkh_texrefs = hiigalaxy[good_zkh].texref
          zkh_texrefs = strtrim(zkh_texrefs[uniq(zkh_texrefs,sort(zkh_texrefs))],2)
          result[k].hii_zkh_texrefs = strjoin(zkh_texrefs,',')

          oh     = hiigalaxy[good_zkh].zstrong_12oh_zkh[0]
          oh_err = hiigalaxy[good_zkh].zstrong_12oh_zkh[1]
          result[k].hii_zkh_log12oh_avg[0] = djs_mean(oh)
          if (ngood_zkh gt 1L) then $
            result[k].hii_zkh_log12oh_avg[1] = djsig(oh) else $
            result[k].hii_zkh_log12oh_avg[1] = oh_err
;            result[k].hii_zkh_log12oh_avg[0] = total(oh/oh_err^2)/total(1.0/oh_err^2)
;            result[k].hii_zkh_log12oh_avg[1] = 1.0/sqrt(total(1.0/oh_err^2))

          rad = where(hiigalaxy[good_zkh].rc3_rr25 gt -900.0,nrad)
          if (nrad ne 0L) then $
            result[k].hii_zkh_rr25_avg = [djs_mean(hiigalaxy[good_zkh[rad]].rc3_rr25),$
              djsig(hiigalaxy[good_zkh[rad]].rc3_rr25)]

       endif

       good_m91 = where((hiigalaxy.zstrong_12oh_m91_upper[0] gt -900.0),ngood_m91) ; upper or lower would work here
       if (ngood_m91 ne 0L) then begin

          result[k].hii_m91_nhii_used = ngood_m91

          up = where((hiigalaxy[good_m91].r23_branch eq 'U') and (hiigalaxy[good_m91].hii_m91_r23_flag eq 0L),nup)
          if (nup ne 0L) then hiigalaxy[good_m91[up]].zstrong_12oh_m91 = hiigalaxy[good_m91[up]].zstrong_12oh_m91_upper
          lo = where((hiigalaxy[good_m91].r23_branch eq 'L') and (hiigalaxy[good_m91].hii_m91_r23_flag eq 0L),nlo)
          if (nlo ne 0L) then hiigalaxy[good_m91[lo]].zstrong_12oh_m91 = hiigalaxy[good_m91[lo]].zstrong_12oh_m91_lower

          m91_avg = where((hiigalaxy[good_m91].zstrong_12oh_m91[0] gt -900.0),nm91_avg)
          if (nm91_avg ne 0L) then begin

             result[k].hii_m91_nhii_used = nm91_avg ; exclude ambigious objects

             m91_texrefs = hiigalaxy[good_m91[m91_avg]].texref
             m91_texrefs = strtrim(m91_texrefs[uniq(m91_texrefs,sort(m91_texrefs))],2)
             result[k].hii_m91_texrefs = strjoin(m91_texrefs,',')

             oh     = hiigalaxy[good_m91[m91_avg]].zstrong_12oh_m91[0]
             oh_err = hiigalaxy[good_m91[m91_avg]].zstrong_12oh_m91[1]
             result[k].hii_m91_log12oh_avg[0] = djs_mean(oh)
             if (ngood_m91 gt 1L) then $
               result[k].hii_m91_log12oh_avg[1] = djsig(oh) else $
               result[k].hii_m91_log12oh_avg[1] = oh_err
;               result[k].hii_m91_log12oh_avg[0] = total(oh/oh_err^2)/total(1.0/oh_err^2)
;               result[k].hii_m91_log12oh_avg[1] = 1.0/sqrt(total(1.0/oh_err^2))

             rad = where(hiigalaxy[good_m91[m91_avg]].rc3_rr25 gt -900.0,nrad)
             if (nrad ne 0L) then $
               result[k].hii_m91_rr25_avg = [djs_mean(hiigalaxy[good_m91[m91_avg[rad]]].rc3_rr25),$
                 djsig(hiigalaxy[good_m91[m91_avg[rad]]].rc3_rr25)]
          endif
          
       endif

; ###########################################################################       
; compute and plot the abundance gradient
; ###########################################################################       

       if result[k].gradient_flag then begin

          result[k].hii_r25frac = (max(hiigalaxy.rc3_rr25)-min(hiigalaxy.rc3_rr25))
          
;         print, galaxy[k], result[k].hii_r25frac
;         splog, 'Computing abundance gradient for '+strtrim(result[k].galaxy,2)+' using '+$
;           string(nindx,format='(I0)')+' HII regions.'
          
          gradient_counter = gradient_counter + 1L
          
          keep = where(hiigalaxy.rc3_radius gt -900.0,nkeep,comp=nodata,ncomp=nnodata)
;         if (nnodata ne 0L) and keyword_set(debug) then begin
;            splog, 'The following HII regions need RC3 radius measurements.'
;            niceprint, strtrim(hiigalaxy[nodata].galaxy,2), strtrim(hiigalaxy[nodata].region,2), $
;              hiigalaxy[nodata].zstrong_12oh_o3n2[0], hiigalaxy[nodata].rc3_radius, $
;              strtrim(hiigalaxy[nodata].reference,2)
;         endif

          ghiigalaxy = hiigalaxy[keep]

; --------------------------------------------------
; measure the gradient: M91
; --------------------------------------------------

          good_m91 = where(ghiigalaxy.zstrong_12oh_m91[0] ge -900.0,ngood_m91)
          if (ngood_m91 ne 0L) then begin

;            if (ngood_m91 ne result[k].hii_m91_nhii_used) then message, 'Problem!'

             rr25_m91 = ghiigalaxy[good_m91].rc3_rr25
             rr25_m91_err = rr25_m91*0.0

             oh_m91 = ghiigalaxy[good_m91].zstrong_12oh_m91[0]
             oh_m91_err = ghiigalaxy[good_m91].zstrong_12oh_m91[1]

             if (ngood_m91 ge nhii_min) and (result[k].hii_r25frac gt r25frac_min) then begin

                res_m91 = sings_fit_gradient(rr25_m91,rr25_m91_err,$
                  oh_m91,oh_m91_err,nmonte=nmonte,xchar=rr25_char,$
                  noslope=1,weightedfit=weightedfit)
                result[k].hii_m91_slope           = res_m91.slope
                result[k].hii_m91_slope_flag      = res_m91.slope_flag
                result[k].hii_m91_log12oh_central = res_m91.int
                result[k].hii_m91_log12oh_char    = res_m91.ychar

                if (result[k].int_obs_log12oh_m91[0] gt -900.0) and (res_m91.slope_flag eq 0L) then $
                  result[k].int_obs_rr25_m91 = interpol(res_m91.xaxis,res_m91.yfit,$
                    result[k].int_obs_log12oh_m91[0])
                
                if (result[k].int_cor_log12oh_m91[0] gt -900.0) and (res_m91.slope_flag eq 0L) then $
                  result[k].int_cor_rr25_m91 = interpol(res_m91.xaxis,res_m91.yfit,$
                    result[k].int_cor_log12oh_m91[0])
                
                if (result[k].int_ew_log12oh_m91[0] gt -900.0) and (res_m91.slope_flag eq 0L) then $
                  result[k].int_ew_rr25_m91 = interpol(res_m91.xaxis,res_m91.yfit,$
                    result[k].int_ew_log12oh_m91[0])

             endif
             
          endif

; --------------------------------------------------
; measure the gradient: ZKH
; --------------------------------------------------

          good_zkh = where((ghiigalaxy.zstrong_12oh_zkh[0] gt -900.0) and (ghiigalaxy.r23_branch ne 'A') and $
            (ghiigalaxy.r23_branch ne '?'),ngood_zkh)
;         good_zkh = where(ghiigalaxy.zstrong_12oh_zkh[0] ge -900.0,ngood_zkh)
          if (ngood_zkh ne 0L) then begin

;            if (ngood_zkh ne result[k].hii_zkh_nhii_used) then message, 'Problem!'

             rr25_zkh = ghiigalaxy[good_zkh].rc3_rr25
             rr25_zkh_err = rr25_zkh*0.0

             oh_zkh = ghiigalaxy[good_zkh].zstrong_12oh_zkh[0]
             oh_zkh_err = ghiigalaxy[good_zkh].zstrong_12oh_zkh[1]

             if (ngood_zkh ge nhii_min) and (result[k].hii_r25frac gt r25frac_min) then begin

                res_zkh = sings_fit_gradient(rr25_zkh,rr25_zkh_err,$
                  oh_zkh,oh_zkh_err,nmonte=nmonte,xchar=rr25_char,$
                  noslope=1,weightedfit=weightedfit)
                result[k].hii_zkh_slope           = res_zkh.slope
                result[k].hii_zkh_slope_flag      = res_zkh.slope_flag
                result[k].hii_zkh_log12oh_central = res_zkh.int
                result[k].hii_zkh_log12oh_char    = res_zkh.ychar

                if (result[k].int_obs_log12oh_zkh[0] gt -900.0) and (res_zkh.slope_flag eq 0L) then $
                  result[k].int_obs_rr25_zkh = interpol(res_zkh.xaxis,res_zkh.yfit,$
                    result[k].int_obs_log12oh_zkh[0])
                
                if (result[k].int_cor_log12oh_zkh[0] gt -900.0) and (res_zkh.slope_flag eq 0L) then $
                  result[k].int_cor_rr25_zkh = interpol(res_zkh.xaxis,res_zkh.yfit,$
                    result[k].int_cor_log12oh_zkh[0])
                
                if (result[k].int_ew_log12oh_zkh[0] gt -900.0) and (res_zkh.slope_flag eq 0L) then $
                  result[k].int_ew_rr25_zkh = interpol(res_zkh.xaxis,res_zkh.yfit,$
                    result[k].int_ew_log12oh_zkh[0])

             endif
             
          endif

; --------------------------------------------------
; measure the gradient: PT05
; --------------------------------------------------
          
          good_pt05 = where((ghiigalaxy.zstrong_12oh_pt05[0] ge -900.0),ngood_pt05)
          if (ngood_pt05 ne 0L) then begin

;            if (ngood_pt05 ne result[k].hii_pt05_nhii_used) then message, 'Problem!'

             rr25_pt05 = ghiigalaxy[good_pt05].rc3_rr25
             rr25_pt05_err = rr25_pt05*0.0

             oh_pt05 = ghiigalaxy[good_pt05].zstrong_12oh_pt05[0]
             oh_pt05_err = ghiigalaxy[good_pt05].zstrong_12oh_pt05[1]

             if (ngood_pt05 ge nhii_min) and (result[k].hii_r25frac gt r25frac_min) then begin

                res_pt05 = sings_fit_gradient(rr25_pt05,rr25_pt05_err,$
                  oh_pt05,oh_pt05_err,nmonte=nmonte,xchar=rr25_char,$
                  noslope=1,weightedfit=weightedfit)
                result[k].hii_pt05_slope           = res_pt05.slope
                result[k].hii_pt05_slope_flag      = res_pt05.slope_flag
                result[k].hii_pt05_log12oh_central = res_pt05.int
                result[k].hii_pt05_log12oh_char    = res_pt05.ychar

                if (result[k].int_obs_log12oh_pt05[0] gt -900.0) and (res_pt05.slope_flag eq 0L) then $
                  result[k].int_obs_rr25_pt05 = interpol(res_pt05.xaxis,res_pt05.yfit,$
                    result[k].int_obs_log12oh_pt05[0])
                
                if (result[k].int_cor_log12oh_pt05[0] gt -900.0) and (res_pt05.slope_flag eq 0L) then $
                  result[k].int_cor_rr25_pt05 = interpol(res_pt05.xaxis,res_pt05.yfit,$
                    result[k].int_cor_log12oh_pt05[0])
                
                if (result[k].int_ew_log12oh_pt05[0] gt -900.0) and (res_pt05.slope_flag eq 0L) then $
                  result[k].int_ew_rr25_pt05 = interpol(res_pt05.xaxis,res_pt05.yfit,$
                    result[k].int_ew_log12oh_pt05[0])

             endif
             
          endif

; --------------------------------------------------
; measure the gradient: KK04
; --------------------------------------------------
          
          good_kk04 = where(ghiigalaxy.zstrong_12oh_kk04[0] ge -900.0,ngood_kk04)
          if (ngood_kk04 ne 0L) then begin

;            if (ngood_kk04 ne result[k].hii_kk04_nhii_used) then message, 'Problem!'

             rr25_kk04 = ghiigalaxy[good_kk04].rc3_rr25
             rr25_kk04_err = rr25_kk04*0.0

             oh_kk04 = ghiigalaxy[good_kk04].zstrong_12oh_kk04[0]
             oh_kk04_err = ghiigalaxy[good_kk04].zstrong_12oh_kk04[1]

             if (ngood_kk04 ge nhii_min) and (result[k].hii_r25frac gt r25frac_min) then begin

                res_kk04 = sings_fit_gradient(rr25_kk04,rr25_kk04_err,$
                  oh_kk04,oh_kk04_err,nmonte=nmonte,xchar=rr25_char,$
                  noslope=1,weightedfit=weightedfit)
                result[k].hii_kk04_slope           = res_kk04.slope
                result[k].hii_kk04_slope_flag      = res_kk04.slope_flag
                result[k].hii_kk04_log12oh_central = res_kk04.int
                result[k].hii_kk04_log12oh_char    = res_kk04.ychar

                if (result[k].int_obs_log12oh_kk04[0] gt -900.0) and (res_kk04.slope_flag eq 0L) then $
                  result[k].int_obs_rr25_kk04 = interpol(res_kk04.xaxis,res_kk04.yfit,$
                    result[k].int_obs_log12oh_kk04[0])
                
                if (result[k].int_cor_log12oh_kk04[0] gt -900.0) and (res_kk04.slope_flag eq 0L) then $
                  result[k].int_cor_rr25_kk04 = interpol(res_kk04.xaxis,res_kk04.yfit,$
                    result[k].int_cor_log12oh_kk04[0])
                
                if (result[k].int_ew_log12oh_kk04[0] gt -900.0) and (res_kk04.slope_flag eq 0L) then $
                  result[k].int_ew_rr25_kk04 = interpol(res_kk04.xaxis,res_kk04.yfit,$
                    result[k].int_ew_log12oh_kk04[0])

             endif
             
          endif

; --------------------------------------------------
; measure the gradient: O3N2
; --------------------------------------------------

          good_o3n2 = where((ghiigalaxy.zstrong_12oh_o3n2[0] ge -900.0) and (ghiigalaxy.r23_branch ne 'A') and $
            (ghiigalaxy.r23_branch ne '?'),ngood_o3n2)
          if (ngood_o3n2 ne 0L) then begin

;            if (ngood_o3n2 ne result[k].hii_o3n2_nhii_used) then message, 'Problem!'

             rr25_o3n2 = ghiigalaxy[good_o3n2].rc3_rr25
             rr25_o3n2_err = rr25_o3n2*0.0

             oh_o3n2 = ghiigalaxy[good_o3n2].zstrong_12oh_o3n2[0]
             oh_o3n2_err = ghiigalaxy[good_o3n2].zstrong_12oh_o3n2[1]

             if (ngood_o3n2 ge nhii_min) and (result[k].hii_r25frac gt r25frac_min) then begin

                res_o3n2 = sings_fit_gradient(rr25_o3n2,rr25_o3n2_err,$
                  oh_o3n2,oh_o3n2_err,nmonte=nmonte,xchar=rr25_char,$
                  noslope=1,weightedfit=weightedfit)
                result[k].hii_o3n2_slope           = res_o3n2.slope
                result[k].hii_o3n2_slope_flag      = res_o3n2.slope_flag
                result[k].hii_o3n2_log12oh_central = res_o3n2.int
                result[k].hii_o3n2_log12oh_char    = res_o3n2.ychar

                if (result[k].int_obs_log12oh_o3n2[0] gt -900.0) and (res_o3n2.slope_flag eq 0L) then $
                  result[k].int_obs_rr25_o3n2 = interpol(res_o3n2.xaxis,res_o3n2.yfit,$
                    result[k].int_obs_log12oh_o3n2[0])
                
             endif
             
          endif 

; --------------------------------------------------
; plot the abundance gradient
; --------------------------------------------------

          if keyword_set(debug) or keyword_set(write) then begin

             if keyword_set(write) then begin

;               psname = strlowcase(strtrim(result[k].galaxy,2))+'.eps'
;               dfpsplot, pspath+psname, xsize=6.6, ysize=6.6, /color, /encap
;               
;               ytitle1 = ytitle
;               noerase = 0L

                if ((gradient_counter-1L) mod ncols) eq 0L then begin
                   delvarx, ytickname
                   ytitle1 = ytitle
                endif else begin
                   ytickname = replicate(' ',10)
                   ytitle1 = ''
                endelse

                if (gradient_counter ge 12L) then begin
                   delvarx, xtickname
                   xtitle1 = xtitle
                endif else begin
                   xtickname = replicate(' ',10)
                   xtitle1 = ''
                endelse
                
                noerase = (gradient_counter gt 1L)
                position = pos[*,gradient_counter-1L]

             endif else begin

                xtitle1 = xtitle
                ytitle1 = ytitle
                noerase = 0L
                
             endelse

;            pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
;              xmargin=[1.5,0.4], ymargin=[0.6,1.1], xpage=8.5, ypage=8.3, $
;              position=position, /normal

;            xrange = [-0.04,(max(ghiigalaxy.hii_rc3_rr25)*1.1)>1.0]
             
             djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=pcharsize, $
               charthick=postthick2, thick=postthick, xtitle=xtitle1, ytitle=ytitle1, $
               xstyle=3, ystyle=1, xrange=xrange, yrange=yrange, position=position, $
               ytickname=ytickname, xtickname=xtickname, noerase=noerase;, /xlog
             legend, strtrim(result[k].nice_galaxy,2), /right, /top, box=0, $
               charsize=lcharsize, charthick=postthick2
;            djs_oplot, 2.1*[1,1], !y.crange, line=2
             
; -------------------------                
; M91
; -------------------------                
             if (ngood_m91 ne 0L) and keyword_set(plot_m91) then begin
                plotsym, m91_psym, psize, fill=1L, thick=postthick2
                if keyword_set(errorbars) then oploterror, rr25_m91, oh_m91, rr25_m91_err, oh_m91_err, $
                  ps=8, color=djs_icolor(m91_color), errcolor=djs_icolor(m91_color), thick=postthick2 else $
                  djs_oplot, rr25_m91, oh_m91, ps=8, color=m91_color
                ambig_m91 = where((ghiigalaxy.r23_branch eq 'A') or (ghiigalaxy.hii_m91_r23_flag),nambig_m91)
                if (nambig_m91 ne 0L) then begin
                   plotsym, m91_psym, psize, fill=0L, thick=postthick2
                   for ii = 0L, nambig_m91-1L do begin
                      xm91 = ghiigalaxy[ambig_m91[ii]].rc3_rr25
                      ym91_up = ghiigalaxy[ambig_m91[ii]].zstrong_12oh_m91_upper[0]
                      ym91_lo = ghiigalaxy[ambig_m91[ii]].zstrong_12oh_m91_lower[0]
                      djs_oplot, xm91*[1,1], [ym91_up,ym91_lo], $
                        ps=-8, line=0, thick=postthick2, color=djs_icolor(m91_color)
                   endfor
                endif
                if (result[k].hii_m91_slope[0] ne -999.0) then $
                  djs_oplot, res_m91.xaxis, res_m91.yfit, line=0, thick=postthick2 ;, color=m91_color
                plotsym, avg_psym, psize2, fill=0, thick=postthick2
                if keyword_set(plotavg) then $
                  oploterror, result[k].hii_m91_rr25_avg[0], result[k].hii_m91_log12oh_avg[0], $
                    result[k].hii_m91_rr25_avg[1], result[k].hii_m91_log12oh_avg[1], $
                    ps=8, color=djs_icolor(avg_color), errcolor=djs_icolor(avg_color), thick=postthick2
;               if (result[k].int_cor_rr25_m91 gt -900.0) and (result[k].int_cor_log12oh_m91[0] gt -900.0) then begin
;                  plots, result[k].int_cor_rr25_m91, result[k].int_cor_log12oh_m91[0], ps=7, sym=2.0, color=djs_icolor('cyan')
;                  plots, 0.4, result[k].int_cor_log12oh_m91[0], ps=4, sym=2.0, color=djs_icolor('cyan')
;               endif
             endif
; -------------------------                
; KK04
; -------------------------                
             if (ngood_kk04 ne 0L) and keyword_set(plot_kk04) then begin
                plotsym, kk04_psym, psize, fill=1L, thick=postthick2
                djs_oplot, rr25_kk04, oh_kk04, ps=8, color=kk04_color
                ambig_kk04 = where((ghiigalaxy.r23_branch eq 'A') or (ghiigalaxy.hii_kk04_r23_flag),nambig_kk04)
                if (nambig_kk04 ne 0L) then begin
                   plotsym, kk04_psym, psize, fill=0L, thick=postthick2
                   for ii = 0L, nambig_kk04-1L do begin
                      xkk04 = ghiigalaxy[ambig_kk04[ii]].rc3_rr25
                      ykk04_up = ghiigalaxy[ambig_kk04[ii]].zstrong_12oh_kk04_upper[0]
                      ykk04_lo = ghiigalaxy[ambig_kk04[ii]].zstrong_12oh_kk04_lower[0]
                      djs_oplot, xkk04*[1,1], [ykk04_up,ykk04_lo], $
                        ps=-8, line=0, thick=postthick2, color=djs_icolor(kk04_color)
                   endfor
                endif
                if (result[k].hii_kk04_slope[0] ne -999.0) then $
                  djs_oplot, res_kk04.xaxis, res_kk04.yfit, line=0, thick=postthick2 ;, color=kk04_color
                plotsym, avg_psym, psize2, fill=0, thick=postthick2
                if keyword_set(plotavg) then $
                  oploterror, result[k].hii_kk04_rr25_avg[0], result[k].hii_kk04_log12oh_avg[0], $
                    result[k].hii_kk04_rr25_avg[1], result[k].hii_kk04_log12oh_avg[1], $
                    ps=8, color=djs_icolor(avg_color), errcolor=djs_icolor(avg_color), thick=postthick2
             endif
; -------------------------                
; PT05
; -------------------------                
             if (ngood_pt05 ne 0L) and keyword_set(plot_pt05) then begin
                plotsym, pt05_psym, psize, fill=1L, thick=postthick2
                if keyword_set(errorbars) then oploterror, rr25_pt05, oh_pt05, rr25_pt05_err, oh_pt05_err, $
                  ps=8, color=djs_icolor(pt05_color), errcolor=djs_icolor(pt05_color), thick=postthick2 else $
                  djs_oplot, rr25_pt05, oh_pt05, ps=8, color=pt05_color
                ambig_pt05 = where((ghiigalaxy.r23_branch eq 'A') or (ghiigalaxy.hii_pt05_r23_flag),nambig_pt05)
                if (nambig_pt05 ne 0L) then begin
                   plotsym, pt05_psym, psize, fill=0L, thick=postthick2
                   for ii = 0L, nambig_pt05-1L do begin
                      xpt05 = ghiigalaxy[ambig_pt05[ii]].rc3_rr25*1.02 ; this factor offsets from the KK04 ambigious line
                      ypt05_up = ghiigalaxy[ambig_pt05[ii]].zstrong_12oh_pt05_upper[0]
                      ypt05_lo = ghiigalaxy[ambig_pt05[ii]].zstrong_12oh_pt05_lower[0]
                      djs_oplot, xpt05*[1,1], [ypt05_up,ypt05_lo], $
                        ps=-8, line=0, thick=postthick2, color=djs_icolor(pt05_color)
                   endfor
                endif
                if (result[k].hii_pt05_slope[0] ne -999.0) then $
                  djs_oplot, res_pt05.xaxis, res_pt05.yfit, line=2, thick=postthick2 ;, color=pt05_color
                plotsym, avg_psym, psize2, fill=0, thick=postthick2
                if keyword_set(plotavg) then $
                  oploterror, result[k].hii_pt05_rr25_avg[0], result[k].hii_pt05_log12oh_avg[0], $
                    result[k].hii_pt05_rr25_avg[1], result[k].hii_pt05_log12oh_avg[1], $
                    ps=8, color=djs_icolor(avg_color), errcolor=djs_icolor(avg_color), thick=postthick2
             endif
; -------------------------                
; ZKH
; -------------------------                
             if (ngood_zkh ne 0L) and keyword_set(plot_zkh) then begin
                plotsym, zkh_psym, psize, fill=1L, thick=postthick2
                if keyword_set(errorbars) then oploterror, rr25_zkh, oh_zkh, rr25_zkh_err, oh_zkh_err, $
                  ps=8, color=djs_icolor(zkh_color), errcolor=djs_icolor(zkh_color), thick=postthick2 else $
                  djs_oplot, rr25_zkh, oh_zkh, ps=8, color=zkh_color
                ambig_zkh = where((ghiigalaxy.r23_branch eq 'A') or (ghiigalaxy.r23_branch eq '?'),nambig_zkh)
                if (nambig_zkh ne 0L) then begin
                   plotsym, zkh_psym, psize, fill=0L, thick=postthick2
;                  if keyword_set(errorbars) then oploterror, ghiigalaxy[ambig_zkh].rc3_rr25*[1,1], $
;                    ghiigalaxy[ambig_zkh].zstrong_12oh_zkh[0], ghiigalaxy[ambig_zkh].rc3_rr25*[0,0], $
;                    ghiigalaxy[ambig_zkh].zstrong_12oh_zkh[1], ps=8, color=djs_icolor(zkh_color), $
;                    errcolor=djs_icolor(zkh_color), thick=postthick2 else $
                     if (nambig_zkh eq 1L) then plots, ghiigalaxy[ambig_zkh].rc3_rr25*[1,1], $
                       ghiigalaxy[ambig_zkh].zstrong_12oh_zkh[0], ps=8, color=djs_icolor(zkh_color) else $
                     djs_oplot, ghiigalaxy[ambig_zkh].rc3_rr25*[1,1], ghiigalaxy[ambig_zkh].zstrong_12oh_zkh[0], $
                       ps=8, color=djs_icolor(zkh_color)
                endif 
                if (result[k].hii_zkh_slope[0] ne -999.0) then $
                  djs_oplot, res_zkh.xaxis, res_zkh.yfit, line=0, thick=postthick2 ;, color=zkh_color
                plotsym, avg_psym, psize2, fill=0, thick=postthick2
                if keyword_set(plotavg) then $
                  oploterror, result[k].hii_zkh_rr25_avg[0], result[k].hii_zkh_log12oh_avg[0], $
                    result[k].hii_zkh_rr25_avg[1], result[k].hii_zkh_log12oh_avg[1], $
                    ps=8, color=djs_icolor(avg_color), errcolor=djs_icolor(avg_color), thick=postthick2
             endif
; -------------------------                
; O3N2
; -------------------------                
             if (ngood_o3n2 ne 0L) and keyword_set(plot_o3n2) then begin
                plotsym, o3n2_psym, psize, fill=1L, thick=postthick2
                djs_oplot, rr25_o3n2, oh_o3n2, ps=8, color=o3n2_color
                ambig_o3n2 = where((ghiigalaxy.r23_branch eq 'A') or (ghiigalaxy.r23_branch eq '?'),nambig_zkh)
                if (nambig_o3n2 ne 0L) then begin
                   plotsym, o3n2_psym, psize, fill=0L, thick=postthick2
;                  if keyword_set(errorbars) then oploterror, ghiigalaxy[ambig_o3n2].rc3_rr25*[1,1], $
;                    ghiigalaxy[ambig_o3n2].zstrong_12oh_o3n2[0], ghiigalaxy[ambig_o3n2].rc3_rr25*[0,0], $
;                    ghiigalaxy[ambig_o3n2].zstrong_12oh_o3n2[1], ps=8, color=djs_icolor(o3n2_color), $
;                    errcolor=djs_icolor(o3n2_color), thick=postthick2 else $
                     if (nambig_o3n2 eq 1L) then plots, ghiigalaxy[ambig_o3n2].rc3_rr25*[1,1], $
                       ghiigalaxy[ambig_o3n2].zstrong_12oh_o3n2[0], ps=8, color=djs_icolor(o3n2_color) else $
                     djs_oplot, ghiigalaxy[ambig_o3n2].rc3_rr25*[1,1], ghiigalaxy[ambig_o3n2].zstrong_12oh_o3n2[0], $
                       ps=8, color=djs_icolor(o3n2_color)
                endif 
                if (result[k].hii_3n2_slope[0] ne -999.0) then $
                  djs_oplot, res_o3n2.xaxis, res_o3n2.yfit, line=2, thick=postthick2 ;, color=o3n2_color
                plotsym, avg_psym, psize2, fill=0, thick=postthick2
                if keyword_set(plotavg) then $
                  oploterror, result[k].hii_o3n2_rr25_avg[0], result[k].hii_o3n2_log12oh_avg[0], $
                    result[k].hii_o3n2_rr25_avg[1], result[k].hii_o3n2_log12oh_avg[1], $
                    ps=8, color=djs_icolor(avg_color), errcolor=djs_icolor(avg_color), thick=postthick2
             endif
; -------------------------                
             
          endif

       endif                    ; close the abunndance gradient logical statement

       if keyword_set(debug) then begin
          splog, 'Waiting for keystroke.'
          cc = get_kbrd(1)
       endif

    endfor

; at what normalized radius does the reddening-corrected integrated
; abundance intercept the gas-phase metallicity?

;   good = where((result.int_cor_rr25_m91 gt -900.0) and (result.int_log12oh_lower_limit eq 0))
;   stats = im_stats(result[good].int_cor_rr25_m91,/verbose)
;   stats = im_stats(result[good].int_cor_rr25_pt05,/verbose)

; ---------------------------------------------------------------------------
; compute some statistics
; ---------------------------------------------------------------------------

    good = where((result.int_log12oh_lower_limit eq 0),ngood)

    ebv   = result[good].int_ebv[0]
    incl  = result[good].incl
    slope = result[good].hii_m91_slope[0]
    hasb = result[good].hasb

; M91/observed integrated abundances
    
    oh_char = result[good].hii_m91_log12oh_char[0]
    oh_int_obs = result[good].int_obs_log12oh_m91[0]
    oh_int_cor = result[good].int_cor_log12oh_m91[0]
    oh_int_ew = result[good].int_ew_log12oh_m91[0]
    
    resid_obs = oh_int_obs-oh_char

    rcor_oh_obs   = r_correlate(oh_char,oh_int_obs,zd=zd_oh_obs,probd=probd_oh_obs)
    rcor_ebv_obs   = r_correlate(ebv,resid_obs,zd=zd_ebv_obs,probd=probd_ebv_obs)
    rcor_slope_obs = r_correlate(slope,resid_obs,zd=zd_slope_obs,probd=probd_slope_obs)
    rcor_incl_obs  = r_correlate(incl,resid_obs,zd=zd_incl_obs,probd=probd_incl_obs)
    rcor_hasb_obs  = r_correlate(hasb,resid_obs,zd=zd_hasb_obs,probd=probd_hasb_obs)

    splog, '12+log(O/H) [M91, Obs vs Char]: ', rcor_oh_obs 
    
    result_stats.mean_int_obs_char_m91[0]     = mean(resid_obs)
    result_stats.mean_int_obs_char_m91[1]     = stddev(resid_obs)
    result_stats.median_int_obs_char_m91      = median(resid_obs)
    result_stats.int_obs_char_m91_ebv_coeff   = rcor_ebv_obs[0]
    result_stats.int_obs_char_m91_ebv_prob    = rcor_ebv_obs[1]
    result_stats.int_obs_char_m91_slope_coeff = rcor_slope_obs[0]
    result_stats.int_obs_char_m91_slope_prob  = rcor_slope_obs[1]
    result_stats.int_obs_char_m91_incl_coeff  = rcor_incl_obs[0]
    result_stats.int_obs_char_m91_incl_prob   = rcor_incl_obs[1]
    result_stats.int_obs_char_m91_hasb_coeff  = rcor_hasb_obs[0]
    result_stats.int_obs_char_m91_hasb_prob   = rcor_hasb_obs[1]

; M91/corrected integrated abundances
    
    resid_cor = oh_int_cor-oh_char

    rcor_oh_cor   = r_correlate(oh_char,oh_int_cor,zd=zd_oh_cor,probd=probd_oh_cor)
    rcor_ebv_cor   = r_correlate(ebv,resid_cor,zd=zd_ebv_cor,probd=probd_ebv_cor)
    rcor_slope_cor = r_correlate(slope,resid_cor,zd=zd_slope_cor,probd=probd_slope_cor)
    rcor_incl_cor  = r_correlate(incl,resid_cor,zd=zd_incl_cor,probd=probd_incl_cor)
    rcor_hasb_cor  = r_correlate(hasb,resid_cor,zd=zd_hasb_cor,probd=probd_hasb_cor)

    splog, '12+log(O/H) [M91, Cor vs Char]: ', rcor_oh_cor 

    result_stats.mean_int_cor_char_m91[0]     = mean(resid_cor)
    result_stats.mean_int_cor_char_m91[1]     = stddev(resid_cor)
    result_stats.median_int_cor_char_m91      = median(resid_cor)
    result_stats.int_cor_char_m91_ebv_coeff   = rcor_ebv_cor[0]
    result_stats.int_cor_char_m91_ebv_prob    = rcor_ebv_cor[1]
    result_stats.int_cor_char_m91_slope_coeff = rcor_slope_cor[0]
    result_stats.int_cor_char_m91_slope_prob  = rcor_slope_cor[1]
    result_stats.int_cor_char_m91_incl_coeff  = rcor_incl_cor[0]
    result_stats.int_cor_char_m91_incl_prob   = rcor_incl_cor[1]
    result_stats.int_cor_char_m91_hasb_coeff  = rcor_hasb_cor[0]
    result_stats.int_cor_char_m91_hasb_prob   = rcor_hasb_cor[1]

; M91/EW integrated abundances
    
    resid_ew = oh_int_ew-oh_char

    rcor_oh_ew   = r_correlate(oh_char,oh_int_ew,zd=zd_oh_ew,probd=probd_oh_ew)
    rcor_ebv_ew   = r_correlate(ebv,resid_ew,zd=zd_ebv_ew,probd=probd_ebv_ew)
    rcor_slope_ew = r_correlate(slope,resid_ew,zd=zd_slope_ew,probd=probd_slope_ew)
    rcor_incl_ew  = r_correlate(incl,resid_ew,zd=zd_incl_ew,probd=probd_incl_ew)
    rcor_hasb_ew  = r_correlate(hasb,resid_ew,zd=zd_hasb_ew,probd=probd_hasb_ew)

    splog, '12+log(O/H) [M91, EW vs Char]: ', rcor_oh_ew 

    result_stats.mean_int_ew_char_m91[0]     = mean(resid_ew)
    result_stats.mean_int_ew_char_m91[1]     = stddev(resid_ew)
    result_stats.median_int_ew_char_m91      = median(resid_ew)
    result_stats.int_ew_char_m91_ebv_coeff   = rcor_ebv_ew[0]
    result_stats.int_ew_char_m91_ebv_prob    = rcor_ebv_ew[1]
    result_stats.int_ew_char_m91_slope_coeff = rcor_slope_ew[0]
    result_stats.int_ew_char_m91_slope_prob  = rcor_slope_ew[1]
    result_stats.int_ew_char_m91_incl_coeff  = rcor_incl_ew[0]
    result_stats.int_ew_char_m91_incl_prob   = rcor_incl_ew[1]
    result_stats.int_ew_char_m91_hasb_coeff  = rcor_hasb_ew[0]
    result_stats.int_ew_char_m91_hasb_prob   = rcor_hasb_ew[1]

; PT05/observed integrated abundances
    
    oh_char = result[good].hii_pt05_log12oh_char[0]
    oh_int_obs = result[good].int_obs_log12oh_pt05[0]
    oh_int_cor = result[good].int_cor_log12oh_pt05[0]
    oh_int_ew = result[good].int_ew_log12oh_pt05[0]
    
    resid_obs = oh_int_obs-oh_char

    rcor_oh_obs   = r_correlate(oh_char,oh_int_obs,zd=zd_oh_obs,probd=probd_oh_obs)
    rcor_ebv_obs   = r_correlate(ebv,resid_obs,zd=zd_ebv_obs,probd=probd_ebv_obs)
    rcor_slope_obs = r_correlate(slope,resid_obs,zd=zd_slope_obs,probd=probd_slope_obs)
    rcor_incl_obs  = r_correlate(incl,resid_obs,zd=zd_incl_obs,probd=probd_incl_obs)
    rcor_hasb_obs  = r_correlate(hasb,resid_obs,zd=zd_hasb_obs,probd=probd_hasb_obs)

    splog, '12+log(O/H) [PT05, Obs vs Char]: ', rcor_oh_obs 

    result_stats.mean_int_obs_char_pt05[0]     = mean(resid_obs)
    result_stats.mean_int_obs_char_pt05[1]     = stddev(resid_obs)
    result_stats.median_int_obs_char_pt05      = median(resid_obs)
    result_stats.int_obs_char_pt05_ebv_coeff   = rcor_ebv_obs[0]
    result_stats.int_obs_char_pt05_ebv_prob    = rcor_ebv_obs[1]
    result_stats.int_obs_char_pt05_slope_coeff = rcor_slope_obs[0]
    result_stats.int_obs_char_pt05_slope_prob  = rcor_slope_obs[1]
    result_stats.int_obs_char_pt05_incl_coeff  = rcor_incl_obs[0]
    result_stats.int_obs_char_pt05_incl_prob   = rcor_incl_obs[1]
    result_stats.int_obs_char_pt05_hasb_coeff  = rcor_hasb_obs[0]
    result_stats.int_obs_char_pt05_hasb_prob   = rcor_hasb_obs[1]

; PT05/corrected integrated abundances
    
    resid_cor = oh_int_cor-oh_char

    rcor_oh_cor   = r_correlate(oh_char,oh_int_cor,zd=zd_oh_cor,probd=probd_oh_cor)
    rcor_ebv_cor   = r_correlate(ebv,resid_cor,zd=zd_ebv_cor,probd=probd_ebv_cor)
    rcor_slope_cor = r_correlate(slope,resid_cor,zd=zd_slope_cor,probd=probd_slope_cor)
    rcor_incl_cor  = r_correlate(incl,resid_cor,zd=zd_incl_cor,probd=probd_incl_cor)
    rcor_hasb_cor  = r_correlate(hasb,resid_cor,zd=zd_hasb_cor,probd=probd_hasb_cor)

    splog, '12+log(O/H) [PT05, Cor vs Char]: ', rcor_oh_cor 

    result_stats.mean_int_cor_char_pt05[0]     = mean(resid_cor)
    result_stats.mean_int_cor_char_pt05[1]     = stddev(resid_cor)
    result_stats.median_int_cor_char_pt05      = median(resid_cor)
    result_stats.int_cor_char_pt05_ebv_coeff   = rcor_ebv_cor[0]
    result_stats.int_cor_char_pt05_ebv_prob    = rcor_ebv_cor[1]
    result_stats.int_cor_char_pt05_slope_coeff = rcor_slope_cor[0]
    result_stats.int_cor_char_pt05_slope_prob  = rcor_slope_cor[1]
    result_stats.int_cor_char_pt05_incl_coeff  = rcor_incl_cor[0]
    result_stats.int_cor_char_pt05_incl_prob   = rcor_incl_cor[1]
    result_stats.int_cor_char_pt05_hasb_coeff  = rcor_hasb_cor[0]
    result_stats.int_cor_char_pt05_hasb_prob   = rcor_hasb_cor[1]

; PT05/EW integrated abundances
    
    resid_ew = oh_int_ew-oh_char

    rcor_oh_ew   = r_correlate(oh_char,oh_int_ew,zd=zd_oh_ew,probd=probd_oh_ew)
    rcor_ebv_ew   = r_correlate(ebv,resid_ew,zd=zd_ebv_ew,probd=probd_ebv_ew)
    rcor_slope_ew = r_correlate(slope,resid_ew,zd=zd_slope_ew,probd=probd_slope_ew)
    rcor_incl_ew  = r_correlate(incl,resid_ew,zd=zd_incl_ew,probd=probd_incl_ew)
    rcor_hasb_ew  = r_correlate(hasb,resid_ew,zd=zd_hasb_ew,probd=probd_hasb_ew)

    splog, '12+log(O/H) [PT05, EW vs Char]: ', rcor_oh_ew 

    result_stats.mean_int_ew_char_pt05[0]     = mean(resid_ew)
    result_stats.mean_int_ew_char_pt05[1]     = stddev(resid_ew)
    result_stats.median_int_ew_char_pt05      = median(resid_ew)
    result_stats.int_ew_char_pt05_ebv_coeff   = rcor_ebv_ew[0]
    result_stats.int_ew_char_pt05_ebv_prob    = rcor_ebv_ew[1]
    result_stats.int_ew_char_pt05_slope_coeff = rcor_slope_ew[0]
    result_stats.int_ew_char_pt05_slope_prob  = rcor_slope_ew[1]
    result_stats.int_ew_char_pt05_incl_coeff  = rcor_incl_ew[0]
    result_stats.int_ew_char_pt05_incl_prob   = rcor_incl_ew[1]
    result_stats.int_ew_char_pt05_hasb_coeff  = rcor_hasb_ew[0]
    result_stats.int_ew_char_pt05_hasb_prob   = rcor_hasb_ew[1]

; close everything and write out     

    if keyword_set(write) then begin
       dfpsclose
       cleanplot, /silent
    endif; else cc = get_kbrd(1)

    if keyword_set(write) then begin
       outfile = 'zintegrated_gradients.fits'
       splog, 'Writing '+outpath+outfile+'.'
       mwrfits, result, outpath+outfile, /create
       mwrfits, result_stats, outpath+outfile
       spawn, ['gzip -f '+outpath+outfile], /sh
    endif

return
end 
