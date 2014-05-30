function get_deep2_completeness_params, field, nbins=nbins, $
  hmin=hmin, hmax=hmax

    params = {$
      field:           field,$
      targ_filter:        '',$
      select_filter:      '',$
      color_filters: strarr(2,2),$
      nice_filter:        '',$
      nice_targ_filter:   '',$
      nice_color_filter: ['',''],$
      nbins:         [0,0,0],$ ; histogram NBINS, min/max
      hmin:    [0.0,0.0,0.0],$
      hmax:    [0.0,0.0,0.0],$
      binsz:   [0.0,0.0,0.0],$
      binsz_midbin: [0.0,0.0,0.0]}
    bright = 18.5
    faint = 24.1

    params.nbins = [26,20,20]
;   params.nbins = [22,16,16]
    
;need a g-r, r-z, etc.
    params.targ_filter = 'capak_cfht_megaprime_sagem_r.par'
    params.select_filter = 'capak_cfht_megaprime_sagem_r.par'
    params.color_filters = 'capak_cfht_megaprime_sagem_'+[['g','r'],['r','z']]+'.par'
    params.nice_filter = 'r'
    params.nice_targ_filter = 'r'
    params.nice_color_filter = ['g-r','r-z']
    params.hmin = [bright,-0.3,-0.5]
    params.hmax = [faint,2.2,2.5]
    
    params.binsz = (params.hmax-params.hmin)/float(params.nbins-1)
    params.binsz_midbin = (params.hmax-params.hmin)/float(params.nbins)
    
return, params
end

