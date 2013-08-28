pro check_mzages

    agesmass = read_mz_sample(/mzhii_mass)
    agesispec = read_mz_sample(/mzhii_ispec)
    agesohdust = read_mz_sample(/mzhii_log12oh)
    agesanc = read_mz_sample(/mzhii_ancillary)

    ww = where(agesmass.mass_50 gt 11.0 and agesohdust.ew_log12oh_t04 gt -900 and $
      agesohdust.ohlimit eq 0 and agesmass.zobj gt 0.55)
    srt = sort(agesohdust[ww].ew_log12oh_t04)
    niceprint, agesmass[ww[srt]].zobj, agesmass[ww[srt]].mass_50, agesohdust[ww[srt]].ew_log12oh_t04, $
      agesohdust[ww[srt]].ew_log12oh_t04_err, agesanc[ww[srt]].final_weight
    qaplot_ages_gandalf_specfit, agesispec[ww[srt]], psfile='junk.ps'    

stop    
    
return
end
