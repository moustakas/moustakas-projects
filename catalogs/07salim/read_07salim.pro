function read_07salim
; jm11jun16ucsd - read the Salim+07 catalog

    catfile = getenv('CATALOGS_DIR')+'/07salim/salim_data_pms3.dat'
    splog, 'Reading '+catfile
    readcol, catfile, plate, mjd, fiber, zobj, chi2, ra, dec, $
      age_avg, age_25, age_975, $
      logsfr_avg, logsfr_25, logsfr_975, $
      logsfrm_avg, logsfrm_25, logsfrm_975, $
      av_avg, av_25, av_975, $
      logb_avg, logb_25, logb_975, $
      fuv, nuv, u, g, r, i, z, fuv_err, nuv_err, u_err, g_err, r_err, i_err, z_err, $
      format='L,L,I,F,F,D,D,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F', /silent
    ngal = n_elements(plate)

    cat = {$
      plate: 0L,$
      mjd: 0L,$
      fiber: 0,$
      z: 0.0,$
      chi2: 0.0,$
      ra: 0D,$
      dec: 0D,$
      age_avg: 0.0,$
      age_25: 0.0,$
      age_975: 0.0,$
      logsfr_avg: 0.0,$
      logsfr_25: 0.0,$
      logsfr_975: 0.0,$
      logsfrm_avg: 0.0,$
      logsfrm_25: 0.0,$
      logsfrm_975: 0.0,$
      av_avg: 0.0,$
      av_25: 0.0,$
      av_975: 0.0,$
      logb_avg: 0.0,$
      logb_25: 0.0,$
      logb_975: 0.0,$
      mab: fltarr(7),$
      mab_err: fltarr(7)}
    cat = replicate(cat,ngal)

    cat.plate = plate
    cat.mjd = mjd
    cat.fiber = fiber
    cat.z = zobj
    cat.chi2 = chi2
    cat.ra = ra
    cat.dec = dec
    cat.age_avg = age_avg
    cat.age_25 = age_25
    cat.age_975 = age_975
    cat.logsfr_avg = logsfr_avg
    cat.logsfr_25 = logsfr_25
    cat.logsfr_975 = logsfr_975
    cat.logsfrm_avg = logsfrm_avg
    cat.logsfrm_25 = logsfrm_25
    cat.logsfrm_975 = logsfrm_975
    cat.av_avg = av_avg
    cat.av_25 = av_25
    cat.av_975 = av_975
    cat.logb_avg = logb_avg
    cat.logb_25 = logb_25
    cat.logb_975 = logb_975
    cat.mab = transpose([[fuv],[nuv],[u],[g],[r],[i],[z]])
    cat.mab_err = transpose([[fuv_err],[nuv_err],[u_err],[g_err],[r_err],[i_err],[z_err]])
    
return, cat
end
    
    
