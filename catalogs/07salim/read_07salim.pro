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
      sfr_avg: 0.0,$
      sfr_25: 0.0,$
      sfr_975: 0.0,$
      sfrm_avg: 0.0,$
      sfrm_25: 0.0,$
      sfrm_975: 0.0,$
      av_avg: 0.0,$
      av_25: 0.0,$
      av_975: 0.0,$
      b_avg: 0.0,$
      b_25: 0.0,$
      b_975: 0.0,$
      mab: fltarr(7),$
      mab_err: fltarr(7),$
      mass_avg: 0.0}
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
    cat.sfr_avg = logsfr_avg
    cat.sfr_25 = logsfr_25
    cat.sfr_975 = logsfr_975
    cat.sfrm_avg = logsfrm_avg
    cat.sfrm_25 = logsfrm_25
    cat.sfrm_975 = logsfrm_975
    cat.av_avg = av_avg
    cat.av_25 = av_25
    cat.av_975 = av_975
    cat.b_avg = logb_avg
    cat.b_25 = logb_25
    cat.b_975 = logb_975
    cat.mab = transpose([[fuv],[nuv],[u],[g],[r],[i],[z]])
    cat.mab_err = transpose([[fuv_err],[nuv_err],[u_err],[g_err],[r_err],[i_err],[z_err]])

    cat.mass_avg = cat.sfr_avg-cat.sfrm_avg
    
return, cat
end
    
    
