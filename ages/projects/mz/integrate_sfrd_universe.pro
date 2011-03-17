pro integrate_sfrd_universe, integrate=integrate

    mzpath = ages_path(/projects)+'mz/'
    gridpath = getenv('PEGASE_HR_SFHGRID_DIR')+'/' ; SFH grid path
    litpath = getenv('PAPERSPATH')+'/literature/data/'

    pegfile = mzpath+'sfrd_universe.fits'

; I need this scale factor because there is a volume factor in the
; SFRD which Pegase can't handle    
    sfrfactor = 1D-9
    yr2Gyr = 1D9
    zsun = 0.0122               ; Asplund et al. (2005)
    omega_baryon = 0.044        ; Pettini (2005, h=0.7)
    
    h100 = mz_h100()
    H0 = (100.0*h100)*3.1556926D7*1D6/3.086D19 ; [Myr^-1]
    bigG = 4.4983D-21 ; [Mpc^3/Myr^2/M_sun]
    rhocrit = 3.0*H0^2.0/(8.0*!pi*bigG) ; [M_sun/Mpc^3, h=0.7]
;   rhocrit = 1.35D11 ; [M_sun/Mpc^3, h=0.7, from Pettini 2006]

; coefficients from plot_redshift_vs_sfr_mass_density

    ztoday = 0.0 & zform = 6.0
    sfrd_zaxis = reverse(findgen((zform-ztoday)/0.01+1)*0.01+ztoday)
    sfrd_time = 1D3*(getage(sfrd_zaxis)-getage(zform)) ; time since zform [Myr]
    
    coeff = [0.010464392D,0.11303816D,3.4534319D,3.7109016D]
    sfrd = (coeff[0]+coeff[1]*sfrd_zaxis)/(1.0+(sfrd_zaxis/coeff[2])^coeff[3]) ; [M_sun/yr/Mpc^3]
    sfrd_Myr = sfrd*1D6 ; [M_sun/Myr/Mpc^3]
    djs_plot, sfrd_zaxis, alog10(sfrd), xsty=3, ysty=3, xtitle='Redshift', $
      ytitle='log (\rho_{SFR}) (M_{\odot} yr^{-1} Mpc^{-3})'
    djs_plot, sfrd_zaxis, alog10(1D9*sfrd/rhocrit), xsty=3, ysty=3, $
      xtitle='Redshift', ytitle='log (\Omega_{SFR}) (Gyr^{-1})'

    if keyword_set(integrate) then begin
       
       sfrdfile = 'sfrd_universe.dat' ; Pegase-HR chokes if there is a path name here!
       openw, lun, sfrdfile, /get_lun
       for jj = 0L, n_elements(sfrd_time)-1L do $
         printf, lun, sfrd_time[jj], sfrd_Myr[jj]*sfrfactor, format='(2F12.5)'
       free_lun, lun

; test
       if keyword_set(test) then begin
;      sfrd_time = findgen((20100.0-100.0)/100.0+1)*100.0
          sfrd_time = 10.0^(findgen((alog10(19.0)-alog10(1.0))/0.01)*0.01+alog10(1.0)) ; [Myr]
          sfrd = sfrd_time*0.0+1E-8                                                    ; [M_sun/Myr]
          openw, lun, sfrdfile, /get_lun
          for jj = 0L, n_elements(sfrd_time)-1L do $
            printf, lun, sfrd_time[jj], sfrd[jj], format='(2F12.7)'
          free_lun, lun
       endif
       
; build the Pegase files; see also WRITE_PEGASE_SFH_GRID

       imf = 'salp'
       library = '1'
       library_str = 'basel'
       
       sspsfile = gridpath+imf+'_'+library_str+'_SSPs.dat'
       scenarios_infile = mzpath+'scenarios_sfrd_universe.input'    ; scenarios.f input
       scenarios_outfile = repstr(scenarios_infile,'.input','.output') ; scenarios.f output
       scenarios_logfile = repstr(scenarios_infile,'.input','.logfile')
       spectra_logfile = repstr(scenarios_logfile,'scenarios','spectra')

       if file_test(pegfile,/regular) then begin
          splog, 'Deleting '+pegfile
          spawn, '/bin/rm '+pegfile, /sh
       endif
       
       openw, lun, scenarios_infile, /get_lun ; scenarios.f input file
       printf, lun, scenarios_outfile         ; scenarios.f output file
       printf, lun, sspsfile                  ; SSPs pre-computed by SSPs.f
       printf, lun, '0.05'                    ; default binary fraction
       printf, lun, library                   ; stellar library
       printf, lun, pegfile                   ; output file
       printf, lun, '0.0'                     ; metallicity of the ISM at t=0
       printf, lun, 'n'                       ; infall - no
       printf, lun, '-1'                      ; star formation scenario (here, an input file)
       printf, lun, sfrdfile                  ; input file name
       printf, lun, 'y'                       ; consistent evolution of the stellar metallicity
       printf, lun, '0.0'                     ; mass fraction of substellar objects formed
       printf, lun, 'n'                       ; no galactic winds
       printf, lun, 'y'                       ; include nebular emission
       printf, lun, '0'                       ; no extinction
       printf, lun, 'end'                     ; exit
       free_lun, lun

       if file_test(scenarios_outfile,/regular) then begin
          splog, 'Deleting '+scenarios_outfile
          spawn, '/bin/rm '+scenarios_outfile, /sh
       endif

; generate the models

       t0 = systime(1)
       splog, 'Integrating the SFRD of the universe'
       spawn, 'scenarios_HR < '+scenarios_infile+' > '+scenarios_logfile, /sh     ; run scenarios.f
       spawn, 'echo "'+scenarios_outfile+'" | '+'spectra_HR > '+spectra_logfile, /sh ; run spectra.f
       splog, format='("Total time = ",G0," minutes.")', (systime(1)-t0)/60.0

    endif 
       
; do it

    jj = mrdfits(pegfile,2)
    peg = im_read_peg(pegfile)
    peg = im_measure_peg(peg,/silent)
    rev = reverse(sort(peg.age)) ; reverse the time array!
    peg = peg[rev]
    peg_zaxis = getredshift(peg.age/1D3+getage(zform))
    peg_lookback = getage(0.0)-getage(peg_zaxis)

    peg_good = where(peg_zaxis gt 0.0)
    peg_zaxis = peg_zaxis[peg_good]
    peg_lookback = peg_lookback[peg_good]
    peg = peg[peg_good]

    wilk = rsex(litpath+'08wilkins.sex')
    nwilk = n_elements(wilk)    
    zz = fltarr(nwilk)
    zzerr = fltarr(nwilk)
    for ii = 0L, nwilk-1L do begin
       zz[ii] = mean([wilk[ii].zmin,wilk[ii].zmax])
       zzerr[ii] = stddev([wilk[ii].zmin,wilk[ii].zmax])
    endfor

; Omega(SFR)
    djs_plot, peg_zaxis, alog10(yr2Gyr*peg.sfr/sfrfactor/rhocrit), $
      xsty=3, ysty=3, xtitle='Redshift', ytitle='log (\Omega_{SFR}) (Gyr^{-1})', $
      yrange=[-2.5,-5.0]
    djs_oplot, sfrd_zaxis, alog10(yr2Gyr*sfrd/rhocrit), color='green'
; Omega(Mass)
    djs_plot, peg_zaxis, alog10(peg.mstar/sfrfactor/rhocrit), $
      xsty=3, ysty=3, xtitle='Redshift', ytitle='log (\Omega_{*})', $
      yrange=[-4.5,-2.0]
    djs_oplot, zz, alog10(wilk.omega), psym=symcat(16)
stop
    notzero = where(peg.zgas gt 0.0)
    djs_plot, peg_zaxis[notzero], alog10(peg[notzero].zgas*peg[notzero].mgas/sfrfactor/rhocrit), xsty=3, ysty=3, $
      xtitle='Redshift', ytitle='log (\Omega_{Z})', $
      yrange=[-6,-3.8]
    
stop    
    
    djs_plot, peg.age, peg.zgas, xsty=3, ysty=3

stop    
    
return
end
    
