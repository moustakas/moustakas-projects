function lir_qpint1d_func, x, wave=wave, flux=flux
return, interpol(flux,wave,x)
end

pro build_cosmicimf_pegase, imfs=imfs, ssps=ssps, sfhgrid=sfhgrid, $
  sfrd=sfrd, properties=properties, qaplots=qaplots, witt=witt, $
  debug=debug
; jm10mar12ucsd - build the Pegase models for various IMF choices
    
    outpath = ages_path(/projects)+'cosmicimf/'
    pegpath = outpath+'pegase/'
    hrpath = getenv('PEGASE_HR_DIR')+'/data/user_defined/'

; define some IMF variables; all the computations are for Salpeter
; plus a grid of high-mass slopes plus (and fixed low-mass slope)
    alpha2 = [2.35,cosmicimf_slope()] ; high-mass slopes
    nimf = n_elements(alpha2)
    imf = ['Salpeter','cosmicimf_'+string(alpha2[1:nimf-1],format='(F4.2)')]

; specify the SFH grid
;   tau = [0.0,2.0,10.0] ; test grid
;   const = [0.1,1.0]
;   time = [0.01,0.1,1.0,15.0] ; test grid
;   istau = fix([const*0])
;   tau = im_array(2.0,5.0,0.5)   ; [Gyr]
    tau = im_array(0.0,10.0,0.5)   ; [Gyr]
    const = [0.01,0.1,1.0,15.0]    ; [Gyr]
    istau = fix([tau*0+1,const*0])

;   sfhfile_root = 'tau_'+tau2string(tau)
;   sfhfile_root = ['const_'+const2string(const)]+'.fits'
    sfhfile_root = ['tau_'+tau2string(tau),'const_'+const2string(const)]+'.fits'

;   time = tau
    time = [tau,const] ; characteristic time for star formation
    time_str = strtrim(string(1E3*time,format='(F12.1)'),2) ; [Myr]
    nsfh = n_elements(time)

; we want the constant SFR models for 1 M_sun/yr but that causes
; numerical problems; so use 1E-4 M_sun/Myr = 1E-10 M_sun/yr when
; building the grid, and then scale the output quantities by 1E4 (see
; cosmicimf_read_peg)
    sfr_const = 1D-4 ; [M_sun/Myr]
    sfr_const_str = strtrim(string(sfr_const,format='(E12.1)'),2)
    
; --------------------------------------------------
; build the IMF files that Pegase needs (see the WRITE_PEGASE_IMF
; documentation) 
    if keyword_set(imfs) then begin
;      alpha01 = [0.3,1.3] ; from Wilkins+08
;      mass = [0.01,0.08,0.5,150.0]
       alpha01 = [1.3]
       mass = [0.1,0.5,120.0]
; build all the IMF files except Salpeter, which already exists
       imffile = hrpath+'IMF_'+imf+'.dat' 
       for ii = 1, nimf-1 do begin 
          slope = 1.0-[alpha01,alpha2[ii]]
          write_pegase_imf, mass, slope, outfile=imffile[ii]
       endfor
    
; be sure the copy of the *original* IMF file exists (list_IMFs.dat)
; because we are going to be overwriting it
       if (file_test(hrpath+'list_IMFs.original.dat') eq 0) then begin
          splog, 'Please manually copy/backup the *original* list_IMF.dat file'
          return
       endif

; write the new set of IMFs plus Salpeter to the list_IMFs file 
       outfile = hrpath+'list_IMFs.dat'
       splog, 'Writing '+outfile
       openw, lun, outfile, /get_lun
       printf, lun, '! File written by BUILD_COSMICIMF_PEGASE on '+im_today()
       for ii = 0, nimf-1 do printf, lun, file_basename(imffile[ii])
       free_lun, lun
    endif

; --------------------------------------------------
; convolve the stellar models with each IMF
    if keyword_set(ssps) then begin
       for ii = 0, nimf-1 do begin
          imfpath = pegpath+imf[ii]+'/'
          if (file_test(imfpath,/dir) eq 0) then $
            spawn, 'mkdir -p '+imfpath, /sh
; clean up old files
          allfiles = file_search(imfpath+imf[ii]+'_*.dat',count=nall)
          if (nall ne 0) then rmfile, allfiles
          sspsfile = imfpath+imf[ii]+'_SSPs.dat'
          infile = repstr(sspsfile,'.dat','.input')
          logfile = repstr(sspsfile,'.dat','.log')
          
          openw, lun, infile, /get_lun            ; SSPs_HR.f input file
          printf, lun, string(ii+1,format='(I0)') ; IMF number
          printf, lun, '0.1'           ; lowest mass M_sun
          printf, lun, '120.0'         ; highest mass M_sun
          printf, lun, 'B'             ; SNII ejecta model of Woosley & Weaver
          printf, lun, 'y'             ; stellar winds
          printf, lun, '1'             ; Basel stellar library
          printf, lun, imf[ii]         ; prefix
          printf, lun, 'end'
          free_lun, lun

          t0 = systime(1)
          splog, 'Building SSPs for IMF: '+imf[ii]
          pushd, imfpath
          spawn, 'SSPs_HR < '+infile+' > '+logfile, /sh ; run SSPs_HR.f
          popd
          splog, format='("Time = ",G0," minutes")', (systime(1)-t0)/60.0
       endfor
    endif

; --------------------------------------------------
; build some SFHs that fit the observed star formation rate density
; (SFRD) evolution, for each IMF
    if keyword_set(sfrd) then begin
; get the observed, minimum, and maximum SFRD evolution; remember the
; normalization is arbitrary; compute the SFRD for two different
; formation redshifts: ZFORM=1 and ZFORM=5; in the latter case assume
; ZCUT=1 
       zcut = [1.2,1.0]
       zform = [1.2,5.0]
;      for iz = 0, 0 do begin
       for iz = 1, 1 do begin
;      for iz = 0, 1 do begin
          sfrd = 10^get_sfrd_evolution(zcut=zcut[iz],zform=zform[iz],zaxis=zaxis)
          sfrdmin = 10^get_sfrd_evolution(zcut=zcut[iz],zform=zform[iz],/sfrdmin)
          sfrdmax = 10^get_sfrd_evolution(zcut=zcut[iz],zform=zform[iz],/sfrdmax)
          taxis1 = (getage(zaxis)-getage(zform[iz]))*1E3 ; [Myr]

;         sfrd = 1D3*sfrd/int_tabulated(taxis1,sfrd) ; [M_sun/Myr???]
;         sfrdmin = 1D3*sfrdmin/int_tabulated(taxis1,sfrdmin)
;         sfrdmax = 1D3*sfrdmax/int_tabulated(taxis1,sfrdmax)

          sfrd = sfrd/qpint1d('lir_qpint1d_func',min(taxis1),$
            max(taxis1),functargs={wave:taxis1,flux:sfrd})
          sfrdmin = sfrdmin/qpint1d('lir_qpint1d_func',min(taxis1),$
            max(taxis1),functargs={wave:taxis1,flux:sfrdmin})
          sfrdmax = sfrdmax/qpint1d('lir_qpint1d_func',min(taxis1),$
            max(taxis1),functargs={wave:taxis1,flux:sfrdmax})
          
;         sfrd_const = 100D/max(sfrd) ; [M_sun/yr; see notes on SFR_CONST, above] 
;         sfrd_const = 1D-4/max(sfrd)  ; [M_sun/Myr; see notes on SFR_CONST, above] 
;         sfrd = sfrd_const*sfrd
;         sfrdmin = sfrd_const*sfrdmin
;         sfrdmax = sfrd_const*sfrdmax
          
; expand the time axis so that it extends to 20 Gyr with the SFR equal
; to zero beyond the maximum age (this is to avoid numerical problems
; in Pegase, since spectra are written out up to 20 Gyr after t=0)
;         taxis = taxis1
;         final_sfrd = sfrd
;         final_sfrdmin = sfrdmin
;         final_sfrdmax = sfrdmax
          
          moretime = range(max(taxis1)*1.05,20E3,50,/log)
          taxis = [reverse(moretime),taxis1]
          final_sfrd = [moretime*0.0,sfrd]
          final_sfrdmin = [moretime*0.0,sfrdmin]
          final_sfrdmax = [moretime*0.0,sfrdmax]
;         final_sfrd = [moretime*0.0+min(sfrd),sfrd]
;         final_sfrdmin = [moretime*0.0+min(sfrdmin),sfrdmin]
;         final_sfrdmax = [moretime*0.0+min(sfrdmax),sfrdmax]
          
; pack everything together          
          allsfrd = [[final_sfrd],[final_sfrdmin],[final_sfrdmax]]
;         djs_plot, taxis/1E3, final_sfrdmax, line=2, xsty=3, ysty=3
;         djs_oplot, taxis/1E3, final_sfrdmin, line=2
;         djs_oplot, taxis/1E3, final_sfrd, line=0

; write out the SFR scenarios for each IMF
;         for ii = 0, 10 do begin
;         for ii = 11, nimf-1 do begin
;         for ii = 0, 0 do begin
          for ii = 0, nimf-1 do begin
             imfpath = pegpath+imf[ii]+'/'
             pushd, imfpath
             sfrdfile = 'sfrd_zform'+string(zform[iz],format='(F3.1)')+$
               '_'+['best','min','max']+'.dat'
             sfrd_sfhfile = repstr(sfrdfile,'.dat','.fits')
             for jj = 0, n_elements(sfrdfile)-1 do begin
                splog, 'Writing '+sfrdfile[jj]
                openw, lun, sfrdfile[jj], /get_lun
; write the file out in reverse order
                for kk = n_elements(taxis)-1, 0, -1 do $ 
                  printf, lun, taxis[kk], allsfrd[kk,jj], $
                  format='(3x,F8.2,3x,E12.5)'
                free_lun, lun
             endfor

; now build the spectra       
             sspsfile = imf[ii]+'_SSPs.dat'
             if (file_test(sspsfile) eq 0) then message, splog, sspsfile+' not found!'

             scenarios_infile = 'scenarios_sfrd_zform'+$
               string(zform[iz],format='(F3.1)')+'_'+imf[ii]+'.input' ; scenarios.f input
             scenarios_outfile = repstr(scenarios_infile,'.input','.output') ; scenarios.f output
             scenarios_logfile = repstr(scenarios_infile,'.input','.log')
             spectra_logfile = repstr(scenarios_logfile,'scenarios','spectra')
             
             openw, lun, scenarios_infile, /get_lun ; scenarios.f input file
             printf, lun, scenarios_outfile         ; scenarios.f output file
             printf, lun, sspsfile                  ; SSPs pre-computed by SSPs.f
             printf, lun, '0.05'                    ; default binary fraction
             printf, lun, '1'                       ; Basel stellar library
             
;            for jj = 0, 0 do begin
             for jj = 0, n_elements(sfrdfile)-1 do begin
                sfhfile = sfrd_sfhfile[jj]
                if file_test(sfhfile) then spawn, '/bin/rm '+sfhfile, /sh
                printf, lun, sfhfile              ; spectra.f output file name
                printf, lun, '0.02'               ; metallicity of the ISM at t=0 (assume "solar")
                printf, lun, 'n'                  ; infall - no
                printf, lun, '-1'                 ; star formation scenario: file giving the SFR
                printf, lun, strtrim(sfrdfile[jj],2) ; file name
                printf, lun, 'n'                     ; consistent evolution of the stellar metallicity - no
                printf, lun, '0.02'                  ; (fixed!) stellar metallicity
                printf, lun, '0.0'                   ; mass fraction of substellar objects formed
                printf, lun, 'n'                     ; no galactic winds
                printf, lun, 'y'                     ; include nebular emission
                printf, lun, '0'                     ; no extinction
             endfor 
             printf, lun, 'end' ; exit
             free_lun, lun

; clean up files    
             if file_test(scenarios_outfile) then spawn, '/bin/rm '+scenarios_outfile, /sh
; generate the models
             t0 = systime(1)
             splog, 'Building the SFRD for IMF '+imf[ii]
             spawn, 'scenarios_HR < '+scenarios_infile+' > '+scenarios_logfile, /sh     ; run scenarios.f
             spawn, 'echo "'+scenarios_outfile+'" | '+'spectra_HR > '+spectra_logfile, /sh ; run spectra.f
             splog, format='("Time = ",G0," minutes.")', (systime(1)-t0)/60.0
             popd
          endfor                ; close IMF loop
       endfor                   ; close ZFORM loop
    endif

stop       

; --------------------------------------------------
; build the SFH grid
    if keyword_set(sfhgrid) then begin
; now loop on each IMF       
;      for ii = 0, 0 do begin
       for ii = 0, nimf-1 do begin
          imfpath = pegpath+imf[ii]+'/'
          sspsfile = imfpath+imf[ii]+'_SSPs.dat'
          if (file_test(sspsfile) eq 0) then message, splog, sspsfile+' not found!'

          scenarios_infile = imfpath+'scenarios_'+imf[ii]+'.input'        ; scenarios.f input
          scenarios_outfile = repstr(scenarios_infile,'.input','.output') ; scenarios.f output
          scenarios_logfile = repstr(scenarios_infile,'.input','.log')
          spectra_logfile = repstr(scenarios_logfile,'scenarios','spectra')
    
          openw, lun, scenarios_infile, /get_lun ; scenarios.f input file
          printf, lun, scenarios_outfile         ; scenarios.f output file
          printf, lun, sspsfile                  ; SSPs pre-computed by SSPs.f
          printf, lun, '0.05'                    ; default binary fraction
          printf, lun, '1'                       ; Basel stellar library
    
          for jj = 0, nsfh-1 do begin
             sfhfile = imfpath+sfhfile_root[jj]
             if file_test(sfhfile) then spawn, '/bin/rm '+sfhfile, /sh
             printf, lun, sfhfile ; spectra.f output file name
; treat the tau- and continuous star formation scenarios differently 
             if (istau[jj] eq 1) then begin
                if (time[jj] eq 0.0) then begin ; tau=0
                   printf, lun, '0.02' ; metallicity of the ISM at t=0 (assume "solar")
                   printf, lun, 'n'    ; infall - no
                   printf, lun, '0'    ; star formation scenario: instantaneous burst
                endif else begin
                   printf, lun, '0.02'       ; metallicity of the ISM at t=0 (for no consistent evolution)
                   printf, lun, 'n'          ; infall - no
                   printf, lun, '2'          ; star formation scenario: e.g., exponentially declining
                   printf, lun, time_str[jj] ; characteristic time for star formation [tau or p1, Myr]
                   printf, lun, '1.0'        ; star formation law pre-factor [p2]
                endelse
             endif 
             if (istau[jj] eq 0) then begin ; continuous star formation
                printf, lun, '0.02'        ; metallicity of the ISM at t=0
                printf, lun, 'n'           ; infall - no
                printf, lun, '1'           ; star formation scenario: continuous star formation
                printf, lun, sfr_const_str ; star formation rate [p1, M_sun/Myr]
                printf, lun, time_str[jj]  ; star formation duration [p2, Myr]
             endif
             printf, lun, 'n'    ; consistent evolution of the stellar metallicity - no
             printf, lun, '0.02' ; (fixed!) stellar metallicity
             printf, lun, '0.0'  ; mass fraction of substellar objects formed
             printf, lun, 'n'    ; no galactic winds
             printf, lun, 'y'    ; include nebular emission
             printf, lun, '0'    ; no extinction
;            printf, lun, '2'    ; inclination-averaged extinction for a disk geometry
          endfor 
          printf, lun, 'end'    ; exit
          free_lun, lun

; clean up files    
          if file_test(scenarios_outfile) then spawn, '/bin/rm '+scenarios_outfile, /sh
; generate the models
          t0 = systime(1)
          splog, 'Building the SFH grid for IMF '+imf[ii]
          spawn, 'scenarios_HR < '+scenarios_infile+' > '+scenarios_logfile, /sh     ; run scenarios.f
          spawn, 'echo "'+scenarios_outfile+'" | '+'spectra_HR > '+spectra_logfile, /sh ; run spectra.f
          splog, format='("Time = ",G0," minutes.")', (systime(1)-t0)/60.0
       endfor
    endif

; --------------------------------------------------
; collect the properties we care about; treat the Salpeter IMF
; separately as a baseline comparison
    if keyword_set(properties) then begin
       out = {$
         imf:               '',$
         slope:            0.0,$

         C_ha:             0.0,$ ; conversion from L(Ha) (erg/s) to SFR(Ha) (M_sun/yr); K98 = 7.9D-42
         C_1500:           0.0,$ ; K98: 7.0D-44 at 1500 A (or 1.4D-28 for erg/s/Hz)
         C_1500_err:       0.0,$ ; scatter due to different burst time scales; 
         C_ir:             0.0,$ ; conversion from L(IR) (erg/s) to SFR(IR) (M_sun/yr); K98: 4.5D-44
         C_ir_err:         0.0,$
         C_24_lo:          0.0,$ ; see eqs. (10) and (11) in Rieke+09
         C_24_lo_err:      0.0,$
         C_24_hi:          0.0,$
         C_24_hi_err:      0.0,$

         dmsalp:           0.0,$ ; M(salpeter)=M(IMF)+DMSALP
         dmsalp_err:       0.0,$
         dmsalp_ssp:       0.0,$ ; limiting case of tau=0
         rfrac:            0.0,$ ; return fraction for a range of tau-models
         rfrac_err:        0.0,$
         rfrac_ssp:        0.0}  ; return fraction for tau=0 at 13 Gyr (SSP limiting case)
       out = replicate(out,nimf)
       out.imf = imf
       out.slope = alpha2
    
; use the <1 Gyr continuous star formation models to get C_Ha, C_IR,
; C24, and C_1500
       iconst = where(strmatch(sfhfile_root,'*const*',/fold) and $
         (time le 1.0),nconst)
       sfhfile_const = sfhfile_root[iconst]
       time_const = time[iconst]

; loop on each IMF       
       for ii = 0, nimf-1 do begin
          splog, 'IMF '+imf[ii]
          imfpath = pegpath+imf[ii]+'/'
          C_1500 = dblarr(nconst)
          C_ir = dblarr(nconst)
          for jj = 0, nconst-1 do begin
             peg = cosmicimf_read_peg(imfpath+sfhfile_const[jj],$
               sfr_const=sfr_const,age=time_const[jj]*1E3)
             l1500 = 1500.0*interpolate(transpose(peg.flux),findex(peg[0].wave,1500.0)) ; [erg/s]
             lir = qpint1d('lir_qpint1d_func',min(peg.wave),10D4,$
               functargs={wave:peg.wave,flux:peg.flux})
; note that the H-alpha conversion factor does not depend on age, so
; just store the constant factor from the first continuous star
; formation model; for the UV take the conversion factor at the
; maximum age of the burst
             if (jj eq 0) then out[ii].c_ha = 1D/peg.lineflux[1]
             C_1500[jj] = 1D/l1500
             C_ir[jj] = 1D/lir
; L(1500) debugging plot             
             if keyword_set(debug) then begin
                if (jj eq 0) then begin
                   djs_plot, [0], [0], /nodata, xrange=[1,1E3], $
                     yr=1D/[1D41,5D44], /xlog, /ylog, xsty=3, ysty=3
                endif
                djs_oplot, peg.age, l1500, psym=-6
                djs_oplot, 1.0*[1,1], 10^!y.crange, color='red'
                djs_oplot, time_const[jj]*1E3*[1,1], 10^!y.crange, color='red'
                djs_oplot, 10^!x.crange, alog10(C_1500[jj])*[1,1], color='orange'
                cc = get_kbrd(1)
;; H-alpha debugging plot             
;                if (jj eq 0) then begin
;                   djs_plot, [0], [0], /nodata, xrange=[1,1E3], $
;                     yr=[1D38,1D42], /xlog, /ylog, xsty=3, ysty=3
;                endif
;                djs_oplot, peg.age, peg.lineflux[1], psym=-6
             endif 
          endfor ; close SFH loop
          
; compute the mean UV and IR conversion factors; the errors are due to
; assuming different star formation timescales (0.01-1 Gyr), and are
; ~0.1 dex
          out[ii].C_1500 = djs_mean(C_1500)
          out[ii].C_1500_err = djsig(C_1500)
          out[ii].C_ir = djs_mean(C_ir)
          out[ii].C_ir_err = djsig(C_ir)
; finally get the L(24) SFR conversion factor by scaling eqs. (10) and
; (11) in Rieke+09, recalling that there is a conversion factor for
; low- (transparent) and high-luminosity (opaque) starburst galaxies
          scale = djs_mean([out[ii].c_ir/4.5D-44,out[ii].c_ha/7.9D-42])
          scale_err = djsig([out[ii].c_ir/4.5D-44,out[ii].c_ha/7.9D-42])
          factor = scale/0.66;/3.826D33
          factor_err = scale_err/0.66;/3.826D33
          out[ii].C_24_lo = 7.8D-10*factor
          out[ii].C_24_hi = 7.76D-11*factor
          out[ii].C_24_lo_err = 7.8D-10*factor_err
          out[ii].C_24_hi_err = 7.76D-11*factor_err

; next, use the tau-models to get the return fraction; consider
; tau-models with tau>2 Gyr and ages between 5-12 Gyr, which
; corresponds to observing galaxies at z=0-0.6 which formed at z=2-3
; these = where((time ge 8.0),nthese)
          these = where(istau and (time gt 2.0),nthese)
;         these = lindgen(nsfh) & nthese = nsfh
          rfrac = fltarr(nthese)
          for jj = 0, nthese-1 do begin
             peg = cosmicimf_read_peg(imfpath+sfhfile_root[these[jj]],$
               sfr_const=sfr_const)
             inrange = where((peg.age/1E3 gt 5.0) and (peg.age/1E3 lt 12.0))
             rfrac[jj] = djs_median(peg[inrange].mgas)
; debugging plot
             if keyword_set(debug) then begin
                if (jj eq 0) then begin
                   djs_plot, [0], [0], /nodata, $
                     xr=[0,20], yr=[0,1.05], xsty=3, ysty=3
                   djs_oplot, 5.0*[1,1], !y.crange, color='red'
                   djs_oplot, 12.0*[1,1], !y.crange, color='red'
                endif
                djs_oplot, peg.age/1E3, peg.mgas, color='orange'
                djs_oplot, peg.age/1E3, peg.mstar, color='blue'
                djs_oplot, !x.crange, rfrac[jj]*[1,1], color='green'
;               if (jj eq nthese-1) then cc = get_kbrd(1)
             endif
          endfor ; SFH
; also compute the limiting return fraction of a 13 Gyr old SSP (tau=0)
          issp = (where(time eq 0.0))[0]
          ssp = cosmicimf_read_peg(imfpath+$
            sfhfile_root[issp],age=13E3)
          out[ii].rfrac = djs_mean(rfrac)
          out[ii].rfrac_err = djsig(rfrac)
          out[ii].rfrac_ssp = ssp.mgas
       endfor ; IMF

; compare my UV/Salpeter conversion factor with K98:
;      print, 10D^(-out[0].c_1500)*im_light(/ang)/1500.0

; build some QAplots (eventually move these below)       
       struct_print, out
;      ploterror, out.slope, out.c_1500, out.c_1500_err, psym=6, xsty=3, ysty=3
       ploterror, out.slope, out.rfrac, out.rfrac_err, psym=6, $
         xsty=3, ysty=3, yr=[0,1], xtitle=textoidl('IMF Slope (0.5-120 M_{\odot})'), $
         ytitle=textoidl('Return Fraction R(\xi)')
       djs_oplot, out.slope, out.rfrac_ssp, psym=4, color='green'
       oploterror, out[1:nimf-1].slope, out[1:nimf-1].rfrac, out[1:nimf-1].rfrac_err, $
         psym=6, sym=2.5, color=djs_icolor('orange'), errcolor=djs_icolor('orange')
       plots, out[0].slope, out[0].rfrac_ssp, psym=6, $
         sym=2.5, color=djs_icolor('blue')

; compute the conversion factor to Salpeter using the tau-models; use
; the M/L ratio, which does not depend on initial mass, etc.; this
; conversion is *relatively* insensitive to age, so just use the mean
; and standard deviations of the values at 13 Gyr across all the SFHs
       these = where(istau,nthese)
       dmsalp = fltarr(nimf,nthese)
       for jj = 0, nthese-1 do begin
          for ii = 0, nimf-1 do begin
             imfpath = pegpath+imf[ii]+'/'
             peg = cosmicimf_read_peg(imfpath+sfhfile_root[these[jj]],$
               sfr_const=sfr_const)
             nage = n_elements(peg.age)
             if (ii eq 0) then begin
; get the Salpeter M/L ratio for this SFH
                ml_age_salp = fltarr(nage)
                for aa = 0L, nage-1 do ml_age_salp[aa] = -2.5*alog10($
                  k_project_filters(k_lambda_to_edges(peg[aa].wave),$
                  peg[aa].flux/(4*!dpi*(10.0*3.085678D18)^2),$
                  filterlist='bessell_V.par'))
             endif
; get the M/L ratio for this IMF and SFH
             ml_age = fltarr(n_elements(peg.age))
             for aa = 0L, n_elements(peg.age)-1 do ml_age[aa] = $
               -2.5*alog10(k_project_filters(k_lambda_to_edges(peg[aa].wave),$
               peg[aa].flux/(4*!dpi*(10.0*3.085678D18)^2),$
               filterlist='bessell_V.par'))
             dmsalp_age = 0.4*(ml_age_salp-ml_age)                 ; M(salpeter)=M(IMF)+DMSALP
             dmsalp[ii,jj] = interpol(dmsalp_age,peg.age/1E3,13.0) ; at 13 Gyr
;            splog, imf[ii], dmsalp[ii,jj]
; debugging plot
             if keyword_set(debug) then begin
                if (ii eq 0) then djs_plot, [0], [0], /nodata, $
                  xrange=[1,15], yrange=[-1,1], xsty=3, ysty=3 else $
                 djs_oplot, peg.age/1E3, dmsalp_age
                cc = get_kbrd(1)
             endif
          endfor ; IMF
       endfor    ; SFH

       out.dmsalp_ssp = reform(dmsalp[*,0]) ; limiting case of an SSP
       for ii = 0, nimf-1 do begin
          out[ii].dmsalp = djs_mean(dmsalp[ii,*])
          out[ii].dmsalp_err = djsig(dmsalp[ii,*])
       endfor

;      struct_print, out
;      ploterror, out.slope, out.dmsalp, out.dmsalp_err, psym=6, $
;        xsty=3, ysty=3, yr=[-0.5,+0.2], xtitle=textoidl('IMF Slope (0.5-120 M_{\odot})'), $
;        ytitle=textoidl('log (\Delta M)')
;      djs_oplot, out.slope, out.dmsalp_ssp, psym=4, color='green'
;      plots, out[0].slope_salp, 0.0, psym=6, sym=2.5, color=djs_icolor('orange')
       
; write out
       outfile = outpath+'pegase_sfrs.fits'
       im_mwrfits, out, outfile, /clobber
       
    endif
       
; --------------------------------------------------
; build the A(lambda) vs lambda look-up table using the Witt & Gordon
; (2000) dust models and the flux-ratio method of Gordon et al. (2000)
    if keyword_set(witt) then begin

; choose three metallicities four ages, and five a_d values 
       read_many_fit, dust_info, sf_type='c' ; continuous star formation
       metal = [-0.4,0.0,0.4]
       ages = [10.0,50.0,100.0,1000.0]
       ad = [0.01,0.25,0.5,0.75,0.95]
       ndust = n_elements(metal)*n_elements(ages)*n_elements(ad)

       irx = range(0.1,1D6,1D3,/log)
       nirx = n_elements(irx)
       witt_dust = {irx: float(irx), a1500: float(irx*0.0), $
         a1500_grid: fltarr(nirx,ndust), params: strarr(ndust)}

       counter = 0
       for im = 0, n_elements(metal)-1 do begin
          for ia = 0, n_elements(ages)-1 do begin
             for id = 0, n_elements(ad)-1 do begin
                witt_dust.params[counter] = 'Z'+repstr(string(metal[im],format='(F4.1)'),' ','+')+$
                  '_Age'+string(ages[ia],format='(I4.4)')+'_ad'+string(ad[id],format='(F4.2)')
                for jj = 0L, nirx-1L do witt_dust.a1500_grid[jj,counter] = fr2atten(irx[jj],$
                  '',1500.0,ages[ia],metal[im],ad[id],dust_info,/silent)
                counter++
             endfor
          endfor
       endfor

; compute the media A(1500)-IRX relation and write out
       witt_dust.a1500 = total(witt_dust.a1500_grid,2)/float(ndust)
;      for jj = 0L, nirx-1L do witt_dust.a1500[jj] = djs_median(witt_dust.a1500_grid[jj,*])

       outfile = outpath+'witt_dust.fits'
       im_mwrfits, witt_dust, outfile, /clobber
       
       if keyword_set(debug) then begin
          djs_plot, 1.0+irx, witt_dust.a1500_grid[*,0], xr=[1,1D3], yr=[0,6], /xlog, ps=3
          for jj = 1, ndust-1 do djs_oplot, 1.0+irx, witt_dust.a1500_grid[*,jj], ps=3
          djs_oplot, 1.0+irx, witt_dust.a1500, color='red', thick=3
       endif
       
;; ###########################################################################       
;; all of the code below works, for the most part; it turns out that
;; the IMF does *not* matter in the computation of A(1500) from IRX, so
;; we just use Karl's results, above
;       
;       geometry = ['dusty','cloudy','shell']
;       dust_type = ['SMC','MW']
;       structure = ['c','h']
;       ng = n_elements(geometry)
;       nd = n_elements(dust_type)
;       ns = n_elements(structure)
;       ndust = ng*nd*ns
;       ntauv = 25 ; see READ_WITT_DUSTMODELS()
;
;; just focus on the young (<1 Gyr) continuous SFH models
;       iconst = where(strmatch(sfhfile_root,'*const*',/fold) and $
;         (time le 1.0),nconst)
;       sfhfile_const = sfhfile_root[iconst]
;       time_const = time[iconst]
;       
;; initialize the output data structure; compute the mean relation
;; between A(1500) and IRX for the full set of dust models, and for all
;; the continuous star formation models
;       irx = range(0.1,1D6,1D4,/log)
;       nirx = n_elements(irx)
;       witt_dust = {imf: '', alpha2: 0.0, coeff: fltarr(9), $
;         a1500: irx*0.0, irx: irx}
;       witt_dust = replicate(witt_dust,nimf)
;       witt_dust.imf = imf
;       witt_dust.alpha2 = alpha2
;       
;       for ii = 0, nimf-1 do begin
;          splog, 'IMF '+imf[ii]
;          imfpath = pegpath+imf[ii]+'/'
;          for jj = 0, nconst-1 do begin
;; read the Pegase model and then iterate on all the dust parameters 
;             peg = cosmicimf_read_peg(imfpath+sfhfile_const[jj],$
;               sfr_const=sfr_const,age=time_const[jj]*1E3)
;             out = replicate({a1500: 0.0D, irx: 0.0D},ng,nd,ns,ntauv)
;             for ig = 0, ng-1 do begin
;                for id = 0, nd-1 do begin
;                   for is = 0, ns-1 do begin
;                      model = read_witt_dustmodels(geometry=geometry[ig],$
;                        dust_type=dust_type[id],structure=structure[is])
;;                     plot, model[2].wave, model[2].alambda, /xlog, ps=-6, $
;;                       /ylog, xr=[800,3E4], xsty=3, ysty=3
;                      for it = 0, ntauv-1 do begin
;; compute the attenuation at 1500 A, for this model
;                         out[ig,id,is,it].a1500 = interpol(model[it].alambda,model[it].wave,1500.0)
;; compute the attenuation L(1500) and integrate the intrinsic minus
;; the attenuated flux from 0.912 to 3 microns as an estimate of L(IR),
;; for this model
;                         linterp, model[it].wave, model[it].alambda, peg.wave, alambda, missing=0.0
;;                        alambda = interpol(model[it].alambda,model[it].wave,peg.wave)
;;                        djs_plot, peg.wave, alambda, psym=-6, xsty=3, ysty=3                            
;                         flux = peg.flux*10D^(-0.4*alambda)
;                         l1500 = interpol(flux,peg.wave,1500.0)*1500D ; [erg/s]
;                         lir = qpint1d('lir_qpint1d_func',912D,3D4,$
;                           functargs={wave:peg.wave,flux:peg.flux-flux})
;; cap the infrared excess at 10^6!
;                         out[ig,id,is,it].irx = lir/l1500<1D6 ; keep *linear*
;; debugging plot                            
;;                        if keyword_set(debug) and (ig eq 1) then begin
;;                           djs_plot, peg.wave, peg.flux, psym=10, xr=[800,1E4], xsty=3, ysty=3, /xlog
;;                           djs_oplot, peg.wave, flux, psym=10, color='red'
;;                           djs_oplot, peg.wave, peg.flux-flux, psym=10, color='yellow'
;;                           print, out[ig,id,is,it].a1500, l1500, lir, lir/l1500
;;                           cc = get_kbrd(1)
;;                        endif
;                      endfor    ; tau_v
;                   endfor       ; structure
;                endfor          ; dust type
;             endfor             ; geometry
;             out = reform(out,ndust*ntauv)
;             keep = where(out.a1500 lt 50.0) ; toss out crazy attenuations
;             out = out[keep]
;             if (jj eq 0) then allout = out else allout = [allout,out]
;             
;;; fit it!  (it's not perfect!)             
;;             coeff = fit_irx_a1500(out.irx,out.a1500,model_irx=witt_dust[ii].irx,$
;;               model_a1500=model_a1500)
;;             witt_dust[ii].a1500_sfh[*,jj] = model_a1500
;;             
;;;            witt_dust[ii].a1500_sfh[*,jj] = interpol(out.a1500,$
;;;              alog10(out.irx),alog10(witt_dust[ii].irx),/spline)
;;             djs_plot, out.irx, out.a1500, psym=6, xsty=3, $
;;               ysty=3, /xlog, xr=[0.1,1D3], yr=[0,6]
;;             djs_oplot, witt_dust[ii].irx, witt_dust[ii].a1500_sfh[*,0], color='red'
;          endfor                ; SFH
;
;; fit the A(1500) vs IRX relation          
;          coeff = fit_irx_a1500(allout.irx,allout.a1500,$
;            model_irx=witt_dust[ii].irx,model_a1500=model_a1500,$
;            debug=0)
;          witt_dust[ii].coeff = coeff
;          witt_dust[ii].a1500 = model_a1500
;          
;          if keyword_set(debug) then begin
;             djs_plot, 1+allout.irx, allout.a1500, psym=6, xsty=1, $
;               ysty=1, xr=[1.0,1D3], yr=[0,6], /xlog
;             djs_oplot, witt_dust[ii].irx, witt_dust[ii].a1500, color='blue'
;          endif
;       endfor                   ; IMF
    endif
    
; --------------------------------------------------
; build a few QAplots comparing the various Pegase models
    if keyword_set(qaplots) then begin

; [1] compare the SFR=constant SEDs as a function of high-mass slope
       psfile = 'test.ps'
       im_plotconfig, 0, pos, psfile=psfile
       for ii = 0, nimf-1 do begin
          imfpath = pegpath+imf[ii]+'/'
          peg = im_read_peg(imfpath+'const_00.10Gyr.fits')
          if (ii eq 0) then begin
             djs_plot, [0], [0], /nodata, /xlog, /ylog, xr=[0.1,1], $
               yr=[5D27,2E31], xsty=3, ysty=3, xtitle='Wavelength (\mu'+'m)', $
               ytitle='Flux (erg s^{-1} \AA^{-1})', position=pos
          endif
          get_element, peg.age, 100.0, indx
          djs_oplot, peg[indx].wave/1D4, peg[indx].flux,$;+peg[indx].lineflux, $
            thick=2.0, color='grey'
       endfor
       im_plotconfig, psfile=psfile, /psclose, /gzip

    endif
       
return
end
