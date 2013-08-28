function tau2string, alltau
    nalltau = n_elements(alltau)
    tau = alltau[uniq(alltau,sort(alltau))]
    ntau = n_elements(tau)
    taustring = strarr(nalltau)

    for ii = 0L, ntau-1L do begin
       these = where(tau[ii] eq alltau)
       if (tau[ii] lt 10.0) then taustring[these] = $
         '0'+string(tau[ii],format='(F3.1)')+'Gyr'
       if (tau[ii] ge 10.0) and (tau[ii] lt 100.0) then $
         taustring[these] = string(tau[ii],format='(F4.1)')+'Gyr'
       if (tau[ii] ge 100.0) then taustring[these] = $
         string(tau[ii],format='(I3.3)')+'Gyr'
    endfor
return, taustring
end

pro mzpegase_build_models, dossps=dossps, doscenarios=doscenarios
; jm10nov04ucsd - generate a set of Pegase models that we will use to
; model the evolution of the LZ and MZ relations; assume the Kroupa+01
; IMF; see also CREATE_PEGASE_SSPS for a complementary piece of code 
    
    hrpath = getenv('PEGASE_HR_DIR')+'/data/user_defined/'
    ssppath = ages_path(/projects)+'mz/pegase/'
    if (file_test(ssppath,/dir) eq 0) then spawn, 'mkdir -p '+ssppath, /sh

; --------------------------------------------------    
; [1] specify the IMF
    imf = 'Kroupa01_100' ; fiducial set
    imfstr = 'kroupa01'
    minmass = 0.01
    maxmass = 100.0
    imffile = hrpath+'IMF_'+imf+'.dat' 
    nimf = n_elements(imf)

; update the list_IMFs.dat file; be sure a copy of the *original* IMF
; file exists because we are going to be overwriting it
    if (file_test(hrpath+'list_IMFs.original.dat') eq 0) then begin
       splog, 'Please manually copy/backup the *original* list_IMF.dat file'
       return
    endif
     
    outfile = hrpath+'list_IMFs.dat'
    splog, 'Writing '+outfile
    openw, lun, outfile, /get_lun
    printf, lun, '! File written by MZPEGASE_BUILD_MODELS on '+im_today()
    for ii = 0, nimf-1 do printf, lun, file_basename(imffile[ii])
    free_lun, lun

; --------------------------------------------------    
; [2] convolve the stellar models with each IMF with a default set of
; optional parameters, if necessary
    if keyword_set(dossps) then begin
       pushd, ssppath
       for ii = 0, nimf-1 do begin
          splog, 'Working on IMF '+imf[ii]
; clean up old files
          allfiles = file_search(ssppath+imf[ii]+'_*.dat',count=nall)
          if (nall ne 0) then rmfile, allfiles
          sspsfile = imf[ii]+'_SSPs.dat'
          infile = repstr(sspsfile,'.dat','.input')
          logfile = repstr(sspsfile,'.dat','.log')
          
          openw, lun, infile, /get_lun         ; SSPs_HR.f input file
          printf, lun, string(ii+1,format='(I0)') ; IMF number
          printf, lun, strtrim(minmass[ii],2)     ; lowest mass M_sun
          printf, lun, strtrim(maxmass[ii],2)     ; highest mass M_sun
          printf, lun, 'B'                        ; SNII ejecta model of Woosley & Weaver
          printf, lun, 'y'                        ; stellar winds
          printf, lun, '1'                        ; Basel stellar library
          printf, lun, imf[ii]                    ; prefix
          printf, lun, 'end'
          free_lun, lun
          
          t0 = systime(1)
          splog, 'Building SSPs for IMF: '+imf[ii]
;         pushd, ssppath
          spawn, 'SSPs_HR < '+infile+' > '+logfile, /sh ; run SSPs_HR.f
;         popd
          splog, format='("Time = ",G0," minutes")', (systime(1)-t0)/60.0
       endfor
       popd
    endif

; --------------------------------------------------    
; [3] call scenarios with a range of tau values 
    if keyword_set(doscenarios) then begin
; define the tau-values between 0 Gyr (SSP) and 13.5 Gyr, the maximum
; age of the universe
;      tau = [range(0.0,10.0,21),range(12.0,100.0,23)]
       tau = [range(0.0,10.0,21),range(12.0,100.0,5)]
       ntau = n_elements(tau)
       tau_str = strtrim(string(1E3*tau,format='(F12.1)'),2) ; [Myr]

       timf = systime(1)
       pushd, ssppath
       for ii = 0, nimf-1 do begin
          sspsfile = imf[ii]+'_SSPs.dat'
          if (file_test(sspsfile) eq 0) then begin
             splog, sspsfile+' not found!'
             splog, 'Rerun with /DOSSPS!'
             continue
          endif

          scenarios_infile = 'scenarios_'+imf[ii]+'.input'        ; scenarios.f input
          scenarios_outfile = repstr(scenarios_infile,'.input','.output') ; scenarios.f output
          scenarios_logfile = repstr(scenarios_infile,'.input','.log')
          spectra_logfile = repstr(scenarios_logfile,'scenarios','spectra')
    
          openw, lun, scenarios_infile, /get_lun ; scenarios.f input file
          printf, lun, scenarios_outfile         ; scenarios.f output file
          printf, lun, sspsfile                  ; SSPs pre-computed by SSPs.f
          printf, lun, '0.05'                    ; default binary fraction
          printf, lun, '1'                       ; Basel stellar library
          
; build the tau-models
          taufile = imfstr[ii]+'_tau_'+tau2string(tau)+'.fits'
          for jj = 0, ntau-1 do begin
             if file_test(taufile[jj]) then spawn, '/bin/rm '+taufile[jj], /sh
             printf, lun, taufile[jj]    ; spectra.f output file name
             if (tau[jj] eq 0.0) then begin ; tau=0
                printf, lun, '0.0001'       ; metallicity of the ISM at t=0
                printf, lun, 'n'            ; infall - no
                printf, lun, '0'            ; star formation scenario: instantaneous burst
             endif else begin
                printf, lun, '0.0001' ; metallicity of the ISM at t=0 (for consistent evolution)
                printf, lun, 'n'      ; infall - no
                printf, lun, '2'      ; star formation scenario: e.g., exponentially declining
                printf, lun, tau_str[jj] ; characteristic time for star formation [tau or p1, Myr]
                printf, lun, '1.0'       ; star formation law pre-factor [p2]
             endelse
             printf, lun, 'y'   ; consistent evolution of the stellar metallicity - yes
             printf, lun, '0.0' ; mass fraction of substellar objects formed
             printf, lun, 'n'   ; no galactic winds
             printf, lun, 'n'   ; include nebular emission - no
             printf, lun, '0'   ; no extinction
;            printf, lun, '2'      ; inclination-averaged extinction for a disk geometry
          endfor 
          printf, lun, 'end'    ; exit
          free_lun, lun

; clean up files    
          if file_test(scenarios_outfile) then spawn, '/bin/rm '+scenarios_outfile, /sh

; generate the models
          t0 = systime(1)
          splog, 'Building the SFH grid for IMF '+imfstr
          spawn, 'scenarios_HR < '+scenarios_infile+' > '+scenarios_logfile, /sh     ; run scenarios.f
          spawn, 'echo "'+scenarios_outfile+'" | '+'spectra_HR > '+spectra_logfile, /sh ; run spectra.f
          splog, format='("Time = ",G0," minutes.")', (systime(1)-t0)/60.0

; now go back through all the models, convert the spectra to an
; isedfit-compatible format, and then overwrite
          for itau = 0, ntau-1 do begin
             peg = im_read_peg(taufile[itau])
             ised = peg2isedfit(peg)
             ised = struct_addtags({imf: imfstr, tau: float(tau[itau])},ised)
             im_mwrfits, ised, ssppath+taufile[itau], /clobber
          endfor
       endfor                    ; close IMF loop
       popd
       splog, format='("Time for all IMFs = '+$
         '",G0," minutes.")', (systime(1)-timf)/60.0
    endif
       
return
end
