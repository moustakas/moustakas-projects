pro burstmodels, basemodels=basemodels
; jm05mar17uofa
; construct the models needed to explore the time-evolution of the
; L(U)/L(Ha) luminosity ratio

    diskage = 10.0 ; disk/galaxy age [Gyr] 
    
    datapath = '/home/ioannis/jobs/burstmodels/'
    bc03path = getenv('bc03_dir')+'/src/'         ; FORTRAN source path
    isedpath = getenv('bc03_dir')+'/models/Padova1994/salpeter/'

    isedfile = 'bc2003_hr_m62_salp_ssp.ised' ; input SSP

; first construct the necessary base models

    if keyword_set(basemodels) then begin
       
       cspinput = 'base_models.input'
       cspoutput = 'base_models.output'

       splog, 'Opening '+datapath+cspinput+'.'
       openw, lun, datapath+cspinput, /get_lun

; truncated burst; burst duration = 250 Myr
       
       printf, lun, isedfile
       printf, lun, 'N'         ; include dust attenuation?
       printf, lun, '2'         ; SFH law [truncated burst]
       printf, lun, '0.25'      ; burst duration [Gyr]
       printf, lun, 'burst_250Myr.ised' ; output file

; truncated burst; burst duration = 500 Myr
       
       printf, lun, isedfile
       printf, lun, 'N'         ; include dust attenuation?
       printf, lun, '2'         ; SFH law [truncated burst]
       printf, lun, '0.50'      ; burst duration [Gyr]
       printf, lun, 'burst_500Myr.ised' ; output file

; SSP (tau = 0)
       
       printf, lun, isedfile
       printf, lun, 'N'         ; include dust attenuation?
       printf, lun, '0'         ; SFH law [instantaneous burst - delta function]
       printf, lun, 'tau0.ised' ; output file
       
; continuous star formation: tau=1
       
       printf, lun, isedfile
       printf, lun, 'N'         ; include dust attenuation?
       printf, lun, '1'         ; SFH law [exponential]
       printf, lun, '1.0'       ; e-folding time
       printf, lun, 'N'         ; include gas recycling?
       printf, lun, '20.0'      ; SFR = 0 at t = 20 Gyr
       printf, lun, 'tau1.ised' ; output file
       
; continuous star formation: tau=2
       
       printf, lun, isedfile
       printf, lun, 'N'         ; include dust attenuation?
       printf, lun, '1'         ; SFH law [exponential]
       printf, lun, '2.0'       ; e-folding time
       printf, lun, 'N'         ; include gas recycling?
       printf, lun, '20.0'      ; SFR = 0 at t = 20 Gyr
       printf, lun, 'tau2.ised' ; output file
       
; continuous star formation: tau=3
       
       printf, lun, isedfile
       printf, lun, 'N'         ; include dust attenuation?
       printf, lun, '1'         ; SFH law [exponential]
       printf, lun, '3.0'       ; e-folding time
       printf, lun, 'N'         ; include gas recycling?
       printf, lun, '20.0'      ; SFR = 0 at t = 20 Gyr
       printf, lun, 'tau3.ised' ; output file
       
; continuous star formation: tau=5
       
       printf, lun, isedfile
       printf, lun, 'N'         ; include dust attenuation?
       printf, lun, '1'         ; SFH law [exponential]
       printf, lun, '5.0'       ; e-folding time
       printf, lun, 'N'         ; include gas recycling?
       printf, lun, '20.0'      ; SFR = 0 at t = 20 Gyr
       printf, lun, 'tau5.ised' ; output file
       
; continuous star formation: tau=8
       
       printf, lun, isedfile
       printf, lun, 'N'         ; include dust attenuation?
       printf, lun, '1'         ; SFH law [exponential]
       printf, lun, '8.0'       ; e-folding time
       printf, lun, 'N'         ; include gas recycling?
       printf, lun, '20.0'      ; SFR = 0 at t = 20 Gyr
       printf, lun, 'tau8.ised' ; output file

; constant star formation (tau=infinity): 1 M_sun/yr 
       
       printf, lun, isedfile
       printf, lun, 'N'         ; include dust attenuation?
       printf, lun, '3'         ; SFH law [constant]
       printf, lun, '1.0'       ; SFR [M_sun/yr]
       printf, lun, '20.0'      ; SFR = 0 at t = 20 Gyr
       printf, lun, 'constant.ised' ; output file

       splog, 'Closing '+datapath+cspinput+'.'
       free_lun, lun

; run csp_galaxev and move the resultant models into the output
; directory 

       splog, 'Generating the base models.'
       t0 = systime(1)
       spawn, [bc03path+'csp_galaxev < '+datapath+cspinput+' > '+datapath+cspoutput], /sh
       splog, format='("Total time = ",G0," minutes.")', (systime(1)-t0)/60.0

    endif

; next construct the composite models.  superpose the 250 Myr burst
; onto all the continuous (tau-based) models

    burstinput = 'add_bursts.input'
    burstoutput = 'add_bursts.output'

    isedcontfile = ['tau0','tau1','tau2','tau3','tau5','tau8']+'.ised'
    ncont = n_elements(isedcontfile)

    tburst = [1.0,2.0,3.0,4.0,6.0]    ; time *since* burst [Gyr]
;   tburst = [0.1,1.0,2.0,5.0]    ; time *since* burst [Gyr]
    ageatburst = diskage - tburst ; age of the galaxy at the time of the burst [Gyr]
    nburst = n_elements(tburst)

    tburst_string = strtrim(string(tburst,format='(F4.1)'),2)
    ageatburst_string = strtrim(string(ageatburst,format='(F4.1)'),2)
    
; loop on each continuous model and on each value of TBURST

    splog, 'Generating the composite models.'
    t0 = systime(1)
    for icont = 0L, ncont-1L do begin
    
       for iburst = 0L, nburst-1L do begin

          print, format='("Base model ",I0,"/",I0,", Burst ",I0,"/",I0,".",A1,$)', $
            icont+1L, ncont, iburst+1L, nburst, string(13b)

          isedoutfile = repstr(isedcontfile[icont],'.ised','_burst_500Myr_'+$
;         isedoutfile = repstr(isedcontfile[icont],'.ised','_burst_250Myr_'+$
            tburst_string[iburst])+'GyrAgo.ised'

;         splog, 'Opening '+datapath+burstinput+'.'
          openw, lun, datapath+burstinput, /get_lun

          printf, lun, isedcontfile[icont]              ; first model
          printf, lun, 'burst_500Myr.ised'              ; second model
;         printf, lun, 'burst_250Myr.ised'              ; second model
          printf, lun, '0.0,1.0'                        ; burst 1: beginning time (Gyr), burst amplitude
          printf, lun, ageatburst_string[iburst]+',0.1' ; burst 2 (10% mass fraction)
          printf, lun, isedoutfile                      ; output file

;         splog, 'Closing '+datapath+burstinput+'.'
          free_lun, lun
          
; launch the script!
          
          spawn, [bc03path+'add_bursts < '+datapath+burstinput+' > '+datapath+burstoutput], /sh

       endfor

    endfor
    splog, format='("Total time = ",G0," minutes.")', (systime(1)-t0)/60.0

stop       

return
end
    
