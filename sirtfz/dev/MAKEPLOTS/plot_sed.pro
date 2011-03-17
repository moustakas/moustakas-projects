; jm00sep30uofa
; plot the sed's in the sed library

pro plot_sed, galaxy

    if not keyword_set(galaxy) then begin
       print, 'Please specify a galaxy name : arp220, mrk273, sbc1987, m82'
       print, 'Syntax : plot_sed, galaxy'
       return
    endif

    case galaxy of

       'arp220' : begin
          sedfile = 'arp220r.txt'
          z = 0.18D
          oplotdata = 1L	; overplot SED measurements
       endcase
       'mrk273' : begin
          sedfile = 'mrk273r.txt'
          z = 0.458D
          oplotdata = 0L
       endcase
       'sbc1987' : begin
          sedfile = 'sbc1987r.txt'
          oplotdata = 0L
       endcase
       'm82' : begin
          sedfile = 'm82r.txt'
          z = 0.000677D
          oplotdata = 0L
       endcase

       else: return        

    endcase
 
    set_cosmology   ; establish the cosmology
    common cosmo, omega_lambda, omega_0, h_100

; read in the sed (mu and W/Hz)

    path = (sirtf_datapath())[0]
    readcol, path+'sed/'+sedfile, wave, lum, format='D,D', /silent

    nwave = n_elements(wave)
    nlum = n_elements(lum)
    if nwave ne nlum then begin
       print, "There's a problem with the SED file."
       return
    endif
        
; rest spectrum (Devriendt et al 1998, Fig. 17)
; ----------------------------------------------------------------------    
    nu = 2.99793D14 / wave 
    lum_nu = nu*lum / 3.826D26  ; units of L_Sun_Bol
; ----------------------------------------------------------------------    
       
; reshift the spectrum
; ----------------------------------------------------------------------    
   wave = wave * (1.+z)
   dlum = dluminosity(z) * 3.086E16 ; luminosity distance (meters)
   flux = lum / (4.*!pi*dlum*dlum)
   flux = flux / 1.E-26 	; flux in Jy
; ----------------------------------------------------------------------    

    colortable1
    plotfaves, charsize=1.5
    s = texstrings()
    
    window, 0, xs=450, ys=450
    plot, wave, flux, /xlog, /ylog, xsty=3, ysty=3, color=1, $
      ytit=textoidl('f_{\nu}')+' (Jy)', xtit=s.lambda+' ('+s.mu+'m)', xr=[1,1E4]
;    plot, wave, alog10(lum_nu), /xlog, xsty=3, ysty=3, color=1, xr=[0.02,9000], yrange=[3,12.2], $
;      ytit='log '+s.nu+' '+textoidl('L_{\nu}')+' ('+textoidl('L_{sol}')+')', xtit=s.lambda+' ('+s.mu+'m)'

    if oplotdata then begin
       readcol, path+'data/arp220_data.dat', wdat, ldat, format='D,D', /silent, skip=7
       plotsym, 0, 0.8, /fill
       oplot, wdat*(1.+z), ldat, ps=8, color=5
    endif

stop
    plotfaves, /restore

return
end
