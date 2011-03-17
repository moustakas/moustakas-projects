pro sc1120_calibrate_wrapper, wfits=wfits
; flux-calibrate all the data in SPEC1D_SIGSPEC using the appropriate
; sensitivity functions in Q1-Q4; this routine should be run after
; SC1120_INFO 

    datapath = sc1120_path(/sigspec)
    outpath = sc1120_path(/spec1d)
    analysis_path = sc1120_path(/analysis)

; read in the ancillary data file from SC1120_INFO

    data = mrdfits(analysis_path+'sc1120_data.fits.gz',1,/silent)
    speclist = strtrim(data.specfile,2)
    quadrant = strtrim(data.quadrant,2)
    
    quad = ['Q1','Q2','Q3','Q4']
    sensnamelist = 'sens_'+quad+'.fits'
    quadpath = [sc1120_path(/Q1),sc1120_path(/Q2),sc1120_path(/Q3),sc1120_path(/Q4)]

;   for iquad = 3L, n_elements(quadpath)-1L do begin
    for iquad = 0L, n_elements(quadpath)-1L do begin

       sensname = sensnamelist[iquad]
       match = where(strmatch(quadrant,'*'+quad[iquad]+'*',/fold) eq 1B,nmatch)

       splog, 'Flux-calibrating '+quad[iquad]+': '+string(nmatch,format='(I0)')+' objects.'
       sc1120_calibrate, speclist[match], sensname, datapath=datapath, $
         senspath=quadpath[iquad], outpath=outpath, wfits=wfits

    endfor
    
stop       
       
return
end
    
