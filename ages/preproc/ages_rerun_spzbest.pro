pro ages_rerun_spzbest
; jm06feb05uofa - re-run Schlegel's redshift classification code on
;                 all those objects whose SPZBEST redshifts differ
;                 from the official AGES redshifts (as determined by
;                 AGES_FLAG_BAD_REDSHIFTS); then, in
;                 AGES_RECLASSIFY_BAD_REDSHIFTS we determine if *any*
;                 of the redshift solutions found by HS_REDUCE1D are
;                 close to the official AGES redshifts

; make sure you are in the proper directory: ~/ages/SDSS/spzbest_reclassify    
    
    analysis_path = ages_path(/analysis)
    info = mrdfits(analysis_path+'ages_flag_bad_redshifts.fits.gz',1,/silent)

    allplatefiles = info[uniq(info.platefile,sort(info.platefile))].platefile
    nplate = n_elements(allplatefiles)

;   for i = 0L, 4L do begin
    for i = 0L, nplate-1L do begin

       match = where(allplatefiles[i] eq info.platefile)
       info1 = info[match]
       
       if (info1[0].pass lt 401) then rerun = '0300' else rerun = '0051'
       datapath = ages_path(/rawdata)+rerun+'/'
       
       platefile = info1[0].platefile
       pass = strmid(platefile,7,1)
       field = strmid(platefile,8,2)
       fiberid = info1.aper

       hs_reduce1d, datapath+platefile, pass=pass, field=field, fiberid=fiberid, rerun=rerun

    endfor

stop    
    
return
end
    
