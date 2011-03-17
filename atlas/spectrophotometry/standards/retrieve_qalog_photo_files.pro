pro retrieve_qalog_photo_files
; jm05jun28uofa

    datapath = '/d1/ioannis/spec2datlas/'
    outpath = atlas_path(/analysis)+'spectrophotometry/standards/qalog_photo/'

    subpath = [$
      '94nov/','95mar/',$
      '95oct/','96apr/',$
      '97apr/',$
      '98mar/','98apr/',$
      '98jun/','98oct/',$
      '99apr/','99may/',$
      '99nov/','00apr/',$
      '01nov/','02feb/',$
      '02apr/','02may/',$
      '03may/']
    npath = n_elements(subpath)

    for i = 0L, npath-1L do begin

       loglist = file_search(datapath+subpath[i]+'qalog_sens_4.5*.log',count=logcount)
       if (logcount ne 0L) then spawn, ['/bin/cp '+strjoin(loglist,' ')+' '+outpath]

       driftlist = file_search(datapath+subpath[i]+'qalog_sens_drift*.log',count=driftlogcount)
       if (driftlogcount ne 0L) then spawn, ['/bin/cp '+strjoin(driftlist,' ')+' '+outpath]
       
    endfor     

return
end
