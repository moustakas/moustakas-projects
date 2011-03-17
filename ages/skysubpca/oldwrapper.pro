pro wrapper
; jm06jan22uofa - wrapper for AGES_SKYSUBPCA

    datapath = './'

    datafiles = file_search(datapath+'spectra_???.fits.gz',count=npass)
    infofiles = file_search(datapath+'target_info_???.fits.gz',count=ninfo)

    for ipass = 0L, npass-1L do begin

       splog, 'Reading passfile '+file_basename(datafiles[ipass])+'.'
       data = mrdfits(datafiles[ipass],1,/silent)
       info = mrdfits(infofiles[ipass],1,/silent)
       nobject = n_elements(info)

       qaplotname = datapath+'qaplot_'+string(info[0].pass,format='(I0)')+'.ps'
       newdata = ages_skysubpca(data,info,qaplotname=qaplotname)

    endfor
    
stop       
    
return
end
    
