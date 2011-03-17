pro wrapper
; jm06jan22uofa - wrapper for AGES_SKYSUBPCA

    datapath = '~/data/AGES/SDSS_before_skysubpca/'
    outpath = '~/data/AGES/SDSS_after_skysubpca/'

    datafiles = file_search(datapath+'spectra_???.fits.gz',count=npass)
    infofiles = file_search(datapath+'target_info_???.fits.gz',count=ninfo)

    splog, filename=outpath+'splog.txt'

    for ipass = 0L, npass-1L do begin

       splog, 'Reading passfile '+file_basename(datafiles[ipass])+'.'
       data = mrdfits(datafiles[ipass],1,/silent)
       info = mrdfits(infofiles[ipass],1,/silent)
       nobject = n_elements(info)

       qaplotname = outpath+'qaplot_'+string(info[0].pass,format='(I0)')+'.ps'
       newdata = ages_skysubpca_v2(data,info,pcainfo,qaplotname=qaplotname)

;       save,newdata,pcainfo,file=outpath+'spectra_'+string(info[0].pass,format='(I0)')+'.sav'

    endfor
    
    splog,/close

    
return
end
    
