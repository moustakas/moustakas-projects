pro ages_skysubpca_wrapper, fluxed=fluxed, cleanfits=cleanfits, clobber=clobber
; jm06jan20uofa - wrapper for AGES_SKYSUBPCA

    rootpath = ages_path(/spec1d)
    if keyword_set(fluxed) then begin
       datapath = rootpath+'fluxed/before_skysubpca/'
       outpath = rootpath+'fluxed/after_skysubpca/'
       scale = 1D-17
       suffix = 'fluxed'
    endif else begin
       datapath = rootpath+'unfluxed/before_skysubpca/'
       outpath = rootpath+'unfluxed/after_skysubpca/'
       scale = 1.0D
       suffix = 'unfluxed'
    endelse

    datafiles = file_search(datapath+'spectra_???.fits.gz',count=npass)
    infofiles = file_search(datapath+'target_info_???.fits.gz',count=ninfo)

    if keyword_set(cleanfits) then begin
       if (file_test(outpath,/directory) ne 0L) then begin
          splog, 'Remove all FITS files from '+outpath+' [Y/N]?'
          cc = get_kbrd(1)
          if (strupcase(cc) eq 'Y') then begin
             flist = file_search(outpath+'*.fits*')
             rmfile, flist
;            spawn, ['/bin/rm '+outpath+'*.fits*'], /sh
          endif
       endif
    endif

    for ipass=0L, npass-1L do if (n_elements(allinfo) eq 0L) then $
      allinfo = mrdfits(infofiles[ipass],1,/silent) else $
      allinfo = struct_append(allinfo,mrdfits(infofiles[ipass],1,/silent))
    help, allinfo

    splog, filename=outpath+'ages_skysubpca_wrapper.log'
    for ipass = 0, npass-1 do begin
       
       splog, 'Reading passfile '+file_basename(datafiles[ipass])+'.'
       data = mrdfits(datafiles[ipass],1,/silent)
       splog, 'Reading infofile '+file_basename(infofiles[ipass])+'.'
       info = mrdfits(infofiles[ipass],1,/silent)
       nobject = n_elements(info)

       qaplotname = datapath+'qaplot_'+string(info[0].pass,format='(I0)')+'.ps'
       newdata = ages_skysubpca(data,info,qaplotname=qaplotname,scale=scale)
       newdata = data

       outinfo = {ages_id: lonarr(nobject), galaxy: strarr(nobject), $
         ra: dblarr(nobject), dec: dblarr(nobject), pass: 0, $
         aper: intarr(nobject), class: strarr(nobject), z: fltarr(nobject), $
         spzbest: fltarr(nobject), $ ; fluxing_problem: intarr(nobject), $
         zmerge_multiple: intarr(nobject), minwave: fltarr(nobject), $
         maxwave: fltarr(nobject)}
       outinfo.ages_id = info.ages_id
       outinfo.galaxy = info.galaxy
       outinfo.ra = info.ra
       outinfo.dec = info.dec
       outinfo.pass = info[0].pass
       outinfo.aper = info.aper
       outinfo.class = info.class
       outinfo.z = info.z
       outinfo.spzbest = info.spzbest
       outinfo.minwave = info.minwave
       outinfo.maxwave = info.maxwave
;      outinfo.fluxing_problem = info.fluxing_problem
       outinfo.zmerge_multiple = info.zmerge_multiple

       final_data = struct_addtags(outinfo,newdata)
;      final_data = struct_addtags(outinfo,struct_trimtags(newdata,except=['SKYFLUX','SKYFERR']))
       outfile = 'ages_'+string(info[0].pass,format='(I0)')+'.fits'
       im_mwrfits, final_data, outpath+outfile, clobber=clobber
    endfor 

    allinfofile = 'ages_after_skysubpca_info_'+suffix+'.fits'
    im_mwrfits, allinfo, outpath+allinfofile, clobber=clobber
    splog, /close

return
end
    
