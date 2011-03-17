pro ediscs_merge_photometry, out
; jm09jun24nyu - parse the EDiscS photometry

    path = ediscs_path(/catalogs)+'photometry/'
    
    all = file_search(path+'*.cat',count=nall)
    for ii = 0, nall-1 do begin
       all1 = file_basename(all[ii])
       splog, 'Parsing '+all1
       fmt = strmid(all1,11,strpos(all1,'.v23.cat')-11)
       fmt = repstr(fmt,'phot','_phot')
       fmtfile = path+fmt+'_cat.v23.fmt'
       readcol, fmtfile, tag, format='A', /silent
       ntag = n_elements(tag)
; read the data
       data = djs_readlines(all[ii])
       ngal = n_elements(data)
;      for jj = 0, 10 do begin
       for jj = 0, ngal-1L do begin
          data1 = strsplit(data[jj],' ',/extract)
          for kk = 0, ntag-1 do begin
             if strmatch(tag[kk],'*ediscsid*',/fold) then $
               value = strtrim(data1[kk],2) else $
                 value = float(data1[kk])
             if (strmatch(tag[kk],'*ra_hours*',/fold) or $
               strmatch(tag[kk],'*dec_deg*',/fold)) then $
                 value = double(data1[kk])
             if strmatch(tag[kk],'*flag*',/fold) then $
                 value = long(data1[kk])
;            case datatype(data1[kk]) of
;               'STR': value = strtrim(data1[kk],2)
;               'INT': value = fix(data1[kk],2)
;               else: value = float(data1[kk],2)
;            endcase
             if (kk eq 0) then $
               out1 = create_struct(tag[kk],value) else $
               out1 = create_struct(out1,tag[kk],value)
          endfor
          if (jj eq 0) then out = out1 else out = [out,out1]
       endfor       
; write out
       outfile = path+strmid(all1,0,11)+'.phot.v23.fits'
       im_mwrfits, out, outfile
    endfor

return
end
    
