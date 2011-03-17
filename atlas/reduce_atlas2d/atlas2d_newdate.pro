pro doupdate, flist, date

    fcount = n_elements(flist)
    for k = 0L, fcount-1L do begin
       h = headfits(flist[k])
       sxaddpar, h, 'DATE-OBS', date
       djs_modfits, flist[k], 0, h
    endfor
    
return
end

pro atlas2d_newdate, update=update
; jm02jan30uofa
; update - specify the directory (by number) to update
    
    datapath = '/home/ioannis/kennicutt/data/'
    path = datapath+['98mar/','98apr/','98jun/','98oct/']
    npath = n_elements(path)

    if (n_elements(update) ne 0L) then update = (0L > long(update)) < npath else $
      update = lindgen(npath)
    nupdate = n_elements(update)

    for k = 0L, nupdate-1L do begin

       i = update[k]

       print, 'Updating '+path[i]+'.'
       pushd, path[i]
       
       case i of

          0L: begin
             newdates = ['1998-03-22']
             flist = findfile('a.01*.fits',count=fcount)
             doupdate, flist, newdates[0]
          end

          1L: begin
             newdates = ['1998-04-28','1998-04-29','1998-04-30','1998-05-01','1998-05-02','1998-05-03']
             ndate = n_elements(newdates)
             offset = 2
             for j = 0L, ndate-1L do begin
                num = offset + j
                if num lt 10L then id = '0'+string(num,format='(I1)') else id = string(num,format='(I2)')
                flist = findfile('a.'+id+'*.fits',count=fcount)
                doupdate, flist, newdates[j]
             endfor
          end

          2L: begin
             newdates = ['1998-06-25','1998-06-26','1998-06-27']
             ndate = n_elements(newdates)
             offset = 8
             for j = 0L, ndate-1L do begin
                num = offset + j
                if num lt 10L then id = '0'+string(num,format='(I1)') else id = string(num,format='(I2)')
                flist = findfile('a.'+id+'*.fits',count=fcount)
                doupdate, flist, newdates[j]
             endfor
          end

          3L: begin
             newdates = ['1998-10-16','1998-10-17','1998-10-18','1998-10-19','1998-10-20']
             ndate = n_elements(newdates)
             offset = 11
             for j = 0L, ndate-1L do begin
                num = offset + j
                if num lt 10L then id = '0'+string(num,format='(I1)') else id = string(num,format='(I2)')
                flist = findfile('a.'+id+'*.fits',count=fcount)
                doupdate, flist, newdates[j]
             endfor
          end
       endcase
       
       popd

    endfor

return
end
