; jm00mar19ucb

; routine to read in the fits images and neighbor-subtracted images
; and write them back out (converts the data type in the header)

pro trgb_convertall
	spawn, ['clear']
	spawn, ['pwd'], datapath
        datapath = datapath[0]

        masterdir = rdtxt(datapath+'/master.dir')
        masterdir = masterdir[0]

; read in the image names

        imnames = rdtxt(datapath+'/image.names')        

; find the neighbor-subtracted images

        spawn, 'find '+masterdir+' -name "*nsub.fits" -print', nsublist
        nold = where(strpos(nsublist,'OLD') eq -1,count)
        if count ne 0 then nsublist = nsublist[nold]

        for k = 0L, n_elements(imnames)-1L do begin
            imfile = masterdir+'/'+imnames[k]+'.fits'
            nsubfile = nsublist[k]

            im = readfits(imfile,him)
            print, 'Writing '+imfile+'.'
            writefits, imfile, im, him

            imnsub = readfits(nsubfile,hnsub)
            print, 'Writing '+nsubfile+'.'
            writefits, nsubfile, imnsub, hnsub
        endfor            

return
end

