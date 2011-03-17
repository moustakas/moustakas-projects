function ediscs_readcat, cluster
; jm04feb13uofa
; read the ediscs photometric catalogs

    if (n_elements(cluster) eq 0L) then begin
       print, 'Syntax - '
       return, -1L
    endif
    
    path = ediscs_path(/catalogs)
    pushd, path
    cllist = file_search('cl*',count=clcount)
    fmtlist = file_search('*.fmt',count=fmtcount)
    popd
    
    match = where(strmatch(cllist,'*'+cluster+'*',/fold_case) eq 1B,nmatch)

    if (nmatch eq 0L) then begin
       splog, 'No matching cluster: '+cluster+'.'
       return, -1L
    endif
    
    if (nmatch gt 1L) then begin
       splog, 'Multiple matching clusters: '+cluster+'.'
       return, -1L
    endif

    catname = cllist[match[0]]
    fmtname = repstr(strmid(catname,11),'phot.v21.cat','')
    
    match = where(strmatch(fmtlist,fmtname+'*',/fold_case) eq 1B,nmatch)
    
    if (nmatch eq 0L) then begin
       splog, 'No matching format file: '+fmtname+'.'
       return, -1L
    endif
    
    if (nmatch gt 1L) then begin
       splog, 'Multiple matching format files: '+fmtname+'.'
       return, -1L
    endif

    fmtfile = fmtlist[match[0]];+'_phot_cat.fmt'

; read the format file and build the data structure

    readcol, path+fmtfile, tags, format='A', /silent
    ntags = n_elements(tags)

    cat = create_struct(tags[0],'',tags[1],'')
    for itag = 2L, ntags-1L do cat = create_struct(cat,tags[itag],0.0)
    
; read the catalog

    catdata = djs_readlines(path+catname)
    nobject = n_elements(catdata)
    cat = replicate(cat,nobject)

    for iobj = 0L, nobject-1L do begin

       line = strsplit(catdata[iobj],' ',/extract)
       for itag = 0L, ntags-1L do cat[iobj].(itag) = line[itag]

    endfor
    
return, cat
end    

