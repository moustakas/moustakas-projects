function read_sfh_comments, moustakas=moustakas

    data = {id: '', memb: 0, specfile: '', snr: 0.0, fixed: '', qual: '', comments: ''}

    qapath = ediscs_path(/projects)+'sfh/'
    if keyword_set(moustakas) then begin
       readcol, qapath+'moustakas_qa_comments.txt', id, memb, specfile, $
         snr, fixed, qual, comments, delim=';', skip=1, /silent, $
         format='A,I,A,F,A,A,A'
    endif else begin
       readcol, qapath+'all_qa_comments.txt', id, memb, specfile, $
         snr, qual, comments, delim=';', skip=1, /silent, $
         format='A,I,A,F,A,A'
    endelse
    ngal = n_elements(id)
    data = replicate(data,ngal)

    data.id = strtrim(id,2)
    data.memb = memb
    data.specfile = strtrim(specfile,2)
    data.snr = snr
    data.qual = strtrim(qual,2)
    data.comments = strtrim(comments,2)

    if keyword_set(moustakas) then data.fixed = strtrim(fixed,2)
    
return, data
end
