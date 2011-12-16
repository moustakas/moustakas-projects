function read_sfh_comments

    qapath = ediscs_path(/projects)+'sfh/qaplots/'
    readcol, qapath+'all_qa_comments.txt', id, memb, specfile, $
      snr, qual, comments, delim=';', skip=1, /silent, $
      format='A,I,A,F,A,A'
    ngal = n_elements(id)
    data = {id: '', memb: 0, specfile: '', snr: 0.0, qual: '', comments: ''}
    data = replicate(data,ngal)

    data.id = strtrim(id,2)
    data.memb = memb
    data.specfile = strtrim(specfile,2)
    data.snr = snr
    data.qual = strtrim(qual,2)
    data.comments = strtrim(comments,2)

return, data
end
