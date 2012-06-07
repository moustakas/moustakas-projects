function irclusters_path, isedfit=isedfit, montegrids=montegrids
; jm12jun08ucsd

    path = getenv('IRCLUSTERS_DATA')+'/'
    if keyword_set(isedfit) then path = path+'isedfit/'
    if keyword_set(montegrids) then path = path+'montegrids/'
    
return, path
end
