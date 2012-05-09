function maskpops_path, isedfit=isedfit, montegrids=montegrids
; jm12may07ucsd

    path = getenv('CLASH_DATA')+'/projects/maskpops/'
    if keyword_set(isedfit) then path = path+'isedfit/'
    if keyword_set(montegrids) then path = path+'montegrids/'
    
return, path
end
