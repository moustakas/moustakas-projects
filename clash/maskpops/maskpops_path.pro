function maskpops_path, isedfit=isedfit, montegrids=montegrids
; jm12may07ucsd

    path = getenv('IM_RESEARCH_DIR')+'/projects/clash/maskpops/'
    if keyword_set(isedfit) then path = path
    if keyword_set(montegrids) then path = path+'montegrids/'
    
return, path
end
