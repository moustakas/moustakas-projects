pro pisces_align

    readfast, 'source_list.txt', fixim, skip=1
    readfast, 'temp_list.txt', rotim, skip=1

; rotations are about the center of the image
    
    x_0 = fixim[0,*]-511.5
    y_0 = fixim[1,*]-511.5
    x_i = rotim[0,*]-511.5
    y_i = rotim[1,*]-511.5
    
    rot = fltarr(2,2) ; rotation matrix
    degree = 1

    polywarp, x_i, y_i, x_0, y_0, degree, k_x, k_y

; translational offset   

    deltax = k_x[0,0]
    deltay = k_y[0,0]

; rotation correction (radians)

    rot[0,0] = k_x[0,1]
    rot[0,1] = k_y[0,1]
    rot[1,0] = k_x[1,0]
    rot[1,1] = k_y[1,0]

; plate scale

    print, 'Inital plate scale guess off by factor of ', determ(rot)

; cross terms

    crossx = k_x[1,1]
    crossy = k_y[1,1]
   
    print, 'Cross terms should be small: ', crossx, crossy

stop    

return
end
