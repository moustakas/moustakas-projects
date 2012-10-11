pro phys130_lab_f12_grades
; jm12oct05siena - parse the grades for this class 

    date = '12oct06' ; update this
    semester = 'Fall 2012'
    class = 'Physics 130 - General Physics I - Lab'
    
; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = '~/Dropbox/teaching/12fall/phys130/grades/phys130_lab_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)
    nstudent = n_elements(data)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = 'Lab'
    weight = 1.0

    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, /lab

stop    

return
end
