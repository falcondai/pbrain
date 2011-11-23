#!/usr/bin/python
import os
import subprocess
import procrustes2reg as p2r
def go():
 print "procrustes moves the MRI into the CT space - written by Chuck Pelizzari"
 print "current directory is %s"%os.getcwd()
 in_file =raw_input('input fid file, including path if not in the current directory: ')
 out_file= raw_input('choose a filename for procrustes to write to: ')
 file(out_file,'wb').write('')
 f=open(out_file,'ab')
 subprocess.call(["/home/vtowle/towlelab/procrustes", in_file],stdout=f,stderr=f)
 f.close()

 out_file2 = raw_input('choose a filename for the reg file to be written to: ')

 p2r.run(out_file,out_file2)
 print "successfully wrote reg file to %s and intermediate file to %s"%(out_file2,out_file1)

go()
