#!/usr/bin/python

import os,sys
if len(sys.argv)<3:
    print "Usage: topdf <npages> <file.ps>"
    exit()
else:
    pages = int(sys.argv[1])
    D = dict(zip(['fname','ext'],sys.argv[2].split('.')))

for i in range(pages):
    D['n']=i+1
    os.system("psselect -p%(n)d %(fname)s.%(ext)s %(fname)s%(n)d.%(ext)s"%D)
    os.system("epstopdf --autorotate=All %(fname)s%(n)d.%(ext)s"%D)
inputs = (' '.join(['%s%d.pdf'%(D['fname'],i+1) for i in range(pages)]))
os.system('pdftk %s cat output %s.pdf'%(inputs,D['fname']))
os.system('rm %s'%inputs)
os.system('rm %s'%inputs.replace('.pdf','.ps'))
