#! /usr/bin/env python

# input parameters
header = raw_input('\nhdr: header (1=Yes, 0=No) ?\n')
pos_diff = raw_input('hdr: distance (clustering condition) ?\n')
cl_size = raw_input('hdr: min cluster size ?\n')
collapse_size = raw_input('hdr: min distance to collapse clusters ?\n\n')
infile = raw_input('hdr: data filename ?\n')
cl_out = raw_input('hdr: cl output filename ?\n')
mcl_out = raw_input('hdr: mcl output filename ?\n')
no_cl_out = raw_input('hdr: no_cl output filename ?\n')
no_mcl_out = raw_input('hdr: no_mcl output file ?\n')


# check for invalid input parameters
e1='\nhdr terminated: input parameter error.\n'
e_num='is not a positive integer.\n'
e_file='is not a valid file name.\n'
#header
try:
    num = int(header)
except ValueError, TypeError:
    print e1,'\"',header.strip(),'\"',e_num
    raise SystemExit
if (int(header) !=0 and int(header) !=1):
	print e1,'header flag can only be 0 or 1\n'
	raise SystemExit
#poss_diff
try:
    num = int(pos_diff)
except ValueError, TypeError:
    print e1,'\"',pos_diff.strip(),'\"',e_num
    raise SystemExit
if (int(pos_diff) < 0):
	print e1,'\"',pos_diff.strip(),'\"',e_num
	raise SystemExit
#cl_size
try:
    num = int(cl_size)
except ValueError, TypeError:
    print e1,'\"',cl_size.strip(),'\"',e_num
    raise SystemExit
if (int(cl_size) < 0):
	print e1,'\"',cl_size.strip(),'\"',e_num
	raise SystemExit
#collapse_size
try:
    num = int(collapse_size)
except ValueError, TypeError:
    print e1,'\"',collapse_size.strip(),'\"',e_num
    raise SystemExit
if (int(collapse_size) < 0):
	print e1,'\"',collapse_size.strip(),'\"',e_num
	raise SystemExit
#infile
if(len(infile.split())!=1):
	print e1,'\"',infile.strip(),'\"',e_file
	raise SystemExit
#cl_out
if(len(cl_out.split())!=1):
	print e1,'\"',cl_out.strip(),'\"',e_file
	raise SystemExit
#mcl_out
if(len(mcl_out.split())!=1):
	print e1,'\"',mcl_out.strip(),'\"',e_file
	raise SystemExit
#no_cl_out
if(len(no_cl_out.split())!=1):
	print e1,'\"',no_cl_out.strip(),'\"',e_file
	raise SystemExit
#no_mcl_out
if(len(no_mcl_out.split())!=1):
	print e1,'\"',no_mcl_out.strip(),'\"',e_file
	raise SystemExit


# read data from file
chrom=[]
pos=[]
all=[]
headerline=[]
with open(infile.strip(),'r') as fin:
    read_line = lambda: fin.readline()
    for line in iter(read_line,''):
    
        if int(header)==1:
            headerline.append(line.strip())
            #print '\nhdr: header is',headerline
            header=2
            continue
        
        if(len(line.split())>1):
            chrom.append(line.split()[0])
            pos.append(line.split()[1])
            all.append(line.strip())
        else:
            print 'hdr warning: line with',len(line.split()),'column(s) ignored'
print '\nhdr:',len(pos),'lines','\n','hdr:',len(all[0].split()),'columns (in first line)'
for i in range(len(pos)-1):
    if len(all[i].split()) != len(all[0].split()):
        print 'hdr warning:',len(all[i].split()),'columns in line',i


# identify clusters
c0=[0]
c1=[]
for i in range(1,len(pos)):
    if chrom[i]!=chrom[i-1] or int(pos[i])-int(pos[i-1])>int(pos_diff):
        c1.append(i-1)
        c0.append(i)
c1.append(len(all)-1)
print 'hdr:',len(c0),'cluster(s)'


# identify clusters of cl_size or larger
cl0=[]
cl1=[]
for i in range(1,len(c0)):
    if int(c1[i])-int(c0[i])+1>=int(cl_size):
        cl1.append(c1[i])
        cl0.append(c0[i])
print 'hdr:',len(cl0),'cluster(s) size',int(cl_size),'or larger'


#output clusters of cl_size or larger
with open(cl_out.strip(),'w') as fclout:
    for i in range(len(cl0)):
        for j in range(int(cl0[i]),int(cl1[i])+1):
            if j==int(cl0[i]):
                fclout.write( '#' + chrom[int(cl0[i])] + '\t' + str(int(cl1[i]) - int(cl0[i]) + 1) + '\t' + pos[int(cl0[i])] + '\t' + pos[int(cl1[i])] +  '\n' )
            fclout.write(all[j]+'\n')    
        fclout.write('\n')


# output no-clusters (variants not in any cluster of cl_size or larger)
with open(no_cl_out.strip(),'w') as fnoclout:
    if cl0: #check cl0 is not empty
        for j in range(0,int(cl0[0])):
            fnoclout.write(all[j]+'\n')
        for i in range(len(cl0)-1):
            for j in range(int(cl1[i])+1,int(cl0[i+1])):
                fnoclout.write(all[j]+'\n')
        for j in range(int(cl1[-1])+1,len(pos)):
            fnoclout.write(all[j]+'\n')


# identify megaclusters:
# merged clusters (of cl_size or larger) that are apart collapse_size or less, including the variants between them
mcl0=[cl0[0]]
mcl1=[]
for i in range(1,len(cl0)):
    if chrom[int(cl0[i])]!=chrom[int(cl1[i-1])] or int(pos[int(cl0[i])])-int(pos[int(cl1[i-1])])>int(collapse_size):
        mcl1.append(cl1[i-1])
        mcl0.append(cl0[i])
mcl1.append(cl1[-1])
print 'hdr:',len(mcl0),'megacluster(s)\n'


# output megaclusters
#with open(mcl_out.strip(),'w') as fmclout:
with open(mcl_out.strip(),'w') as fmclout:
    for i in range(len(mcl0)):
        for j in range(int(mcl0[i]),int(mcl1[i])+1):
            if j==int(mcl0[i]):
                fmclout.write( '\n#' + chrom[int(mcl0[i])] + '\t' + str(int(mcl1[i]) - int(mcl0[i]) + 1) + '\t' + pos[int(mcl0[i])] + '\t' + pos[int(mcl1[i])] +  '\n' )
            fmclout.write(all[j]+'\n')


# output no-megaclusters (variants not in any megacluster)
with open(no_mcl_out.strip(),'w') as fnomclout:
    if mcl0: #check mcl0 is not empty
        for j in range(0,int(mcl0[0])):
            fnomclout.write(all[j]+'\n')
        for i in range(len(mcl0)-1):
            for j in range(int(mcl1[i])+1,int(mcl0[i+1])):
                fnomclout.write(all[j]+'\n')
        for j in range(int(mcl1[-1])+1,len(pos)):
                fnomclout.write(all[j]+'\n')
