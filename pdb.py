from datetime import date
from functools import reduce

def header(data, io, rev_date="", united=False):
        # Write header
        # Refer to http://www.wwpdb.org/documentation/format32/sect2.html 
        print('HEADER'+4*' '+'%-40s' %('UNCLASSIFIED')+\
                date.strftime(date.today(),'%d-%b-%g'), file=io)
        if united:
            strtype = 'UNITED'
        else:
            strtype = 'ALL'
        print('TITLE '+4*' '+\
                '%-70s' %(strtype + ' ATOM STRUCTURE FOR MOLECULE ' + data[1]['group']), file=io)
        print('AUTHOR'+4*' '+'GROMOS AUTOMATIC TOPOLOGY BUILDER REVISION ' + rev_date, file=io)
        print('AUTHOR'+3*' '+'2'+'  http://compbio.biosci.uq.edu.au/atb', file=io)

def atoms(data, io, atoms, united=False, optimized=True, use_rnme=True):
        # Write structure
        for i in atoms:
            if united:
                if 'uindex' not in i: continue
                index = i['uindex']
            else:
                index = i['index']
            if optimized and data.completed('has_ocoord'):
                coord = [j*10. for j in i['ocoord']]
            else:
                coord = [j*10. for j in i['coord']]
            if use_rnme:
                rnme = data.var['rnme']
            else:
                rnme = i['group']
            print(i['pdb'][:6] + '%5d' %index + i['pdb'][11:17] + \
                    '%-4s' %rnme + "    0    " + \
                    '%8.3f' %coord[0] + '%8.3f' %coord[1] + '%8.3f' %coord[2] +\
                    i['pdb'][54:], file=io)

def connectivity(data, io, atoms, united=False):
        # Write connectivity
        for i in atoms:
            if united:
                if 'uindex' not in i: continue
                nbr = sorted([data[n]['uindex'] \
                        for n in i['conn'] if 'uindex' in data[n]])
                if len(nbr) == 0: continue
                neighborstr = reduce(lambda x,y:x+y, ['%5s' %n for n in nbr], '')
                print("CONECT" + '%5d' %i['uindex'] + neighborstr, file=io)
            else:
                #sort connectivity
                nbr = sorted([data[n]['index'] for n in i['conn']])
                neighborstr = reduce(lambda x,y:x+y, ['%5s' %n for n in nbr], '')
                print("CONECT" + '%5d' %i['index'] + neighborstr, file=io)

def footer(data, io):
        print('END', file=io)
