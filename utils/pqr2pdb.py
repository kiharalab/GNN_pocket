def pqr2pdb(ifname,ofname):
    with  open(ifname,"r") as ifhandle:
        with open(ofname,"w") as ofhandle:
            for line in ifhandle:
                if line.startswith("ATOM"):
                    writePQRAtom(line,ofhandle)
                else:
                    ofhandle.write(line)

def writePQRAtom(line,ofh):
    spl = line.split()
    
    aname=spl[2]
    if not aname[0].isdigit():
        #print "change aname\"%s\""%aname
        aname = " "+aname
        #print "to \"%s\""%aname
    
    if len(spl)==11:
        outline="ATOM  %5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s%s\n" % (int(spl[1]),aname," ",spl[3],spl[4],int(spl[5])," ",float(spl[6]),float(spl[7]),float(spl[8]),float(spl[9]),float(spl[10]),"    ","  ","  ","")
    elif len(spl)==10:
        outline="ATOM  %5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s%s\n" % (int(spl[1]),aname," ",spl[3]," ",int(spl[4])," ",float(spl[5]),float(spl[6]),float(spl[7]),float(spl[8]),float(spl[9]),"    ","  ","  ","")
    else:
        print(spl)
        raise ValueError("This pqr was not written by editconf...")

        
    ofh.write(outline)
    