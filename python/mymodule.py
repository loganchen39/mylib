#!/usr/bin/python
# #!/usr/common/usg/dynamic_libs/usr/bin/python

# Filename: mymodule.py

#--------------------------------------------------------------------
# functions

import sys
import os
import time
import datetime


def is_leap(year):
    if (year%400 == 0 or (year%100 != 0 and year%4 == 0)):
        return True
    else:
        return False


def days_in_month(year, month):
    dinm = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if (is_leap(year)):
        dinm[1] = 29
    return dinm[month-1]


def juliand(datestr) : #datestr YYYYMMDDHH....
    year=int(datestr[0:4])
    mon =int(datestr[4:6])
    day =int(datestr[6:8])
    t =time.mktime((year,mon,day,0,0,0,0,0,0))
    return str(time.gmtime(t)[7])


def timesplit(time_begin,time_end):
    timelist=[]
    hourstep=[]
    for i in range(0,365*24*50):
        d=datetime.datetime(int(time_begin[0:4]),int(time_begin[4:6]),int(time_begin[6:8]), \
            int(time_begin[8:10]) )+datetime.timedelta(hours=i)
        dstr = d.strftime("%Y%m%d%H")
        if ( int(time_begin) < int(dstr) < int(time_end) ) :
            if ((d.day==1)and(d.hour==0)) : #first day,hour of month
                timelist.append(dstr)
                hourstep.append(i)
        elif ( int(dstr) == int(time_begin) ):
            timelist.append(dstr)
            hourstep.append(i)
        elif ( int(dstr) == int(time_end) ):
            timelist.append(dstr)
            hourstep.append(i)
            break
    return (timelist,hourstep,juliand(time_begin))


def datelist(date_begin,date_end):
    dlist=[]
    for i in range(-1,365*50):
        d=datetime.datetime(int(date_begin[0:4]),int(date_begin[4:6]),int(date_begin[6:8])) \
          +datetime.timedelta(days=i)
        dstr = d.strftime("%Y%m%d")
        if ( int(dstr) <= int(date_end[0:8]) ):
            dlist.append(dstr)
        elif ( int(dstr) > int(date_end[0:8]) ):
            break
    return (dlist)


def symlink(src,dst):
    if ( os.path.exists(dst) ):
        os.remove(dst)
    os.symlink(src,dst)


def sysecho(s):
    os.system("echo \""+s+"\"")
