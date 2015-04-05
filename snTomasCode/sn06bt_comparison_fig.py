__author__ = 'rodney'

def getdata():
    """ read in the 06bt and tomas rest-frame UBVri light curves
    :return:
    """
    import sys
    import os
    from astropy.io import ascii
    from astropy.table import Table, Column
    import sncosmo
    import numpy as np

    thisfile = sys.argv[0]
    if 'ipython' in thisfile :
        thisfile = __file__
    thisdir = os.path.abspath( os.path.dirname(thisfile))
    datadir = os.path.join( thisdir, "data/sn06bt/lightcurve")

    # read in the 06bt and tomas rest-frame UBVri light curves
    sn06bt = ascii.read(os.path.join(datadir,'sn06bt_lightcurve.dat'))
    tomas = ascii.read(os.path.join(datadir,'tomas_UBVri.dat'))

    # simulate a normal SN Ia
    snIa = sncosmo.Model( 'salt2' )
    snIa.set(c=-0.15)
    snIa.set(z=0.001)
    trest=[]
    mag = []
    magerr = []
    bandlist=[]
    t = np.arange(-5,25,0.1)
    zeros = np.zeros(len(t))
    for band in ['bessellux','bessellb','bessellv','sdssr']:
        m = snIa.bandmag( band, 'Vega', t )
        trest = np.append( trest, t )
        mag = np.append( mag, m )
        magerr = np.append( magerr, zeros )
        bandlist = np.append( bandlist, [band for i in range(len(m))])
    snIaUBVr = Table({'trest':trest,'mag':mag,'magerr':magerr,'filter':bandlist})

    return( sn06bt, tomas, snIaUBVr )

def plotcolorcurve( band1='bessellb', band2='bessellv', **plotargs):
    """  plot color vs time
    :return:
    """
    from matplotlib import pyplot as pl
    import numpy as np

    # read in the 06bt and tomas rest-frame UBVri light curves
    sn06bt, tomas, snIa = getdata()

    ax = pl.gca()

    # for band1,band2 in (['bessellb','bessellv'],['bessellb','sdssr'],['sdssr','sdssi']):
    for sn, clr, marker, ls in zip([sn06bt,tomas, snIa],
                                     ['darkorange','0.5','b'],
                                     [' ','d',' '],
                                     ['-',' ','--']):
        time, terr = [], []
        color, colorerr = [], []
        if 'trest' not in sn.colnames :
            sn['trest'] = sn['mjd'] - 53856.4
        ib1list = np.where(sn['filter']==band1)[0]
        ib2list = np.where(sn['filter']==band2)[0]
        if len(ib1list)==0 : continue
        if len(ib2list)==0 : continue

        t1list = sn['trest'][ib1list]
        t2list = sn['trest'][ib2list]
        m1list = sn['mag'][ib1list]
        m2list = sn['mag'][ib2list]
        merr1list = sn['magerr'][ib1list]
        merr2list = sn['magerr'][ib2list]

        for i1,t1,m1,merr1 in zip(ib1list,t1list,m1list,merr1list):
            i2 = np.argmin( np.abs(t2list-t1) )
            t2 = t2list[i2]
            m2 = m2list[i2]
            merr2 = merr2list[i2]
            if np.abs(t2-t1)>2 : continue
            t = (t1+t2)/2.
            if t in time : continue
            time.append(t)
            color.append( m1-m2 )
            colorerr.append( np.sqrt(merr1**2 + merr2**2))
            terr.append(np.std([t1,t2]))

        if marker!=' ':
            ax.errorbar( time, color, colorerr, terr,
                         mec=clr, capsize=0,
                         ls=ls, marker=marker, color=clr, **plotargs )
        else:
            ax.plot( time, color, ls=ls, marker=marker, color=clr, **plotargs )

def colorfig( **plotargs ):
    """  make a figure comparing the 06bt colors to SN Tomas colors
    :return:
    """
    from pytools import plotsetup
    from matplotlib import pyplot as pl


    fig = plotsetup.halfpaperfig()
    fig.clf()

    ax1 = fig.add_subplot(2,2,1)
    plotcolorcurve( 'bessellux', 'bessellb',  **plotargs )
    ax1.xaxis.set_ticks_position('top')
    ax1.xaxis.set_ticks_position('both')
    pl.ylabel('U-B',labelpad=-5)
    ax1.set_ylim(-0.49,0.39)

    ax2 = fig.add_subplot(2,2,2, sharex=ax1)
    plotcolorcurve( 'bessellux', 'bessellv',  **plotargs )
    ax2.yaxis.set_ticks_position('right')
    ax2.yaxis.set_ticks_position('both')
    ax2.xaxis.set_ticks_position('top')
    ax2.xaxis.set_ticks_position('both')
    ax2.yaxis.set_label_position('right')
    pl.ylabel('U-V', rotation=-90)
    ax2.set_ylim(-0.45,1.19)

    ax3 = fig.add_subplot(2,2,3, sharex=ax1)
    plotcolorcurve( 'bessellb', 'bessellv',  **plotargs )
    pl.ylabel('B-V')
    ax3.set_ylim(0.01,0.99)

    ax4 = fig.add_subplot(2,2,4, sharex=ax1)
    plotcolorcurve( 'bessellb', 'sdssr',  **plotargs )
    ax4.yaxis.set_ticks_position('right')
    ax4.yaxis.set_ticks_position('both')
    ax4.yaxis.set_label_position('right')
    pl.ylabel('B-r',rotation=-90, labelpad=10)
    ax4.set_ylim(-0.14,1.09)
    ax4.text(5, 0.8,'SN 2006bt',rotation=45, color='darkorange')
    ax4.text(6.8, 0.37,'\\noindent normal Ia\\\\ $c=-0.12$',rotation=40, color='b', ha='left')

    fig.subplots_adjust(left=0.15,right=0.87,bottom=0.12,top=0.92,hspace=0,wspace=0)
    suplabel( 'x', 'Rest-frame time (days)', labelpad=8, ha='center', va='bottom')

    ax1.set_xlim(-1,18)

    pl.draw()



def plotlc(band='bessellb', offset2=0, **plotargs):
    from matplotlib import pyplot as pl
    import numpy as np

    ax = pl.gca()
    sn06bt, tomas = getdata()
    for sn, offset, marker, ls in zip([sn06bt,tomas],[0,offset2],[' ','d'],['-',' ']):
        if 'trest' not in sn.colnames :
            sn['trest'] = sn['mjd'] - 53856.4
        i1list = np.where(sn['filter']==band)[0]
        if ls=='-':
            ax.plot( sn['trest'][i1list], sn['mag'][i1list],
                     marker=marker, ls=ls, **plotargs)
        else :
            ax.errorbar( sn['trest'][i1list], sn['mag'][i1list]+offset2,
                         sn['magerr'][i1list],
                         marker=marker, ls=ls, **plotargs)
    ax.invert_yaxis()

def lcfig( offset2=-0.76, **plotargs ):
    """  make a figure comparing the 06bt light curve to tomas
    :return:
    """
    from pytools import plotsetup

    fig = plotsetup.halfpaperfig()

    ax1 = fig.add_subplot(2,2,1)
    plotlc( 'bessellux', offset2=offset2, color='b', **plotargs )

    ax2 = fig.add_subplot(2,2,2)
    plotlc( 'bessellb', offset2=offset2, color='g', **plotargs )

    ax3 = fig.add_subplot(2,2,3)
    plotlc( 'bessellv', offset2=offset2, color='darkorange', **plotargs )

    ax4 = fig.add_subplot(2,2,4)
    plotlc( 'sdssr', offset2=offset2, color='r', **plotargs )


def suplabel(axis,label,label_prop=None,
             labelpad=5,
             ha='center',va='center'):
    ''' Add super ylabel or xlabel to the figure
    Similar to matplotlib.suptitle
    axis       - string: "x" or "y"
    label      - string
    label_prop - keyword dictionary for Text
    labelpad   - padding from the axis (default: 5)
    ha         - horizontal alignment (default: "center")
    va         - vertical alignment (default: "center")
    '''
    from matplotlib import pyplot as pl

    fig = pl.gcf()
    xmin = []
    ymin = []
    for ax in fig.axes:
        xmin.append(ax.get_position().xmin)
        ymin.append(ax.get_position().ymin)
    xmin,ymin = min(xmin),min(ymin)
    dpi = fig.dpi
    if axis.lower() == "y":
        rotation=90.
        x = xmin-float(labelpad)/dpi
        y = 0.5
    elif axis.lower() == 'x':
        rotation = 0.
        x = 0.5
        y = ymin - float(labelpad)/dpi
    else:
        raise Exception("Unexpected axis: x or y")
    if label_prop is None:
        label_prop = dict()
    pl.text(x,y,label,rotation=rotation,
               transform=fig.transFigure,
               ha=ha,va=va,
               **label_prop)
