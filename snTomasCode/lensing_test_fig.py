
#------------------------------------------------------------
# The lensing comparison figure
import numpy as np

def mkLensingTestFig( presfig=False ):
    from matplotlib import pyplot as pl
    from matplotlib.patches import FancyArrowPatch
    from pytools import plotsetup, colorpalette as cp
    from astropy.io import ascii
    import os
    import sys

    if presfig :
        fig = plotsetup.presfig( figsize=[12,8])
        labelvalues=True
        pl.clf()
        ax1 = pl.axes( [0.01,0.14,0.67,0.83] )
        ms=15
    else :
        labelvalues=False
        fig = plotsetup.fullpaperfig( figsize=[4,5])
        ms=8
        pl.clf()
        ax1 = pl.axes( [0.01,0.08,0.7,0.86] )

    thisfile = sys.argv[0]
    if 'ipython' in thisfile :
        thisfile = __file__
    thisdir = os.path.dirname( os.path.abspath(thisfile))

    lensdatfile = os.path.join(thisdir,'data/lensing/lensing_medians.dat')
    lensingdat = ascii.read( lensdatfile, format='commented_header',
                             header_start=-1, data_start=0 )

    ytop=15.2
    y = ytop-0.7
    for row in lensingdat:
        y-=1
        color = row['parametric'] and 'k' or 'g'
        mfc = row['postFF'] and color or 'w'
        model = row['model']
        best = row['best']
        med = row['med']
        errplus = row['med+']-med
        errminus = row['med-']-med
        if row['strong+weak']:
            marker='D'
        else:
            marker='o'

        #ax1.plot( best, y, marker='x', ms=ms, color=color,
        #          zorder=100, label='_nolegend_' )
        if med+errplus>3.6 :
            # make an arrow for the really long williams error
            errplus = 3.5-med
            ax1.errorbar( med, y, yerr=None, xerr=[[0],[errplus]],
                          marker=None, capsize=0.6*ms, color=color,
                          xuplims=[True], zorder=10, label='_nolegend_' )
        ax1.errorbar( med, y, yerr=None,
                             xerr=[[abs(errminus)],[errplus]],
                             marker=marker, mfc=mfc, mec=color, ms=ms,
                             capsize=1, color=color, zorder=10,
                             label='_nolegend_' )
        if labelvalues:
            label = '%s: $%.2f^{%+.2f}_{%+.2f}$'%( model, med, errplus, errminus )
            xlabel=1.35
            ha='right'
        else :
            label = '%s'%( model )
            xlabel=1.01
            ha='left'

        ax1.text( xlabel, float(y)/ytop, label, ha=ha, va='center',
                  color=color, transform=ax1.transAxes)
    fig.text( 0.98, 0.7, 'Parametric', ha='right', va='center', size='large',
              color='k', transform=fig.transFigure, rotation=-90)
    fig.text( 0.98, 0.32, 'Free-Form', ha='right', va='center',size='large',
              color='g', transform=fig.transFigure, rotation=-90)

    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['top'].set_visible(False)

    np.std(lensingdat['med'])

    ipar = np.where(lensingdat['parametric'])[0]
    iff = np.where(lensingdat['parametric']==0)[0]
    np.median(lensingdat['med'])
    for muSN, muSNerr, ymin, ymax, color, textcolor, ls, fitter in zip(
            [2.00,1.88,np.median(lensingdat['med'][ipar]),np.median(lensingdat['med'][iff])],
            [0.28,0.35,np.std(lensingdat['med'][ipar]),np.std(lensingdat['med'][iff])],
            [0.,0.,0.4,0.],
            [0.88,0.8,0.92,0.4],
            [cp.lightblue,cp.darkred,cp.lightgrey,cp.darkgreen],
            [cp.darkblue,cp.darkred,cp.black,cp.darkgreen],
            ['-','-','--','--'],
            ['MLCS2k2','SALT2','Lens Models',''] ):
        ax1.axvline( muSN, ls=ls, color=color, lw=2, ymin=ymin, ymax=ymax )
        if not fitter.lower().startswith('foo'):
            ax1.axvspan( muSN-muSNerr, muSN+muSNerr, ymin=ymin, ymax=ymax,
                         color=color, alpha=0.3,zorder=-100 )
        ax1.text( muSN-0.03,  ymax*ytop+0.05, fitter, color=textcolor,
                  ha='center', va='bottom' )
        print( "%s : mu=%.2f +- %.2f"%(fitter,muSN,muSNerr))
    ax1.text( 0.12, 0.4, 'SN HFF14tom', color='k', rotation=90, ha='right', va='bottom', transform=ax1.transAxes )
    #SLpre = ax1.plot(0,0,marker='o',mfc='w',mec='0.5',ls=' ',ms=ms,label='strong')
    #SLWLpre = ax1.plot(0,0,marker='D',mfc='w',mec='0.5',ls=' ',ms=ms,label='str+wk')
    #SLpost = ax1.plot(0,0,marker='o',mfc='k',mec='0.5',ls=' ',ms=ms,label='strong')
    #SLWLpost = ax1.plot(0,0,marker='D',mfc='k',mec='0.5',ls=' ',ms=ms,label='str+wk')
    #best = ax1.plot(0,0,marker='x',mfc='w',mec='0.5',ls=' ',ms=ms,label='best')
    #ax1.legend(loc='lower right',
    #           bbox_to_anchor=[0.98,0], bbox_transform=fig.transFigure,
    #           numpoints=1, frameon=False, handletextpad=-0.4,labelspacing=0.15)

    axleg = pl.axes([0.5,0.87,0.3,0.1], frameon=True )
    SLpre    = axleg.plot(0,1,marker='o',mfc='w',mec='0.5',ls=' ',ms=ms,label='strong')
    SLWLpre  = axleg.plot(0,0,marker='D',mfc='w',mec='0.5',ls=' ',ms=ms,label='str+wk')
    SLpost   = axleg.plot(1,1,marker='o',mfc='k',mec='0.5',ls=' ',ms=ms,label='strong')
    SLWLpost = axleg.plot(1,0,marker='D',mfc='k',mec='0.5',ls=' ',ms=ms,label='str+wk')
    axleg.text( -0.5, 1.5, 'pre/post-HFF', ha='left',va='bottom' )
    axleg.text( 1.5, 1.0, 'strong', va='center' )
    axleg.text( 1.5, 0.0, 'str+wk', va='center' )
    axleg.set_xlim(-0.8,3.5)
    axleg.set_ylim(-0.6,2.5)
    axleg.xaxis.set_ticks([])
    axleg.yaxis.set_ticks([])

    ax1.set_xlim( 1.57, 3.6 )
    ax1.set_ylim( 0, ytop )

    ax1.set_yticklabels( [] )
    ax1.set_yticks( [] )
    ax1.xaxis.set_ticks_position('bottom')
    ax1.set_xlabel(r'Lensing Magnification, $\mu$')
    pl.draw()
