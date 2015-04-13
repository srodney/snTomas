
#------------------------------------------------------------
# The lensing comparison figure
import numpy as np

def mkLensingTestFig( show2snfits=False, showlegend=False, presfig=False ):
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
        if showlegend:
            ax1 = pl.axes( [0.01,0.14,0.7,0.84] )
        else:
            ax1 = pl.axes( [0.01,0.14,0.7,0.84] )


    thisfile = sys.argv[0]
    if 'ipython' in thisfile :
        thisfile = __file__
    thisdir = os.path.dirname( os.path.abspath(thisfile))

    lensdatfile = os.path.join(thisdir,'data/lensing/lensing_medians.dat')
    lensingdat = ascii.read( lensdatfile, format='commented_header',
                             header_start=-1, data_start=0 )

    ytop=len(lensingdat['med'])+1
    y = ytop
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
        if med+errplus>3.65 :
            # make an arrow for the really long williams error
            errplus = 3.6-med
            ax1.errorbar( med, y, yerr=None, xerr=[[0],[errplus]],
                          marker=None, capsize=0.6*ms, color=color,
                          xlolims=[True], zorder=10, label='_nolegend_' )
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
    fig.text( 0.98, 0.75, 'Parametric', ha='right', va='center', size='large',
              color='k', transform=fig.transFigure, rotation=-90)
    fig.text( 0.98, 0.42, 'Free-Form', ha='right', va='center',size='large',
              color='g', transform=fig.transFigure, rotation=-90)

    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['top'].set_visible(False)

    np.std(lensingdat['med'])

    ipar = np.where(lensingdat['parametric'])[0]
    iff = np.where(lensingdat['parametric']==0)[0]
    np.median(lensingdat['med'])
    for muSN, muSNerr, ymin, ymax, color, textcolor, ls, fitter in zip(
            [1.81,1.79,np.median(lensingdat['med'][ipar]),np.median(lensingdat['med'][iff])],
            [0.26,0.34,np.std(lensingdat['med'][ipar]),np.std(lensingdat['med'][iff])],
            [0.,0.,0.4,0.],
            [1.,1.,1.,0.4],
            [cp.lightblue,cp.darkred,cp.lightgrey,cp.darkgreen],
            [cp.darkblue,cp.darkred,cp.black,cp.darkgreen],
            ['-','-','--','--'],
            ['MLCS2k2','SALT2','Parametric','Free-Form'] ):
        if not show2snfits :
            if fitter=='MLCS2k2': continue
        if not ls=='--':
            ax1.axvline( muSN, ls=ls, color=color, lw=2, ymin=ymin, ymax=ymax )
            ax1.axvspan( muSN-muSNerr, muSN+muSNerr, ymin=ymin, ymax=ymax,
                         color=color, alpha=0.3,zorder=-100 )
        #ax1.text( muSN-0.03,  ymax*ytop+0.05, fitter, color=textcolor,
        #          ha='center', va='bottom' )
        print( "%s : mu=%.2f +- %.2f"%(fitter,muSN,muSNerr))
    ax1.text( 1.75, 2*ytop/3, 'SN HFF14tom', rotation=90,
              ha='right', va='center', color=cp.darkred )

    if showlegend:
        if showlegend=='top':
            axleg = pl.axes([0.5,0.87,0.3,0.1], frameon=True )
        else :
            axleg = pl.axes([0.68,0.02,0.29,0.11], frameon=True )

        SLpre    = axleg.plot(0,1,marker='o',mfc='w',mec='0.5',ls=' ',ms=ms,label='strong')
        SLWLpre  = axleg.plot(0,0,marker='D',mfc='w',mec='0.5',ls=' ',ms=ms,label='str+wk')
        SLpost   = axleg.plot(1,1,marker='o',mfc='k',mec='0.5',ls=' ',ms=ms,label='strong')
        SLWLpost = axleg.plot(1,0,marker='D',mfc='k',mec='0.5',ls=' ',ms=ms,label='str+wk')
        axleg.text( -0.5, 1.5, 'pre/post-HFF', ha='left',va='bottom' )
        axleg.text( 1.5, 1.0, 'strong', va='center' )
        axleg.text( 1.5, 0.0, 'str+wk', va='center' )
        axleg.set_xlim(-0.7,3.4)
        axleg.set_ylim(-0.6,2.7)
        axleg.xaxis.set_ticks([])
        axleg.yaxis.set_ticks([])

    ax1.set_xlim( 1.36, 3.7 )
    ax1.set_ylim( 0, ytop )

    ax1.set_yticklabels( [] )
    ax1.set_yticks( [] )
    ax1.xaxis.set_ticks_position('bottom')
    ax1.set_xlabel(r'Lensing Magnification, $\mu$')
    pl.draw()
