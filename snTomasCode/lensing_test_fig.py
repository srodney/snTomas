
#------------------------------------------------------------
# The lensing comparison figure
import numpy as np

def mkLensingTestFig( show2snfits=False, showlegend=True, presfig=False,
                      labelvalues=False ):
    from matplotlib import pyplot as pl
    from matplotlib.patches import FancyArrowPatch
    from pytools import plotsetup, colorpalette as cp
    from astropy.io import ascii
    import os
    import sys

    muSNmlcs = 2.03
    muSNmlcserr = 0.29
    muSNsalt = 1.99
    muSNsalterr = 0.38

    if presfig :
        fig = plotsetup.presfig( figsize=[8,12])
        pl.clf()
        ax1 = pl.axes( [0.01,0.14,0.67,0.83] )
        ms=15
    else :
        fig = plotsetup.fullpaperfig( figsize=[4,6])
        ms=8
        pl.clf()
        if showlegend:
            ax1 = pl.axes( [0.01,0.14,0.68,0.84] )
        else:
            ax1 = pl.axes( [0.01,0.14,0.68,0.84] )


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
        # color = row['parametric'] and 'k' or 'g'
        #mfc = row['postFF'] and color or 'w'
        color = row['postFF'] and 'g' or 'k'
        mfc = row['parametric'] and color or 'w'
        model = row['model']
        best = row['best']
        med = row['med']
        errplus = row['med+']-med
        errminus = row['med-']-med
        if row['strong+weak']:
            marker='D'
        else:
            marker='o'
        mec=color
        mew=ms/8.

        #ax1.plot( best, y, marker='x', ms=ms, color=color,
        #          zorder=100, label='_nolegend_' )
        if med+errplus>3.65 :
            # make an arrow for the really long williams error
            errplus = 3.6-med
            ax1.errorbar( med, y, yerr=None, xerr=[[0],[errplus]],
                          marker=None, capsize=0.6*ms, color=color, mec=mec,
                          ms=ms, mew=mew, xlolims=[True], zorder=10,
                          label='_nolegend_' )
        ax1.errorbar( med, y, yerr=None,
                      xerr=[[abs(errminus)],[errplus]],
                      marker=marker, mfc=mfc, mec=mec, mew=mew, ms=ms,
                      capsize=1, color=color, zorder=10,
                      label='_nolegend_' )
        if row['postSN']:
            ax1.plot( med, y, marker='s', mfc='w', mec='darkorange',
                      mew=ms/7, ms=1.3*ms, zorder=-10,
                      label='_nolegend_' )

        print( "%s discrepancy = %.2f sigma"%(model,(med-muSNmlcs)/np.sqrt(errminus**2+muSNmlcserr**2)))
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
    fig.text( 0.98, 0.75, 'Pre-HFF', ha='right', va='center', size='large',
              color='k', transform=fig.transFigure, rotation=-90)
    fig.text( 0.98, 0.34, 'Post-HFF', ha='right', va='center',size='large',
              color='g', transform=fig.transFigure, rotation=-90)
    #fig.text( 0.98, 0.26, 'Post-SN', ha='right', va='center',size='large',
    #          color='g', transform=fig.transFigure, rotation=-90)

    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['top'].set_visible(False)

    np.std(lensingdat['med'])

    ipar = np.where(lensingdat['parametric'])[0]
    iff = np.where(lensingdat['parametric']==0)[0]

    ipostFF = np.where(lensingdat['postFF'])[0]
    ipreFF = np.where(lensingdat['postFF']==0)[0]
    ipostSN = np.where(lensingdat['postSN'])[0]

    np.median(lensingdat['med'])
    for muSN, muSNerr, ymin, ymax, color, textcolor, ls, fitter in zip(
            [muSNmlcs,muSNsalt,np.median(lensingdat['med'][ipreFF]),np.median(lensingdat['med'][ipostFF])],
            [muSNmlcserr,muSNsalterr,np.std(lensingdat['med'][ipreFF]),np.std(lensingdat['med'][ipostFF])],
            [0.,0.,0.4,0.],
            [1.,1.,1.,0.4],
            [cp.lightblue,cp.darkred,cp.lightgrey,cp.darkgreen],
            [cp.darkblue,cp.darkred,cp.black,cp.darkgreen],
            ['-','-','--','--'],
            ['MLCS2k2','SALT2','Pre-HFF','Post-HFF'] ):
        if not show2snfits :
            if fitter=='SALT2': continue
        if not ls=='--':
            ax1.axvline( muSN, ls=ls, color=textcolor, lw=2, ymin=ymin, ymax=ymax )
            ax1.axvspan( muSN-muSNerr, muSN+muSNerr, ymin=ymin, ymax=ymax,
                         color=color, alpha=0.3,zorder=-100 )
        #ax1.text( muSN-0.03,  ymax*ytop+0.05, fitter, color=textcolor,
        #          ha='center', va='bottom' )
        print( "%s : mu=%.2f +- %.2f"%(fitter,muSN,muSNerr))
    ax1.text( 1.75, 2*ytop/3, 'SN HFF14tom', rotation=90,
              ha='right', va='center', color=cp.darkblue )

    ax1.axvline( np.mean(lensingdat['med']), ls='--', color='0.5', lw=1 )
    print( "All models : mu=%.2f +- %.2f"%(np.mean(lensingdat['med']),np.std(lensingdat['med'])))

    if showlegend:
        if showlegend=='top':
            axleg = pl.axes([0.5,0.87,0.3,0.1], frameon=True )
        else :
            axleg = pl.axes([0.68,0.02,0.3,0.14], frameon=True )

        SLpre    = axleg.plot(0,2,marker='o',mfc='w',mec='0.5',ls=' ',ms=ms,label='strong')
        SLWLpre  = axleg.plot(0,1,marker='D',mfc='w',mec='0.5',ls=' ',ms=ms,label='str+wk')
        SLpost   = axleg.plot(1,2,marker='o',mfc='k',mec='0.5',ls=' ',ms=ms,label='strong')
        SLWLpost = axleg.plot(1,1,marker='D',mfc='k',mec='0.5',ls=' ',ms=ms,label='str+wk')
        unblind  = axleg.plot(0.5,0,marker='s',mfc='w',mec='darkorange',ls=' ',ms=ms,label='unblind')
        #axleg.text( -0.5, 2.5, 'pre/post-HFF', ha='left',va='bottom' )
        axleg.text( -0.45, 2.5, 'Fr. \ Par.', ha='left',va='bottom' )
        axleg.text( 1.5, 2.0, 'strong', va='center' )
        axleg.text( 1.5, 1.0, 'str+wk', va='center' )
        axleg.text( 1.4, 0.0, 'unblind', va='center' )
        axleg.set_xlim(-0.7,3.4)
        axleg.set_ylim(-0.6,3.7)
        axleg.xaxis.set_ticks([])
        axleg.yaxis.set_ticks([])

    ax1.set_xlim( 1.36, 3.7 )
    ax1.set_ylim( 0, ytop )

    ax1.set_yticklabels( [] )
    ax1.set_yticks( [] )
    ax1.xaxis.set_ticks_position('bottom')
    ax1.set_xlabel(r'Lensing Magnification, $\mu$')
    pl.draw()


def mkTensionFig(presfig=False, showlegend=True, showlines=True):
    from matplotlib import ticker,pyplot as pl, rcParams
    from pytools import plotsetup
    from astropy.io import ascii

    muSNmlcs = 2.03
    muSNmlcserr = 0.29
    muSNsalt = 1.99
    muSNsalterr = 0.38

    lensdatfile = 'snTomasCode/data/lensing/lensing_medians.dat'
    lensingdat = ascii.read( lensdatfile, format='commented_header',
                             header_start=-1, data_start=0 )

    if presfig :
        fig = plotsetup.presfig( wide=True )
        fig.subplots_adjust( wspace=0, bottom=0.17, left=0.1, right=0.95, top=0.95)
        ms=15
        rcParams['axes.labelsize']=26
    else :
        fig = plotsetup.fullpaperfig()
        fig.subplots_adjust( wspace=0, bottom=0.17, left=0.08, right=0.95, top=0.95)
        ms=8
    fig.clf()
    ax1 = pl.subplot(1,2,1)
    ax2 = pl.subplot(1,2,2, sharey=ax1)

    for row in lensingdat:
        # color = row['parametric'] and 'k' or 'g'
        #mfc = row['postFF'] and color or 'w'
        color = row['postFF'] and 'g' or 'k'
        mfc = row['parametric'] and color or 'w'

        model = row['model']
        best = row['best']
        med = row['med']
        errplus = row['med+']-med
        errminus = row['med-']-med
        if row['strong+weak']:
            marker='D'
        else:
            marker='o'
        mec=color
        mew=ms/8.
        specfrac = row['nzSpec']/float(row['nSys'])

        nim = row['nIm']/float(row['nSys'])

        nsys = row['nSys']
        tension = (med-muSNmlcs)/np.sqrt(errminus**2+muSNmlcserr**2)

        ax1.plot( nsys, tension, ls=' ', alpha=0.8,
                 marker=marker, mfc=mfc, mec=mec, mew=mew, ms=ms,
                 color=color )
        ax2.plot( specfrac, tension, ls=' ', alpha=0.8,
                 marker=marker, mfc=mfc, mec=mec, mew=mew, ms=ms,
                 color=color )
        if row['postSN']:
            ax1.plot( nsys, tension, ls=' ', alpha=0.8,
                      marker='s', mfc='w', mec='darkorange',
                      mew=ms/7, ms=1.3*ms, zorder=-10, color=color )
            ax2.plot( specfrac, tension, ls=' ', alpha=0.8,
                      marker='s', mfc='w', mec='darkorange',
                      mew=ms/7, ms=1.3*ms, zorder=-10, color=color )

    ax1.axhline(0,ls='--',color='0.5')
    ax2.axhline(0,ls='--',color='0.5')
    ax1.set_ylabel('Tension with SN ($\sigma$)')
    ax1.set_xlabel('\# multiply imaged systems')
    ax2.set_xlabel('Fraction with spec-z')
    ax1.set_ylim(-0.9,5.1)
    ax1.set_xlim(1,58)
    ax2.set_xlim(0.01,0.65)
    ax1.yaxis.set_major_locator( ticker.MultipleLocator( 1 ) )
    ax1.yaxis.set_minor_locator( ticker.MultipleLocator( 0.25 ) )
    ax1.xaxis.set_major_locator( ticker.MultipleLocator( 10 ) )
    ax1.xaxis.set_minor_locator( ticker.MultipleLocator( 2 ) )
    ax2.xaxis.set_major_locator( ticker.MultipleLocator( 0.2 ) )
    ax2.xaxis.set_minor_locator( ticker.MultipleLocator( 0.05 ) )
    ax2.yaxis.set_ticks_position('right')
    # pl.setp(ax2.get_yticklabels(), visible=False)

    if showlegend:
        axleg = pl.axes([0.75,0.68,0.18,0.22], frameon=True )
        SLpre    = axleg.plot(0,2,marker='o',mfc='w',mec='0.5',ls=' ',ms=ms,label='strong')
        SLWLpre  = axleg.plot(0,1,marker='D',mfc='w',mec='0.5',ls=' ',ms=ms,label='str+wk')
        SLpost   = axleg.plot(1,2,marker='o',mfc='k',mec='0.5',ls=' ',ms=ms,label='strong')
        SLWLpost = axleg.plot(1,1,marker='D',mfc='k',mec='0.5',ls=' ',ms=ms,label='str+wk')
        unblind  = axleg.plot(0.5,0,marker='s',mfc='w',mec='darkorange',ls=' ',ms=ms,label='unblind')
        axleg.text( -0.5, 2.5, 'Free \ Par.', ha='left',va='bottom' )
        axleg.text( 1.5, 2.0, 'strong', va='center' )
        axleg.text( 1.5, 1.0, 'str+wk', va='center' )
        axleg.text( 1.5, 0.0, 'unblind', va='center' )
        axleg.set_xlim(-0.8,3.8)
        axleg.set_ylim(-0.6,3.7)
        axleg.xaxis.set_ticks([])
        axleg.yaxis.set_ticks([])


    if showlines:
        for modname in ['CATS(v2','Bradac','Lam']:
            ithismod = np.where([mod.startswith(modname)
                                 for mod in lensingdat['model']])[0]
            med = lensingdat['med'][ithismod]
            errminus = lensingdat['med-'][ithismod] - med
            tension = (med-muSNmlcs)/np.sqrt(errminus**2+muSNmlcserr**2)

            nsys = lensingdat['nSys'][ithismod]
            specfrac = lensingdat['nzSpec'][ithismod]/nsys.astype(float)
            ax1.plot( nsys, tension, ls='-', marker=' ', color='0.5' )
            ax2.plot( specfrac, tension, ls='-', marker=' ', color='0.5' )
        ax1.text(48,4.1,'CATS(v2)',color='0.5',fontsize='small',ha='right',va='bottom')
        ax1.text(55,2.8,'2.1',color='0.5',fontsize='small',ha='center',va='top')
        ax1.text(25,2.8,'2.2',color='0.5',fontsize='small',ha='center',va='top')
        ax1.text(27,0.6,'v2',color='0.5',fontsize='small',ha='left',va='bottom')
        ax1.text(14,1.2,'Bradac(v1)',color='0.5',fontsize='small',ha='right',va='bottom')
        ax1.text(22,1.7,'Lam(v1)',color='0.5',fontsize='small',ha='left',va='bottom')
        ax1.text(8,-0.1,'v2',color='0.5',fontsize='small',ha='right',va='top')

        ax2.text(0.15,4.4,'CATS(v2)',color='0.5',fontsize='small',ha='right',va='bottom')
        ax2.text(0.13,2.8,'2.1',color='0.5',fontsize='small',ha='center',va='top')
        ax2.text(0.3,2.8,'2.2',color='0.5',fontsize='small',ha='center',va='top')
        ax2.text(0.1,1.3,'Bradac',color='0.5',fontsize='small',ha='right',va='bottom')
        ax2.text(0.08,1.25,'(v1)',color='0.5',fontsize='small',ha='right',va='top')
        ax2.text(0.28,0.3,'v2',color='0.5',fontsize='small',ha='left',va='bottom')
        ax2.text(0.23,1.6,'Lam(v1)',color='0.5',fontsize='small',ha='left',va='bottom')
        ax2.text(0.53,-0.2,'v2',color='0.5',fontsize='small',ha='left',va='center')

    pl.draw()



def mkLensingTestFig_simplified():
    """ Make a simplified version of the lensing test figure,
    for use in proposals, presentations and the like
    :param show2snfits:
    :param showlegend:
    :param presfig:
    :param labelvalues:
    :return:
    """
    from matplotlib import pyplot as pl
    from pytools import plotsetup, colorpalette as cp
    from astropy.io import ascii
    import os
    import sys

    muSNmlcs = 2.03
    muSNmlcserr = 0.29
    muSNsalt = 1.99
    muSNsalterr = 0.38

    fig = plotsetup.fullpaperfig( figsize=[4,4])
    ms=8
    pl.clf()
    ax1 = pl.axes( [0.03,0.14,0.92,0.80] )

    thisfile = sys.argv[0]
    if 'ipython' in thisfile :
        thisfile = __file__
    thisdir = os.path.dirname( os.path.abspath(thisfile))

    lensdatfile = os.path.join(thisdir,'data/lensing/lensing_medians_simplified.dat')
    lensingdat = ascii.read( lensdatfile, format='commented_header',
                             header_start=-1, data_start=0 )

    ytop=len(lensingdat['med'])+1
    y = ytop

    for row in lensingdat:
        y-=1
        model = row['model']
        med = row['med']
        errplus = row['med+']-med
        errminus = row['med-']-med
        mec='k'
        mew=ms/8.

        disc = (med-muSNmlcs)/np.sqrt(errminus**2+muSNmlcserr**2)
        if disc<1:
            mfc='forestgreen'
            marker='D'
            color='forestgreen'
        else:
            mfc = 'w'
            marker='o'
            color='k'

        if med+errplus>3.65 :
            # make an arrow for the really long williams error
            errplus = 3.6-med
            ax1.errorbar( med, y, yerr=None, xerr=[[0],[errplus]],
                          marker=None, capsize=0.6*ms, color=color, mec=mec,
                          ms=ms, mew=mew, xlolims=[True], zorder=10,
                          label='_nolegend_' )
        ax1.errorbar( med, y, yerr=None,
                      xerr=[[abs(errminus)],[errplus]],
                      marker=marker, mfc=mfc, mec=mec, mew=mew, ms=ms,
                      capsize=1, color=color, zorder=10,
                      label='_nolegend_' )

    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.axvline( muSNmlcs, ls='-', color=cp.darkblue, lw=2 )
    ax1.axvspan( muSNmlcs-muSNmlcserr, muSNmlcs+muSNmlcserr,
                 color=cp.lightblue, alpha=0.3, zorder=-100 )
    ax1.text( 1.99, 0.8*ytop,"""Measurement
    from Type Ia SN
    HFF14tom""",
              rotation=90, ha='right', va='center', color=cp.darkblue,
              linespacing=1.5)
    ax1.text(3.6, 0.4*ytop,
             """Lens Model
             Magnification
             Predictions
             within 1$\sigma$ of
             SN Measurement""",
             ha='right', va='top', color='forestgreen', linespacing=1.5)
    ax1.set_xlim( 1.36, 3.7 )
    ax1.set_ylim( 0, ytop )
    ax1.set_yticklabels( [] )
    ax1.set_yticks( [] )
    ax1.xaxis.set_ticks_position('bottom')
    ax1.set_xlabel(r'Lensing Magnification, $\mu$')
    pl.draw()

