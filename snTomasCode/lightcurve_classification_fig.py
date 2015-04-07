
# Results of the classification run using sncosmo.classify()
# sncosmo.classify.classify( snTomas, zhost = 1.30, zhosterr=0.5,
#                       t0_range=[56800,56840], zminmax=[0.9,1.8], verbose=3 )

_bestfitIa = {
    'model':'salt2-extended', 'class':'Ia',
    'z': 1.35126, 'zerr':0.01608,
    'c':-0.09321, 'cerr':0.01946,
    'x1': 0.03530,'x1err':0.18911,
    'mB':24.66244,'mBerr':0.01881,
    't0':56816.28,'t0err':0.33732,
    'x0':2.171e-6,'x0err':0.038e-6,
}

_bestfitII = {
    'model':'snana-2007pg', 'class':'II',
    'z':1.79865, 'zerr':0.00134,
    't0':56808.6, 't0err':0.25432,
    'amplitude':2.98840e-18, 'amplitudeerr':9.81812e-20,
    'hostebv':0.00870956, 'hostebverr':0.00694,
    'hostr_v':2.94742, 'hostr_verr':0.25432,
    }


_bestfitIbc = {
    'model':'snana-sdss014475', 'class':'Ib/c',
    'z':0.695, 'zerr':0.016,
    't0':56805.4, 't0err':0.34,
    'amplitude':2.263e-17, 'amplitudeerr':3.289e-18,
    'hostebv':0.349, 'hostebverr':0.044,
    'hostr_v':3.72, 'hostr_verr':0.24,
    }


def mk_figure(yunit='flux'):
    """ make a figure showing the SN Tomas light curve classification
    best-fit models for each class
    :return:
    """
    import numpy as np
    from matplotlib import pyplot as pl, ticker, gridspec, rcParams
    import sys
    import os
    from pytools import colorpalette as cp, plotsetup
    from astropy.io import ascii
    import sncosmo

    bandlist = ['f814w','f105w','f125w','f140w','f160w']
    thisfile = sys.argv[0]
    if 'ipython' in thisfile :
        thisfile = __file__
    thisdir = os.path.abspath( os.path.dirname(thisfile))

    obsdat = ascii.read(os.path.join(thisdir,'data/HST_FFSN_tomas.dat'),
                        format='commented_header',
                        header_start=-1, data_start=0)
    if yunit=='mag':
        yobs = obsdat['mag']
        yerrobs = obsdat['magerr']
    else :
        yobs = obsdat['flux'] * 10**(-0.4*(obsdat['zpt']-25))
        yerrobs = obsdat['fluxerr']  * 10**(-0.4*(obsdat['zpt']-25))

    mjdmod = np.arange(obsdat['mjd'].min()-5, obsdat['mjd'].max()+20,1)

    #fig = plotsetup.fullpaperfig([8,3.5])
    fig = plotsetup.halfpaperfig(figsize=[4,7.5])
    rcParams['xtick.major.pad'] = 4
    pl.clf()
    gs = gridspec.GridSpec( 4,2)

    axdict = {
        'F814W':pl.subplot(gs[0, :]),
        'F105W':pl.subplot(gs[1, :]),
        'F125W':pl.subplot(gs[2, 0]),
        'F140W':pl.subplot(gs[2, 1]),
        'F160W':pl.subplot(gs[3, :]),
        }

    # plot the observations as grey points, one band per axis
    for band in bandlist :
        ax = axdict[band.upper()]
        ifilt = np.where(obsdat['filter']==band.upper())
        ax.errorbar( obsdat['mjd'][ifilt], yobs[ifilt], yerrobs[ifilt],
                     color='k',mfc='k', mec='k',
                     ms=8, marker='o', capsize=0,
                     ls=' ', alpha=0.5, zorder=-100)
        ax.text(0.95,0.9,band.upper(),color='k',ha='right',va='top',
                transform=ax.transAxes,fontsize='large')

    # overlay the models as colored lines
    for dashes,color,bestfitdict in zip(
            [ None,[8,3],[8,3,3,3]],
            [cp.darkorange,cp.darkgreen,cp.bluegrey],
            [_bestfitIa,_bestfitIbc,_bestfitII] ):
        if bestfitdict['class']!='Ia':
            dust = sncosmo.CCM89Dust( )
            model = sncosmo.Model(bestfitdict['model'], effects=[dust],
                                  effect_names=['host'],
                                  effect_frames=['rest'])
        else :
            model = sncosmo.Model(bestfitdict['model'])

        model.parameters = [bestfitdict[parname]
                            for parname in model.param_names]
        # if bestfitdict['class']=='Ic':
        #    import pdb; pdb.set_trace()
        for band in bandlist :
            ax = axdict[band.upper()]
            if yunit=='mag':
                ymod = model.bandmag(band,'ab',mjdmod)
            else:
                ymod = model.bandflux(band,mjdmod,zp=25,zpsys='ab')
            if dashes:
                ax.plot( mjdmod, ymod, dashes=dashes, color=color, lw=1.5,
                         label=bestfitdict['class'] )
            else :
                ax.plot( mjdmod, ymod, ls='-', color=color, lw=2.0,
                         label=bestfitdict['class'] )


    for ax in axdict.values():
        ax.set_ylim(-0.2,ax.get_ylim()[1])

    ax814 = axdict['F814W']
    # ax814.set_xlim(56790,56849)
    ax814.set_xlim(56790,56969)
    ax814.xaxis.set_major_locator( ticker.MultipleLocator(50) )
    ax814.xaxis.set_minor_locator( ticker.MultipleLocator(10) )
    ax814.axvline(56813.5,ls='-',color='0.5',lw=0.75)
    ax814.legend(loc='right', frameon=False, labelspacing=0.08,
                 fontsize=11.5, handletextpad=0.08)

    ax105 = axdict['F105W']
    ax105.set_xlim(56790,56969)
    ax105.xaxis.set_major_locator( ticker.MultipleLocator(50) )
    ax105.xaxis.set_minor_locator( ticker.MultipleLocator(10) )
    ax105.axvline(56813.5,ls='-',color='0.5',lw=0.75)
    ax105.text(56812,0.2,'spectrum',rotation=90,ha='right',va='bottom',
               color='0.5', fontsize='small')

    ax125 = axdict['F125W']
    ax125.set_xlim(56790,56877)
    ax125.xaxis.set_major_locator( ticker.MultipleLocator(50) )
    ax125.xaxis.set_minor_locator( ticker.MultipleLocator(10) )
    ax125.axvline(56813.5,ls='-',color='0.5',lw=0.75)

    ax140 = axdict['F140W']
    pl.setp(ax140.yaxis.get_ticklabels(), visible=False)
    ax140.set_ylim(ax125.get_ylim())
    ax140.set_ylabel('')
    ax140.set_xlim(56883,56969)
    ax140.xaxis.set_major_locator( ticker.MultipleLocator(50) )
    ax140.xaxis.set_minor_locator( ticker.MultipleLocator(10) )
    ax140.axvline(56813.5,ls='-',color='0.5',lw=0.75)

    ax160 = axdict['F160W']
    ax160.set_xlim(56790,56969)
    ax160.xaxis.set_major_locator( ticker.MultipleLocator(50) )
    ax160.xaxis.set_minor_locator( ticker.MultipleLocator(10) )
    ax160.set_xlabel('Observer-frame Days (MJD)')
    ax160.axvline(56813.5,ls='-',color='0.5',lw=0.75)

    for ax in axdict.values():
        ax.yaxis.set_major_locator( ticker.MultipleLocator(1) )
        ax.yaxis.set_minor_locator( ticker.MultipleLocator(0.5) )
        ax.ticklabel_format(axis='x',style='plain',useOffset=False)
        # ax.xaxis.get_major_formatter().set_powerlimits((0, 1))

    if yunit=='mag':
        ax105.set_ylabel('AB mag',va='center',ha='center')
    else :
        ax105.set_ylabel('Flux (zp$_{\\rm AB}$=25)',va='center',ha='center')
    ax105.yaxis.set_label_coords(0.05, 0.5, transform=fig.transFigure)

    fig.subplots_adjust(left=0.12,right=0.97,bottom=0.08,top=0.97,wspace=0.05)
    pl.draw()


