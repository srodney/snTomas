
def getchi2(templatedatfile, obsdatfile='snTomas_acs_g800l_1d.dat'):
    import numpy as np
    from matplotlib import pyplot as pl, ticker
    import os
    import sys

    thisfile = sys.argv[0]
    if 'ipython' in thisfile:
        thisfile = __file__
    thisdir = os.path.dirname( thisfile )
    obsdatfile = thisdir+'/data/spectrum/'+obsdatfile

    wave, flux, err, sens = np.loadtxt(obsdatfile, unpack=True)#, usecols=[0,1])
    wavemod,fluxmod= np.loadtxt(templatedatfile, unpack=True)
    chi2 = np.sum( (fluxmod-flux)**2 / err**2 )
    return( chi2 )


def mkspecfig3panel(plotflam=False, showerr=False):
    """ original function for making a 3-panel spectrum
    figure, showing the z=1.31 redshift solution
    :param plotflam:
    :param showerr:
    :return:
    """

    from astropy.io import ascii
    import numpy as np
    from matplotlib import pyplot as pl, ticker
    from pytools import plotsetup
    import sys
    import os

    thisfile = sys.argv[0]
    if 'ipython' in thisfile:
        thisfile = __file__
    thisdir = os.path.dirname( thisfile )

    plotsetup.halfpaperfig( figsize=[4,6])
    specfile = thisdir+'/data/spectrum/snTomas_acs_g800l_1d.dat'
    wave, flux, err, sens = np.loadtxt(specfile, unpack=True)#, usecols=[0,1])

    flambda = flux / sens * 1e-17
    flambdaErr = err / sens * 1e-17

    templist = [thisdir+'/data/spectrum/sm-cont-1.31000-sn2012cg-visit1-hst.flm',
                thisdir+'/data/spectrum/sm-1.31000-sn2011fe-visit4-hst.flm',
                thisdir+'/data/spectrum/sm-0.980000-sn2011iv-hst.flm',
                ]
    zlist = [1.31, 1.31, 0.98]
    labellist = ['SN 2012cg (warped)\nz=1.31\nage=+3',
                 'SN 2011fe\nz=1.31\nage=-3',
                 'SN 2011iv\nz=0.98\nage=0',
                 ]
    chi2list = [ getchi2(tempfile) for tempfile in templist ]
    ndoflist = [ len(wave)-6, len(wave)-2, len(wave)-2 ]

    fig = pl.gcf()

    fig.subplots_adjust( left=0.03, right=0.97, top=0.9, bottom=0.1, hspace=0.5)
    pl.clf()
    nax = len( templist )
    iax = 0
    colorlist = ['r','g','b']
    for z, template, color, label, chi2, ndof in zip( zlist, templist, colorlist, labellist, chi2list, ndoflist) :
        iax += 1
        wavemod,fluxmod= np.loadtxt(template, unpack=True)
        flambdamod = fluxmod / sens * 1e-17

        if iax==1 :
            ax = fig.add_subplot(nax,1,iax)
            ax1 = ax
        else :
            ax = fig.add_subplot(nax,1,iax, sharey=ax1)

        #if iax==2 :
        #ax.set_ylabel('Flux')
        ax.yaxis.set_ticklabels([])

        ax.set_xlim(5800, 9550)

        if plotflam :
            pl.plot( wavemod, flambdamod, color=color, marker=' ' )
            pl.plot( wave, flambda, color='k', marker=' ', linestyle='steps-mid' )
            if showerr:
                pl.fill_between( wave, flambda-flambdaerr, flambda+flambdaerr, color='k', alpha=0.3 )
        else :
            pl.plot( wavemod, fluxmod, color=color, marker=' ' )
            pl.plot( wave, flux, color='k', marker=' ', linestyle='steps-mid' )
            if showerr:
                pl.fill_between( wave, flux-err, flux+err, color='k', alpha=0.3 )
        axtop = ax.twiny()
        axtop.set_xlim( ax.get_xlim()[0]/(1+z), ax.get_xlim()[1]/(1+z) )
        #axtop.set_xlabel('Rest Wavelength at z=%.2f'%z)
        axtop.set_xlabel('Rest $\lambda$:')
        axtop.xaxis.set_label_coords( 0.1, 1.05)
        #ax.set_xlabel('Obs:')
        #axtop.set_ylabel('Flux')
        axtop.xaxis.set_major_locator( ticker.MultipleLocator( 500 ) )
        if iax==3:
            axtop.xaxis.set_ticks( [3500,4000,4500] )
        axtop.xaxis.set_minor_locator( ticker.MultipleLocator( 100 ) )
        ax.xaxis.set_major_locator( ticker.MultipleLocator( 500 ) )
        ax.xaxis.set_minor_locator( ticker.MultipleLocator( 100 ) )

        ax.text(0.05,0.88,label, color=color, ha='left', va='top', transform=ax.transAxes)
        ax.text(0.98,0.88,r'$\chi^2/\nu=$%.1f/%i'%(chi2,ndof), color=color, ha='right', va='top', transform=ax.transAxes)

        for tick in axtop.get_xaxis().get_major_ticks():
            tick.set_pad(2.)
            tick.label1 = tick._get_text1()

    ax1.yaxis.set_major_locator( ticker.MultipleLocator( 0.02 ) )
    ax1.yaxis.set_minor_locator( ticker.MultipleLocator( 0.01 ) )
    ax1.set_ylim(-0.009, 0.059)

    ax.set_xlabel('Observed Wavelength (\AA)')
    pl.draw()

def mkspecfig5panel(plotflam=False, showerr=False):
    """ make a 4-panel spec figure, showing the best-fit
    spectra with/without warping and with/without floating
    the redshift
    :param plotflam:
    :param showerr:
    :return:
    """
    from astropy.io import ascii
    import numpy as np
    from matplotlib import pyplot as pl, ticker
    from pytools import plotsetup
    import sys
    import os

    thisfile = sys.argv[0]
    if 'ipython' in thisfile:
        thisfile = __file__
    thisdir = os.path.dirname( thisfile )

    plotsetup.halfpaperfig( figsize=[4,8.5])
    specfile = thisdir+'/data/spectrum/snTomas_acs_g800l_1d.dat'
    wave, flux, err, sens = np.loadtxt(specfile, unpack=True)#, usecols=[0,1])

    flambda = flux / sens * 1e-17
    flambdaErr = err / sens * 1e-17

    templist = [thisdir+'/data/spectrum/sm-1.34570-sn2014j-visit3-hst.flm',
                thisdir+'/data/spectrum/sm-cont-1.31000-sn2012cg-visit1-hst.flm',
                thisdir+'/data/spectrum/sm-1.34570-sn2011fe-visit4-hst.flm',
                thisdir+'/data/spectrum/sm-1.31000-sn2011fe-visit4-hst.flm',
                thisdir+'/data/spectrum/sm-0.980000-sn2011iv-hst.flm',
                ]
    zlist = [ 1.3457, 1.31,  1.3457, 1.31, 0.98 ]
    labellist = ['SN 2014J (warped)\nz=1.3457 (fixed)\nage=-3',
                 'SN 2012cg (warped)\nz=1.31 (free)\nage=+3',
                 'SN 2011fe (unwarped)\nz=1.3457 (fixed)\nage=-3',
                 'SN 2011fe (unwarped)\nz=1.31 (free)\nage=-3',
                 'SN 2011iv (unwarped)\nz=0.98 (free)\nage=0',
                 ]
    chi2list = [ getchi2(tempfile) for tempfile in templist ]
    ndoflist = [ len(wave)-4, len(wave)-5, len(wave)-2, len(wave)-3, len(wave)-3 ]

    fig = pl.gcf()

    fig.subplots_adjust( left=0.03, right=0.97, top=0.95, bottom=0.08, hspace=0.5)
    pl.clf()
    nax = len( templist )
    iax = 0
    colorlist = ['r','g','r','g','b']
    for z, template, color, label, chi2, ndof in zip( zlist, templist, colorlist, labellist, chi2list, ndoflist) :
        iax += 1
        wavemod,fluxmod= np.loadtxt(template, unpack=True)
        flambdamod = fluxmod / sens * 1e-17

        if iax==1 :
            ax = fig.add_subplot(nax,1,iax)
            ax1 = ax
        else :
            ax = fig.add_subplot(nax,1,iax, sharey=ax1)

        #if iax==2 :
        #ax.set_ylabel('Flux')
        ax.yaxis.set_ticklabels([])

        ax.set_xlim(5800, 9550)

        if plotflam :
            pl.plot( wavemod, flambdamod, color=color, marker=' ' )
            pl.plot( wave, flambda, color='k', marker=' ', linestyle='steps-mid' )
            if showerr:
                pl.fill_between( wave, flambda-flambdaerr, flambda+flambdaerr, color='k', alpha=0.3 )
        else :
            pl.plot( wavemod, fluxmod, color=color, marker=' ' )
            pl.plot( wave, flux, color='k', marker=' ', linestyle='steps-mid' )
            if showerr:
                pl.fill_between( wave, flux-err, flux+err, color='k', alpha=0.3 )
        axtop = ax.twiny()
        axtop.set_xlim( ax.get_xlim()[0]/(1+z), ax.get_xlim()[1]/(1+z) )
        #axtop.set_xlabel('Rest Wavelength at z=%.2f'%z)
        axtop.set_xlabel('Rest $\lambda$:')
        axtop.xaxis.set_label_coords( 0.1, 1.05)
        #ax.set_xlabel('Obs:')
        #axtop.set_ylabel('Flux')
        axtop.xaxis.set_major_locator( ticker.MultipleLocator( 500 ) )
        axtop.xaxis.set_ticks( [3000,3500,4000] )
        if iax==4:
            axtop.xaxis.set_ticks( [3500,4000,4500] )
        axtop.xaxis.set_minor_locator( ticker.MultipleLocator( 100 ) )
        ax.xaxis.set_major_locator( ticker.MultipleLocator( 500 ) )
        ax.xaxis.set_minor_locator( ticker.MultipleLocator( 100 ) )

        ax.text(0.05,0.88,label, color=color, ha='left', va='top', transform=ax.transAxes)
        ax.text(0.98,0.88,r'$\chi^2/\nu=$%.1f/%i'%(chi2,ndof), color=color, ha='right', va='top', transform=ax.transAxes)

        for tick in axtop.get_xaxis().get_major_ticks():
            tick.set_pad(2.)
            tick.label1 = tick._get_text1()

    ax1.yaxis.set_major_locator( ticker.MultipleLocator( 0.02 ) )
    ax1.yaxis.set_minor_locator( ticker.MultipleLocator( 0.01 ) )
    ax1.set_ylim(-0.009, 0.059)

    ax.set_xlabel('Observed Wavelength (\AA)')
    pl.draw()





