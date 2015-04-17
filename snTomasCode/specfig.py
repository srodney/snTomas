
def getchi2(templatedatfile, obsdatfile='SPECTRUM/snTomas_acs_g800l_1d.dat'):
    import numpy as np
    from matplotlib import pyplot as pl, ticker

    wave, flux, err, sens = np.loadtxt(obsdatfile, unpack=True)#, usecols=[0,1])
    wavemod,fluxmod= np.loadtxt(templatedatfile, unpack=True)
    chi2 = np.sum( (fluxmod-flux)**2 / err**2 )
    return( chi2 )

def isigclip( valarray, sigclip, igood=[], maxiter=10, thisiter=0 ) :
    """ find the indices of valarray that
    survive after a clipping all values more than
    sigclip from the mean.  recursively iterative """
    import numpy as np
    if not type(valarray)==np.ndarray :
        valarray = np.array( valarray )
    if not len(igood) : igood = range(len(valarray))

    Ngood = len(igood)
    mnval = np.mean( valarray[igood] )
    sigma = np.std( valarray[igood] )
    igood = np.where( (np.abs(valarray-mnval)<(sigclip*sigma))  )[0]

    # import pdb; pdb.set_trace()
    if len(igood) == Ngood : return( igood )
    if thisiter>=maxiter :
        print("WARNING : Stopping after %i recursions"%maxiter)
        return( igood )
    thisiter+=1
    igood = isigclip( valarray, sigclip, igood=igood, maxiter=maxiter, thisiter=thisiter )
    return( igood )


def binspecdat( wavelength, flux, fluxerr=[], binwidth=10, sigclip=0, sumerrs=False,
                wstart=0, wend=0 ):
    """  bin up the given wavelength and flux arrays
    and return the binned values.
    binwidth is in the wavelength units of the wavelength
    array  (typically Angstroms)
    """
    import numpy as np

    w,f = wavelength, flux
    wbinned, fbinned = [], []
    wbin,fbin,dfbin = np.array([]), np.array([]), np.array([])
    dw, df = [], []
    if wstart : istart = np.where( w>wstart )[0][0]
    else : istart = 0
    if wend : iend = np.where( w<wend )[0][-1]
    else : iend = len(w)
    w0 = w[istart]
    for i in range(istart,iend):
        fullbin = False
        if wend and w[i]>wend : break
        if w[i]>w0+binwidth :
            # determine the mean value in this bin
            w0 = w[i]
            igoodval = []
            if sigclip :
                # use sigma clipping to reject outliers
                igoodval = isigclip( fbin, sigclip )
                if len(igoodval) :
                    wbinval = np.mean( wbin[igoodval] )
                    fbinval = np.mean( fbin[igoodval] )
                    dwbinval = (wbin[igoodval].max() - wbin[igoodval].min())/2.
                    #dwbinval = (wbin.max() - wbin.min())/2.
                    if sumerrs :
                        # flux uncertainty is the quadratic sum of the mean flux error
                        # and the error of the mean
                        dfbinval1 = np.std( fbin[igoodval] ) / np.sqrt(len(igoodval)-2)
                        dfbinval2 = np.mean( dfbin[igoodval] ) / np.sqrt(len(igoodval)-2)
                        dfbinval = np.sqrt( dfbinval1**2 + dfbinval2**2 )
                    else :
                        # flux uncertainty is the std error of the mean
                        dfbinval = np.std( fbin[igoodval] ) / np.sqrt(len(igoodval)-2)

                    fullbin = True
                # note: if the binning is not successful, we continue building the bin
            else :
                # use a straight median
                wbinval = np.median( wbin )
                fbinval = np.median( fbin )
                dwbinval = (wbin[-1]-wbin[0])/2.
                if sumerrs :
                    # flux uncertainty is the quadratic sum of the mean flux error
                    # and the error of the mean
                    dfbinval1 = np.std( fbin )/np.sqrt(len(fbin)-2)
                    dfbinval2 = np.mean( dfbin )
                    dfbinval = np.sqrt( dfbinval1**2 + dfbinval2**2 )
                else :
                    # flux uncertainty is the std error of the mean
                    dfbinval = np.std( fbin ) / np.sqrt(max(1,len(fbin)))
                fullbin = True

            if fullbin :
                wbinned.append( wbinval )
                fbinned.append( fbinval )
                dw.append( dwbinval )
                df.append( dfbinval )

                # start a new bin
                wbin,fbin,dfbin = np.array([]), np.array([]), np.array([])

        # add a new data point to the bin
        wbin = np.append( wbin, w[i] )
        fbin = np.append( fbin, f[i] )
        if len(fluxerr):
            dfbin = np.append( dfbin, fluxerr[i] )
        else : dfbin = np.append( dfbin, 0 )

    return( np.array( wbinned ), np.array(dw), np.array(fbinned), np.array(df) )



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

def mkspecfig4panel(plotflam=False, showerr=False):
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

    plotsetup.halfpaperfig( figsize=[4,7])
    specfile = thisdir+'/data/spectrum/snTomas_acs_g800l_1d.dat'
    wave, flux, err, sens = np.loadtxt(specfile, unpack=True)#, usecols=[0,1])

    flambda = flux / sens * 1e-17
    flambdaErr = err / sens * 1e-17

    templist = [thisdir+'/data/spectrum/sm-1.34570-sn2014j-visit3-hst.flm',
                thisdir+'/data/spectrum/sm-cont-1.31000-sn2012cg-visit1-hst.flm',
                thisdir+'/data/spectrum/sm-1.31000-sn2011fe-visit4-hst.flm',
                thisdir+'/data/spectrum/sm-0.980000-sn2011iv-hst.flm',

                ]
    zlist = [1.31, 1.3457, 1.31, 0.98 ]
    labellist = ['SN 2014J (warped) \nz=1.3457 (fixed)\nage=-3',
                 'SN 2012cg (warped)\nz=1.31 (free)\nage=+3',
                 'SN 2011fe\nz=1.31 (free)\nage=-3',
                 'SN 2011iv\nz=0.98 (free)\nage=0',
                 ]
    chi2list = [ getchi2(tempfile) for tempfile in templist ]
    ndoflist = [ len(wave)-5, len(wave)-6, len(wave)-3, len(wave)-3 ]

    fig = pl.gcf()

    fig.subplots_adjust( left=0.03, right=0.97, top=0.9, bottom=0.1, hspace=0.5)
    pl.clf()
    nax = len( templist )
    iax = 0
    colorlist = ['r','darkorange','g','b']
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





