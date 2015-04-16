



def plot_light_curve_fit( fitter='mlcs2k2', plotmags=True,
                          showblue=False, compare=False ) :
    """ make the SN Tomas light curve fit figure"""
    from pytools import plotsetup, colorpalette as cp
    from astropy.io import ascii
    import numpy as np
    from matplotlib import pyplot as pl
    from matplotlib import ticker, gridspec
    import os
    import sys
    from hstphot import hstzpt


    thisfile = sys.argv[0]
    if 'ipython' in thisfile :
        thisfile = __file__
    thisdir = os.path.abspath( os.path.dirname(thisfile))

    sn = ascii.read(os.path.join(thisdir,'data/HST_FFSN_tomas.dat'),
                    format='commented_header', header_start=-1, data_start=0)

    alpha2filter = { 'H':'F160W','N':'F140W','J':'F125W',
                     'Y':'F105W','I':'F814W','V':'F606W','B':'F435W' }
    colordict = { 'H':'k', 'N':cp.darkred,
                  'J':cp.darkgold, 'Y':cp.cadetblue,
                  'I':cp.darkgreen, 'V':cp.coral,
                  'B':cp.purple}

    if fitter=='both' :
        plotsetup.fullpaperfig([8,4])
        gs = gridspec.GridSpec( 2,5 )
        rowdict = {'salt2':1,'mlcs2k2':0}
        fitterlist=['mlcs2k2','salt2']
    else :
        plotsetup.fullpaperfig([8,4])
        gs = gridspec.GridSpec(1,5)
        rowdict = {fitter:0}
        fitterlist = [fitter]
    fig = pl.gcf()
    pl.clf()

    if showblue:
        bandlist = ['B','V','I','Y','J','N','H']
    else:
        bandlist = ['V','I','Y','J','N','H']
    for fitter in fitterlist :
        axdict = {
            'B':fig.add_subplot(gs[rowdict[fitter],0]),
            'Y':fig.add_subplot(gs[rowdict[fitter],1]),
            'J':fig.add_subplot(gs[rowdict[fitter],2]),
            'N':fig.add_subplot(gs[rowdict[fitter],3]),
            'H':fig.add_subplot(gs[rowdict[fitter],4]),
            }
        axdict['V'] = axdict['B']
        axdict['I'] = axdict['B']

        for band in bandlist:
            ax = axdict[band.upper()]
            filtname = alpha2filter[band]
            color = colordict[band]

            if band=='B':
                ax.text(0.85,0.45,'%s'%filtname.upper(), color=color,
                        fontsize='small',ha='right',va='top',
                        transform=ax.transAxes)
            if band=='V':
                ax.text(0.85,0.55,'%s'%filtname.upper(), color=color,
                        fontsize='small',ha='right',va='top',
                        transform=ax.transAxes)


            if plotmags :
                ifilt = np.where((sn['filter']==filtname) &
                                 (sn['magerr']>0))[0]
                y,yerr = sn['mag'][ifilt], sn['magerr'][ifilt]
            else :
                ifilt = np.where((sn['filter']==filtname))[0]
                f,ferr = sn['flux'][ifilt], sn['fluxerr'][ifilt]
                y = f*10**(-0.4*(sn['zpt'][ifilt]-25))
                yerr = ferr*10**(-0.4*(sn['zpt'][ifilt]-25))

            ax.errorbar( sn['mjd'][ifilt], y, yerr, color=color,
                         marker='D', ls=' ', zorder=-100)
            # ax.errorbar( mjd[iobs], y[iobs], yerr[iobs], color=color,
            #                     marker='D', ls=' ', zorder=-100 )


            fluxfile = os.path.join(
                thisdir,'data/%s_fit/HST_FFSN_tomas-1101-%s.flux'%(fitter,band) )
            if not os.path.isfile(fluxfile): continue

            fluxdat = ascii.read(fluxfile)

            mjd = fluxdat['col2']
            tobs = fluxdat['col3']
            fluxcal = fluxdat['col4']
            fluxcalerr = fluxdat['col5']

            if filtname.lower().startswith('f1'):
                zptvega = hstzpt.getzptWFC3IR(filtname,'vega')
                zptab = hstzpt.getzptWFC3IR(filtname,'ab')
            else :
                zptvega = hstzpt.getzptACS(filtname,'vega')
                zptab = hstzpt.getzptACS(filtname,'ab')
            flux        = fluxcal * 10**(-0.4*(27.5-zptvega))
            fluxab25    = flux * 10**(-0.4*(zptab-25))
            fluxerr     = fluxcalerr * 10**(-0.4*(27.5-zptvega))
            fluxerrab25 = fluxerr * 10**(-0.4*(zptab-25))

            modkey = fluxdat['col6']
            mag = np.where( fluxab25>0, -2.5*np.log10( fluxab25 ) + 25, 35 )
            magerr = np.where( fluxab25>0, 1.0857 * fluxerrab25 / fluxab25, 0.2 )

            if plotmags :
                y,yerr = mag, magerr
            else :
                y,yerr = fluxab25, fluxerrab25

            iobs = np.where( modkey>0 )[0]
            imod = np.where( modkey==0 )[0]

            ax.plot( mjd[imod], y[imod], color=color, marker=' ', ls='-', zorder=-100 )
            ax.fill_between( mjd[imod], y[imod]+yerr[imod], y[imod]-yerr[imod],
                             color=color, alpha=0.3, zorder=-1000 )

            ax.set_xlim( 56780,56950 )
            ax.xaxis.set_major_locator( ticker.MultipleLocator( 100 ) )
            ax.xaxis.set_minor_locator( ticker.MultipleLocator( 20 ) )
            if plotmags :
                ax.set_ylim( 27.9, 22.25 )
            else :
                ax.set_ylim( -0.2, 2.9 )
                ax.yaxis.set_major_locator( ticker.MultipleLocator( 1 ) )
                ax.yaxis.set_minor_locator( ticker.MultipleLocator( 0.5 ) )

            ax.text(0.85,0.65,'%s'%filtname.upper(), color=color,
                    fontsize='small',ha='right',va='top',
                    transform=ax.transAxes)
            if band=='H':
                ax.text(0.9,0.9,'%s'%fitter.upper(), color='k',
                        fontsize='large',ha='right',va='top',
                        transform=ax.transAxes)

            if rowdict[fitter]==0:
                axtop = ax.twiny()
                z=1.31
                mjdpk=56816.3
                axtop.set_xlim( (ax.get_xlim()[0]-mjdpk)/(1+z), (ax.get_xlim()[1]-mjdpk)/(1+z) )
                axtop.xaxis.set_major_locator( ticker.MultipleLocator( 20 ) )
                axtop.xaxis.set_minor_locator( ticker.MultipleLocator( 10 ) )
                axtop.xaxis.set_label_coords( -0.03, 1.01)
                axtop.tick_params( pad=2)
                if band=='I':
                    axtop.set_xlabel( "t$_{\\rm rest}$ : ")
            if ax != axdict['B']:
                pl.setp( ax.get_yticklabels(),visible=False)


        if plotmags :
            axdict['B'].set_ylabel( r"AB mag")#, fontsize='x-large')
        else :
            axdict['B'].set_ylabel( r"Flux (zp$_{\rm AB}$=25)")#, fontsize='x-large')

    # fig.suptitle( 'SN HFF14tom at $z=1.33\pm0.02$' )
    axdict['B'].set_xlabel( "MJD : ")#, fontsize='x-large')
    axdict['B'].xaxis.set_label_coords( -0.22, -0.06)
    #ax.yaxis.set_major_locator( ticker.MultipleLocator( 1 ) )
    #ax.yaxis.set_minor_locator( ticker.MultipleLocator( 0.5 ) )



    fig.subplots_adjust( left=0.08, bottom=0.12, right=0.98, top=0.91,
                         hspace=0, wspace=0)
    pl.draw()



def do_minchi2_fit( datfilename='HST_FFSN_tomas.dat', modelname='salt2', verbose=True ) :
    import sncosmo
    from astropy.io import ascii
    import numpy as np
    from pytools import cosmo
    import os
    import sys
    dm = cosmo.mu
    MBmodel = -19.223 # AB mags ( = -19.12 Vega )


    thisfile = sys.argv[0]
    if 'ipython' in thisfile :
        thisfile = __file__
    thisdir = os.path.dirname( os.path.abspath(thisfile))
    datfilename = os.path.join( thisdir, 'data/'+datfilename )
    tomas = ascii.read( datfilename, format='commented_header', header_start=-1, data_start=0 )

    if modelname.startswith('salt2'):
        model = sncosmo.Model( source=modelname)
        model.source.set_peakmag( 0., 'bessellb', 'ab' )
        x0_AB0 = model.get('x0')
        model.set( z=1.343, t0=56814.6 )
        res, fit = sncosmo.fit_lc( tomas, model, ['z','t0','x0','x1','c',],
                                   bounds={'z':(1.30,1.36),'t0':(56804,56824),
                                           'x1':(-5.,5.), 'c':(-0.3,1.5) })
        x0 = fit.get( 'x0' )
        z = fit.get( 'z' )
        MBobs = -2.5*np.log10(  x0 / x0_AB0 ) - dm(z)
        deltaM = MBobs - MBmodel
        distmod = 10**(-0.4*(deltaM))
        print( "MB_obs = %.2f"%MBobs )
        print( "delta(M) = %.2f"%deltaM )
        print( "dist.mod. = %.2f"%distmod )

    # any other model (2011fe, nugent-*, ccsn, etc)
    else :
        dust = sncosmo.CCM89Dust( )
        model = sncosmo.Model( source=modelname, effects=[dust], effect_names=['host'], effect_frames=['rest'] )
        model.set( z=1.343, t0=56814.6, hostr_v=3.1, hostebv=0.1 )

        res, fit = sncosmo.fit_lc( tomas, model, ['z','t0','amplitude','hostebv','hostr_v'],
                                   bounds={'z':(1.30,1.36),'t0':(56804,56824),
                                           'hostebv':(0.0,1.5), 'hostr_v':(1.0,5.0) })
    print( 'chi2/ndof = %.1f / %i = %.3f'%(res.chisq, res.ndof, res.chisq/res.ndof) )
    for parname in res.param_names :
        ipar = res.param_names.index( parname )
        parval = res.parameters[ ipar ]
        err = res.errors[parname]
        if verbose :
            if np.abs(err)>=0.1:
                print( '  %s =  %.2f +- %.2f'%( parname, np.round(parval,2), np.round(err,2))  )
            elif np.abs(err)>=0.01:
                print( '  %s =  %.3f +- %.3f'%( parname, np.round(parval,3), np.round(err,3)) )
            elif np.abs(err)>=0.001:
                print( '  %s =  %.4f +- %.4f'%( parname, np.round(parval,4), np.round(err,4)) )
            else :
                print( '  %s = %.3e +- %.3e'%( parname, parval, err) )


        if parname == 'x0' :
            salt2 = sncosmo.Model( source='salt2')
            salt2.source.set_peakmag( 0., 'bessellb', 'ab' )
            x0_AB0 = salt2.get('x0')
            mB = -2.5*np.log10(  parval / x0_AB0 )
            mBerr = 2.5*np.log10( np.e ) *  err / parval
            print( '  %s =  %.3f +- %.3f'%( 'mB', np.round(mB,3), np.round(mBerr,3)) )



    return( tomas, fit, res )


def sncosmoplot_minchi2fit( sn, fit, res ):
    import sncosmo
    sncosmo.plot_lc( sn, model=fit, errors=res.errors )

def plot_phot_compare():
    """ make a comparison plot showing aperture and psf photometry """
    from matplotlib import pyplot as pl
    import sncosmo
    import numpy as np
    from .import _ALPHA2FILTER

    datadir='/Users/rodney/Dropbox/Papers/snTomas/LCFIT/'
    ap = sncosmo.read_snana_ascii( datadir+'HST_FFSN_tomas.snana_apphot.dat')
    psf = sncosmo.read_snana_ascii(datadir+'HST_FFSN_tomas.snana_psfphot.dat')

    apdat = ap[1]['OBS']
    psfdat = psf[1]['OBS']

    fig = pl.gcf()
    pl.clf()
    iax = 0
    for band, color in zip( ['Y','J','N','H'],
                            ['c','r','g','m'] ):
        iax += 1
        fig.add_subplot( 2, 2, iax )
        for dat,marker,c in zip( [apdat,psfdat], ['o','D'], ['k',color]):
            iband = np.where( dat['FLT']==band )
            pl.errorbar( dat['MJD'][iband], dat['FLUXCAL'][iband],
                         dat['FLUXCALERR'][iband],
                         marker=marker, color=c, ms=10, alpha=0.5 )
            ax = pl.gca()
            pl.text(0.95,0.95,_ALPHA2FILTER[band],color=c,ha='right',va='top',transform=ax.transAxes)

    pl.draw()

def plot_lightcurve_snana():
    from matplotlib import pyplot as pl
    import sncosmo
    import numpy as np
    from .import _ALPHA2FILTER

    datadir='/Users/rodney/Dropbox/Papers/snTomas/LCFIT/'
    dat = sncosmo.read_snana_ascii( datadir+'HST_FFSN_tomas.snana.dat')
    obsdat = dat[1]['OBS']

    fig = pl.gcf()
    pl.clf()
    iax = 0
    marker='o'
    for band, color in zip( ['B','V','I','Y','J','N','H'],
                            ['c','g','r','c','g','r','k'] ):
        iax += 1
        fig.add_subplot( 2, 4, iax )
        iband = np.where( obsdat['FLT']==band )
        pl.errorbar( obsdat['MJD'][iband], obsdat['FLUXCAL'][iband],
                     obsdat['FLUXCALERR'][iband],
                     marker=marker, color=color, ms=10, alpha=0.5 )
        ax = pl.gca()
        pl.text(0.95,0.95,_ALPHA2FILTER[band],color=color,ha='right',va='top',transform=ax.transAxes)

    pl.draw()


def convert_light_curve_file( ):
    """ read in the SNANA-style light curve file (from D.Scolnic)
    and convert it into normal fluxes (not SNANA fluxcal) and
    AB and Vega mags
    :return:
    """
    import sncosmo
    import numpy as np
    from .import _ALPHA2FILTER
    from hstphot import hstzpt
    from astropy import table
    from astropy.utils import OrderedDict
    import sys
    import os

    thisfile = sys.argv[0]
    if 'ipython' in thisfile :
        thisfile = __file__
    thisdir = os.path.dirname( os.path.abspath(thisfile))
    datadir = os.path.join(thisdir,'data/')

    dat = sncosmo.read_snana_ascii( datadir+'HST_FFSN_tomas.snana.dat')
    obsdat = dat[1]['OBS']

    filtname= np.array([_ALPHA2FILTER[obsdat['FLT'][i]].lower()
                         for i in range(len(obsdat))])

    iIR = np.array([ i for i in range(len(obsdat))
                     if filtname[i].startswith('f1') ] )
    iACS = np.array([ i for i in range(len(obsdat))
                     if not filtname[i].startswith('f1') ] )
    zpt = np.ones(len(obsdat))
    zpt[iIR] = hstzpt.getzptWFC3IR(filtname[iIR],'AB')
    zpt[iACS] = hstzpt.getzptACS(filtname[iACS],'AB')
    zptsys = ['AB' for i in range(len(obsdat))]

    zptVega = np.ones(len(obsdat))
    zptVega[iIR] = hstzpt.getzptWFC3IR(filtname[iIR],'Vega')
    zptVega[iACS] = hstzpt.getzptACS(filtname[iACS],'Vega')

    fAB = obsdat['FLUXCAL'] * 10**(0.4*(27.5-zpt))
    ferrAB = obsdat['FLUXCALERR'] * 10**(0.4*(27.5-zpt))

    mAB = np.where( fAB>0, -2.5*np.log10(fAB)+zpt,
                    -2.5*np.log10(3*np.abs(ferrAB))+zpt,)
    merr = np.where( fAB>0, 1.0857*(ferrAB/fAB), -9*np.ones(len(fAB)))

    fcalVega = obsdat['FLUXCAL'] * 10**(0.4*(zpt-zptVega))
    fcalerrVega = obsdat['FLUXCALERR'] * 10**(0.4*(zpt-zptVega))
    mVega = np.where( fcalVega>0, -2.5*np.log10(fcalVega)+27.5,
                          -2.5*np.log10(3*np.abs(fcalerrVega))+27.5,)

    data=OrderedDict(
        [ ['mjd',obsdat['MJD']],
          ['filter',filtname],
          ['flux',fAB],
          ['fluxerr',ferrAB],
          ['mag',mAB],
          ['magerr',merr],
          ['zpt',zpt],
          ['zptsys',zptsys],
          ['fluxcalVega',fcalVega],
          ['fluxcalerrVega',fcalerrVega],
          ['magVega',mVega],
          ['zptVega',zptVega],
          ])
    outdat = table.Table(data=data)
    outdat['mjd'].format='%.1f'
    for colname in ['flux','fluxerr','fluxcalVega','fluxcalerrVega']:
        outdat[colname].format='%8.5f'
    for colname in ['mag','magerr','magVega']:
        outdat[colname].format='%8.3f'

    outdat.write( datadir+'HST_FFSN_tomas.dat',
                  format='ascii.fixed_width')
                  # format='ascii.commented_header')
    return(outdat)

