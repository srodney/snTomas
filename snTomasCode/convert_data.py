def snana_to_sncosmo( snanafile='HST_FFSN_tomas_snana_ab.dat') :
    """ make the SN Tomas light curve fit figure"""
    import sncosmo
    import numpy as np
    from . import _ALPHA2FILTER

    snanadat = sncosmo.read_snana_ascii(snanafile)
    mjd = snanadat[1]['OBS']['MJD']
    band = snanadat[1]['OBS']['FLT']
    fluxcal = snanadat[1]['OBS']['FLUXCAL']
    fluxcalerr = snanadat[1]['OBS']['FLUXCALERR']
    mag  = snanadat[1]['OBS']['MAG']
    magerr = snanadat[1]['OBS']['MAGERR']
    zptab  = snanadat[1]['OBS']['ZPTAB']

    flux = fluxcal * 10**(-0.4*(27.5-zptab))
    fluxerr = fluxcalerr * 10**(-0.4*(27.5-zptab))

    magab = -2.5*np.log10( flux ) + zptab
    magaberr = 1.0857 * fluxerr / flux


    filter = [ _ALPHA2FILTER[b] for b in band ]

    fout = open('HST_FFSN_tomas.sncosmo.dat','w')
    print >> fout, "# mjd     filter    flux     fluxerr    mag     magerr  zpt    magsys"

    for i in range(len(mjd)) :
        if flux[i]>0 :
            print >> fout, "%8.2f  %6s  %8.3f %8.3f  %8.2f %6.2f  %7.3f   AB"%(
                mjd[i], filter[i], flux[i], fluxerr[i], magab[i], magaberr[i], zptab[i]
            )
        else :
            print >> fout, "%8.2f  %6s  %8.3f %8.3f  %8.2f %6.2f  %7.3f   AB"%(
                mjd[i], filter[i], flux[i], fluxerr[i],
                -2.5*np.log10(3*fluxerr[i]) + zptab[i], -9.0, zptab[i]
            )

    fout.close()

