measurements = ['asymmetry', 'fraction']
partitions = ['full', 'hiM', 'loM', 'hiY', 'loY']


def measurement_pars(measure='asymmetry', partition='full'):
    fields = ('R0_', 'signal', 'profile')
    base = 2
    hemicycle = 4
    asymmetry = (base, 'fitTopQueuedBin7TridiscriminantWTopQCD', ('falphaL', 'falphaT'))
    fraction = (base+hemicycle, 'fitTopTanhRapiditySum_triD', ('d_qq',))

    mtype = dict([(m, dict(zip(fields, eval(m)))) for m in measurements])
    pars = mtype[measure]
    N = pars['R0_']
    cycle = 2*hemicycle
    pDirs = [N,
             N + 1*cycle, (N, N + 1*cycle),
             N + 2*cycle, (N, N + 2*cycle)]
    pars.update({'R0_': dict(zip(partitions,pDirs))[partition]})
    if measure=='fraction' and partition=='hiM': pars.update({'hackZeroBins':True})
    return pars


def central():
    return {'d_lumi': 0,
            'd_xs_dy': 0,
            'd_xs_st': 0,
            'tag': 'ph_sn_jn_20',
            'genPre': '',
            'sigPre': '',
            'dirIncrement': 0,
            'label': 'central'
            }


def systematics():
    return ([
            {'label': 'upLumi', "d_lumi": +0.04},
            {'label': 'dnLumi', "d_lumi": -0.04},

            {'label': 'upDY', "d_xs_dy": +0.04},
            {'label': 'dnDY', "d_xs_dy": -0.04},

            {'label': 'upST', "d_xs_st": +0.04},
            {'label': 'dnST', "d_xs_st": -0.04},

            #{'label': 'upRFS', 'tag': 'up_sn_jn_20'},
            #{'label': 'dnRFS', 'tag': 'dn_sn_jn_20'},

            #{'label': 'upJER', 'tag': 'ph_su_jn_20'},
            #{'label': 'dnJER', 'tag': 'ph_sd_jn_20'},

            #{'label': 'upJES', 'tag': 'ph_sn_ju_20'},
            #{'label': 'dnJES', 'tag': 'ph_sn_jd_20'},

            {'label': 'upPU', 'dirIncrement': 1, 'sigPre': '001_'},
            {'label': 'dnPU', 'dirIncrement': 1, 'sigPre': '000_'},

            ] +
            [{'label': '%del' % i,
              'dirIncrement': 2,
              'sigPre': '%03d_' % i,} for i in range(4)
            ]+
            [{'label': '%dmu' % i,
              'dirIncrement': 3,
              'sigPre': '%03d_' % i,} for i in range(4)
            ]+
            [{'label': '%02dPDF' % i,
              'genPre': '%03d_' % i,
              'sigPre': '%03d_' % i} for i in range(1, 53)]
            )


if __name__ == '__main__':
    sys = systematics()
    for s in sys:
        cen = central()
        cen.update(s)
        print cen
    print
    print central()
