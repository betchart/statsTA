measurements = ['asymmetry', 'fraction']
partitions = ['full', 'hiM', 'loM', 'hiY', 'loY']


def measurement_pars(measure='asymmetry', partition='full'):
    fields = ('R0_', 'signal', 'profile')
    base = 3
    hemicycle = 5
    asymmetry = (base, 'fitTopQueuedBin5_TridiscriminantWTopQCD', ('falphaL', 'falphaT'))
    fraction = (base+hemicycle, 'fitTopTanhRapiditySum_triD', ('d_qq',))

    mtype = dict([(m, dict(zip(fields, eval(m)))) for m in measurements])
    pars = mtype[measure]
    N = pars['R0_']
    cycle = 2*hemicycle
    pDirs = [N,
             N + 1*cycle, (N, N + 1*cycle),
             N + 2*cycle, (N, N + 2*cycle), (N, N + 2*cycle), (N, N + 1*cycle)]
    pars.update({'R0_': dict(zip(partitions,pDirs))[partition]})
    if measure=='fraction': pars.update({'hackZeroBins':True})
    pars.update({'label':'_'.join([measure,partition])})
    return pars


def central():
    return {'d_lumi': 0,
            'd_xs_dy': 0,
            'd_xs_st': 0,
            'tag': 'ph_sn_jn_20',
            'genDirPre': 'R01',
            'genPre': '',
            'sigPre': '',
            'dirIncrement': 0,
            'd_wbb':0,
            'label': 'central'
            }


def systematics():
    sys =  ([
            {'label': 'JER_up', 'tag': 'ph_su_jn_20'},
            {'label': 'JER_dn', 'tag': 'ph_sd_jn_20'},

            {'label': 'JES_up', 'tag': 'ph_sn_ju_20'},
            {'label': 'JES_dn', 'tag': 'ph_sn_jd_20'},

            {'label': 'PU_up', 'dirIncrement': 1, 'sigPre': '001_'},
            {'label': 'PU_dn', 'dirIncrement': 1, 'sigPre': '000_'},

            {'label': 'PT', 'dirIncrement': 4, 'sigPre': '001_', 'genPre': '001_', 'genDirPre':'R02'},

            {'label': 'WBB_up', 'd_wbb':+0.2},
            {'label': 'WBB_dn', 'd_wbb':-0.2},

            {'label': 'lumi_up', "d_lumi": +0.044},
            {'label': 'lumi_dn', "d_lumi": -0.044},
            
            {'label': 'DY_up', "d_xs_dy": +0.20},
            {'label': 'DY_dn', "d_xs_dy": -0.20},
            
            {'label': 'ST_up', "d_xs_st": +0.20},
            {'label': 'ST_dn', "d_xs_st": -0.20},

            {'label': 'Q_dn', 'genPre': '053_', 'sigPre': '053_'},
            {'label': 'Q_up', 'genPre': '054_', 'sigPre': '054_'},

            {'label': 'as_dn', 'genPre': '055_', 'sigPre': '055_'},
            {'label': 'as_up', 'genPre': '056_', 'sigPre': '056_'},

            {'label': 'thr30', 'tag': 'ph_sn_jn_30'}
            ] +
            [{'label': 'el%d' % i,
              'dirIncrement': 2,
              'sigPre': '%03d_' % i,} for i in range(4)
            ]+
            [{'label': 'mu%d' % i,
              'dirIncrement': 3,
              'sigPre': '%03d_' % i,} for i in range(4)
            ]+
            [{'label': 'PD_%02d' % i,
              'genPre': '%03d_' % i,
              'sigPre': '%03d_' % i} for i in range(1, 53)
            ]
        )
    return sys


if __name__ == '__main__':
    sys = systematics()
    for s in sys:
        cen = central()
        cen.update(s)
        print cen
    print
    print central()
