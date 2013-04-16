measurements = ['asymmetry', 'fraction']
partitions = ['full', 'hiM', 'loM', 'hiY', 'loY','loYalt','loMalt']


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
             N + 2*cycle, (N, N + 2*cycle), (N, N + 2*cycle), (N, N + 1*cycle)]
    pars.update({'R0_': dict(zip(partitions,pDirs))[partition]})
    if measure=='fraction' and partition=='hiM': pars.update({'hackZeroBins':True})
    if measure=='asymmetry' and 'alt'==partition[-3:]: pars.update({'alternateModel':True})
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
    sys =  ([
            {'label': 'lumi_up', "d_lumi": +0.04},
            {'label': 'lumi_dn', "d_lumi": -0.04},

            {'label': 'DY_up', "d_xs_dy": +0.04},
            {'label': 'DY_dn', "d_xs_dy": -0.04},

            {'label': 'ST_up', "d_xs_st": +0.04},
            {'label': 'ST_dn', "d_xs_st": -0.04},

            #{'label': 'RFS_up', 'tag': 'up_sn_jn_20'},
            #{'label': 'RFS_dn', 'tag': 'dn_sn_jn_20'},

            #{'label': 'JER_up', 'tag': 'ph_su_jn_20'},
            #{'label': 'JER_dn', 'tag': 'ph_sd_jn_20'},

            #{'label': 'JES_up', 'tag': 'ph_sn_ju_20'},
            #{'label': 'JES_dn', 'tag': 'ph_sn_jd_20'},

            {'label': 'PU_up', 'dirIncrement': 1, 'sigPre': '001_'},
            {'label': 'PU_dn', 'dirIncrement': 1, 'sigPre': '000_'},

            ] +
            [{'label': 'el%d' % i,
              'dirIncrement': 2,
              'sigPre': '%03d_' % i,} for i in range(4)
            ]+
            [{'label': 'mu%d' % i,
              'dirIncrement': 3,
              'sigPre': '%03d_' % i,} for i in range(4)
            ]+
            [{'label': 'PDF%02d' % i,
              'genPre': '%03d_' % i,
              'sigPre': '%03d_' % i} for i in range(1, 53)]
            )
    #return [s for s in sys if s['label']=='PDF33']
    return sys


if __name__ == '__main__':
    sys = systematics()
    for s in sys:
        cen = central()
        cen.update(s)
        print cen
    print
    print central()
