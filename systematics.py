

def central():
    return {'d_lumi': 0,
            'd_xs_dy': 0,
            'd_xs_st': 0,
            'tag': 'ph_sn_jn_20',
            'genPre': '',
            'sigPre': '',
            'dirPreInc': 0,
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

            {'label': 'upRFS', 'tag': 'up_sn_jn_20'},
            {'label': 'dnRFS', 'tag': 'dn_sn_jn_20'},

            {'label': 'upJER', 'tag': 'ph_su_jn_20'},
            {'label': 'dnJER', 'tag': 'ph_sd_jn_20'},

            {'label': 'upJES', 'tag': 'ph_sn_ju_20'},
            {'label': 'dnJES', 'tag': 'ph_sn_jd_20'},

            {'label': 'upPU', 'dirPreInc': 1, 'sigPre': '001'},
            {'label': 'dnPU', 'dirPreInc': 1, 'sigPre': '000'},

            ] +
            [{'label': '%02dPDF' % i,
              'genPre': '%03d' % i,
              'sigPre': '%03d' % i} for i in range(1, 53)]
            )


if __name__ == '__main__':
    sys = systematics()
    for s in sys:
        cen = central()
        cen.update(s)
        print cen
    print
    print central()
