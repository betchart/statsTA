

def ensemble_specs():
    base = {'lumiFactor':1.0, 'Nens':1000}

    ups = [{'label':'A','letter':'A'},
           {'label':'B','letter':'B'},
           {'label':'C','letter':'C'},
           {'label':'D','letter':'D'},
           {'label':'A4','letter':'A','lumiFactor':0.25, 'Nens':100},
           {'label':'A2','letter':'A','lumiFactor':0.50, 'Nens':100},
           {'label':'2A','letter':'A','lumiFactor':2.00, 'Nens':100},
           {'label':'4A','letter':'A','lumiFactor':4.00, 'Nens':100},
           {'label':'8A','letter':'A','lumiFactor':8.00, 'Nens':100},
    ]
    specs = [dict(base) for u in ups]
    for s,u in zip(specs,ups): s.update(u)
    return specs


def calibration_specs():
    return [{'which':'mg', 'sample':'calib_mg.pu.sf'},
            {'which':'mn', 'sample':'calib_mn.pu.sf'},
            {'which':'ZP', 'sample':'calib_ZP.pu.sf'},
            {'which':'A2K', 'sample':'calib_A2K.pu.sf'},
            {'which':'R2K', 'sample':'calib_R2K.pu.sf'},
            {'which':'A.2K', 'sample':'calib_A200.pu.sf'},
            {'which':'R.2K', 'sample':'calib_R200.pu.sf'},
            {'which':'L.2K', 'sample':'calib_L200.pu.sf'}]
